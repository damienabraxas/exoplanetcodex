"""Reads one element's run report + product feed + litscan reference; emits a verdict and a problem-line list; decides no abundance.

RYA-1234 -- governing process steps 12 (check against known literature values) and 13
(report problem lines). It changes no A value and no line selection (RYA-161 firewall),
never combines engines (RYA-712), and holds no tolerance of its own: every band comes
from the element's litscan YAML through `pipeline.litscan`.

STEP 0 TRUST AUDIT (RYA-1234), and what it decided here
-------------------------------------------------------
  litscan.py           TRUST    pure reader of data/reference/litscan/<El>.yaml; the band,
                                its scale and the DEVIATE bound all come from the YAML.
                                Gap: the path carries no STAR, so a litscan is the solar
                                literature -- used as a pass band only for a benchmark star.
  validate_element.py  ADAPT    its judging rule (scale mismatch is flagged, never silently
                                compared; VIS validates, frontier bands report) is reused
                                through `VALIDATING_BANDS` and the same scale test; its
                                BandProduct has no instrument / holding / tier, so it cannot
                                take the current cell key and is not called directly.
  problem_children.py  REPLACE  keys a problem on (species, lambda-or-scope, class): no EP,
                                band or instrument (RYA-1037: lambda alone is never a line
                                key), and it aggregates the pre-1000 *_ew_integrity.csv.
                                Problem lines go to <star>_problem_lines.csv instead.
  rejection_ledger.py  REPLACE  accounts for the pre-1000 lines_fit EW fitter only; the
                                current per-line products carry status + reason directly.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

from pipeline import run_matrix as rm

ROOT = rm.ROOT

# ── statuses (RYA-1234 spec) ─────────────────────────────────────────────────
IN_BAND, OUT_OF_BAND, NO_REFERENCE, REPORTED = "IN_BAND", "OUT_OF_BAND", "NO_REFERENCE", "REPORTED"
PASS, REVIEW, INCOMPLETE = "PASS", "REVIEW", "INCOMPLETE"

#: RYA-1234: problem_children is REPLACE, so the problem-line reasons are this small
#: vocabulary, mapped from the per-line products' own exclusion codes (the text before
#: the first ':' of `reason_note` / `excluded_reason`).
REASON_CODES = {
    "FIT-NOT-PHYSICAL": "FIT_INVALID", "NON-MINIMUM": "FIT_INVALID",
    "CHI2-GATE": "FIT_INVALID", "SYNTHESIS": "FIT_INVALID", "EDGE-PINNED": "FIT_INVALID",
    "OVER-PHYSICAL-WIDTH": "FIT_INVALID",
    "QUARANTINED-TELLURIC": "TELLURIC_SATURATED", "TELLURIC-SATURATED": "TELLURIC_SATURATED",
    "ENGINE-A-NOT-SERVED": "NOT_SERVED", "NO-NLTE-LABEL": "NOT_SERVED",
}
OTHER = "EXCLUDED_OTHER"


def reason_code(text: str) -> str:
    head = (text or "").split(":", 1)[0].strip().upper()
    if head.startswith("REW") and "SATURATION" in (text or "").upper():
        return "EW_SATURATED"
    return REASON_CODES.get(head, OTHER)


def is_benchmark(star: str) -> tuple[bool, str]:
    """RYA-872 benchmark fork: `is_gbs` true, or the calibration anchor (the Sun)."""
    from pipeline.system_catalog import load_catalog
    for row in load_catalog():
        if row.get("star_params_key") == star:
            gbs = str(row.get("is_gbs", "")).strip().lower() == "true"
            anchor = str(row.get("role", "")).strip() == "calibration_anchor"
            return gbs or anchor, (f"system_catalog: is_gbs={row.get('is_gbs')} "
                                   f"role={row.get('role')}")
    return False, f"{star} is not in system_catalog.csv -- treated as non-benchmark"


def literature_for(symbol: str, ion: str):
    """(LiteratureRange | None, why). The litscan is per element AND ion (Fe II is not
    folded into the Fe I range -- the YAML says so)."""
    from pipeline.litscan import load_litscan, literature_range, litscan_path
    doc = load_litscan(symbol)
    if doc is None:
        return None, f"no {litscan_path(symbol).relative_to(ROOT)}"
    ions = [str(i).strip() for i in (doc.get("ions") or [doc.get("ion", "I")])]
    if ion not in ions:
        return None, f"the {symbol} litscan is for {symbol} {'/'.join(ions)}, not {symbol} {ion}"
    lit = literature_range(symbol)
    return lit, ("" if lit else f"{symbol} litscan declares no `range`")


def _scale(product: dict) -> str:
    from pipeline.treatment_axes import axes_for
    try:
        return axes_for(product["treatment"], route_token=product.get("route")).scale
    except Exception as exc:                                   # noqa: BLE001
        print(f"WARNING: no scale for treatment {product.get('treatment')!r}: "
              f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return ""


def _cell_products(feed: dict, symbol: str, cell: dict) -> list[dict]:
    emits = set(rm.DECK_EMITS.get(cell["engine"], ()))
    return [p for p in feed.get("products", [])
            if str(p.get("element")) == symbol and str(p.get("ion")) == cell["ion"]
            and str(p.get("band")) == cell["band"]
            and str(p.get("instrument")) == cell["instrument"]
            and str(p.get("holding")) == cell["holding"]
            and str(p.get("route")) == cell.get("route", "")
            and str(p.get("selector")) == cell.get("pool", "")
            and str(p.get("treatment")) in emits]


def judge(product: dict, lit, why_no_lit: str, benchmark: bool) -> dict:
    from pipeline.validate_element import VALIDATING_BANDS
    A = product.get("A")
    out = {"treatment": product.get("treatment"), "A": A, "scale": _scale(product),
           "tier": "validate" if product.get("band") in VALIDATING_BANDS else "report"}
    if not benchmark:
        out.update(status=REPORTED, reason="non-benchmark star: reported, no pass/fail",
                   literature=(None if lit is None else
                               {"central": lit.central, "min": lit.min, "max": lit.max}))
        return out
    if lit is None or not isinstance(A, (int, float)):
        out.update(status=NO_REFERENCE,
                   reason=why_no_lit or "the product carries no numeric A")
        return out
    out.update(literature={"central": lit.central, "min": lit.min, "max": lit.max,
                           "scale": lit.scale},
               offset_vs_central=round(lit.offset(A), 4),
               deviate=lit.is_deviate(A))
    # validate_element's own rule, kept: a scale mismatch is FLAGGED, never hidden.
    if lit.scale and out["scale"] and out["scale"] != lit.scale:
        out["scale_mismatch"] = (f"product on {out['scale']}, literature on {lit.scale}: "
                                 f"not scale-consistent until moved onto {lit.scale}")
    if lit.contains(A):
        out.update(status=IN_BAND, delta=0.0)
    else:
        edge = lit.max if A > lit.max else lit.min
        out.update(status=OUT_OF_BAND, delta=round(A - edge, 4))
    return out


def problem_lines(star: str, symbol: str, judged: list[dict], cells: list[dict]) -> list[dict]:
    """Step 13. Excluded lines of every judged product (from the published per-line
    product) and of every cell this run built (from its recorded *_lines.csv). Identity
    is lambda + EP, never lambda alone."""
    rows: list[dict] = []
    keys = {(j["band"], j["instrument"], j["holding"], j["route"], j["pool"], j["treatment"])
            for j in judged}
    perline = ROOT / "data" / "products" / star / f"{symbol}_perline.csv"
    if perline.exists() and keys:
        with perline.open(encoding="utf-8") as fh:
            for r in csv.DictReader(l for l in fh if not l.startswith("#")):
                if r.get("status") != "excluded":
                    continue
                k = (r.get("band"), r.get("instrument"), r.get("holding"), r.get("route"),
                     r.get("selector"), r.get("treatment"))
                if k not in keys:
                    continue
                note = r.get("reason_note") or r.get("reason_code") or ""
                rows.append({"element": symbol, "ion": r.get("ion"),
                             "wavelength_air_A": r.get("wavelength_air_A"),
                             "ep_eV": r.get("excitation_potential_eV"),
                             "band": r.get("band"), "instrument": r.get("instrument"),
                             "holding": r.get("holding"), "engine": r.get("treatment"),
                             "reason_code": reason_code(note), "reason": note[:200],
                             "source": perline.name})
    led = rm.load_ledger()
    for key, entry in led.items():
        parts = key.split("|")
        if len(parts) < 4 or parts[0] != star or parts[1] != symbol:
            continue
        for a in entry.get("artifacts", []):
            if not a["path"].endswith("_lines.csv"):
                continue
            path = rm.ARTIFACT_ROOT / a["path"]
            if not path.exists():
                continue
            with path.open(encoding="utf-8") as fh:
                for r in csv.DictReader(fh):
                    if str(r.get("in_aggregate")) == "True":
                        continue
                    note = r.get("excluded_reason") or ""
                    rows.append({"element": symbol, "ion": r.get("ion"),
                                 "wavelength_air_A": r.get("wavelength_air_A"),
                                 "ep_eV": r.get("ep_eV"), "band": parts[3],
                                 "instrument": r.get("instrument"), "holding": parts[5]
                                 if len(parts) > 5 else "",
                                 "engine": r.get("treatment"),
                                 "reason_code": reason_code(note), "reason": note[:200],
                                 "source": a["path"]})
    seen, out = set(), []
    for r in rows:
        k = (star, r["element"], r["wavelength_air_A"], r["ep_eV"], r["band"],
             r["instrument"], r["engine"], r["reason_code"])
        if k not in seen:
            seen.add(k)
            out.append(r)
    return out


def verdict(star: str, element: str, *, report_dir: Path | None = None,
            write: bool = True) -> dict:
    """The element's verdict + problem lines, written INTO its run report (latest and
    the stamped twin) and appended to <star>_problem_lines.csv. Decides no abundance."""
    dest = report_dir or rm.REPORT_DIR
    safe = element.replace(" ", "")
    latest = dest / f"{star}_{safe}_latest.json"
    doc = json.loads(latest.read_text(encoding="utf-8"))
    symbol, _ = rm.split_symbol(element)
    feed = rm.load_feed(star, symbol)
    bench, bench_why = is_benchmark(star)

    judged: list[dict] = []
    for c in doc["cells"]:
        prods = _cell_products(feed, symbol, c)
        if not prods:
            continue
        lit, why = literature_for(symbol, c["ion"])
        for p in prods:
            j = judge(p, lit, why, bench)
            j.update(band=c["band"], instrument=c["instrument"], holding=c["holding"],
                     ion=c["ion"], engine=c["engine"], route=c.get("route", ""),
                     pool=c.get("pool", ""), cell_status=c["status"])
            judged.append(j)

    n_failed = sum(1 for c in doc["cells"] if c["status"] == rm.FAILED)
    if not bench:
        ev = REPORTED
    elif not judged:
        ev = INCOMPLETE
    elif n_failed or any(j["status"] == OUT_OF_BAND for j in judged):
        ev = REVIEW
    elif not any(j["status"] == IN_BAND for j in judged):
        # Process step 12 is a CHECK against literature. Every product NO_REFERENCE means
        # nothing was checked -- that is not a pass (Fe II, whose litscan does not exist,
        # read PASS on zero comparisons before this).
        ev = INCOMPLETE
    else:
        ev = PASS
    v = {"element_verdict": ev, "benchmark": bench, "benchmark_basis": bench_why,
         "n_products_judged": len(judged), "n_failed_cells": n_failed,
         "counts": {s: sum(j["status"] == s for j in judged)
                    for s in (IN_BAND, OUT_OF_BAND, NO_REFERENCE, REPORTED)},
         "n_scale_mismatch": sum(1 for j in judged if j.get("scale_mismatch")),
         "cells": judged}
    probs = problem_lines(star, symbol, judged, doc["cells"])
    if write:
        doc["verdict"], doc["problem_lines"] = v, probs
        text = json.dumps(doc, indent=2) + "\n"
        latest.write_text(text, encoding="utf-8")
        stamp = doc["generated_at"].replace(":", "").replace("-", "")
        twin = dest / f"{star}_{safe}_{stamp}.json"
        if twin.exists():
            twin.write_text(text, encoding="utf-8")
        _append_problem_csv(dest / f"{star}_problem_lines.csv", star, probs)
    return {"verdict": v, "problem_lines": probs}


_PROBLEM_COLS = ("star", "element", "ion", "wavelength_air_A", "ep_eV", "band",
                 "instrument", "holding", "engine", "reason_code", "reason", "source")


def _append_problem_csv(path: Path, star: str, rows: list[dict]) -> None:
    """Dedupe on (star, element, lambda, EP, band, instrument, engine, reason_code)."""
    have: dict = {}
    if path.exists():
        with path.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                have[tuple(r[k] for k in _PROBLEM_COLS[:7]) + (r["engine"], r["reason_code"])] = r
    for r in rows:
        r = {"star": star, **{k: r.get(k, "") for k in _PROBLEM_COLS if k != "star"}}
        have[tuple(str(r[k]) for k in _PROBLEM_COLS[:7]) + (str(r["engine"]),
                                                             str(r["reason_code"]))] = r
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_PROBLEM_COLS)
        w.writeheader()
        for r in have.values():
            w.writerow({k: r.get(k, "") for k in _PROBLEM_COLS})
