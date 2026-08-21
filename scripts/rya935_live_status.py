#!/usr/bin/env python3
"""Derive the live tracker's status from products on disk (RYA-935).

The tracker's JSON was HAND-TYPED. Two things follow from that, and both had
already happened by the time this was written:

  * it violated RYA-686 -- a result artifact landed with no generating harness,
    which is the RYA-559 hole that convention exists to close, and CI said so;
  * it went stale silently. The committed copy listed RYA-931 as BACKLOG on the
    day RYA-931 merged, and showed HARPS as "correction owed" after the
    corrected holding existed. A dashboard that is typed is a dashboard that
    lies eventually.

So this derives every cell it can from artifacts that already exist, and is
element-agnostic on purpose: Fe appears the moment its products land, with no
edit here.

What it deliberately does NOT emit: the ticket pipeline. That is Linear state,
not repository state, and a copy of it in a committed file is guaranteed to
drift -- which is exactly the failure above. The dashboard shows what the repo
can prove.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

#: `<El><Ion>_<lo>_<hi>_<instrument>_<holding>_<HANDLER>_products.csv`. The holding
#: is in the stem because RYA-933/934 put it there -- before that, two holdings of
#: one instrument wrote the same filename and the second overwrote the first.
STEM = re.compile(r"^(?P<el>[A-Z][a-z]?)(?P<ion>I+|IV|VI*)_(?P<lo>\d+)_(?P<hi>\d+)_"
                  r"(?P<rest>.+?)_(?P<handler>PROFILEFIT|SYNTH)_products\.csv$")


def parse_stem(name: str, instruments: set[str], holdings: set[str]) -> dict | None:
    m = STEM.match(name)
    if not m:
        return None
    rest = m.group("rest")
    instrument = next((i for i in sorted(instruments, key=len, reverse=True)
                       if rest == i or rest.startswith(i + "_")), None)
    if instrument is None:
        return None
    tail = rest[len(instrument):].lstrip("_")
    # No holding in the stem = a product from BEFORE RYA-933/934, when the stem
    # keyed on instrument alone. Say so; do not guess. Every product committed
    # before that change is attributable to an instrument and NOT to a holding,
    # and an instrument can serve several -- including a corrected and an
    # uncorrected one. Inferring which would be exactly the collapse the stem
    # change was made to prevent.
    if tail and tail in holdings:
        holding, source = tail, "filename"
    elif tail:
        holding, source = tail, "filename (unregistered holding)"
    else:
        holding, source = None, "absent -- product predates RYA-933/934"
    return {"element": m.group("el"), "ion": m.group("ion"),
            "lo_A": float(m.group("lo")), "hi_A": float(m.group("hi")),
            "instrument": instrument, "holding": holding,
            "holding_source": source, "handler": m.group("handler")}


def collect_products(roots: list[Path], instruments: set[str],
                     holdings: set[str]) -> list[dict]:
    rows: list[dict] = []
    for root in roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*_products.csv")):
            meta = parse_stem(path.name, instruments, holdings)
            if meta is None:
                continue
            try:
                frame = pd.read_csv(path)
            except Exception:                                   # noqa: BLE001
                continue
            for _, r in frame.iterrows():
                rows.append({
                    **meta,
                    "band": r.get("band"), "treatment": r.get("treatment"),
                    "display": f"{r.get('route', '')} · {r.get('scale', '')}".strip(" ·"),
                    "A": None if pd.isna(r.get("A")) else float(r["A"]),
                    "sigma_stat": None if pd.isna(r.get("stat_dex")) else float(r["stat_dex"]),
                    "sigma_syst": None if pd.isna(r.get("syst_dex")) else float(r["syst_dex"]),
                    "n_lines": None if pd.isna(r.get("n_lines")) else int(r["n_lines"]),
                    "source": str(path.relative_to(ROOT)),
                })
    return rows


def collect_instruments() -> list[dict]:
    """Coverage and telluric state, read from the registries -- never retyped."""
    import measure_band_ew as M
    from pipeline.telluric_policy import applied_state
    out = []
    for instrument, specs in M._INSTRUMENT_HOLDINGS.items():
        for spec in specs:
            try:
                telluric = applied_state(spec.holding_id)
            except KeyError:
                telluric = "unregistered"
            out.append({
                "instrument": instrument, "holding": spec.holding_id,
                "telluric_applied": telluric,
                "pre_normalised": spec.pre_normalised,
                "coverage_A": list(spec.span_A) if spec.span_A else None,
                "coverage_note": ("declared" if spec.span_A else
                                  "discovered by the reader; not declarable as one "
                                  "interval"),
            })
    return out


def collect_reference(root: Path) -> dict:
    """Literature anchor per species, from the FROZEN gold reference.

    Read, not typed. The previous hand-written page carried Al 6.43 as a literal;
    it is right, but a literal cannot follow the reference when it is re-frozen,
    and this project's whole discipline is that a value cites its source. The
    pointer file `data/reference/solar/CURRENT` names the live version.
    """
    current = (root / "data" / "reference" / "solar" / "CURRENT")
    version = current.read_text().strip() if current.exists() else "v5"
    table = root / "data" / "reference" / "solar" / f"solar_abundances_{version}.csv"
    if not table.exists():
        return {}
    frame = pd.read_csv(table, comment="#")
    out = {}
    for _, r in frame.iterrows():
        if pd.isna(r.get("asplund2021")):
            continue
        out[f"{r['element']}{r['ion']}"] = {
            "asplund2021": float(r["asplund2021"]),
            "codex_A_X": None if pd.isna(r.get("A_X")) else float(r["A_X"]),
            "verdict": str(r.get("verdict")),
            "source": f"data/reference/solar/solar_abundances_{version}.csv",
            "sigma_external": None, "band": None, "best_external": None,
            "scale": None, "deviate_beyond": None,
        }

    # The gold table carries the literature VALUE but no uncertainty. litscan does,
    # and it is the ratified comparator: best-external +/- sigma_external, with the
    # source named. Take the band from there wherever an element has a litscan.
    #
    # NOTE the band is AGREEMENT WITH THE LITERATURE, not a pass/fail gate. litscan's
    # own basis text warns against conflating it with the FE_GATE policy window
    # ([7.41, 7.51], RYA-166) -- they answer different questions.
    try:
        from pipeline import litscan
        for element in litscan.available_elements():
            rng = litscan.literature_range(element)
            if rng is None or rng.sigma_external is None:
                continue
            for key in [k for k in out if k.startswith(element)
                        and k[len(element):].strip("IV") == ""]:
                out[key].update({
                    "asplund2021": rng.central,
                    "sigma_external": rng.sigma_external,
                    "band": [rng.min, rng.max],
                    "deviate_beyond": rng.deviate_beyond,
                    "best_external": rng.best_external,
                    "scale": rng.scale,
                    "band_meaning": ("agreement with the literature (best external +/- "
                                     "sigma_ext) -- NOT a pass/fail gate"),
                    "source": f"pipeline/litscan.py :: {element}.yaml",
                })
    except Exception:                                           # noqa: BLE001
        pass
    return out


def why_no_product(element: str, ion: str, holding: str, instrument: str,
                   lo: float, hi: float) -> str:
    """Why this cell is empty -- from the SAME resolver that plans the runs.

    An empty cell with no explanation is the RYA-833 shape: "we do not hold this"
    becomes indistinguishable from "nobody looked". Every blank on this grid has
    to say which it is.
    """
    from pipeline.run_descriptor import RunDescriptor, resolve
    d = RunDescriptor(element, ion, instrument, holding, lo, hi)
    # Interpreter and engine dir are supplied so that only REAL blockers -- coverage,
    # wiring, the telluric gate -- surface here. Whether a run host happens to have
    # the right numpy is not a fact about the science and does not belong on the grid.
    r = resolve(d, interpreter="(host)", ispec_dir="(host)")
    if r.blocked_reason:
        return r.blocked_reason
    return "no run yet"


def collect_telluric(audit_root: Path) -> list[dict]:
    """Before/after residuals from the correction tickets' own evidence."""
    out = []
    for manifest in sorted(audit_root.glob("rya940_kp1984_telluric/*/fit_manifest.json")):
        d = json.loads(manifest.read_text())
        c = d.get("correction")
        if not c:
            out.append({"window": f"{int(d['band_A'][0])}-{int(d['band_A'][1])} A",
                        "product": "solar_kpno_molecfit_corrected",
                        "before_pct_below_0.5": None, "after_pct_below_0.5": None,
                        "note": "no admissible fit; band NOT corrected"})
            continue
        out.append({
            "window": f"{int(d['band_A'][0])}-{int(d['band_A'][1])} A",
            "product": "solar_kpno_molecfit_corrected",
            "before_pct_below_0.5": c["before"]["pct_below_0.5"],
            "after_pct_below_0.5": c["after"]["pct_below_0.5"],
            "externally_validated": d.get("externally_validated"),
        })
    verification = audit_root / "rya931_molecfit_runtime" / "verification.json"
    if verification.exists():
        v = json.loads(verification.read_text())["o2b_gate"]
        out.append({"window": "6867-6884 A", "product": "solar_harps_molecfit_corrected",
                    "before_pct_below_0.5": v["o2b_before"]["pct_below_0.5"],
                    "after_pct_below_0.5": v["o2b_after"]["pct_below_0.5"],
                    "externally_validated": True})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--products-root", type=Path, action="append", default=None)
    ap.add_argument("--audit-root", type=Path, default=ROOT / "data" / "audit")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "results" / "rya935"
                    / "live_status.json")
    ap.add_argument("--refresh-seconds", type=int, default=5)
    args = ap.parse_args()

    import measure_band_ew as M
    instruments = set(M._INSTRUMENT_HOLDINGS)
    holdings = {s.holding_id for specs in M._INSTRUMENT_HOLDINGS.values() for s in specs}
    roots = args.products_root or [ROOT / "data" / "results"]

    products = collect_products(roots, instruments, holdings)
    status = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "scripts/rya935_live_status.py",
        "refresh_seconds": args.refresh_seconds,
        "derivation_note": ("Every value here is read from a product or a registry. "
                            "Nothing is typed. The ticket pipeline is deliberately "
                            "absent: it is Linear state, and a committed copy drifts."),
        # EVERY species in the frozen gold reference, not merely those with products.
        # This page is the progress framework for the whole solar calibration and
        # then for every star after it, so early on the absences ARE the content --
        # listing only what is finished would hide the work that remains.
        "elements": None,   # filled below, once the reference is read
        "elements_with_products": sorted({p["element"] + p["ion"] for p in products}),
        "bands": ["near-UV", "VIS", "red-optical", "NIR"],
        "instruments": collect_instruments(),
        "products": products,
        "telluric": collect_telluric(args.audit_root),
        "reference": collect_reference(ROOT),
    }
    status["elements"] = sorted(status["reference"])
    status["system"] = "solar"   # the framework is per-star; this build is the Sun

    # WHY a cell is empty. Computed once per (holding, band) and NOT per element,
    # because every blocker the resolver reports -- wiring, coverage, the telluric
    # gate -- is a property of the holding and the band. Recomputing it per element
    # would be 26x the work for identical answers, and would invite someone to read
    # element-specific meaning into a reason that has none.
    from pipeline import band_policy
    reachability = {}
    for inst in status["instruments"]:
        for policy in band_policy.POLICIES:
            try:
                reason = why_no_product("Fe", "I", inst["holding"], inst["instrument"],
                                        policy.lo_A, policy.hi_A)
            except Exception as exc:                            # noqa: BLE001
                reason = f"could not resolve: {exc}"
            reachability[f"{inst['holding']}|{policy.name}"] = reason
    status["reachability"] = reachability
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(status, indent=2) + "\n")
    have = ", ".join(status["elements_with_products"]) or "(none)"
    print(f"{len(products)} product rows across {have}; "
          f"{len(status['elements'])} species tracked -> {args.out}")
    unattributed = sum(1 for p in products if p["holding"] is None)
    status["unattributed_products"] = unattributed
    args.out.write_text(json.dumps(status, indent=2) + "\n")
    for i in status["instruments"]:
        have = {p["band"] for p in products if p["holding"] == i["holding"]}
        print(f"  {i['holding']:<34} telluric={i['telluric_applied']:<12} "
              f"bands: {sorted(have) or '—'}")
    if unattributed:
        print(f"\n  {unattributed} product rows predate RYA-933/934 and name no "
              f"holding in their filename.\n  They are attributable to an INSTRUMENT "
              f"only. Re-run to attribute them.")


if __name__ == "__main__":
    main()
