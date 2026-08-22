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


#: Sub-paths that mark a DIAGNOSTIC variant rather than a headline product --
#: RYA-877's before/after control pair, RYA-847's gated sweep. They are real
#: products and must not be deleted from the page; they are also not what the
#: matrix is asking about, so they are labelled and hidden behind a toggle.
VARIANT_MARKERS = ("control", "gated")


#: The dates the telluric-corrected holdings first existed in the repository.
#: A product committed before its instrument had a corrected holding CANNOT have
#: used one -- this is provenance, not inference, and it is the first question
#: anyone asks of a number on this page.
CORRECTED_HOLDING_BORN = {
    "harps": "2026-08-20",              # RYA-931, commit 4d8abf8
    "kpno_solar_atlas": "2026-08-21",   # RYA-940, commit c0465b1
}


def telluric_state_of(row: dict, committed: str | None) -> dict:
    """Was this product made before or after its instrument had a corrected holding?

    Where the row names a holding, the holding answers it outright. Where it does
    not -- every product predating RYA-933/934 -- the date still answers it: the
    corrected holdings did not exist, so nothing could have used them.
    """
    if row.get("holding"):
        return {"telluric_basis": "named holding", "telluric_epoch": None}
    born = CORRECTED_HOLDING_BORN.get(row["instrument"])
    if born and committed and committed < born:
        return {"telluric_basis": "PRE-correction (provenance)",
                "telluric_epoch": f"committed {committed}; {row['instrument']} had no "
                                  f"corrected holding until {born}"}
    return {"telluric_basis": "unknown", "telluric_epoch": None}


def run_context(path: Path) -> dict:
    """WHICH RUN produced this row. Part of a product's identity, not decoration.

    Six identities appear more than once across ticket output directories, and
    they are NOT duplicates: rya877/control/before gives Fe II 7.568 where
    rya877/control/after gives 7.542. That pair is the whole point of a control.
    Deduplicating on (species, holding, band, engine) would silently keep one and
    discard the other -- picking a winner between two numbers that were produced
    to be compared.
    """
    import subprocess
    rel = path.relative_to(ROOT / "data" / "results").parent
    parts = rel.parts
    try:
        committed = subprocess.run(
            ["git", "log", "--format=%ad", "--date=short", "-1", "--", str(path)],
            cwd=ROOT, capture_output=True, text=True, timeout=20).stdout.strip() or None
    except Exception:                                           # noqa: BLE001
        committed = None
    return {"run_context": str(rel) if parts else "(root)",
            "ticket_dir": parts[0] if parts else None,
            "committed": committed,
            "is_variant": any(m in parts for m in VARIANT_MARKERS)}


def display_name(row) -> str:
    """The physics-axis name, DERIVED from the stored axes (RYA-906).

    Not `treatment`. "ENGINE-A" / "ENGINE-B" are letters, not physics: they say
    nothing about route, scale or model, and the same letter means different
    things on different rows. RYA-906 stored the five axes precisely so the name
    could be derived, and the tracker was still showing the legacy labels.

    Route comes from the HANDLER where the row carries one -- never from the
    label. RYA-906 measured this over 2153 committed rows: `1D-LTE` and
    `1D-LTE-LABGF` each pair with BOTH ProfileFitHandler and SynthesisHandler, so
    on those labels the legacy string is not merely lossy, it is FALSE.
    """
    from pipeline.treatment_axes import Axes
    stored = {k: row.get(k) for k in ("route", "scale", "model", "atmos", "gf")}
    if all(v is not None and not pd.isna(v) for v in stored.values()):
        return Axes(route=str(stored["route"]), scale=str(stored["scale"]),
                    model=str(stored["model"]), atmos=str(stored["atmos"]),
                    gf=str(stored["gf"]),
                    route_basis=str(row.get("route_basis") or "stored")).display
    # No stored axes = a row that predates RYA-906. Derive from the legacy label
    # plus whatever route evidence the row carries, and fail loudly on an unknown
    # label rather than defaulting -- a silent default is how RYA-869 published
    # four wrong systematics.
    from pipeline.treatment_axes import display_for, UnknownTreatment
    try:
        return display_for(str(row.get("treatment")), handler=row.get("handler"))
    except (UnknownTreatment, ValueError):
        return f"{row.get('treatment')} (unresolved axes)"


#: `gf rung: gf rung N (term): reason` -- the reason text states the graded count in one
#: of two shapes, and both are parsed here rather than recomputed. Recomputing would mean
#: re-grading every line against the line list inside the tracker, which is a SECOND
#: implementation of membership (the RYA-845 two-homes shape); the budget file is the
#: artifact the product was actually CHARGED on, so it is the honest source.
_RUNG_MIXED = re.compile(r"MIXED POOL:\s*(\d+)\s+of\s+(\d+)\s")
_RUNG_ALL = re.compile(r"every one of the\s+(\d+)\s")
_RUNG_HEAD = re.compile(r"gf rung:\s*gf rung\s*(\d+)\s*\(([^)]*)\)")


def graded_counts(products_csv: Path) -> dict:
    """(n_graded, n_pool, rung) per treatment, read from the sibling *_budgets.txt.

    Returns {} when there is no budget file -- an older artifact predates the gf rung and
    must read as UNKNOWN, never as zero. Zero graded lines is a real, different statement
    from "this product was written before we recorded the rung".
    """
    b = products_csv.with_name(products_csv.name.replace("_products.csv", "_budgets.txt"))
    if not b.exists():
        return {}
    out, treatment = {}, None
    for line in b.read_text(errors="replace").splitlines():
        t = line.strip()
        # Budget blocks open with the cell header, e.g. "Fe . VIS . n=6"; the gf rung line
        # belongs to the block above it. Track the most recent non-indented header.
        if t and not line.startswith(" ") and "gf rung" not in t:
            treatment = t
        if "gf rung:" not in t:
            continue
        head = _RUNG_HEAD.search(t)
        rung = int(head.group(1)) if head else None
        m = _RUNG_MIXED.search(t)
        if m:
            n_graded, n_pool = int(m.group(1)), int(m.group(2))
        else:
            m2 = _RUNG_ALL.search(t)
            if not m2:
                continue
            n_graded = n_pool = int(m2.group(1))
        # KEYED BY POOL SIZE, not by treatment name. The budget block header is
        # "Fe · VIS · n=148" -- it carries the cell and the pool size but NOT the
        # treatment string the products table uses, and the two vocabularies do not
        # match (RYA-906: `display` is derived, `treatment` is the legacy label). The
        # header's n IS the pool the budget was charged on, so it joins to the row's
        # n_lines exactly.
        out[int(n_pool)] = {"n_graded": n_graded, "n_pool": n_pool, "gf_rung": rung,
                            "budget_cell": treatment}
    return out


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
            ctx = run_context(path)
            graded = graded_counts(path)
            for _, r in frame.iterrows():
                rows.append({
                    **meta,
                    "band": r.get("band"), "treatment": r.get("treatment"),
                    "display": display_name(r),
                    "A": None if pd.isna(r.get("A")) else float(r["A"]),
                    "sigma_stat": None if pd.isna(r.get("stat_dex")) else float(r["stat_dex"]),
                    "sigma_syst": None if pd.isna(r.get("syst_dex")) else float(r["syst_dex"]),
                    "n_lines": None if pd.isna(r.get("n_lines")) else int(r["n_lines"]),
                    # RYA-946 — how many of THIS engine's lines carry a primary-lab gf.
                    # None means the artifact predates the gf rung, which is not 0.
                    **((graded.get(int(r["n_lines"]))
                        if not pd.isna(r.get("n_lines")) else None)
                       or {"n_graded": None, "n_pool": None, "gf_rung": None}),
                    "source": str(path.relative_to(ROOT)),
                    **ctx,
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

    # A SECOND comparator, kept separate rather than averaged: Lodders, Bergemann &
    # Palme 2025 Table 6. Its PRESENT column is the photospheric-era value; the
    # proto-solar column runs ~0.09 dex higher and quoting it against a photospheric
    # measurement would manufacture a discrepancy.
    #
    # Comparators are a LIST, and each may be scoped to particular bands. A source
    # whose determination only covers the infrared must not be drawn across the
    # optical as if it applied there.
    lodders = root / "data" / "reference" / "solar" / "lodders2025_table6.csv"
    if lodders.exists():
        frame = pd.read_csv(lodders, comment="#")
        by_element = {r["element"]: r for _, r in frame.iterrows()}
        for key, entry in out.items():
            element = "".join(c for c in key if not c.isupper() or c == key[0])
            element = key.rstrip("IV") or key
            row = by_element.get(element)
            if row is None:
                continue
            entry.setdefault("comparators", []).append({
                "name": "Lodders+ 2025", "value": float(row["A_present"]),
                "sigma": float(row["sigma_present"]), "colour": "gold",
                "bands": None,          # applies everywhere
                "note": "Table 6, present-day Sun (proto-solar is "
                        f"{float(row['A_protosolar']):.2f}, ~0.09 dex higher)",
                "source": "data/reference/solar/lodders2025_table6.csv",
            })
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


def collect_graded(root: Path) -> list[dict]:
    """The GRADED lab-gf cells -- RYA-850's primary reported value.

    These were missing from the tracker entirely, which is why the page appeared to
    show error bars that had grown. They had not: the page was showing only the
    MIXED pools, which carry the blanket 0.17 dex ungraded-gf placeholder. A fully
    graded cell carries its own CITED pool sigma instead, and the systematic drops
    to 0.061-0.113.

    They live in `rya850_summary.json`, not in a `*_products.csv`, so the product
    glob never saw them.

    RYA-851's reporting contract: GRADED is primary, UNGRADED is secondary, and the
    ungraded value is never the headline. Both are shown -- more lines buy a wider
    gf floor, which is a trade, not a defect.
    """
    summary = root / "data" / "results" / "rya850" / "rya850_summary.json"
    if not summary.exists():
        return []
    doc = json.loads(summary.read_text())
    out = []
    for cell in doc.get("graded_cells", []):
        out.append({
            "element": "Fe", "ion": cell.get("ion", "I"), "band": cell["band"],
            "engine": cell.get("engine"), "A": cell["A"], "n_lines": cell["n_lines"],
            "sigma_stat": cell["stat_dex"], "sigma_syst": cell["syst_dex"],
            "total_dex": cell.get("total_dex"),
            "ungraded_total_dex": cell.get("ungraded_total_dex"),
            "graded_beats_ungraded": cell.get("graded_beats_ungraded"),
            "gf_term": doc.get("gf_term_published"),
            "role": "PRIMARY (graded lab-gf)",
            "source": "data/results/rya850/rya850_summary.json",
        })
    return out


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
        "graded": collect_graded(ROOT),
        "reporting_contract": {
            "primary": "graded lab-gf pool, on its own CITED pool sigma (RYA-850)",
            "secondary": "ungraded all-lines pool, on the 0.17 dex gf placeholder",
            "headline_rule": "the ungraded value is NEVER the headline (RYA-851)",
            "bars": "statistical SOLID, systematic WIREFRAME -- never summed; "
                    "error_budget.py deliberately provides no combined()",
        },
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
    status["variant_products"] = sum(1 for p in products if p["is_variant"])
    status["run_contexts"] = sorted({p["run_context"] for p in products})
    for row in products:
        row.update(telluric_state_of(row, row.get("committed")))
    status["pre_correction_products"] = sum(
        1 for p in products if p["telluric_basis"].startswith("PRE"))
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
