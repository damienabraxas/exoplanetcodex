#!/usr/bin/env python3
"""RYA-1187 — the (holding x band x engine x grade) applicability matrix. READ-ONLY.

    python3 scripts/rya1187_applicability_matrix.py [--check]

Ryan's completeness principle (2026-09-03): every cell must hold a product WHERE
APPLICABLE, and every empty cell must carry a REASON. A silent empty is the defect.

🔴 THE HEADLINE QUESTION IS ANSWERED "NO LINES", NOT "OWED PRODUCT". The ticket flags
zero Deep Grade in red-optical/NIR as the gap. It is not a gap: there is no
primary-lab-gf Fe line deep enough to be Deep-gradeable anywhere in either band.

    Fe I red-optical   70 lab lines   66 shallow   0 deep   deepest feature 0.583
    Fe I NIR           29 lab lines   29 shallow   0 deep   deepest feature 0.409
    Fe II red-optical    0 lab lines                        no lab gf at all
    Fe II NIR            0 lab lines                        no lab gf at all

The gate is 0.60, so red-optical misses it by 0.017 and NIR by 0.191. Emitting a Deep
product there would require either a gf we do not have or a grade reassignment, and
RYA-161 forbids both. The cell is an HONEST EMPTY and is published as one.

WHERE EVERY AXIS COMES FROM — none of it from the product feed (the ticket's instruction):

  depth grade      `derive_band_products._feature_depth` — the MEASURED `central_depth` of
                   the blended FEATURE a line belongs to, grouped by `GROUP_A` exactly as
                   `line_accounting_rya709.features()` does it. This is the quantity the
                   harness actually gates on; the accounting table's `predicted_depth` is
                   NOT it (RYA-1052 got 0 deep NIR lines from that column and it was the
                   wrong column, not the wrong answer).
  gf provenance    `pipeline.gf_grades.canonical_species` -> `gf_tier` containing LAB.
  instrument reach `data/catalog/instrument_catalog.csv` wavelength_min/max_nm, intersected
                   with the band. ⚠️ This is why harps has NO red-optical or NIR cell at
                   all: harps stops at 691.0 nm, which is the VIS red edge.
  engine reach     `pipeline.treatment_axes` + the RYA-1015 availability matrix.
  live products    `data/products/solar/Fe.json` — used ONLY to mark a cell LIVE, never to
                   decide whether it is applicable.

⚠️ DEPTH IS A LINE PROPERTY, NOT A PER-SPECTRUM MEASUREMENT. `central_depth` comes from
`data/linelists/linelist_solar.csv`, one list for all holdings, so a cell's depth census
varies between holdings ONLY through wavelength coverage. The ticket asks for per-holding
counts and they are reported per holding — but two holdings covering the same span report
the same census BY CONSTRUCTION, and that is a property of the quantity, not a bug. A
genuinely per-spectrum depth would be a different measurement and a different ticket.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from derive_band_products import _feature_depth              # noqa: E402
from line_accounting_rya709 import DEPTH_HI, DEPTH_LO        # noqa: E402
from pipeline.gf_grades import canonical_species             # noqa: E402

FEED = ROOT / "data/products/solar/Fe.json"
INSTRUMENTS = ROOT / "data/catalog/instrument_catalog.csv"
HOLDINGS = ROOT / "data/catalog/holdings_manifest_registry.csv"
ASPLUND = ROOT / "data/reference/asplund2021_fe/asplund2021_fe_lines.csv"
OUT = ROOT / "data" / "audit" / "rya1187_applicability_matrix"

#: The four DEFINED bands, from `pipeline.band_policy`. H/K are proposed, not defined, and
#: are reported by RYA-1052 rather than here.
BANDS = [("near-UV", 3000.0, 3780.0), ("VIS", 3780.0, 6910.0),
         ("red-optical", 6910.0, 9199.0), ("NIR", 9199.0, 13000.0)]

#: The solar Fe holdings actually in play — the ones the live feed measures on.
HOLDINGS_IN_PLAY = [
    ("solar_kpno_kurucz2005_corrected", "kpno_solar_atlas"),
    ("solar_kpno_molecfit_corrected", "kpno_solar_atlas"),
    ("solar_harps_molecfit_corrected", "harps"),
    ("solar_iag", "iag_fts_solar_atlas"),
    ("solar_crires_plus_y_wide_rya1054", "crires_plus"),
]

#: 🔴 THE SPAN IS THE HOLDING'S, NOT THE INSTRUMENT'S — AND THIS MATRIX GOT IT WRONG ONCE.
#: RYA-1109 names the trap in as many words: "A holding can cover less than its instrument's
#: full range (a single setting, a trimmed product), so these counts are an UPPER BOUND."
#: Measured here: the instrument catalog gives `kpno_solar_atlas` 296-1300 nm, but the
#: holdings registry's own note for `solar_kpno_kurucz2005_corrected` says the Kurucz-2005
#: irradiance product is **300-1000 nm**. On the instrument span this matrix reported a NIR
#: Codex GAP for that holding — 29 lab lines "owed a product" in a band the holding only
#: reaches to 10000 A, where the live feed has never put a product and molecfit (which does
#: reach 12976 A) has one. Declared holding spans win; anything without one falls back to
#: the instrument and SAYS SO in `span_source`.
HOLDING_SPAN_A = {
    "solar_kpno_kurucz2005_corrected": (
        (3000.0, 10000.0),
        "holdings_manifest_registry note: 'Kurucz R.L. 2005 high-resolution irradiance "
        "product ... 300-1000 nm'"),
}

#: ⚠️ CRIRES+ coverage is MEASURED from real pixels (RYA-805), not read off WAVELMIN/MAX:
#: the catalogue span 950-5300 nm is the instrument's, while the actual comb of settings
#: starts at 9479.3 A. Using the catalogue value would claim red-optical reach it does not
#: have (its floor is 280.3 A above that band's red edge — RYA-1052, Ryan 2026-08-26).
MEASURED_REACH_A = {"crires_plus": (9479.3, 24855.0)}

GRADES = ("Codex Grade", "Deep Grade", "Reference Grade")


def coverage_span_A(inst: pd.DataFrame, holding: str,
                    instrument_id: str) -> tuple[float, float, str]:
    """(lo, hi, source) — the HOLDING's span where one is declared, else the instrument's."""
    if holding in HOLDING_SPAN_A:
        (lo, hi), why = HOLDING_SPAN_A[holding]
        return lo, hi, f"HOLDING — {why}"
    if instrument_id in MEASURED_REACH_A:
        lo, hi = MEASURED_REACH_A[instrument_id]
        return lo, hi, "INSTRUMENT (measured from real pixels, RYA-805)"
    r = inst[inst.instrument_id.astype(str) == instrument_id]
    if r.empty:
        raise SystemExit(f"{instrument_id} is not in the instrument catalog")
    r = r.iloc[0]
    return (float(r.wavelength_min_nm) * 10.0, float(r.wavelength_max_nm) * 10.0,
            "INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own")


def depth_census(species: str, lo: float, hi: float) -> dict:
    """Lab-gf lines in [lo, hi] split on the depth gate. The product path's own depth."""
    cg = canonical_species(species)
    lab = cg[cg.gf_tier.astype(str).str.contains("LAB", na=False)]
    inband = lab[lab.wavelength_air_A.between(lo, hi)]
    n = len(inband)
    if not n:
        return {"lab_lines": 0, "shallow": 0, "deep": 0, "below_floor": 0,
                "depth_max": None, "depth_median": None}
    d = _feature_depth(inband.wavelength_air_A.values.astype(float))
    return {"lab_lines": n,
            "shallow": int(((d >= DEPTH_LO) & (d <= DEPTH_HI)).sum()),
            "deep": int((d > DEPTH_HI).sum()),
            "below_floor": int((d < DEPTH_LO).sum()),
            "depth_max": round(float(d.max()), 3),
            "depth_median": round(float(np.median(d)), 3)}


def asplund_in(lo: float, hi: float) -> int:
    if not ASPLUND.exists():
        return 0
    a = pd.read_csv(ASPLUND)
    return int(a.wavelength_air_A.between(lo, hi).sum())


#: 🔴 ENGINE REACH IS GATED BY THE BAND'S LINELIST, NOT BY THE DECK. `gerber_nlte.DECKS`
#: holds Fe and Fe@mean3D and neither declares a wavelength limit — RYA-713's finding in
#: its own title, "Engine A is level-indexed not wavelength-limited". What actually stops
#: an NLTE engine in a band is whether that band's linelist carries NLTE LABELS for Fe:
#: without them iSpec drops the lines into `nlte_ignored` and synthesises LTE **without
#: raising** (RYA-764), which is a silent LTE product wearing an NLTE name.
#:
#: MEASURED HERE, from the committed artifact:
#:   NIR (`ispec_ir_9200_13000/atomic_lines.tsv`)  237 Fe I rows, `nlte` == 'F' on ALL 237.
#:
#: ⚠️ NOT DETERMINABLE FROM THIS CHECKOUT, and said so rather than assumed:
#:   near-UV     `ispec_nearuv_3000_3780/atomic_lines.tsv` is DELIBERATELY not committed
#:               (its README: regenerable, 12 MB, built on Sirius).
#:   VIS + red-optical  bind `ispec_ges_v6`, an iSpec-BUNDLED list resolved at runtime and
#:               not in the repo at all.
#: A cell that cannot be decided here is published as UNDETERMINED-IN-CHECKOUT with the
#: reason. Guessing it would be the RYA-833 shape: a well-formed answer to a question that
#: was never asked.
NLTE_LABEL_EVIDENCE = {
    "NIR": ("NO_NLTE_LABELS",
            "data/linelists/ispec_ir_9200_13000/atomic_lines.tsv carries 237 Fe I rows and "
            "`nlte` reads 'F' on every one. An NLTE engine here would be silently "
            "synthesised in LTE (RYA-764), so the NLTE treatments are NOT applicable."),
    "near-UV": ("UNDETERMINED_IN_CHECKOUT",
                "ispec_nearuv_3000_3780/atomic_lines.tsv is deliberately not committed "
                "(regenerable on Sirius, 12 MB). NLTE labelling cannot be read here."),
    "VIS": ("LABELS_PRESENT",
            "the live feed carries ENGINE-B-NLTE and synth-mean3D-NLTE products measured "
            "in VIS, which is direct evidence the bound list is NLTE-labelled for Fe."),
    "red-optical": ("UNDETERMINED_IN_CHECKOUT",
                    "binds `ispec_ges_v6`, an iSpec-bundled list resolved at runtime and "
                    "not present in the repo. Same list as VIS, where labels ARE present, "
                    "but the red-optical SLICE cannot be confirmed from this checkout."),
}


def build() -> dict:
    inst = pd.read_csv(INSTRUMENTS, comment="#")
    feed = json.loads(FEED.read_text())
    live = {(p["holding"], p["band"], p["grade"], p["treatment"]) for p in feed["products"]}
    live_cell = {(p["holding"], p["band"], p["grade"]) for p in feed["products"]}

    rows = []
    for holding, instrument in HOLDINGS_IN_PLAY:
        ilo, ihi, span_src = coverage_span_A(inst, holding, instrument)
        for band, blo, bhi in BANDS:
            lo, hi = max(ilo, blo), min(ihi, bhi)
            reaches = hi > lo
            cen = {sp: depth_census(sp, lo, hi) for sp in ("Fe I", "Fe II")} if reaches else {}
            n_asp = asplund_in(lo, hi) if reaches else 0
            for grade in GRADES:
                if not reaches:
                    verdict, why = "ABSENT_NO_COVERAGE", (
                        f"this holding covers {ilo:.0f}-{ihi:.0f} A ({span_src}) and the "
                        f"{band} band is {blo:.0f}-{bhi:.0f} A — they do not overlap. No "
                        f"spectrum exists to measure.")
                    n = 0
                elif grade == "Reference Grade":
                    n = n_asp
                    verdict, why = (("LIVE" if (holding, band, grade) in live_cell else "GAP"),
                                    f"{n} AGSS21 Table A.2 line(s) in {band} within this "
                                    f"holding's coverage") if n else (
                        "ABSENT_NO_REFERENCE_LINE",
                        f"AGSS21 Table A.2 publishes no Fe line inside {band} ∩ this "
                        f"holding's span ({lo:.0f}-{hi:.0f} A). The reference set is "
                        f"optical; this is a published SELECTION, not a coverage hole.")
                else:
                    key = "shallow" if grade == "Codex Grade" else "deep"
                    n = sum(cen[sp][key] for sp in cen)
                    if n:
                        verdict = "LIVE" if (holding, band, grade) in live_cell else "GAP"
                        why = (f"{n} primary-lab-gf Fe line(s) at this depth grade "
                               f"({'<=' if key == 'shallow' else '>'}{DEPTH_HI}) in "
                               f"{lo:.0f}-{hi:.0f} A")
                    else:
                        verdict = "ABSENT_NO_LAB_LINE_AT_THIS_DEPTH"
                        dmax = max([c["depth_max"] for c in cen.values()
                                    if c["depth_max"] is not None], default=None)
                        why = (f"no primary-lab-gf Fe line with feature depth "
                               f"{'<=' if key == 'shallow' else '>'}{DEPTH_HI} in "
                               f"{lo:.0f}-{hi:.0f} A "
                               + (f"(deepest lab feature in band = {dmax})" if dmax is not None
                                  else "(no primary-lab-gf Fe line in this band at all)"))
                rows.append({
                    "holding": holding, "instrument": instrument, "band": band,
                    "grade": grade, "band_lo_A": blo, "band_hi_A": bhi,
                    "covered_lo_A": round(lo, 1) if reaches else None,
                    "covered_hi_A": round(hi, 1) if reaches else None,
                    "instrument_reaches_band": bool(reaches),
                    "span_source": span_src,
                    "n_lines_at_grade": int(n),
                    "fe1_lab": cen.get("Fe I", {}).get("lab_lines", 0) if reaches else 0,
                    "fe1_shallow": cen.get("Fe I", {}).get("shallow", 0) if reaches else 0,
                    "fe1_deep": cen.get("Fe I", {}).get("deep", 0) if reaches else 0,
                    "fe1_depth_max": cen.get("Fe I", {}).get("depth_max") if reaches else None,
                    "fe2_lab": cen.get("Fe II", {}).get("lab_lines", 0) if reaches else 0,
                    "fe2_deep": cen.get("Fe II", {}).get("deep", 0) if reaches else 0,
                    "verdict": verdict, "reason": why,
                })
    return {"cells": pd.DataFrame(rows), "live": live, "feed_version": feed["version"]}


def engine_rows(cells: pd.DataFrame, live: set) -> pd.DataFrame:
    """The ENGINE axis: for each (band, treatment), is it applicable and is it live?"""
    feed = json.loads(FEED.read_text())
    treatments = sorted({p["treatment"] for p in feed["products"]})
    live_bt = {(p["band"], p["treatment"]) for p in feed["products"]}
    nlte = tuple(t for t in treatments if "NLTE" in t.upper())
    out = []
    for band, _, _ in BANDS:
        any_line = bool(cells[(cells.band == band) & cells.instrument_reaches_band
                              & (cells.n_lines_at_grade > 0)].shape[0])
        state, ev = NLTE_LABEL_EVIDENCE[band]
        for t in treatments:
            is_nlte = t in nlte
            if not any_line:
                verdict, why = "ABSENT_NO_LINES_IN_BAND", f"no gradeable Fe line in {band}"
            elif (band, t) in live_bt:
                verdict, why = "LIVE", "product exists in the feed"
            elif t == "ENGINE-A-3DNLTE":
                verdict, why = ("ABSENT_ENGINE_BOXED", (
                    "the Amarsi MLP is the 3D-NLTE route and every product it has ever "
                    "emitted is VIS; RYA-1106 ran it on AGSS21's optical set. Whether its "
                    "grid reaches this band is NOT decidable from this checkout — the "
                    "network's own domain check is on STELLAR parameters, not wavelength "
                    "(pipeline/amarsi3d.py). Reported as owed-verification, not as reach."))
            elif is_nlte and state == "NO_NLTE_LABELS":
                verdict, why = "ABSENT_NO_NLTE_LABELS_IN_BAND_LINELIST", ev
            elif is_nlte and state == "UNDETERMINED_IN_CHECKOUT":
                verdict, why = "UNDETERMINED_IN_CHECKOUT", ev
            else:
                verdict, why = "GAP", (
                    f"{band} has gradeable Fe lines and this treatment is not NLTE-gated "
                    f"here, so a product is applicable and missing")
            out.append({"band": band, "treatment": t, "is_nlte_treatment": is_nlte,
                        "verdict": verdict, "reason": why})
    return pd.DataFrame(out)


def render(cells: pd.DataFrame, engines: pd.DataFrame, feed_version: str) -> str:
    n_silent = 0
    lines = [f"# RYA-1187 — (holding x band x engine x grade) applicability matrix",
             "",
             f"Built from LINE-LEVEL data against `Fe.json` v{feed_version}. Every cell is a "
             f"live product or a stated reason; a silent empty is the defect this exists to "
             f"remove.", "",
             "## The headline question: Deep Grade in red-optical and NIR", "",
             "**Answered NO LINES, not owed product.** Per-band primary-lab-gf Fe census "
             "(the product path's own feature depth, gate 0.60):", "",
             "| band | species | lab lines | shallow (Codex) | deep (Deep) | deepest feature |",
             "|---|---|---|---|---|---|"]
    for band, lo, hi in BANDS:
        for sp in ("Fe I", "Fe II"):
            c = depth_census(sp, lo, hi)
            lines.append(f"| {band} | {sp} | {c['lab_lines']} | {c['shallow']} | "
                         f"**{c['deep']}** | {c['depth_max']} |")
    lines += ["",
              "Red-optical misses the 0.60 gate by **0.017** (deepest lab feature 0.583) and "
              "NIR by **0.191** (0.409). Fe II has no primary-lab gf in either band at all. "
              "Emitting a Deep product there would need a gf we do not hold or a grade "
              "reassignment — RYA-161 forbids both.", "",
              "## Cell matrix — holding x band x grade", "",
              "| holding | band | grade | verdict | n lines | reason |", "|---|---|---|---|---|---|"]
    for _, r in cells.iterrows():
        if r.verdict in ("LIVE", "GAP"):
            n_silent += 0
        lines.append(f"| {r.holding} | {r.band} | {r.grade} | **{r.verdict}** | "
                     f"{r.n_lines_at_grade} | {r.reason} |")
    lines += ["", "## Engine axis — band x treatment", "",
              "| band | treatment | verdict | reason |", "|---|---|---|---|"]
    for _, r in engines.iterrows():
        lines.append(f"| {r.band} | {r.treatment} | **{r.verdict}** | {r.reason} |")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    b = build()
    cells, live = b["cells"], b["live"]
    engines = engine_rows(cells, live)
    md = render(cells, engines, b["feed_version"])

    #: 🔴 THE GAPS ARE IDENTIFIED HERE AND CANNOT BE EMITTED HERE. Emitting a product means
    #: RUNNING a synthesis, and this checkout cannot: `import ispec` fails on the Mac, the
    #: spectra live at /mnt/codex-data on Sirius, and the near-UV/GES line lists are not
    #: committed. Writing rows without a run would be the RYA-1034 laundering this repo
    #: exists to prevent — a product enters through `publish_product.py` from a real
    #: artifact or it does not exist. So the owed products are NAMED with what each needs.
    owed = cells[cells.verdict.eq("GAP")]
    silent = cells[cells.verdict.eq("")]
    verdict = {
        "ticket": "RYA-1187", "feed_version": b["feed_version"],
        "n_cells": int(len(cells)), "n_engine_cells": int(len(engines)),
        "cells_by_verdict": {k: int(v) for k, v in cells.verdict.value_counts().items()},
        "engine_by_verdict": {k: int(v) for k, v in engines.verdict.value_counts().items()},
        "silent_empties": int(len(silent)),
        "owed_products": [
            {"holding": r.holding, "band": r.band, "grade": r.grade,
             "n_lines_at_grade": int(r.n_lines_at_grade), "covered_A": [r.covered_lo_A, r.covered_hi_A]}
            for _, r in owed.iterrows()],
        "🔴 emit_blocked_here": {
            "blocked": True,
            "why": ("this checkout cannot run a synthesis: `import ispec` fails on the Mac, "
                    "the spectra are staged at /mnt/codex-data on Sirius, and the near-UV "
                    "linelist + iSpec's GES v6 are not in the repo. A product enters the "
                    "feed through publish_product.py from a REAL artifact or it does not "
                    "exist (RYA-1034); rows written without a run would be laundering."),
            "runs_on": "Sirius",
            "what_each_needs": (
                "Reference Grade red-optical/NIR: measure the AGSS21 line set on that "
                "holding via `scripts/measure_reference_lineset.py --line-set asplund`, "
                "then publish with `--line-set asplund`. NIR Codex on "
                "solar_kpno_kurucz2005_corrected: a derive_band_products run over "
                "9199-10000 A on the GRADED selector."),
        },
        "headline": {
            "question": "is there a Deep Grade gap in red-optical / NIR?",
            "answer": "NO — there is no primary-lab-gf Fe line deep enough in either band",
            "red_optical": depth_census("Fe I", 6910.0, 9199.0),
            "NIR": depth_census("Fe I", 9199.0, 13000.0),
            "gate": DEPTH_HI,
        },
    }
    if a.check:
        bad = []
        for name, txt in (("applicability_matrix.md", md),
                          ("verdict.json", json.dumps(verdict, indent=2) + "\n")):
            p = OUT / name
            if not p.exists() or p.read_text() != txt:
                bad.append(name)
        for name, df in (("cells.csv", cells), ("engine_axis.csv", engines)):
            import io
            buf = io.StringIO(); df.to_csv(buf, index=False)
            p = OUT / name
            if not p.exists() or p.read_text() != buf.getvalue():
                bad.append(name)
        print("DRIFT: " + ", ".join(bad) if bad else "matrix reproduces")
        return 1 if bad else 0

    OUT.mkdir(parents=True, exist_ok=True)
    cells.to_csv(OUT / "cells.csv", index=False)
    engines.to_csv(OUT / "engine_axis.csv", index=False)
    (OUT / "applicability_matrix.md").write_text(md)
    (OUT / "verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(f"cells: {verdict['cells_by_verdict']}")
    print(f"engine: {verdict['engine_by_verdict']}")
    print(f"silent empties: {verdict['silent_empties']}  (must be 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
