#!/usr/bin/env python3
"""RYA-1213 — the Reference Grade (band x holding x ion x engine) matrix. READ-ONLY.

    python3 scripts/rya1213_reference_matrix.py [--check]

Ryan's completeness principle, applied to the tier RYA-946 defines by gf PEDIGREE: every
cell holds a Reference product where applicable, and every empty cell carries a REASON.

WHAT THIS ANSWERS, AND WHERE EACH ANSWER COMES FROM — none of it from the ticket's prose:

  the Reference pool   canonical_gf rows whose `gf_tier` contains LAB, for the species,
                       inside the product's OWN window. No depth term: that is the whole
                       definition of the tier (`derive_band_products._cand_reference`).
  the depth split      `derive_band_products._feature_depth` + `line_accounting_rya709`'s
                       gate. REPORTED, never applied — it is what a reader needs to see
                       how this pool relates to its Codex and Deep siblings.
  engine N/A           `data/results/rya1208/rya1208_applicability.json` and
                       `rya1208_nlte_refusals.json` — MEASURED refusals from real runs,
                       not asserted here. The Amarsi domain verdict comes from RYA-1187's
                       `amarsi_domain_check.csv`, measured on Sirius.
  live products        `data/products/solar/Fe.json` — used ONLY to mark a cell LIVE.

🔴 THE HEADLINE IS THAT "429 LAB LINES HAVE NO REFERENCE PRODUCT" IS TRUE OF THE LABEL
AND NOT OF THE MEASUREMENT. Outside VIS the lab pool does not straddle the depth gate at
all, so the Reference pool is not a new population — it is an existing one under a
different name:

    near-UV      58 of 59 Fe I lab lines are DEEP, 0 shallow   -> Reference == Deep + 1
    red-optical   0 of 70 are deep                             -> Reference == Codex + 1
    NIR           0 of 29 are deep                             -> Reference == Codex
    H             0 of 25 are deep                             -> Reference == Codex
    VIS          67 shallow / 109 deep                         -> genuinely a third pool

So in NIR and H the Reference product measures the SAME LINES as its Codex sibling, and
the new statement is the GRADE rather than the number.

🔴 AND IT REPRODUCES THE SAME NUMBER ONLY WHEN THE COMPARISON IS PAIRED. Against the
STORED Codex values the CRIRES+ Y Reference product looked different -- 7.546 vs 7.551
(1D-LTE), 7.486 vs 7.492 (ENGINE-A) -- on pools that are five identical wavelengths. The
GRADED leg re-run on the same commit gives 7.546 and 7.486: the selector difference is
ZERO and all of the apparent difference was CODE DRIFT between a 2026-08-26 artifact and
today (RYA-1204). Reading a Reference product against a stored Codex value measures drift
as well as the selector. Measured: nir_paired_selector_control.json.

⚠️ TWO POPULATIONS OF LAB LINES ARE UNREACHABLE AND NEITHER IS A MISSING RUN.
  * 68 lab lines (67 Fe I, 1 Fe II) sit in 3780-4200 A, between the near-UV list's red
    edge and GES v6's blue edge. NO synthesis line list covers that span, so the
    harness refuses there rather than emit a product naming a range it did not measure
    (`derive_band_products` line-coverage guard, RYA-911/913/967).
  * 2 Ruffoni-2013 lab lines (14679.831, 14826.408 A) sit between NIR's red edge
    (12976 A) and the CRIRES+ H arm's blue edge (15007 A). No holding reaches them.
Both are counted here so the ticket's "450 lab lines" reconciles to the last line.
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

FEED = ROOT / "data/products/solar/Fe.json"
CANON = ROOT / "data/linelists/canonical_gf.csv"
RYA1208 = ROOT / "data/results/rya1208/rya1208_applicability.json"
RYA1208_REF = ROOT / "data/results/rya1208/rya1208_nlte_refusals.json"
AMARSI_DOMAIN = ROOT / "data/audit/rya1187_applicability_matrix/amarsi_domain_check.csv"
OUT = ROOT / "data" / "audit" / "rya1213_reference_matrix"

#: 🔴 THE WINDOW IS THE PRODUCT'S, NOT THE BAND'S, AND THE DIFFERENCE IS 68 LINES. Every
#: live VIS product runs 4200-6910 A, not the band's 3780-6910, because the VIS synthesis
#: list (iSpec-vendored GES v6) starts at 4200 A and `derive_band_products` refuses a run
#: it cannot cover. Censusing the Reference pool over the BAND would report 253 VIS lab
#: lines "owed a product" when 68 of them are in a span no synthesis list reaches. Each
#: entry is the window the Reference run actually used, so the count here is the count the
#: product can contain.
#: ⚠️ NIR's red edge is 12976 A, not the band's 13000: the NIR list ends at 12976.14 and
#: a 13000 A request is refused by the same guard (measured, this ticket).
#: 🔴 A HOLDING CAN COVER LESS THAN ITS BAND, AND CENSUSING THE BAND MANUFACTURES A GAP.
#: RYA-1187 names the trap and MEASURED both entries below; they are imported in spirit
#: rather than re-derived, and the windows in CELLS are intersected with them. Without
#: this, CRIRES+ Y would be censused over 9199-12976 A and report 29 NIR lab lines "owed
#: a product" when its real comb of settings starts at 9479.3 A and the holding is a
#: 9800-10796 A window — the Reference run measured 5 lines there, not 29.
HOLDING_SPAN_A = {
    # holding -> ((lo, hi), source)
    "solar_kpno_kurucz2005_corrected": (
        (3000.0, 10000.0),
        "holdings_manifest_registry note: the Kurucz-2005 irradiance product is "
        "300-1000 nm (RYA-1187)"),
    "solar_crires_plus_y_wide_rya1054": (
        (9800.0, 10796.0),
        "the RYA-1054 wide Y window, as the holdings registry declares it"),
    "solar_crires_plus_h_rya1094": (
        (15007.11, 17493.69),
        "the RYA-1094 H arm span, as the holdings registry declares it"),
    "solar_iag": (
        (3000.0, 11083.0),
        "MEASURED from the live NIR products' own red edge (9199-11083 A); the IAG FTS "
        "atlas stops there and a 12976 A census would claim reach it does not have"),
    "solar_harps_molecfit_corrected": (
        (3782.6, 6910.0),
        "MEASURED, RYA-1208: HARPS begins 2.6 A above the near-UV product window's red "
        "edge and ends at the VIS red edge"),
}


def _clip(holding: str, lo: float, hi: float) -> tuple[float, float, str]:
    """The window this holding can actually serve, and where the bound came from."""
    if holding not in HOLDING_SPAN_A:
        return lo, hi, "BAND — this holding declares no span of its own"
    (hlo, hhi), why = HOLDING_SPAN_A[holding]
    return max(lo, hlo), min(hi, hhi), f"HOLDING — {why}"


CELLS = [
    # (band, ion, holding, instrument, lo_A, hi_A)
    ("near-UV", "I",  "solar_kpno_kurucz2005_corrected", "kpno_solar_atlas", 3000.0, 3780.0),
    ("near-UV", "I",  "solar_kpno_molecfit_corrected",   "kpno_solar_atlas", 3000.0, 3780.0),
    ("near-UV", "II", "solar_kpno_kurucz2005_corrected", "kpno_solar_atlas", 3000.0, 3780.0),
    ("near-UV", "II", "solar_kpno_molecfit_corrected",   "kpno_solar_atlas", 3000.0, 3780.0),
    ("VIS", "I",  "solar_kpno_kurucz2005_corrected", "kpno_solar_atlas", 4200.0, 6910.0),
    ("VIS", "I",  "solar_kpno_molecfit_corrected",   "kpno_solar_atlas", 4200.0, 6910.0),
    ("VIS", "I",  "solar_harps_molecfit_corrected",  "harps",            4200.0, 6910.0),
    ("VIS", "I",  "solar_iag",                       "iag_fts_solar_atlas", 4200.0, 6910.0),
    ("VIS", "II", "solar_kpno_kurucz2005_corrected", "kpno_solar_atlas", 4200.0, 6910.0),
    ("VIS", "II", "solar_kpno_molecfit_corrected",   "kpno_solar_atlas", 4200.0, 6910.0),
    ("VIS", "II", "solar_harps_molecfit_corrected",  "harps",            4200.0, 6910.0),
    ("red-optical", "I", "solar_kpno_kurucz2005_corrected", "kpno_solar_atlas", 6910.0, 9199.0),
    ("red-optical", "I", "solar_kpno_molecfit_corrected",   "kpno_solar_atlas", 6910.0, 9199.0),
    ("red-optical", "I", "solar_iag",                       "iag_fts_solar_atlas", 6910.0, 9199.0),
    ("NIR", "I", "solar_kpno_kurucz2005_corrected",    "kpno_solar_atlas", 9199.0, 12976.0),
    ("NIR", "I", "solar_kpno_molecfit_corrected",      "kpno_solar_atlas", 9199.0, 12976.0),
    ("NIR", "I", "solar_iag",                          "iag_fts_solar_atlas", 9199.0, 12976.0),
    ("NIR", "I", "solar_crires_plus_y_wide_rya1054",   "crires_plus",      9199.0, 12976.0),
    ("H", "I", "solar_crires_plus_h_rya1094", "crires_plus", 15007.0, 17494.0),
]

#: The five treatments a Reference run can emit, and what decides each. Every N/A string
#: below is either quoted from a MEASURED refusal artifact or names the measurement that
#: establishes it — none is an assertion made here.
ENGINES = ("1D-LTE", "ENGINE-A", "ENGINE-B", "ENGINE-B-NLTE", "synth-1D-LTE-gerber",
           "synth-mean3D-LTE-gerber-stagger", "synth-mean3D-NLTE-gerber-stagger",
           "ENGINE-A-3DNLTE")

#: 🔴 QUOTED, NOT RESTATED. The Fe II Gerber limit and the two zero-label bands are
#: RYA-1055/1208 measurements; reproducing the wording here rather than the verdict would
#: let this file drift from the artifact that owns it, so the reasons are READ from the
#: artifacts at run time and only the KEYS live here.
NLTE_LABEL_BANDS = {"near-UV": "near-UV", "NIR": "NIR", "H": "H"}



def _canon() -> pd.DataFrame:
    return pd.read_csv(CANON, low_memory=False)


def reference_pool(cg: pd.DataFrame, ion: str, lo: float, hi: float) -> dict:
    """The Reference pool for one cell, plus the depth split it does NOT apply."""
    species = f"Fe {ion}"
    lab = cg[(cg.species == species)
             & cg.gf_tier.astype(str).str.contains("LAB", na=False)
             & cg.wavelength_air_A.between(lo, hi)]
    n = len(lab)
    if not n:
        return {"reference_pool": 0, "would_be_codex": 0, "would_be_deep": 0,
                "no_known_depth": 0, "depth_max": None, "lab_sources": {}}
    d = _feature_depth(lab.wavelength_air_A.values.astype(float))
    return {
        "reference_pool": n,
        # `would_be_*` and not `codex`/`deep`: these counts are a COMPARISON, and naming
        # them as though the Reference product applied them would be the defect this
        # selector exists to avoid.
        "would_be_codex": int((d <= DEPTH_HI).sum()),
        "would_be_deep": int((d > DEPTH_HI).sum()),
        "below_depth_floor": int((d < DEPTH_LO).sum()),
        "no_known_depth": int(np.isnan(d).sum()),
        "depth_max": (None if np.all(np.isnan(d)) else round(float(np.nanmax(d)), 3)),
        "lab_sources": {str(k): int(v)
                        for k, v in lab.lab_source_tag.value_counts().items()},
    }


def _amarsi_out_of_domain(band: str) -> tuple[bool, int]:
    if not AMARSI_DOMAIN.exists():
        return False, 0
    d = pd.read_csv(AMARSI_DOMAIN)
    s = d[d.synth_bands_band == band]
    return (bool(len(s)) and not bool(s.in_amarsi_domain.any())), int((~s.in_amarsi_domain).sum())


def engine_verdicts(band: str, ion: str, live: set, holding: str,
                    app: dict, refus: dict) -> dict:
    """Per-engine LIVE / GAP / N/A-with-reason for one cell."""
    out = {}
    labels = app.get("nlte_label_coverage_by_band", {})
    fe2_limit = (app.get("nlte_ion_capability", {}).get("Fe II", {}).get("limit") or "")
    for eng in ENGINES:
        na = None
        if eng in ("ENGINE-B-NLTE",):
            if ion == "II":
                na = ("Fe II NLTE is unavailable on the Gerber deck — MEASURED, "
                      "data/results/rya1055/atom_ion_reach.json: " + fe2_limit)
            elif band in labels and labels[band].get("labelled") == 0:
                r = refus.get(band, {}).get("reason", "")
                na = (f"zero NLTE level labels for Fe in this band's synthesis list "
                      f"({labels[band]['linelist']}: {labels[band]['fe_lines']} Fe lines, "
                      f"{labels[band]['labelled']} labelled). {r}")
        elif eng in ("synth-mean3D-LTE-gerber-stagger", "synth-mean3D-NLTE-gerber-stagger"):
            # Ryan, 2026-09-12: applicable, just not being run now. The mean-3D legs wait
            # on the Bride. Not N/A and not a gap — DEFERRED.
            out[eng] = {"verdict": "DEFERRED",
                        "reason": ("the Gerber mean-3D legs are applicable but are not "
                                   "being run now; they wait on the Bride (Ryan, "
                                   "2026-09-12).")}
            continue
        elif eng == "ENGINE-A":
            if band == "H":
                na = ("the Bergemann MPIA grid does not reach the H arm. MEASURED by the "
                      "Reference run itself, which queried it and reported '20 of 20 "
                      "line(s) unserved' rather than emitting a product: the service "
                      "returns Fe I through 12648.742 A and nothing from 15051.700 A. "
                      "⚠️ 'N of N unserved' is ALSO what a DOWN service says, so this is "
                      "recorded as a reach limit only because the same harness served "
                      "the NIR cells in this ticket's own pool minutes earlier.")
        elif eng == "ENGINE-B":
            na = ("`ENGINE-B` is a RETIRED SPELLING of `1D-LTE`, not a second engine -- "
                  "`treatment_axes.DEPRECATED_ALIASES` maps it to ('1D-LTE', 'synth'), "
                  "and `plot_grid._pick` refuses to let it win the 1D-LTE slot from the "
                  "product it aliases. A Reference cell here would duplicate this run's "
                  "own 1D-LTE product under a name that can never render. Not run.")
        elif eng == "ENGINE-A-3DNLTE":
            if band != "VIS":
                allout, n_out = _amarsi_out_of_domain(band)
                na = ("the Amarsi 2022 MLP is VIS-boxed: every non-VIS Fe I line fails "
                      "`transition energy Eup-Elo outside training [1.8190, 2.5898] eV`, "
                      "and dE IS the wavelength, so a redder line is below the box by "
                      "construction. `amarsi3d.classify_line` refuses to extrapolate. "
                      "MEASURED: rya1187_applicability_matrix/amarsi_domain_check.csv"
                      + (f" ({n_out} line(s) out of domain in this band)" if n_out else ""))
            elif ion == "II":
                na = ("the reactivated Amarsi network we hold is the Fe I MLP; no Fe II "
                      "3D-NLTE leg has been run in any grade (RYA-817/1106 scope).")
        if na:
            out[eng] = {"verdict": "N/A", "reason": na}
        else:
            out[eng] = {"verdict": ("LIVE" if (band, holding, ion, eng) in live else "GAP"),
                        "reason": ""}
    return out


def nearuv_codex_check(cg: pd.DataFrame) -> dict:
    """RYA-1213 Step 3 — do shallow lab-gf Fe lines exist in 3000-3780 A?

    The Codex selector is `depth <= DEPTH_HI`, which INCLUDES lines below the 0.05 floor,
    so this reproduces that comparison exactly rather than the [0.05, 0.60] window a
    reader might assume from the accounting triage.
    """
    out = {}
    for ion in ("I", "II"):
        lab = cg[(cg.species == f"Fe {ion}")
                 & cg.gf_tier.astype(str).str.contains("LAB", na=False)
                 & cg.wavelength_air_A.between(3000.0, 3780.0)]
        if lab.empty:
            out[f"Fe {ion}"] = {"lab_lines": 0, "at_or_below_gate": 0,
                                "verdict": "N/A — no LAB-tier line in the band"}
            continue
        d = _feature_depth(lab.wavelength_air_A.values.astype(float))
        sel = lab[d <= DEPTH_HI]
        n = len(sel)
        out[f"Fe {ion}"] = {
            "lab_lines": int(len(lab)),
            "at_or_below_gate": n,
            "above_gate": int((d > DEPTH_HI).sum()),
            "lines_A": [round(float(x), 3) for x in sel.wavelength_air_A.values],
            "depths": [round(float(x), 3) for x in d[d <= DEPTH_HI]],
            "verdict": (
                "NO CODEX GRADE PRODUCT IS BUILDABLE — "
                f"{n} LAB-tier Fe {ion} line(s) sit at or below the {DEPTH_HI} depth "
                f"gate, and `_cand_graded` refuses below 2: line-to-line scatter has "
                f"n-1 degrees of freedom, so at n<2 the statistical term cannot be "
                f"computed from the data at all and would be invented (RYA-1031). "
                f"This is a POPULATION fact about the band, not a missing run — the "
                f"near-UV lab pool is saturated almost to a line, which is why the "
                f"band's live products are all Deep Grade."
                if n < 2 else
                f"BUILDABLE — {n} LAB-tier Fe {ion} line(s) at or below the gate.")}
    return out


def unreachable_lab_lines(cg: pd.DataFrame) -> dict:
    """The lab lines no cell above can contain, and WHY. Reconciles the 450."""
    lab = cg[cg.species.isin(["Fe I", "Fe II"])
             & cg.gf_tier.astype(str).str.contains("LAB", na=False)]
    blue = lab[lab.wavelength_air_A.between(3780.0, 4200.0, inclusive="left")]
    gap = lab[lab.wavelength_air_A.between(12976.0, 15007.0, inclusive="neither")]
    # 🔴 THE UNION OF EVERY CELL'S LINES, NOT THE SUM OF THEIR COUNTS. Four NIR holdings
    # cover overlapping slices of one band, so adding their pools would count the same
    # physical line up to four times and the reconciliation against 450 would be
    # meaningless. This is the set of distinct wavelengths at least one cell can contain.
    reach: set = set()
    for band, ion, holding, instrument, lo, hi in CELLS:
        clo, chi, _ = _clip(holding, lo, hi)
        sel = lab[(lab.species == f"Fe {ion}") & lab.wavelength_air_A.between(clo, chi)]
        reach |= set(sel.wavelength_air_A.astype(float))
    covered = len(reach)
    unreached = lab[~lab.wavelength_air_A.astype(float).isin(reach)]
    blue_or_gap = set(blue.wavelength_air_A.astype(float)) | set(gap.wavelength_air_A.astype(float))
    other = unreached[~unreached.wavelength_air_A.astype(float).isin(blue_or_gap)]
    return {
        "total_lab_fe_lines": int(len(lab)),
        "no_synthesis_list_3780_4200_A": {
            "n": int(len(blue)),
            "by_species": {str(k): int(v) for k, v in blue.species.value_counts().items()},
            "reason": ("no synthesis line list covers 3780-4200 A. The near-UV list "
                       "(ispec_nearuv_3000_3780) stops at 3780 A and the VIS list "
                       "(iSpec-vendored GES v6) starts at 4200 A. `derive_band_products` "
                       "refuses a run it cannot cover rather than emit a product naming "
                       "a range it was not synthesised over (RYA-911/913/967). This is a "
                       "LINE LIST gap; it is not a Reference-tier gap and cannot be "
                       "closed by running anything."),
            "lines_A": [round(float(x), 3) for x in sorted(blue.wavelength_air_A.values)],
        },
        "no_holding_reaches_12976_15007_A": {
            "n": int(len(gap)),
            "reason": ("between the NIR list's red edge (12976.14 A) and the CRIRES+ H "
                       "arm holding's blue edge (15007.11 A). No solar holding we own "
                       "covers this span, so these two Ruffoni-2013 lines are not "
                       "measurable in any grade. The ticket counts them in its H-band "
                       "total of 27; the H holding contains 25."),
            "lines_A": [round(float(x), 3) for x in sorted(gap.wavelength_air_A.values)],
            "lab_source": sorted(set(gap.lab_source_tag.astype(str))),
        },
        "distinct_lab_lines_inside_a_reference_window": int(covered),
        # ⚠️ ANY REMAINDER IS A CELL THIS MATRIX DOES NOT COVER, AND IT IS LISTED RATHER
        # THAN LET TO BALANCE SILENTLY. 450 = covered + the two named populations + this;
        # a non-empty `other` means a lab line exists that no cell above reaches and that
        # neither documented reason explains — which is a GAP, not an N/A.
        "unreached_and_unexplained": {
            "n": int(len(other)),
            "lines_A": [round(float(x), 3)
                        for x in sorted(other.wavelength_air_A.values)][:60],
            "by_species": {str(k): int(v)
                           for k, v in other.species.value_counts().items()},
        },
    }


def build() -> dict:
    cg = _canon()
    feed = json.loads(FEED.read_text())
    live = {(p["band"], p["holding"], p["ion"], p["treatment"])
            for p in feed["products"] if p.get("tier") == "REFERENCE"}
    app = json.loads(RYA1208.read_text()) if RYA1208.exists() else {}
    refus = json.loads(RYA1208_REF.read_text()) if RYA1208_REF.exists() else {}

    rows = []
    for band, ion, holding, instrument, lo, hi in CELLS:
        clo, chi, span_src = _clip(holding, lo, hi)
        pool = reference_pool(cg, ion, clo, chi)
        engines = engine_verdicts(band, ion, live, holding, app, refus)
        rows.append({"band": band, "ion": f"Fe {ion}", "holding": holding,
                     "instrument": instrument,
                     "band_window_A": [lo, hi], "window_A": [clo, chi],
                     "window_source": span_src, **pool, "engines": engines})
    return {
        "ticket": "RYA-1213",
        "note": ("The Reference Grade applicability matrix. DIAGNOSTIC — no abundance "
                 "here, and nothing is calibrated (RYA-161)."),
        "definition": ("Reference Grade = every canonical_gf row whose gf_tier contains "
                       "LAB, for the species, inside the product's window, with the "
                       f"{DEPTH_HI} feature-depth gate NOT applied. RYA-946 defines the "
                       "tier by gf pedigree; depth is orthogonal to pedigree."),
        "nearuv_codex_check": nearuv_codex_check(cg),
        "lab_line_reconciliation": unreachable_lab_lines(cg),
        "cells": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="print only; write nothing")
    a = ap.parse_args()
    doc = build()

    print(f"\nRYA-1213 Reference Grade matrix  ({len(doc['cells'])} cells)\n")
    print(f"{'band':12s} {'ion':6s} {'holding':34s} {'REF':>4s} "
          f"{'(codex':>7s}{'/deep':>6s}{'/nodepth)':>10s}  live/gap/na/deferred")
    for r in doc["cells"]:
        v = [e["verdict"] for e in r["engines"].values()]
        print(f"{r['band']:12s} {r['ion']:6s} {r['holding']:34s} "
              f"{r['reference_pool']:4d} {r['would_be_codex']:7d}{r['would_be_deep']:6d}"
              f"{r['no_known_depth']:10d}  "
              f"{v.count('LIVE')}/{v.count('GAP')}/{v.count('N/A')}/"
              f"{v.count('DEFERRED')}")
    print("\nnear-UV Codex Grade check (RYA-1213 Step 3):")
    for k, s in doc["nearuv_codex_check"].items():
        print(f"  {k}: {s['lab_lines']} lab lines, {s['at_or_below_gate']} at/below the "
              f"gate -> {s['verdict'].split(' — ')[0]}")
    rec = doc["lab_line_reconciliation"]
    print(f"\nlab-line reconciliation: {rec['total_lab_fe_lines']} LAB Fe lines total; "
          f"{rec['distinct_lab_lines_inside_a_reference_window']} inside a Reference "
          f"window; {rec['no_synthesis_list_3780_4200_A']['n']} in 3780-4200 A with no "
          f"synthesis list; {rec['no_holding_reaches_12976_15007_A']['n']} in "
          f"12976-15007 A with no holding; "
          f"{rec['unreached_and_unexplained']['n']} unreached and UNEXPLAINED.")

    if not a.check:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "rya1213_reference_matrix.json").write_text(json.dumps(doc, indent=1) + "\n")
        flat = []
        for r in doc["cells"]:
            for eng, e in r["engines"].items():
                flat.append({"band": r["band"], "ion": r["ion"], "holding": r["holding"],
                             "instrument": r["instrument"], "engine": eng,
                             "reference_pool": r["reference_pool"],
                             "verdict": e["verdict"], "reason": e["reason"]})
        pd.DataFrame(flat).to_csv(OUT / "rya1213_engine_cells.csv", index=False)
        print(f"\nwrote {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
