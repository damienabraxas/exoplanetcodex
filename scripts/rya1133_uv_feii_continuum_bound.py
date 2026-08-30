#!/usr/bin/env python3
"""RYA-1133 — why near-UV Fe II breaches the 0.10 pseudo-continuum bound, and what is owed.

    python3 scripts/rya1133_uv_feii_continuum_bound.py [--check]

READ-ONLY (RYA-161). Reads the feed, canonical_gf, RYA-846's artifact and the pipeline's own
constants. Changes no value, line, gf or threshold. The firewall is a SYSTEMATIC-BOUND
question only: the bar is too small, and the fix is an honest LARGER bar — never chasing the
two holdings together (RYA-777, frontier band).

THE BREACH (RYA-1113): the same 12 Fe II lines return A = 7.957 on kurucz2005 and 7.697 on
molecfit — 0.260 dex apart, half-spread 0.130, against a band-flat 0.10 dex bound meant to
cover exactly this. Fe I does not breach (0.046 apart, half-spread 0.023).

Q1 — IS RYA-846'S MEASURED TERM APPLIED TO FE II?
-------------------------------------------------
🔴 NO, AND NOT TO ANYTHING. `error_budget.py:603` charges `Term("pseudo-continuum", 0.10)`
unconditionally wherever the BandPolicy says the continuum is pseudo. RYA-846's measurement
is referenced by nothing outside RYA-841/846's own audit scripts — asserted here by grepping
the tree, not by reading a comment. 846's own artifact records `assumed_term_dex: 0.1` and a
`revised_cells` list containing ONLY Fe I: Fe II near-UV was a zero-product band at the time
(RYA-908), so it was never in scope to revise.

Q2 — DOES 846'S TERM COVER THE kurucz2005-vs-molecfit SPREAD?
--------------------------------------------------------------
🔴 NO — IT MEASURES A DIFFERENT PAIR. 846's method line says it compared

    the 1984 Kitt Peak atlas (Kurucz/Brault `lm` files)   vs   Wallace+2011 (NSO)

while the breach is between

    the Kurucz 2005 irradiance distribution              vs   the 1984 Kitt Peak atlas

⚠️ Kurucz 2005 appears in NEITHER of 846's spectra. Its term is a KP-vs-NSO *reduction*
disagreement; the breach is a *source-spectrum* disagreement between two different products.
They are not the same quantity, and 846 never measured the second one.

⚠️ AND APPLYING 846 WOULD MOSTLY SHRINK THE BAR, NOT GROW IT. Every candidate term it
publishes sits at or below the observed 0.130 half-spread, and three of five sit below the
0.10 already charged — so "apply 846" is not a route to an honest larger bar.

WHY FE II BREACHES AND FE I DOES NOT — MEASURED, NOT ARGUED
-----------------------------------------------------------
The two ions' lab-gf pools do not sample the same band. Fe II's twelve lines are a narrow
clump ~300 A bluer than Fe I's, entirely inside the worst part of the band, where the flux is
lowest and the blanketing heaviest. Fe I averages across the band; Fe II samples only the
zone where two continuum realisations have the most room to disagree.

🔴 AND THE WHOLE BAND SITS BELOW THE REPO'S OWN CONTINUUM-LOCATABLE EDGE.
`normalization_intake.CONTINUUM_LIMITED_BLUE_EDGE_A = 4500 A` is where a continuum can be
located AT ALL — cross-derived twice (RYA-451/460's blue-edge no-true-continuum class, and
`diagnostics_abundance.BLUE_EDGE`). `detect()` REFUSES to classify any spectrum lying
entirely below it, returning UNKNOWN because "near-UV line blanketing leaves no true
continuum in any window ... the known-normalised KP1984 atlas fails this test on itself below
4200 A". This script demonstrates that refusal against a VIS control, so it is the
instrument's behaviour and not a reading of its docstring.

⚠️ THAT REFINES RYA-1113's F4. F4 reported that neither holding's normalisation was probed
inside 3000-3780 A and read it as a gap. It is not an oversight: the probe would REFUSE
there by construction. The honest statement is stronger — in this band "is it normalised?"
is unanswerable by the standard instrument, which is precisely the vacuum a band-flat
placeholder filled.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import normalization_intake as NI     # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
RYA846 = ROOT / "data" / "results" / "rya846" / "rya846_sigma_delta.json"
OUT = ROOT / "data" / "results" / "rya1133" / "uv_feii_continuum_bound.json"

BAND_LO_A, BAND_HI_A = 3000.0, 3780.0
PSEUDO_CONTINUUM_DEX = 0.10


def term_application() -> dict:
    """Q1 — does RYA-846's measured term reach any product? Grepped, not assumed."""
    def grep(pattern: str) -> list[str]:
        r = subprocess.run(["git", "grep", "-l", pattern, "--", "pipeline", "scripts"],
                           cwd=ROOT, capture_output=True, text=True)
        return sorted(x for x in r.stdout.split() if x)

    readers = grep("rya846")
    # 🔴 THE INSTRUMENT MUST NOT ENTER ITS OWN MEASUREMENT. This script greps for
    # "rya846" and CONTAINS that string, so the moment it is tracked it appears in its own
    # result set — and, not being matched by the rya846/rya841 filter below, it would land
    # in `consumers` and flip the Q1 verdict from NOT APPLIED to APPLIED. That is a false
    # conclusion produced by the act of measuring. It excludes itself BY NAME, never by a
    # pattern like "rya11*", which would hide a genuine future consumer.
    #
    # ⚠️ It only becomes true at `git add` — `git grep` sees TRACKED files — so this was
    # fixed before staging rather than discovered by a red suite afterwards.
    SELF = "scripts/rya1133_uv_feii_continuum_bound.py"
    readers = [p for p in readers if p != SELF]
    # the audit scripts that MEASURED it are the instrument, not the product path
    own = [p for p in readers if "rya846" in p or "rya841" in p]
    consumers = [p for p in readers if p not in own]
    charged = [p for p in grep('Term("pseudo-continuum"') if p != SELF]
    return {
        "files_referencing_rya846": readers,
        "its_own_audit_scripts": own,
        "product_path_consumers": consumers,
        "verdict": ("NOT APPLIED — RYA-846's measured term is read by nothing outside its "
                    "own audit scripts" if not consumers else
                    f"APPLIED via {consumers}"),
        "what_is_actually_charged": {
            "site": charged,
            "value_dex": PSEUDO_CONTINUUM_DEX,
            "shape": "band-flat, unconditional wherever BandPolicy.continuum_treatment "
                     "is pseudo — identical for both ions, both holdings, every window",
        },
    }


def rya846_scope() -> dict:
    """Q2 — what did 846 actually compare, and how big are its terms vs the breach?"""
    d = json.loads(RYA846.read_text())
    terms = {k: round(v["term_dex"], 4) for k, v in d["terms_dex"].items()}
    terms["net_at_product_lines"] = round(d["decomposition"]["term_net_dex"], 4)
    terms["raw_at_product_lines"] = round(d["decomposition"]["term_raw_dex"], 4)
    revised = [f"Fe {c['ion']} {c['treatment']}" for c in d["revised_cells"]]
    return {
        "method": d["method"],
        "compared_pair": ["1984 Kitt Peak atlas (Kurucz/Brault lm files)",
                          "Wallace+2011 (NSO)"],
        "breach_pair": ["Kurucz 2005 irradiance distribution",
                        "1984 Kitt Peak atlas (molecfit-corrected)"],
        "kurucz2005_in_846s_comparison": False,
        "assumed_term_dex": d["assumed_term_dex"],
        "cells_846_revised": revised,
        "revised_any_Fe_II": any("Fe II" in c for c in revised),
        "candidate_terms_dex": terms,
        "lever_dex_per_unit_delta": d["lever"]["value"],
    }


def pool_geometry() -> dict:
    """Why Fe II breaches and Fe I does not: the two pools sample different sub-bands."""
    c = pd.read_csv(CANON, low_memory=False)
    uv = c[(c.wavelength_air_A >= BAND_LO_A) & (c.wavelength_air_A < BAND_HI_A)]
    out = {}
    for ion in ("I", "II"):
        w = np.sort(uv[(uv.species.astype(str) == f"Fe {ion}")
                       & (uv.gf_tier.astype(str) == "LAB")].wavelength_air_A.values)
        out[f"Fe {ion}"] = {
            "n_LAB": int(len(w)),
            "span_A": [round(float(w.min()), 2), round(float(w.max()), 2)],
            "median_A": round(float(np.median(w)), 1),
            "p25_A": round(float(np.percentile(w, 25)), 1),
            "p75_A": round(float(np.percentile(w, 75)), 1),
            "frac_below_3300A": round(float((w < 3300).mean()), 3),
        }
    out["contrast"] = ("Fe II is a narrow clump ~300 A bluer than Fe I and lies ENTIRELY "
                       "below 3300 A, where the flux is lowest and the blanketing "
                       "heaviest. Fe I averages across the band. That is why the same "
                       "two holdings disagree by 0.260 dex on Fe II and 0.046 on Fe I.")
    return out


def continuum_locatability() -> dict:
    """🔴 Demonstrate the instrument's refusal below 4500 A, against a VIS control.

    The control matters: without it, 'UNKNOWN in the near-UV' could be a property of the
    synthetic flux rather than of the wavelength regime. Same flux, two windows.
    """
    def probe(lo, hi):
        w = np.linspace(lo, hi, 4000)
        f = np.ones_like(w) * 0.9 + 0.1 * np.sin(w)
        e = NI.detect(w, f)
        return {"window_A": [lo, hi], "value": e.value, "reason": (e.reason or "")[:220]}

    feii = probe(3000.0, 3280.0)
    fei = probe(3026.0, 3780.0)
    vis = probe(5000.0, 5300.0)
    return {
        "blue_edge_A": NI.CONTINUUM_LIMITED_BLUE_EDGE_A,
        "band_entirely_below_edge": BAND_HI_A < NI.CONTINUUM_LIMITED_BLUE_EDGE_A,
        "probe_FeII_pool_window": feii,
        "probe_FeI_span": fei,
        "control_VIS_window": vis,
        "control_is_valid": bool(feii["value"] == NI.UNKNOWN
                                 and vis["value"] == NI.NORMALISED),
        "refines_RYA1113_F4": (
            "F4 read 'never probed in 3000-3780 A' as an oversight. It is not: the probe "
            "REFUSES there by construction. In this band 'is it normalised?' is "
            "unanswerable by the standard instrument — which is the vacuum the band-flat "
            "placeholder filled."),
    }


def source_census() -> dict:
    """🔴 WHICH SOURCES ARE ACTUALLY IN PLAY — and the answer is: Kitt Peak, all of them.

    Both near-UV holdings sit on `kpno_solar_atlas`. No HARPS, IAG or CRIRES+ product
    reaches this band in the live feed — HARPS is VIS-only there. So the "two holdings"
    are two SOURCE PRODUCTS on one instrument:

        solar_kpno_kurucz2005_corrected   Kurucz 2005 irradiance distribution
        solar_kpno_molecfit_corrected     the 1984 KP atlas, molecfit-corrected

    ⚠️ AND SO IS THE THIRD SPECTRUM. RYA-846's "independent reduction", Wallace, Hinkle,
    Livingston & Davis 2011, is itself a KPNO solar flux atlas. Every spectrum anyone has
    compared in this band is Kitt Peak, so there is no genuinely independent observatory
    to triangulate against: the disagreement can be MEASURED but not adjudicated, and no
    holding can be shown to be the right one. That is what makes the observed spread a
    floor rather than an error estimate.
    """
    feed = json.loads(FEED.read_text())
    uv = [p for p in feed["products"] if p["band"] == "near-UV"]
    by_inst = sorted({p["instrument"] for p in uv})
    bands = {}
    for p in feed["products"]:
        bands.setdefault(p["instrument"], set()).add(p["band"])
    return {
        "near_UV_instruments": by_inst,
        "near_UV_holdings": sorted({p["holding"] for p in uv}),
        "is_single_instrument": len(by_inst) == 1,
        "bands_per_instrument": {k: sorted(v) for k, v in sorted(bands.items())},
        "harps_reaches_near_UV": "near-UV" in bands.get("harps", set()),
        "harps_blue_limit_A": 3780.0,
        "why_harps_cannot_help": (
            "HARPS is VIS-only BY DESIGN, not by omission: config/constants.py PIPELINE"
            "['wav_min_A'] = 3780.0, described there as the 'full HARPS optical range' — "
            "which is EXACTLY the near-UV/VIS boundary. The band is defined where HARPS "
            "stops, so no HARPS cross-check of this band can ever exist."),
        "the_two_KP_sources": {
            "solar_kpno_molecfit_corrected": "the 1984 KP atlas, ESO molecfit-corrected",
            "solar_kpno_kurucz2005_corrected": "the Kurucz 2005 irradiance distribution",
        },
        "note": ("all three spectra ever compared in this band — the 1984 KP atlas, the "
                 "Kurucz 2005 irradiance product and RYA-846's Wallace+2011 NSO atlas — "
                 "are Kitt Peak. KP is the ONLY instrument that reaches the band at all. "
                 "So there is no independent observatory in band and the disagreement can "
                 "be measured but never adjudicated: neither source can be shown to be "
                 "the right one, which is what makes the observed spread a floor rather "
                 "than an error estimate."),
    }


def breach() -> dict:
    feed = json.loads(FEED.read_text())
    uv = [p for p in feed["products"] if p["band"] == "near-UV"]
    out = {}
    for ion in ("I", "II"):
        for t in sorted({p["treatment"] for p in uv if p["ion"] == ion}):
            sel = {("kurucz2005" if "kurucz2005" in p["holding"] else "molecfit"): p
                   for p in uv if p["ion"] == ion and p["treatment"] == t}
            if len(sel) != 2:
                continue
            s = float(sel["kurucz2005"]["A"]) - float(sel["molecfit"]["A"])
            out[f"Fe {ion} {t}"] = {
                "A_kurucz2005": sel["kurucz2005"]["A"], "A_molecfit": sel["molecfit"]["A"],
                "n_lines": sel["kurucz2005"]["n_lines"],
                "spread_dex": round(s, 4), "half_spread_dex": round(abs(s) / 2, 4),
                "breaches_bound": bool(abs(s) / 2 > PSEUDO_CONTINUUM_DEX),
            }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    q1, q2 = term_application(), rya846_scope()
    geom, loc, br = pool_geometry(), continuum_locatability(), breach()
    src = source_census()
    worst = max(v["half_spread_dex"] for v in br.values())
    over = {k: v for k, v in q2["candidate_terms_dex"].items() if v >= worst}

    doc = {
        "ticket": "RYA-1133",
        "read_only": ("no value, line, gf or threshold changed (RYA-161); this reads the "
                      "feed, canonical_gf, RYA-846's artifact and pipeline constants"),
        "breach": br,
        "Q1_is_846s_term_applied": q1,
        "Q2_does_846_cover_this_pair": q2,
        "source_census": src,
        "pool_geometry": geom,
        "continuum_locatability": loc,
        "verdict": {
            "answer": "(b) — a distinct term is owed, with a (c) floor",
            "why_not_a": ("846's term is applied to NOTHING, so 'apply it' is available — "
                          "but it would not fix this: it measures 1984KP-vs-Wallace, a "
                          "different pair, and most of its variants are SMALLER than the "
                          "0.10 already charged, so applying it shrinks the bar."),
            "why_b": ("the Kurucz2005-vs-1984KP source-spectrum disagreement has never "
                      "been measured by anyone. It is the direct analogue of what 846 did, "
                      "for the pair that actually underlies these two holdings."),
            "why_c_floor": (f"until that term exists, the honest floor is the OBSERVED "
                            f"holding spread — half-spread {worst} dex on Fe II — and that "
                            f"is itself a LOWER bound: it is two realisations, not the "
                            f"space of them. Reporting it is report-and-note (RYA-777), "
                            f"not chasing the holdings together."),
            "846_terms_at_or_above_the_observed_breach": over or
                "NONE — every RYA-846 candidate term is below the observed half-spread",
        },
        "owed_measurement": {
            "what": ("sigma_delta between the Kurucz 2005 irradiance distribution and the "
                     "1984 Kitt Peak atlas across 3000-3780 A, evaluated AT THE PRODUCT "
                     "LINES, per ion — mirroring RYA-846's method for the correct pair"),
            "why_per_ion": ("the two lab-gf pools sample different sub-bands (Fe II 100% "
                            "below 3300 A, Fe I median 3506.6 A), so a band-flat number "
                            "cannot serve both — which is the defect that produced this "
                            "breach in the first place"),
            "caveat": ("both spectra sit entirely below the CONTINUUM_LIMITED blue edge, "
                       "so the result is a REDUCTION-DISAGREEMENT lower bound, never a "
                       "distance from a true continuum — the same honesty RYA-846 "
                       "attached to its own number"),
            "blocked_on": ("Fe II near-UV ships no *_lines.csv (RYA-1113 F6), so the "
                           "at-product-lines evaluation cannot be done for Fe II from the "
                           "committed tree"),
        },
    }
    print(json.dumps({k: doc[k] for k in ("breach", "verdict")}, indent=2)[:2000])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {a.out.relative_to(ROOT)}")

    if a.check:
        bad = []
        if q1["product_path_consumers"]:
            bad.append("RYA-846's term now HAS a product-path consumer — Q1 has changed")
        if q2["kurucz2005_in_846s_comparison"]:
            bad.append("846 now compares Kurucz2005 — Q2 has changed")
        if not loc["control_is_valid"]:
            bad.append("the locatability control is vacuous: the near-UV probe and the "
                       "VIS control no longer differ, so UNKNOWN proves nothing")
        if not src["is_single_instrument"] or src["harps_reaches_near_UV"]:
            bad.append("a second instrument now reaches the near-UV — the 'no independent "
                       "observatory in band' limitation has changed")
        if not br["Fe II 1D-LTE"]["breaches_bound"]:
            bad.append("Fe II no longer breaches the bound — the premise has moved")
        if bad:
            print("\nFAIL:")
            for b in bad:
                print("  -", b)
            return 1
        print("\nOK — Q1/Q2 hold, the locatability control is non-vacuous, breach stands.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
