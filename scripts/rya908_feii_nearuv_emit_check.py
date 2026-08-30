#!/usr/bin/env python3
"""RYA-908 — the per-line layer of the near-UV Fe II Deep Graded cell, and its two gates.

    python3 scripts/rya908_feii_nearuv_emit_check.py [--check]

This cell had only an aggregate row: a number with no showable derivation. The per-line,
budget and provenance artifacts were retained in the Mac worktree the feed's provenance
NAMES (`~/codex/rya1015`), so this was an EMIT of verified-identical files, never a re-run.

🔴 GATE 1 — THE AGGREGATE REPRODUCES BY IDENTITY, NOT BY RE-DERIVATION. The retained
`_products.csv` sha256 matches the sha256 the feed recorded, on both holdings, so the
emitted per-line rows are provably the derivation of the PUBLISHED number rather than of a
lookalike re-run. Nothing was re-fitted and no abundance moved (RYA-161).

🔴 GATE 2 — LTE-EQUIVALENCE, ON RYA-908's CORRECTED CRITERION. The first criterion asked for
"departure = 1", which this schema cannot express: the column is `nlte_delta_dex`, an
ADDITIVE dex correction where no departure is 0.0 and 1 would be a 10x correction. Ryan
restated it (2026-08-30) as a two-part test with a SOURCED bar:

  1. MAGNITUDE  |nlte_delta_dex| <= the literature's own solar Fe II correction
  2. SHAPE      the correction is UNSTRUCTURED -- uniform across served lines, no outlier

⚠️ THE TOLERANCE DOES NOT DEFANG THE PHANTOM DETECTOR, and that is the point of the shape
half. A phantom -- a Fe II line mis-identified against a Fe I level -- would be LARGE and/or
STRUCTURED, so it breaches (1) or (2) loudly. A strict `== 0` would instead have REJECTED
the true near-LTE physics as a phantom: a threshold too tight to admit reality.

MEASURED: -0.001 dex on all 7 served lines, both holdings, perfectly uniform. Five lines
are NOT SERVED by the departure source, and the SAME five wavelengths on both holdings --
a coverage property of the source, not of either spectrum.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BAND = ROOT / "data" / "results" / "band_products"
FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
OUT = ROOT / "data" / "results" / "rya908" / "feii_nearuv_emit.json"

STEM = "FeII_3000_3780_kpno_solar_atlas_solar_kpno_{h}_corrected_SYNTH_DEEPGRADED"
HOLDINGS = {"kurucz2005": "solar_kpno_kurucz2005_corrected",
            "molecfit": "solar_kpno_molecfit_corrected"}

#: 🔴 SOURCED, NOT CHOSEN. Amarsi's own solar Fe II corrections are -0.0010 / -0.0008 dex
#: and Lind et al. 2017 calls solar Fe II departures insignificant. A departure at or below
#: the literature's own solar Fe II value IS LTE-equivalent by definition, so the bar is a
#: citation rather than an epsilon someone picked to make a test pass (RYA-161).
LTE_EQUIVALENT_TOLERANCE_DEX = 0.0010
LTE_TOLERANCE_SOURCE = ("Amarsi's published solar Fe II corrections (-0.0010 / -0.0008 dex) "
                        "and Lind et al. 2017 ('insignificant' for solar Fe II)")
#: The shape half: served-line departures must not spread, or a single line is driving it.
SHAPE_MAX_SPREAD_DEX = 0.0005

#: RYA-1133's measured floor: the two holdings return 7.957 vs 7.697 on the SAME 12 lines,
#: so half the spread is 0.130 dex -- ABOVE the 0.10 band-flat pseudo-continuum term the
#: run charged. Interim and a LOWER bound: two realisations are not the space of them, and
#: all three near-UV atlases are Kitt Peak so the disagreement cannot be adjudicated.
HONEST_FLOOR_DEX = 0.130
RUN_BUDGET_SYST_DEX = 0.1095


def gate1_identity() -> dict:
    """The emitted per-line rows belong to the PUBLISHED aggregate, proven by sha256."""
    feed = json.loads(FEED.read_text())
    out = {}
    for h, holding in HOLDINGS.items():
        p = BAND / (STEM.format(h=h) + "_products.csv")
        rec = next(r for r in feed["products"]
                   if r["band"] == "near-UV" and r["ion"] == "II"
                   and r["treatment"] == "1D-LTE" and r["holding"] == holding)
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        out[holding] = {"artifact_sha256": got,
                        "feed_recorded_sha256": rec["provenance"]["sha256"],
                        "identical": got == rec["provenance"]["sha256"],
                        "published_A": rec["A"]}
    return out


def gate2_lte_equivalence() -> dict:
    """Ryan's corrected criterion: sourced magnitude AND unstructured shape."""
    out = {}
    for h, holding in HOLDINGS.items():
        d = pd.read_csv(BAND / (STEM.format(h=h) + "_ENGINE-A_lines.csv"))
        served = d[d.in_aggregate.astype(bool)]
        dep = served.nlte_delta_dex.astype(float)
        unserved = d[~d.in_aggregate.astype(bool)]
        mag_ok = bool((dep.abs() <= LTE_EQUIVALENT_TOLERANCE_DEX + 1e-12).all())
        spread = float(dep.max() - dep.min())
        out[holding] = {
            "n_rows": int(len(d)), "n_served": int(len(served)),
            "departures": sorted({round(float(x), 6) for x in dep}),
            "max_abs_departure_dex": round(float(dep.abs().max()), 6),
            "spread_dex": round(spread, 6),
            "magnitude_within_sourced_tolerance": mag_ok,
            "shape_unstructured": bool(spread <= SHAPE_MAX_SPREAD_DEX),
            "not_served_wavelengths": sorted(unserved.wavelength_air_A.round(3).tolist()),
            "verdict": "LTE-EQUIVALENT" if mag_ok and spread <= SHAPE_MAX_SPREAD_DEX
                       else "🔴 BREACH — phantom check fires",
        }
    ks = list(out.values())
    out["_same_lines_unserved_on_both_holdings"] = (
        ks[0]["not_served_wavelengths"] == ks[1]["not_served_wavelengths"])
    return out


def pool_vs_nlte() -> dict:
    """🔴 The decomposition the per-line layer exists to make possible (RYA-1055/1113-F5).

    The aggregate is the MEDIAN -- verified against the published value before anything is
    built on it, because a decomposition computed with the wrong estimator would be
    confidently wrong.
    """
    feed = json.loads(FEED.read_text())
    out = {}
    for h, holding in HOLDINGS.items():
        L = pd.read_csv(BAND / (STEM.format(h=h) + "_1D-LTE_lines.csv"))
        N = pd.read_csv(BAND / (STEM.format(h=h) + "_ENGINE-A_lines.csv"))
        L, N = L[L.in_aggregate.astype(bool)], N[N.in_aggregate.astype(bool)]
        sub = L[L.wavelength_air_A.isin(set(N.wavelength_air_A))]
        m_all, m_sub, m_nl = (float(np.median(L.abundance)), float(np.median(sub.abundance)),
                              float(np.median(N.abundance)))
        pub_l = next(r["A"] for r in feed["products"] if r["band"] == "near-UV"
                     and r["ion"] == "II" and r["treatment"] == "1D-LTE"
                     and r["holding"] == holding)
        pub_n = next(r["A"] for r in feed["products"] if r["band"] == "near-UV"
                     and r["ion"] == "II" and r["treatment"] == "ENGINE-A"
                     and r["holding"] == holding)
        out[holding] = {
            "median_1D_LTE_all12": round(m_all, 4), "published_1D_LTE": pub_l,
            "median_ENGINE_A_on7": round(m_nl, 4), "published_ENGINE_A": pub_n,
            "median_1D_LTE_on_the_same_7": round(m_sub, 4),
            "pool_effect_dex": round(m_sub - m_all, 4),
            "nlte_effect_paired_dex": round(m_nl - m_sub, 4),
            "published_delta_dex": round(pub_n - pub_l, 4),
            "estimator_reproduces_published": bool(round(m_all, 3) == round(pub_l, 3)
                                                   and round(m_nl, 3) == round(pub_n, 3)),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    g1, g2, dec = gate1_identity(), gate2_lte_equivalence(), pool_vs_nlte()
    doc = {
        "ticket": "RYA-908",
        "what": ("emit of the retained per-line/budget/provenance layer for the near-UV "
                 "Fe II Deep Graded synthesis cell; no re-run, no abundance move"),
        "gate1_aggregate_identity": g1,
        "gate2_lte_equivalence": {
            "criterion": ("|nlte_delta_dex| <= sourced tolerance AND unstructured shape "
                          "(RYA-908 corrected criterion, 2026-08-30)"),
            "tolerance_dex": LTE_EQUIVALENT_TOLERANCE_DEX,
            "tolerance_source": LTE_TOLERANCE_SOURCE,
            "shape_max_spread_dex": SHAPE_MAX_SPREAD_DEX,
            "per_holding": g2,
        },
        "pool_vs_nlte_decomposition": dec,
        "bar": {
            "run_budget_syst_dex": RUN_BUDGET_SYST_DEX,
            "published_floor_syst_dex": HONEST_FLOOR_DEX,
            "basis": ("RYA-1133: the two holdings disagree by 0.260 dex on the SAME 12 "
                      "lines, so half-spread 0.130 EXCEEDS the 0.10 band-flat "
                      "pseudo-continuum term the run charged"),
            "status": ("INTERIM LOWER BOUND — two realisations are not the space of them, "
                       "and all three near-UV atlases are Kitt Peak, so the disagreement "
                       "can be measured but never adjudicated. Pending RYA-1133's per-ion "
                       "term."),
        },
    }
    print(json.dumps(doc, indent=2)[:2600])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")

    if a.check:
        bad = []
        for k, v in g1.items():
            if not v["identical"]:
                bad.append(f"{k}: emitted per-line layer does not belong to the published "
                           f"aggregate (sha256 mismatch)")
        for k, v in g2.items():
            if k.startswith("_"):
                continue
            if not v["magnitude_within_sourced_tolerance"]:
                bad.append(f"{k}: departure exceeds the sourced solar Fe II tolerance — "
                           f"PHANTOM check fires")
            if not v["shape_unstructured"]:
                bad.append(f"{k}: departures are STRUCTURED (spread {v['spread_dex']}) — "
                           f"a single line may be driving it")
        if not g2["_same_lines_unserved_on_both_holdings"]:
            bad.append("the unserved set differs between holdings — it is no longer a "
                       "pure coverage property of the departure source")
        for k, v in dec.items():
            if not v["estimator_reproduces_published"]:
                bad.append(f"{k}: the median no longer reproduces the published value — "
                           f"the decomposition is built on the wrong estimator")
        if bad:
            print("\nFAIL:")
            for b in bad:
                print("  -", b)
            return 1
        print("\nOK — aggregate identity proven, LTE-equivalence passes on the sourced "
              "criterion, decomposition reproduces the published values.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
