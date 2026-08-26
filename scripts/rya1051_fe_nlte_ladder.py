#!/usr/bin/env python3
"""RYA-1051: the solar Fe I NLTE ladder as measured, and exactly how it was measured.

🔴 WHY THIS EXISTS. Four separate treatments of solar Fe I now exist, they are only
meaningful as DIFFERENCES of specific pairs, and the pairing rules are not obvious from
the product names. This script derives the whole ladder from the committed per-line
products so the numbers cannot drift from the files that carry them, and states the
process beside every value.

THE LADDER (all on ONE line set: 67 graded lab-gf Fe I lines, 4200-6910 A)

    1D-LTE                    ATLAS9.Castelli, this route's own atmosphere
      |
      |  ATMOSPHERE (1D -> mean-3D)          <3D>-LTE minus 1D-LTE
      v
    <3D>-LTE                  STAGGER mean-3D, departures WITHHELD
      |
      |  NLTE on <3D>                        <3D>-NLTE minus <3D>-LTE
      v
    <3D>-NLTE                 same atmosphere, departures ON

    synth-1D-LTE-gerber       MARCS.GES, departures WITHHELD   <-- the 1D rung's comparand
      |
      |  NLTE on 1D                          ENGINE-B-NLTE minus synth-1D-LTE-gerber
      v
    ENGINE-B-NLTE             MARCS.GES, Gerber 1D deck, departures ON

🔴 EVERY DIFFERENCE IS A PER-LINE PAIRED DIFFERENTIAL, NEVER A DIFFERENCE OF THE PUBLISHED
PRODUCT VALUES. Ryan ratified this on RYA-1040 after the graded run produced two products
that BOTH report A = 7.552 while their per-line differential is +0.032 on 66 of 67 lines --
the medians collided by coincidence (62 unique values, min gap 0.001; not a quantiser).
Subtracting the two published numbers gives a clean, well-formed, entirely false ZERO.

🔴 AND EVERY DIFFERENCE PAIRS PRODUCTS ON ONE ATMOSPHERE. Differencing <3D>-NLTE against
1D-LTE reports the 1D->mean-3D atmosphere shift as non-LTE physics (RYA-542). The 1D rung
had no same-atmosphere comparand at all until RYA-1045 added `synth-1D-LTE-gerber`; before
that, ENGINE-B-NLTE could only be differenced against an ATLAS9 product, mixing a
+0.004 atmosphere change into the answer.

⚠️ THE POOL IS PART OF THE MEASUREMENT. `--lines-tier graded` selects the 67 LAB-tier
Fe I lines AT OR BELOW the 0.6 depth gate, out of 176 in band. It is not a quality filter
bolted onto an existing pool -- it IS a different pool, and the anchor match rate goes
12% -> 64% because lab-graded lines and Amarsi's `Clean` lines are largely the same
population. `--lines-deep-graded` is the OPPOSITE selection (above the gate, the saturated
population) and must not be substituted for it.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

#: The invocations that produced the inputs. Recorded so the ladder states its own process.
INVOCATIONS = {
    "<3D> pair": (
        "scripts/derive_band_products.py --star solar --element Fe --ion I "
        "--lo 4200 --hi 6910 --force-synthesis --instrument kpno_solar_atlas "
        "--holding solar_kpno_molecfit_corrected --lines-tier graded "
        "--engine-b-deck gerber-mean3d   (and --engine-b-deck gerber-mean3d-lte)"),
    "1D pair": (
        "scripts/derive_band_products.py --star solar --element Fe --ion I "
        "--lo 4200 --hi 6910 --force-synthesis --instrument kpno_solar_atlas "
        "--holding solar_kpno_molecfit_corrected --lines-tier graded "
        "--engine-b-deck gerber-nlte     (and --engine-b-deck gerber-1d-lte)"),
}

#: Which product is differenced against which, and what the difference MEANS.
#: Keyed so a reader cannot pair them wrongly by accident.
PAIRS = {
    "atmosphere_1D_to_mean3D": ("mean3d_lte", "lte_1d_atlas9",
                                "<3D>-LTE minus 1D-LTE: the 1D->mean-3D ATMOSPHERE shift"),
    "nlte_on_mean3D": ("mean3d_nlte", "mean3d_lte",
                       "<3D>-NLTE minus <3D>-LTE on ONE atmosphere: the NLTE effect"),
    "nlte_on_1D": ("nlte_1d", "lte_1d_marcs",
                   "ENGINE-B-NLTE minus synth-1D-LTE-gerber on ONE atmosphere (MARCS.GES)"),
}


def _lines(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path)
    if "in_aggregate" in d.columns:
        d = d[d.in_aggregate.astype(str).str.lower().isin(("true", "1"))]
    return d[["wavelength_air_A", "abundance", "ep_eV"]].copy()


def paired(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    m = a.merge(b, on="wavelength_air_A", suffixes=("_a", "_b"))
    m["delta"] = m.abundance_a - m.abundance_b
    return m


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("lte-1d-atlas9", "lte-1d-marcs", "nlte-1d", "mean3d-lte", "mean3d-nlte"):
        ap.add_argument(f"--{k}", required=True, help=f"per-line CSV for {k}")
    ap.add_argument("--anchor", default=str(
        ROOT / "data" / "nlte_grids" / "amarsi2016_fe" / "amarsi2016_mean3d_solar_anchor.csv"))
    ap.add_argument("--anchor-1d", default=None,
                    help="Amarsi nmarcs_lmarcs solar slice (1D NLTE - 1D LTE); optional")
    ap.add_argument("--out", default=str(ROOT / "data" / "results" / "rya1051" /
                                         "fe_nlte_ladder.json"))
    a = ap.parse_args()

    P = {"lte_1d_atlas9": _lines(Path(a.lte_1d_atlas9)),
         "lte_1d_marcs": _lines(Path(a.lte_1d_marcs)),
         "nlte_1d": _lines(Path(a.nlte_1d)),
         "mean3d_lte": _lines(Path(a.mean3d_lte)),
         "mean3d_nlte": _lines(Path(a.mean3d_nlte))}

    doc = {"ticket": "RYA-1051",
           "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "star": "solar", "species": "Fe I", "band_A": [4200, 6910],
           "pool": "--lines-tier graded: the 67 LAB-tier Fe I lines AT OR BELOW the 0.6 "
                   "depth gate, of 176 in band. NOT --lines-deep-graded, which is the "
                   "saturated population above the gate (RYA-984/RYA-954).",
           "holding": "kpno_solar_atlas / solar_kpno_molecfit_corrected",
           "invocations": INVOCATIONS,
           "differencing_rule": "PER-LINE PAIRED, then aggregate. Never the difference of "
                                "the published product values (RYA-1040 median collision), "
                                "and never across two atmospheres (RYA-542).",
           "products": {}, "differences": {}}

    print("PRODUCTS (median of the per-line abundances, on the graded 67)")
    for k, d in P.items():
        doc["products"][k] = {"n": int(len(d)),
                              "median": round(float(d.abundance.median()), 4),
                              "mean": round(float(d.abundance.mean()), 4),
                              "std": round(float(d.abundance.std(ddof=1)), 4)}
        print(f"  {k:<16} n={len(d):>3}  A = {d.abundance.median():.4f}")

    print("\nDIFFERENCES (per-line paired)")
    for name, (hi, lo, why) in PAIRS.items():
        m = paired(P[hi], P[lo])
        v = m.delta
        doc["differences"][name] = {
            "minuend": hi, "subtrahend": lo, "meaning": why,
            "n_paired": int(len(m)),
            "median": round(float(v.median()), 4), "mean": round(float(v.mean()), 4),
            "min": round(float(v.min()), 4), "max": round(float(v.max()), 4),
            "n_nonzero": int((v.abs() > 1e-6).sum()),
            # 🔴 the collision check, stated for every pair rather than only where it bit
            "difference_of_published_medians": round(
                float(P[hi].abundance.median() - P[lo].abundance.median()), 4)}
        e = doc["differences"][name]
        flag = ("   🔴 MEDIAN COLLISION: differencing the published values gives "
                f"{e['difference_of_published_medians']:+.4f}"
                if abs(e["median"] - e["difference_of_published_medians"]) > 0.005 else "")
        print(f"  {name:<26} n={len(m):>3}  {v.median():+.4f}{flag}")

    # the increment that exonerated the <3D> deck (RYA-1045)
    m3 = paired(P["mean3d_nlte"], P["mean3d_lte"])[["wavelength_air_A", "delta"]]
    m1 = paired(P["nlte_1d"], P["lte_1d_marcs"])[["wavelength_air_A", "delta"]]
    inc = m3.merge(m1, on="wavelength_air_A", suffixes=("_3", "_1"))
    inc["inc"] = inc.delta_3 - inc.delta_1
    doc["differences"]["increment_mean3D_minus_1D"] = {
        "meaning": "how much the <3D> deck ADDS over the 1D deck, per line. RYA-1045: this "
                   "is the quantity that agrees with Amarsi and exonerates the <3D> deck.",
        "n_paired": int(len(inc)), "median": round(float(inc.inc.median()), 4)}
    print(f"  {'increment_mean3D_minus_1D':<26} n={len(inc):>3}  {inc.inc.median():+.4f}")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
