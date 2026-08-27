#!/usr/bin/env python3
"""RYA-1071: what is the ROUTE systematic for solar Fe I? Both routes, one pool.

🔴 WHY THIS COULD NOT BE READ OFF EXISTING PRODUCTS. RYA-1044's RYA-783 reconciliation
found the EW-route and synthesis-route Fe I pools shared TWO candidate lines out of 159 and
67, so their 0.075 dex gap decomposed into pool, route, gf rung and chi2 gate at once. The
only way to isolate the route is to run both on ONE pool, which is what this measures.

    route A  PROFILEFIT  measured equivalent widths, inverted to abundance
    route B  SYNTH       the observed flux, fitted directly

Same star, same instrument, same holding, same 67 graded lab-gf lines, same gf rung 3.

🔴 THE ANSWER DEPENDS ON HOW YOU ASK, AND ALL THREE NUMBERS ARE REPORTED.

    difference of published medians   +0.0160   <- 7.463 (n=13) vs 7.447 (n=67): DIFFERENT
                                                  POOLS, so this one is not a route
                                                  difference at all
    difference of medians, shared     +0.0679
    PER-LINE PAIRED median            +0.0540   <- the one that means what it says

The first is the trap RYA-1040 named: two products whose aggregates are differenced when
their pools are not the same set. It is quoted here only to show what it would have said.

⚠️ AND THE MEDIAN HIDES THE SPREAD, WHICH IS THE REAL FINDING.
sd 0.1755, range -0.246 .. +0.420 across 13 lines. The routes agree on AVERAGE to about
0.05 dex and disagree LINE BY LINE by up to 0.42. Reporting +0.054 alone would say the
routes are nearly equivalent, and they are not.

⚠️ ONLY 13 OF 67 LINES SURVIVE BOTH ROUTES -- AN 81% ATTRITION, AND IT IS NOT SYMMETRIC.
The synthesis route accepted all 67. The EW route:

    12  absent from the line-accounting table -- it cannot even consider them
    24  REW above the saturation ceiling (RYA-711) -- the EW->abundance inversion runs
        along the flat curve of growth and is ill-conditioned in BOTH directions
    10  OVER-PHYSICAL-WIDTH -- fitted Voigt sigma beyond physical
     8  BLEND-DOMINATED -- deepest feature offset from the target

Those are the route's own physics constraints, not noise, and they are why a flux fit and
an EW inversion are different instruments rather than two ways of computing one number.
⚠️ Note what that does to the comparison: the 13 survivors are the EW route's EASIEST
lines, so this delta is measured where the EW route is at its most reliable. That is a
selection effect INSIDE the comparison and it bounds how far the number generalises.

FIREWALL (RYA-161). Delta_route is a MEASUREMENT OF A SYSTEMATIC, not a correction.
Neither route is "right". This number belongs in the reported uncertainty layer (RYA-282);
it must never be used to nudge either route toward the other, toward Amarsi, or toward 7.46.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MATCH_TOL_A = 0.05      # the two routes round wavelengths differently: 4372.9817 vs .9820


def _agg(d: pd.DataFrame) -> pd.DataFrame:
    if "in_aggregate" in d.columns:
        return d[d.in_aggregate.astype(str).str.lower().isin(("true", "1"))]
    return d


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ew-lines", required=True, help="PROFILEFIT 1D-LTE per-line CSV")
    ap.add_argument("--synth-lines", required=True, help="SYNTH 1D-LTE per-line CSV")
    ap.add_argument("--ew-measured", default=None,
                    help="the measure_band_profilefit EW file, for the attrition census")
    ap.add_argument("--pool-size", type=int, default=67)
    ap.add_argument("--out", default=str(ROOT / "data" / "results" / "rya1071" /
                                         "route_systematic.json"))
    a = ap.parse_args()

    ew, sy = _agg(pd.read_csv(a.ew_lines)), _agg(pd.read_csv(a.synth_lines))
    sw, sa = sy.wavelength_air_A.values, sy.abundance.values
    rows = []
    for _, r in ew.iterrows():
        d = np.abs(sw - r.wavelength_air_A)
        if len(d) and d.min() <= MATCH_TOL_A:
            i = int(np.argmin(d))
            rows.append({"wave_A": float(r.wavelength_air_A),
                         "A_ew": float(r.abundance), "A_synth": float(sa[i]),
                         "delta": float(r.abundance - sa[i])})
    if not rows:
        raise SystemExit("no line is in BOTH routes' aggregates -- no route delta exists")
    m = pd.DataFrame(rows)

    doc = {"ticket": "RYA-1071",
           "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "quantity": "Delta_route = A(PROFILEFIT/EW) - A(SYNTH/flux fit), same lines",
           "pool": "the graded 67 (--lines-tier graded, 4200-6910 A, corrected KPNO)",
           "n_pool": a.pool_size, "n_both_routes": len(m),
           "delta_route": {
               "per_line_paired_median": round(float(m.delta.median()), 4),
               "mean": round(float(m.delta.mean()), 4),
               "sd": round(float(m.delta.std(ddof=1)), 4),
               "min": round(float(m.delta.min()), 4),
               "max": round(float(m.delta.max()), 4)},
           "other_ways_of_asking": {
               "difference_of_medians_on_shared_lines":
                   round(float(m.A_ew.median() - m.A_synth.median()), 4),
               "note": "a difference of PUBLISHED product medians is NOT a route "
                       "difference here -- the two products aggregate different pools "
                       "(13 vs 67). RYA-1040."},
           "lines": rows}
    if a.ew_measured:
        meas = pd.read_csv(a.ew_measured)
        q = meas[~meas.in_aggregate.astype(bool)] if "in_aggregate" in meas else meas.iloc[0:0]
        reasons = (q["excluded_reason"].astype(str)
                   .str.replace(r"[\d.\-]+", "N", regex=True).str.slice(0, 46)
                   .value_counts().to_dict()) if "excluded_reason" in q else {}
        doc["attrition"] = {
            "requested": a.pool_size, "measurable_by_ew": int(len(meas)),
            "ew_aggregate": int(meas.in_aggregate.astype(bool).sum())
            if "in_aggregate" in meas else None,
            "synth_aggregate": int(len(sy)),
            "quarantine_reasons": reasons,
            "note": "asymmetric on purpose -- these are the EW route's own physics "
                    "constraints, and they are why the 13 survivors are that route's "
                    "EASIEST lines. A selection effect inside the comparison."}

    d = doc["delta_route"]
    print(f"Delta_route, per-line paired, n={len(m)} of {a.pool_size}")
    print(f"   median {d['per_line_paired_median']:+.4f}   mean {d['mean']:+.4f}   "
          f"sd {d['sd']:.4f}   range [{d['min']:+.4f}, {d['max']:+.4f}]")
    print(f"   (difference of medians on the same lines: "
          f"{doc['other_ways_of_asking']['difference_of_medians_on_shared_lines']:+.4f})")
    if "attrition" in doc:
        at = doc["attrition"]
        print(f"\nattrition: {at['requested']} requested -> {at['measurable_by_ew']} "
              f"EW-measurable -> {at['ew_aggregate']} EW aggregate "
              f"(synth kept {at['synth_aggregate']})")
    print(f"\n⚠️ the median is small and the SPREAD is not: sd {d['sd']:.4f} over a "
          f"{d['max'] - d['min']:.3f} dex range. Reporting the median alone would call "
          f"the routes nearly equivalent.")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
