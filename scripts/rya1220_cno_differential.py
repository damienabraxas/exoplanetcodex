#!/usr/bin/env python3
"""
RYA-1220 Work Package D -- the target/Sun CNO differential, done line by line.

WHAT THE PAIR AUDIT SAID. cno_target_sun_pair_audit.json records a 2.404 dex CN_red
difference between Procyon and the Sun and holds the whole differential pending
Jacobians and cross-star covariance.

WHAT IT ACTUALLY IS. The solar CN_red fit is broken, not the star. Its sigma_fit is
0.835 and its frac_rise_weaker is 4.8e-05 -- the flattest chi2 surface of any route in
this campaign -- and it returns A(N) = 7.384, which is 0.45 dex BELOW the solar
reference. Procyon returns 9.788. Two ends of one unusable diagnostic, the same CN_red
the ratified policy already rejected as an abundance route. No Jacobian promotes that,
and none is needed to see it.

[O I] 6300 is weaker evidence for the same reason: sigma_fit 0.421 on the SOLAR side.

WHAT IS ACTUALLY DELIVERABLE. Carbon. Four independent matched indicators -- CH G-band,
C I 5052, C I 5380 and C2 Swan -- all on the same instrument, same windows, same
synthesis setup, with solar sigma_fit 0.005-0.098 and Procyon 0.002-0.008.

METHOD. MEDIAN OF DIFFERENCES, never difference of medians: the pairing is what makes
shared terms cancel, and aggregating first throws that away (pipeline/paired_differential
states the same rule for the per-line Fe case). Both are reported so the gap is visible.

Matching is on the indicator KEY, an exact identity, rather than a wavelength tolerance
-- these are bands, not single lines, so a tolerance would be the wrong instrument.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import statistics
import sys

import pandas as pd

SOLAR = "data/audit/cno_synthesis/solar_vis_cno_per_band.csv"
TARGET = "data/audit/cno_synthesis/procyon_vis_cno_per_band.csv"

#: A diagnostic the ratified policy rejects cannot carry a differential either.
REJECTED = {"CN_red": "rejected as an abundance route (cno_method_policy.json)"}
#: Solar-side fit quality below which a pair is evidence about the FIT, not the star.
SOLAR_SIGMA_FLAG = 0.2


def build(solar: pd.DataFrame, target: pd.DataFrame) -> dict:
    s = solar.set_index("key")
    t = target.set_index("key")
    pairs = []
    for key in s.index:
        if key not in t.index:
            continue
        srow, trow = s.loc[key], t.loc[key]
        delta = float(trow["A_X"]) - float(srow["A_X"])
        blockers = []
        if key in REJECTED:
            blockers.append(REJECTED[key])
        if float(srow["sigma_fit"]) > SOLAR_SIGMA_FLAG:
            blockers.append(
                f"solar sigma_fit {float(srow['sigma_fit']):.3f} > {SOLAR_SIGMA_FLAG}: the "
                f"SOLAR end is unconstrained, so this pair measures the fit, not the star")
        fr = srow.get("frac_rise_weaker")
        pairs.append({
            "key": key, "element": str(srow["element"]), "role": str(srow["role"]),
            "solar_A": float(srow["A_X"]), "solar_sigma_fit": float(srow["sigma_fit"]),
            "solar_frac_rise_weaker": float(fr) if pd.notna(fr) else None,
            "target_A": float(trow["A_X"]), "target_sigma_fit": float(trow["sigma_fit"]),
            "delta_dex": round(delta, 4),
            "usable": not blockers, "blockers": blockers,
        })

    out = {"schema": "rya1220.cno_target_sun_differential.v1", "ticket": "RYA-1220",
           "target": "Procyon", "instrument": "HARPS", "band": "VIS",
           "method": ("median of per-indicator DIFFERENCES, matched on indicator key. "
                      "Difference of medians is reported alongside only to show the gap; "
                      "it is not the differential."),
           "pairs": pairs, "by_element": {}}

    for element in sorted({p["element"] for p in pairs}):
        usable = [p for p in pairs if p["element"] == element and p["usable"]]
        blocked = [p for p in pairs if p["element"] == element and not p["usable"]]
        if not usable:
            out["by_element"][element] = {
                "state": "BLOCKED", "n_usable": 0,
                "why": [b for p in blocked for b in p["blockers"]],
                "blocked_indicators": [p["key"] for p in blocked]}
            continue
        d = [p["delta_dex"] for p in usable]
        naive = (statistics.median([p["target_A"] for p in usable])
                 - statistics.median([p["solar_A"] for p in usable]))
        out["by_element"][element] = {
            "state": "MEASURED", "n_usable": len(d),
            "indicators": [p["key"] for p in usable],
            "median_of_differences_dex": round(statistics.median(d), 4),
            "mean_of_differences_dex": round(statistics.mean(d), 4),
            "sd_dex": round(statistics.stdev(d), 4) if len(d) > 1 else None,
            "min_dex": round(min(d), 4), "max_dex": round(max(d), 4),
            "spread_dex": round(max(d) - min(d), 4),
            "difference_of_medians_dex": round(naive, 4),
            "blocked_indicators": [p["key"] for p in blocked],
            "note": ("indicator-to-indicator spread is the honest uncertainty here; the "
                     "parameter Jacobians that would let it be decomposed are still owed "
                     "on both the target and the solar side."),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solar", default=SOLAR)
    ap.add_argument("--target", default=TARGET)
    ap.add_argument("--out", default="data/output/rya1220/cno_target_sun_differential.json")
    args = ap.parse_args()
    for p in (args.solar, args.target):
        if not pathlib.Path(p).exists():
            print(f"absent: {p}", file=sys.stderr)
            return 2
    doc = build(pd.read_csv(args.solar), pd.read_csv(args.target))
    pathlib.Path(args.out).write_text(json.dumps(doc, indent=2) + "\n")

    print("RYA-1220 CNO target/Sun differential -- Procyon, HARPS VIS\n")
    print(f"  {'key':10s} {'el':2s} {'Sun':>7s} {'sig':>6s} {'Target':>7s} {'sig':>6s} "
          f"{'delta':>8s}  blockers")
    for p in doc["pairs"]:
        print(f"  {p['key']:10s} {p['element']:2s} {p['solar_A']:7.3f} "
              f"{p['solar_sigma_fit']:6.3f} {p['target_A']:7.3f} {p['target_sigma_fit']:6.3f} "
              f"{p['delta_dex']:+8.3f}  {'; '.join(p['blockers']) if p['blockers'] else '-'}")
    print()
    for element, g in doc["by_element"].items():
        if g["state"] != "MEASURED":
            print(f"  [{element}/H]  BLOCKED -- {g['blocked_indicators']}")
            continue
        print(f"  [{element}/H](Procyon - Sun) = {g['median_of_differences_dex']:+.3f} dex "
              f"from {g['n_usable']} matched indicators {g['indicators']}")
        print(f"      spread {g['spread_dex']:.3f}, difference-of-medians would be "
              f"{g['difference_of_medians_dex']:+.3f}")
    print(f"\n  -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
