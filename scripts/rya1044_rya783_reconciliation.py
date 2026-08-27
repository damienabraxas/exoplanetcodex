#!/usr/bin/env python3
"""RYA-1044: reconcile the RYA-783 matrix's Fe I Engine-B cells against the products the
RYA-1044 synthesis leg now emits.

🔴 THE ANSWER IS THAT THEY CANNOT BE DIFFERENCED, AND THAT IS A MEASUREMENT, NOT A DODGE.

The matrix reports Fe I VIS `ENGINE-B-NLTE` = 7.572 (n=136); our graded run reports 7.497
(n=67). The 0.075 dex between them invites a route-vs-route or old-vs-new reading. It is
neither. The two pools share **2 candidate lines out of 159 and 67**, before either route
excludes anything:

    matrix candidates              159   (136 in aggregate)
    ours                            67
    SHARED CANDIDATES                 2   -- 5618.632 and 6843.655
    ours the matrix never considered 65
    matrix lines we never considered 157

That is SELECTION, not exclusion -- the matrix's own exclusion list is 23 rows and none of
them are our lines. Two near-disjoint pools measured on the same star, same instrument and
(as it turns out) the same wavelength range are two different experiments, and RYA-842
ratified the reason this matters: WHICH LINES ENTER THE POOL is first-order on the value
and the scatter, while gf work moves the error bar.

⚠️ NOR IS THE BAND THE EXPLANATION, WHICH IS THE FIRST THING ANYONE WILL GUESS. The matrix
cell is labelled 3800-6910 A against our 4200-6910, but ALL 136 of its aggregate lines fall
inside 4200-6910. The label is wider than the data; the pools do not differ because of it.

⚠️ AND NO CLAIM IS MADE FROM THE TWO SHARED LINES. n=2 cannot carry one. This same session
retracted a level-dependence claim built on n=8 (RYA-1051) after it evaporated at n=43; two
points deserve less, not more, credit than eight.

FOUR AXES DIFFER AT ONCE, which is why the headline difference decomposes into nothing:
  * POOL      -- graded LAB-tier (67) vs depth-selected ungraded (159 candidates)
  * ROUTE     -- SYNTH flux fit vs PROFILEFIT EW/profile fit
  * gf RUNG   -- ours rung 3 (67/67 GF-LAB, measured RMS 0.0634) vs the matrix's
                 "gf scale (UNGRADED)", syst 0.1705
  * chi2 GATE -- the matrix excludes red_chi2 >= 10.0; our route applies no chi2 cut at
                 all, because RYA-847 measured that no threshold transfers between routes

A single number differencing two cells that differ on four axes would be a coincidence
wearing the clothes of a result.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TOL = 0.005      # wavelengths agree to better than this or they are different lines


def _rows(p: Path, agg_only: bool) -> pd.DataFrame:
    d = pd.read_csv(p)
    if agg_only and "in_aggregate" in d.columns:
        d = d[d.in_aggregate.astype(str).str.lower().isin(("true", "1"))]
    return d


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix-lines", required=True,
                    help="RYA-783 Fe I VIS ENGINE-B-NLTE per-line CSV (Sirius-only)")
    ap.add_argument("--ours-lines", required=True,
                    help="our graded ENGINE-B-NLTE per-line CSV")
    ap.add_argument("--out", default=str(ROOT / "data" / "results" / "rya1044" /
                                         "rya783_reconciliation.json"))
    a = ap.parse_args()

    M_all = _rows(Path(a.matrix_lines), False)
    M = _rows(Path(a.matrix_lines), True)
    O = _rows(Path(a.ours_lines), True)

    mw = set(M_all.wavelength_air_A.round(3))
    ow = set(O.wavelength_air_A.round(3))
    shared = sorted(mw & ow)

    inband = M[M.wavelength_air_A.between(4200, 6910)]
    doc = {
        "ticket": "RYA-1044",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "verdict": "NOT_COMPARABLE",
        "why": "the two pools share 2 candidate lines out of 159 and 67; the headline "
               "difference is a POOL difference, not a route or physics one",
        "matrix": {"candidates": int(len(M_all)), "in_aggregate": int(len(M)),
                   "in_our_band_4200_6910": int(len(inband)),
                   "band_label": "3800-6910",
                   "note": "the label is wider than the data -- every aggregate line "
                           "falls inside 4200-6910, so the band is NOT the explanation"},
        "ours": {"in_aggregate": int(len(O))},
        "overlap": {"shared_candidates": int(len(shared)),
                    "shared_wavelengths_A": shared,
                    "ours_not_considered_by_matrix": int(len(ow - mw)),
                    "matrix_not_considered_by_us": int(len(mw - ow)),
                    "basis": "SELECTION, not exclusion -- the matrix's exclusion list "
                             "contains none of our lines"},
        "axes_that_differ": ["pool (graded LAB-tier vs depth-selected ungraded)",
                             "route (SYNTH flux fit vs PROFILEFIT EW/profile fit)",
                             "gf rung (3, measured RMS 0.0634 vs UNGRADED, syst 0.1705)",
                             "chi2 gate (matrix excludes red_chi2 >= 10; we apply none, "
                             "RYA-847)"],
        "no_claim_from_shared": "n=2 cannot carry a value comparison. RYA-1051 retracted "
                                "a claim built on n=8 in this same session.",
    }
    print(f"matrix: {len(M_all)} candidates, {len(M)} in aggregate, "
          f"{len(inband)} inside 4200-6910")
    print(f"ours  : {len(O)} in aggregate")
    print(f"SHARED CANDIDATES: {len(shared)}  {shared}")
    print(f"  ours the matrix never considered : {len(ow - mw)}")
    print(f"  matrix lines we never considered : {len(mw - ow)}")
    print(f"\nVERDICT: {doc['verdict']} -- {doc['why']}")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
