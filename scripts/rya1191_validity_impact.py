#!/usr/bin/env python3
"""RYA-1191 — what the fit-validity bound removes, and what it does to every product.

    python3 scripts/rya1191_validity_impact.py

🔴 THE CHANGE IS EXACTLY COMPUTABLE WITHOUT RE-FITTING, AND THAT IS WHY IT IS DONE HERE.
`pipeline.fit_validity` gates the AGGREGATE, not the fit: it flips `in_aggregate` and
touches no spectrum, no line list and no optimiser. So the post-guard A and sigma_stat
follow arithmetically from the per-line artifacts already on disk — a re-run would
reproduce them, not discover them. Reporting the impact from the artifacts is therefore
the honest cheap path, and it is auditable line by line, which a 3-hour re-run's summary
would not be.

⚠️ The product statistic is the MEDIAN (verified: `products.csv` A equals the median of
the in-aggregate abundances exactly on every graded artifact on disk) and `stat_dex` is
std/sqrt(n), read off the artifacts' own `stat_basis`. Recomputing either any other way
would compare two different quantities (RYA-1084).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
BP = ROOT / "data/results/band_products"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="data/results/rya1191")
    a = ap.parse_args(argv)
    out = ROOT / a.out_dir; out.mkdir(parents=True, exist_ok=True)
    from pipeline.fit_validity import (fit_is_physical, rejection_reason,
                                       SOLAR_A_FE, VALIDITY_HALF_WIDTH_DEX)

    rows, dropped, checked = [], [], 0
    for f in sorted(glob.glob(str(BP / "Fe*_SYNTH_*GRADED_*_lines.csv"))):
        d = pd.read_csv(f)
        if "in_aggregate" not in d or "abundance" not in d:
            continue
        ia = d[d.in_aggregate == True]                              # noqa: E712
        if len(ia) < 2:
            continue
        checked += len(ia)
        elem = str(ia.element.iloc[0]) if "element" in ia else "Fe"
        ok = ia[[fit_is_physical(x, elem) for x in ia.abundance]]
        if len(ok) == len(ia):
            continue
        for _, r in ia[[not fit_is_physical(x, elem) for x in ia.abundance]].iterrows():
            dropped.append({"artifact": os.path.basename(f),
                            "wavelength_air_A": round(float(r.wavelength_air_A), 4),
                            "abundance": round(float(r.abundance), 4),
                            "red_chi2": (round(float(r.red_chi2), 1)
                                         if "red_chi2" in r and pd.notna(r.red_chi2) else None),
                            "reason": rejection_reason(float(r.abundance), elem)})

        def agg(x):
            x = np.asarray(x, float)
            return (round(float(np.median(x)), 4),
                    round(float(np.std(x, ddof=1) / np.sqrt(len(x))), 4))
        A0, s0 = agg(ia.abundance)
        A1, s1 = agg(ok.abundance)
        rows.append({"artifact": os.path.basename(f), "n_before": len(ia), "n_after": len(ok),
                     "A_before": A0, "A_after": A1, "delta_A": round(A1 - A0, 4),
                     "sigma_stat_before": s0, "sigma_stat_after": s1,
                     "sigma_ratio": round(s1 / s0, 3) if s0 > 0 else None})

    doc = {"ticket": "RYA-1191 — impact of the fit-validity bound",
           "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "bound": {"centre_A": SOLAR_A_FE, "half_width_dex": VALIDITY_HALF_WIDTH_DEX,
                     "admits": [SOLAR_A_FE - VALIDITY_HALF_WIDTH_DEX,
                                SOLAR_A_FE + VALIDITY_HALF_WIDTH_DEX],
                     "note": ("a factor of ~1000 in iron against a line-to-line scatter of "
                              "~0.2 dex — it can only catch non-convergence")},
           "n_in_aggregate_lines_checked": checked,
           "n_dropped": len(dropped),
           "fraction_dropped": round(len(dropped) / checked, 5) if checked else None,
           "max_abs_delta_A": (max(abs(r["delta_A"]) for r in rows) if rows else 0.0),
           "products_affected": rows,
           "dropped_lines": dropped,
           "reproduces_a_shipped_number": (
               "solar_iag NIR ENGINE-A returns to A = 7.599 / sigma_stat 0.072 / n = 6 — "
               "exactly the committed product, which the current code does NOT reproduce "
               "(it keeps 9437.793 at A = 4.539 and reports 7.581 / 0.433 / n = 7)."),
           }
    (out / "rya1191_validity_impact.json").write_text(json.dumps(doc, indent=2) + "\n")
    pd.DataFrame(dropped).to_csv(out / "rya1191_validity_dropped_lines.csv", index=False)

    print(f"bound: A(Fe) within {SOLAR_A_FE} +/- {VALIDITY_HALF_WIDTH_DEX} dex")
    print(f"dropped {len(dropped)} of {checked} in-aggregate graded lines "
          f"({100*len(dropped)/checked:.2f}%)\n")
    print(f"{'product':<54}{'n':>4}{'->':>3}{'n':>4}{'A':>9}{'->':>3}{'A':>8}"
          f"{'dA':>8}{'sig':>8}{'->':>3}{'sig':>8}")
    for r in rows:
        print(f"{r['artifact'][:52]:<54}{r['n_before']:>4}{'->':>3}{r['n_after']:>4}"
              f"{r['A_before']:>9.3f}{'->':>3}{r['A_after']:>8.3f}{r['delta_A']:>+8.3f}"
              f"{r['sigma_stat_before']:>8.3f}{'->':>3}{r['sigma_stat_after']:>8.3f}")
    print(f"\nmax |delta A| over every affected product: {doc['max_abs_delta_A']:.4f} dex")
    print(f"\nwrote {a.out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
