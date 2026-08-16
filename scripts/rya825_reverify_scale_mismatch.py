#!/usr/bin/env python3
"""Re-verify RYA-799's SCALE-MISMATCH set against the corrected gf column — RYA-825.

    python3 scripts/rya825_reverify_scale_mismatch.py --pool <FeI_..._PASS.csv>

RYA-799 graded the 271 Fe I IR lines by comparing each line's reference gf against the
pool's `log_gf` column, and classified 48 as SCALE-MISMATCH — "a better-referenced gf
exists in-repo but is not the value the pool was measured with".

RYA-824 then re-inverted and showed that column was **stale intake metadata**: the
production path resolves gf out of `canonical_gf.csv` before inverting, so the column
described an input nobody used. RYA-825 corrected the column.

So the SCALE-MISMATCH count has to be re-derived, and the question this script answers is
the one the ticket asks: **how many of the 48 are genuinely on a different gf than
canonical provides, and how many were pure metadata artifacts?**

It re-grades the same pool twice — once on the stale column, once on the corrected one —
through the SAME `pipeline.gf_grades` code. Nothing about the grading logic changes; only
the gf it is told the pool used. Any difference between the two runs is attributable to
the column and to nothing else.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_grades import (  # noqa: E402
    GRADE_LAB, GRADE_MISMATCH, GRADE_SYSTEMATIC, grade_pool)

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"
OUT_DEFAULT = ROOT / "data" / "results" / "rya825"
WAVE_TOL_A, EP_TOL_EV = 0.01, 0.01


def corrected_gf(pool: pd.DataFrame) -> pd.DataFrame:
    """Attach the corrected (used) gf from the accounting table, per line."""
    acc = pd.read_csv(ACCOUNTING)
    acc = acc[acc["element"] == "Fe"]
    out = pool.copy()
    used, intake, canon = [], [], []
    for r in pool.itertuples():
        m = acc[(acc["ion"] == r.ion)
                & (np.abs(acc["wave_air_A"] - r.wavelength_air_A) <= WAVE_TOL_A)
                & (np.abs(acc["ep_eV"] - r.ep_eV) <= EP_TOL_EV)]
        if m.empty:
            used.append(np.nan); intake.append(np.nan); canon.append(False)
            continue
        used.append(float(m.log_gf.iloc[0]))
        intake.append(float(m.log_gf_intake.iloc[0]))
        canon.append(bool(m.gf_canonical.iloc[0]))
    out["log_gf_used"] = used
    out["log_gf_intake"] = intake
    out["gf_canonical"] = canon
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pool", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args(argv)

    pool = pd.read_csv(args.pool)
    aug = corrected_gf(pool)
    n_tied = int(np.isfinite(aug.log_gf_used).sum())
    n_moved = int((np.abs(aug.log_gf_used - aug.log_gf_intake) > 1e-6).sum())

    print("=" * 78)
    print("RYA-825 — RYA-799's SCALE-MISMATCH, re-derived on the CORRECTED gf column")
    print("=" * 78)
    print(f"\npool {len(pool)} lines | tied to the accounting table {n_tied} | "
          f"corrected gf differs from intake on {n_moved}")
    print(f"resolver-sourced (gf_canonical): {int(aug.gf_canonical.sum())}")

    # BEFORE — exactly what RYA-799 did: grade against the stale column.
    before = grade_pool(pool)
    # AFTER — the same grading code, told the gf the inversion actually used.
    after_in = aug[np.isfinite(aug.log_gf_used)].copy()
    after_in["log_gf"] = after_in["log_gf_used"]
    after = grade_pool(after_in)

    def tally(df, label):
        return {"label": label, "n": len(df),
                GRADE_LAB: int((df.gf_grade == GRADE_LAB).sum()),
                GRADE_SYSTEMATIC: int((df.gf_grade == GRADE_SYSTEMATIC).sum()),
                GRADE_MISMATCH: int((df.gf_grade == GRADE_MISMATCH).sum())}

    tb, ta = tally(before, "stale column (RYA-799 as published)"), tally(after, "corrected column")
    print(f"\n{'verdict':<34}{'BEFORE':>10}{'AFTER':>10}")
    for k in (GRADE_LAB, GRADE_SYSTEMATIC, GRADE_MISMATCH):
        print(f"  {k:<32}{tb[k]:10d}{ta[k]:10d}")
    print(f"  {'(rows graded)':<32}{tb['n']:10d}{ta['n']:10d}")

    # THE ANSWER: of RYA-799's mismatches, which survive on the corrected column?
    b_mis = set(np.round(before.loc[before.gf_grade == GRADE_MISMATCH,
                                    "wavelength_air_A"], 4))
    a_mis = set(np.round(after.loc[after.gf_grade == GRADE_MISMATCH,
                                   "wavelength_air_A"], 4))
    genuine = sorted(b_mis & a_mis)
    artifact = sorted(b_mis - a_mis)
    print(f"\nOF RYA-799's {len(b_mis)} SCALE-MISMATCH LINES:")
    print(f"  {len(genuine):3d} GENUINE — still on a gf that differs from what the "
          f"reference provides")
    print(f"  {len(artifact):3d} PURE METADATA ARTIFACT — the mismatch was the stale "
          f"column, not the science")
    if a_mis - b_mis:
        print(f"  {len(a_mis - b_mis):3d} NEWLY mismatched (visible only once the column "
              f"was honest)")

    print("\nNO ABUNDANCE MOVED — this ticket changes a reported column, not a value.")
    print("  The accounting `log_gf` is written and carried; no consumer computes with it")
    print("  (grep over scripts/ and pipeline/ finds only the writer). The one consumer "
          "that\n  does read it is pipeline/gf_grades.py, whose verdicts are exactly what "
          "this\n  script re-derives above.")

    args.out.mkdir(parents=True, exist_ok=True)
    aug.to_csv(args.out / "rya825_pool_gf_corrected.csv", index=False)
    (args.out / "rya825_reverification.json").write_text(json.dumps({
        "ticket": "RYA-825",
        "pool": str(args.pool), "n_pool": len(pool),
        "n_tied_to_accounting": n_tied,
        "n_gf_corrected": n_moved,
        "before": tb, "after": ta,
        "scale_mismatch_before": len(b_mis),
        "scale_mismatch_genuine": len(genuine),
        "scale_mismatch_metadata_artifact": len(artifact),
        "scale_mismatch_new": len(a_mis - b_mis),
        "genuine_wavelengths": [float(w) for w in genuine],
        "abundance_changed": False,
        "run_date": str(date.today()),
    }, indent=2) + "\n")
    print(f"\nwrote {args.out}/rya825_pool_gf_corrected.csv")
    print(f"wrote {args.out}/rya825_reverification.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
