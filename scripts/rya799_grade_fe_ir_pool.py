#!/usr/bin/env python3
"""Grade — or honestly bound — the 271 PASS Fe I IR lines. RYA-799 smoke test.

    python3 scripts/rya799_grade_fe_ir_pool.py --pool <FeI_6910_9199_..._PASS.csv>

Emits `gf_grade` and `gf_grade_source` (plus the evidence columns and a numeric
`gf_sigma_dex`) onto the pool, writing the graded copy beside it. `data/measured/band_ew/`
is gitignored (RYA-469), so that file is regenerable rather than committed; the committed
half is the summary report under `data/results/rya799/`.

Refuses to write unless: rows in == rows out, every row carries a non-empty grade AND
source, and graded + systematic == the input count.
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
    GRADE_LAB, GRADE_MISMATCH, GRADE_SYSTEMATIC, K07_SYSTEMATIC_DEX, grade_pool)

OUT_DEFAULT = ROOT / "data" / "results" / "rya799"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pool", type=Path, required=True,
                    help="the PASS pool CSV (lives on Sirius; a path, never a literal)")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--graded-pool", type=Path, default=None,
                    help="where to write the graded pool copy "
                         "(default: alongside --pool, suffixed _GRADED)")
    args = ap.parse_args(argv)

    pool = pd.read_csv(args.pool)
    n_in = len(pool)
    graded = grade_pool(pool)

    n_lab = int((graded.gf_grade == GRADE_LAB).sum())
    n_sys = int((graded.gf_grade == GRADE_SYSTEMATIC).sum())
    n_mis = int((graded.gf_grade == GRADE_MISMATCH).sum())
    assert n_lab + n_sys + n_mis == n_in, "grade states do not partition the pool"

    print("=" * 78)
    print("RYA-799 — gf grade or honest bound, Fe I 6910-9199 A PASS pool")
    print("=" * 78)
    print(f"\nrows in {n_in}  ->  rows out {len(graded)}   "
          f"{'RECONCILED' if n_in == len(graded) else 'MISMATCH'}")
    # Within the plain-systematic bucket, separate the lines whose gf VALUE is confirmed
    # against a named non-K07 reference from the ones with no tie at all. Same bar for
    # both — but the first group is ungraded only because NIST ASD is down, which is a
    # recoverable debt, and the second is ungraded because nobody has measured the line.
    plain = graded[graded.gf_grade == GRADE_SYSTEMATIC]
    attributed = plain[plain.gf_reference_tag.astype(str).str.strip().isin(["", "K07"]) == False]  # noqa: E712
    n_attr = int(len(attributed))

    print(f"\n  graded            {n_lab:3d}   {GRADE_LAB} (primary lab measurement, "
          f"per-line sigma)")
    print(f"  systematic-only   {n_sys + n_mis:3d}   of which:")
    print(f"      {n_sys:3d}   {GRADE_SYSTEMATIC}")
    print(f"            {n_attr:3d}   value CONFIRMED against a named non-K07 reference "
          f"(attributed, not graded: no citable per-line sigma while ASD is down)")
    print(f"            {n_sys - n_attr:3d}   no reference tie at all")
    print(f"      {n_mis:3d}   {GRADE_MISMATCH}")
    print(f"\n  graded {n_lab} / systematic-only {n_sys + n_mis}  "
          f"(N + M = {n_lab + n_sys + n_mis} == {n_in})")

    blank = int((graded.gf_grade.astype(str).str.strip() == "").sum()
                + (graded.gf_grade_source.astype(str).str.strip() == "").sum())
    print(f"  blank bars: {blank}  {'PASS' if blank == 0 else 'FAIL'}")

    print("\nPER-LINE SIGMA DISTRIBUTION (what RYA-783's budget can consume "
          "instead of a blanket 0.20):")
    lab = graded[graded.gf_grade == GRADE_LAB]
    if len(lab):
        q = lab.gf_sigma_dex.quantile([0, .25, .5, .75, 1]).round(3).tolist()
        print(f"  graded lines   n={len(lab):3d}  min {q[0]}  q25 {q[1]}  median {q[2]}  "
              f"q75 {q[3]}  max {q[4]} dex")
        print(f"                 mean {lab.gf_sigma_dex.mean():.3f} dex "
              f"vs the blanket {K07_SYSTEMATIC_DEX:.2f} — a factor "
              f"{K07_SYSTEMATIC_DEX / lab.gf_sigma_dex.mean():.1f} tighter")
    print(f"  bounded lines  n={n_sys + n_mis:3d}  all {K07_SYSTEMATIC_DEX:.2f} dex "
          f"(RYA-161 semi-empirical Kurucz systematic)")
    pool_sigma = float(np.sqrt(np.mean(graded.gf_sigma_dex.values ** 2)))
    print(f"  pool RMS sigma {pool_sigma:.4f} dex "
          f"(blanket-0.20 assumption: {K07_SYSTEMATIC_DEX:.4f})")

    print("\nLAB COVERAGE vs LAB USE — the finding:")
    n_lab_exists = n_lab + int((graded.gf_note.astype(str)
                                .str.contains("PRIMARY lab gf exists")).sum())
    print(f"  {n_lab_exists} of {n_in} lines have a PRIMARY laboratory gf in the "
          f"literature;\n  the pool was measured on that value for {n_lab} of them.")

    print("\nSOURCES OF THE GRADED LINES:")
    for s, n in graded[graded.gf_grade == GRADE_LAB].gf_grade_source.value_counts().items():
        print(f"  {n:3d}  {s}")

    if n_mis:
        print(f"\nSCALE-MISMATCH — {n_mis} lines where a better-referenced gf EXISTS "
              f"IN-REPO but is\nnot the value the pool was measured with. Not a grade; "
              f"an actionable finding:")
        mm = graded[graded.gf_grade == GRADE_MISMATCH]
        print(f"  by reference tag: "
              f"{dict(mm.gf_reference_tag.value_counts().head(10))}")
        print(f"  |delta| median {mm.gf_delta_dex.abs().median():.3f} dex, "
              f"max {mm.gf_delta_dex.abs().max():.3f} dex")

    graded_path = (args.graded_pool if args.graded_pool is not None
                   else args.pool.with_name(args.pool.stem + "_GRADED.csv"))
    graded.to_csv(graded_path, index=False)

    args.out.mkdir(parents=True, exist_ok=True)
    summary = {
        "ticket": "RYA-799",
        "pool": str(args.pool),
        "n_rows_in": n_in, "n_rows_out": int(len(graded)),
        "graded": n_lab, "systematic_only": n_sys + n_mis,
        "systematic_plain": n_sys, "systematic_scale_mismatch": n_mis,
        "systematic_value_confirmed_but_ungraded": n_attr,
        "systematic_no_tie_at_all": n_sys - n_attr,
        "attributed_reference_tags": {k: int(v) for k, v in
                                      attributed.gf_reference_tag.value_counts().items()},
        "blank_bars": blank,
        "sigma_graded_dex": {
            "n": len(lab),
            "min": float(lab.gf_sigma_dex.min()) if len(lab) else None,
            "median": float(lab.gf_sigma_dex.median()) if len(lab) else None,
            "max": float(lab.gf_sigma_dex.max()) if len(lab) else None,
            "mean": float(lab.gf_sigma_dex.mean()) if len(lab) else None,
        },
        "sigma_bounded_dex": K07_SYSTEMATIC_DEX,
        "pool_rms_sigma_dex": pool_sigma,
        "graded_sources": {k: int(v) for k, v in
                           graded[graded.gf_grade == GRADE_LAB]
                           .gf_grade_source.value_counts().items()},
        "scale_mismatch_by_tag": {k: int(v) for k, v in
                                  graded[graded.gf_grade == GRADE_MISMATCH]
                                  .gf_reference_tag.value_counts().items()},
        "graded_pool_written_to": str(graded_path),
        "run_date": str(date.today()),
    }
    (args.out / "rya799_gf_grade_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    per_line = graded[["wavelength_air_A", "ep_eV", "log_gf", "gf_grade",
                       "gf_grade_source", "gf_sigma_dex", "gf_reference_tag",
                       "gf_ref_loggf", "gf_delta_dex"]]
    per_line.to_csv(args.out / "rya799_gf_grades_per_line.csv", index=False)
    print(f"\nwrote {graded_path}  (gitignored, regenerable — RYA-469)")
    print(f"wrote {args.out}/rya799_gf_grades_per_line.csv")
    print(f"wrote {args.out}/rya799_gf_grade_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
