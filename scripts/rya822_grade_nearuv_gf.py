#!/usr/bin/env python3
"""
RYA-822 — attach a gf grade-or-bound to every near-UV Fe I line, through RYA-799's scale.

Single-sourced through `pipeline.gf_grades`: the `GF-`/`systematic:` vocabulary is a claim
about the ATOMIC DATUM and is deliberately disjoint from `MQ-` (measurement quality) and
from NIST's bare letters. This ticket adds no vocabulary of its own beyond `GF-NIST`, which
lives in that module with the rest.

WHAT THIS IS EXPECTED TO SHOW, stated before it is run so the result cannot be read as a
discovery after the fact. RYA-799 proved the Fe pool is uniformly `gf:K14` in the optical
AND the IR, and predicted the near-UV would be Kurucz-floored too. Recovering the per-line
source from VALD's fourth physical line confirms it for this band:

    Fe I 3000-3780 A : 84.3% gf:K14, 15.7% no gf tag at all, ZERO primary lab sources.

So the honest deliverable is a bounded count of anything with real accuracy, not a recovery
of grades that were never there.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import gf_grades as gg                                   # noqa: E402

ADJ = ROOT / 'data' / 'audit' / 'nearuv_gf_adjudication'
OUT = ROOT / 'data' / 'results' / 'rya822'


def main() -> None:
    src = ADJ / 'nearuv_fe1_gf_adjudication.csv'
    d = pd.read_csv(src)
    d = d[(d.key_z == 26) & (d.ion == 1)].reset_index(drop=True)

    # The value the synthesis ACTUALLY used: canonical resolves at load, so an adjudicated
    # line was measured on the NIST value and an unadjudicated one on VALD's. Grading the
    # wrong one is exactly the fabricated pedigree 799's SCALE-MISMATCH exists to catch.
    d['log_gf_used'] = np.where(d.adjudicated, d.nist_log_gf, d.gf_linelist_vald)

    graded = gg.grade_pool(d, wave_col='wavelength_air_A',
                           ep_col='excitation_potential_eV', loggf_col='log_gf_used')
    assert len(graded) == len(d), "grade_pool changed the row count"
    assert graded.gf_grade.notna().all(), "a line came back with a blank grade"
    assert (graded.gf_sigma_dex > 0).all(), (
        "a non-positive gf sigma — the RYA-799 positivity guard: a mis-parsed table once "
        "put log gf into the uncertainty column and produced -2.62 dex")

    counts = graded.gf_grade.value_counts()
    print(f"near-UV Fe I lines graded : {len(graded)}")
    for k, v in counts.items():
        print(f"   {k:<34} {v:5d}  ({100 * v / len(graded):5.1f}%)")

    cited = graded[graded.gf_grade.isin([gg.GRADE_LAB, gg.GRADE_NIST])]
    print(f"\n  with a CITABLE per-line sigma : {len(cited)}")
    if len(cited):
        print(f"    sigma  min {cited.gf_sigma_dex.min():.4f}  median "
              f"{cited.gf_sigma_dex.median():.4f}  max {cited.gf_sigma_dex.max():.4f} dex")
        print(f"    vs the blanket K07 systematic {gg.K07_SYSTEMATIC_DEX:.2f} dex — "
              f"a factor {gg.K07_SYSTEMATIC_DEX / cited.gf_sigma_dex.median():.1f} tighter")
    mism = graded[graded.gf_grade == gg.GRADE_MISMATCH]
    print(f"  SCALE-MISMATCH                : {len(mism)}")
    if len(mism):
        print(f"    |delta| median {mism.gf_delta_dex.abs().median():.3f}  "
              f"max {mism.gf_delta_dex.abs().max():.3f} dex")

    OUT.mkdir(parents=True, exist_ok=True)
    graded.to_csv(OUT / 'rya822_nearuv_gf_grades_per_line.csv', index=False)
    summary = {
        'ticket': 'RYA-822', 'band_A': [3000.0, 3780.0], 'species': 'Fe I',
        'n_lines': int(len(graded)),
        'grade_counts': {str(k): int(v) for k, v in counts.items()},
        'n_with_cited_sigma': int(len(cited)),
        'sigma_dex': ({'min': float(cited.gf_sigma_dex.min()),
                       'median': float(cited.gf_sigma_dex.median()),
                       'max': float(cited.gf_sigma_dex.max())} if len(cited) else None),
        'blanket_systematic_dex': gg.K07_SYSTEMATIC_DEX,
        'vald_gf_source_recovered': {'gf:K14': 3680, 'no_gf_tag': 684,
                                     'primary_lab': 0,
                                     'note': "recovered from VALD's 4th physical line; "
                                             "the built list had flattened it to 'VALD3'"},
        'caveat': ('GF-NIST is a COMPILATION grade, not a primary laboratory measurement '
                   '(RYA-760: FMW *is* NIST and VALD copies it). It carries a citable '
                   'per-line sigma but is deliberately outside `is_graded`.'),
    }
    (OUT / 'rya822_nearuv_gf_grade_summary.json').write_text(json.dumps(summary, indent=2))
    print(f"\n[out] {OUT}/rya822_nearuv_gf_grade_summary.json")


if __name__ == '__main__':
    main()
