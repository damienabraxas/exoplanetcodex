#!/usr/bin/env python3
"""
scripts/add_namg_recovery_to_pool_rya465.py
===========================================
RYA-465 — add the measured UNSATURATED recovery lines (from
measure_unsaturated_namg_rya465.py) to the canonical solar EW pool
(data/measured/sol_ew_results_v1.csv). Only lines that land in the LINEAR COG
regime are added (the whole point: give the element lines it can measure on).
Saturated candidates are NOT added (they would only be SAT-culled) — they are the
documented finding, kept in the audit table.

Idempotent: a wavelength already in the pool is skipped. Existing rows (incl. the
Fe anchor) are byte-identical — this only appends. No EW is tuned.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
POOL = ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv'
AUDIT = ROOT / 'data' / 'audit' / 'line_recovery' / 'namg_unsaturated_rya465.csv'
SCHEMA = ['element', 'ion', 'wavelength_air_A', 'ew_mA', 'ew_err_mA', 'profile_type',
          'chi2', 'blend_flag', 'notes']


def main():
    audit = pd.read_csv(AUDIT)
    add = audit[(audit['status'] == 'ok') & (audit['unsaturated'])].copy()

    pool = pd.read_csv(POOL)
    existing = set(zip(pool['element'], pool['ion'], pool['wavelength_air_A'].round(3)))

    raw = POOL.read_text()
    if not raw.endswith('\n'):
        raw += '\n'
    new_lines, added = [], []
    for _, r in add.iterrows():
        key = (r['element'], r['ion'], round(float(r['wavelength_air_A']), 3))
        if key in existing:
            continue
        note = (f"RYA-465 unsaturated recovery; REW {r['rew']}; NIST grade {r['nist_grade']}; "
                f"loggf {r['log_gf']} ({str(r['loggf_reference']).split(';')[0]})")
        row = [r['element'], r['ion'], f"{float(r['wavelength_air_A']):.3f}",
               f"{float(r['ew_mA']):.2f}", f"{float(r['ew_err_mA']):.2f}",
               r['profile'], f"{float(r['chi2']):.4f}", 'False', note]
        assert not any(',' in str(x) for x in row), "field has a comma — would break CSV"
        new_lines.append(','.join(str(x) for x in row))
        added.append(key)

    if new_lines:
        POOL.write_text(raw + '\n'.join(new_lines) + '\n')
    print(f"RYA-465 pool append: {len(added)} unsaturated line(s) added "
          f"({len(add) - len(added)} already present)")
    for k in added:
        print(f"  + {k[0]} {k[1]} {k[2]}")
    # report the saturated candidates NOT added (the finding)
    sat = audit[(audit['status'] == 'ok') & (~audit['unsaturated'])]
    print(f"  {len(sat)} saturated candidate(s) NOT added (finding): "
          f"{', '.join(f'{r.element} {r.wavelength_air_A:.1f}' for r in sat.itertuples())}")


if __name__ == '__main__':
    main()
