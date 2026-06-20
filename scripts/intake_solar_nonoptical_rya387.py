#!/usr/bin/env python3
"""
scripts/intake_solar_nonoptical_rya387.py
==========================================
RYA-387 — RE-INTAKE of the solar non-optical VALD wings, re-extracted at
central-depth **0.001** (synthesis-grade) to supersede the RYA-381 0.05 batch.

Reuses the RYA-381 intake scan verbatim (completeness-by-span, truncation,
coverage, HFS-on content signal, solar [M/H]) and adds the gate that is the
whole point of this ticket: the **RYA-389 extraction-threshold consistency
check** (`vald_parse.verify_extraction_threshold`). For the 0.05 batch the
threshold was a soft FLAG ("don't silently mix"); here the 0.001 depth IS the
acceptance criterion, so a delivery is ACCEPT only when the RYA-381 hard gates
pass AND the RYA-389 gate returns ACCEPT (effective depth == 0.001). A wing that
came back at 0.05 would FLAG here — a real catch that the re-extraction slipped,
not a nuisance (RYA-389 first-real-data test, RYA-387 comment).

Chunk plan (the delivered 0.001 re-extraction, finer in the IR than the 0.05
batch's 9500–14000 / 14000–18000 / 18000–25000):
    1150–2000 / 2000–3780 / [optical core 3780–6910 untouched] /
    6910–9500 / 9500–17000 / 17000–25000

Read-only audit (writes a CSV report only). ACCEPT / REJECT per delivery.

Usage:  python scripts/intake_solar_nonoptical_rya387.py
Out:    data/audit/vald_inventory/rya387_intake.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LL = ROOT / 'data' / 'linelists'
OUT = ROOT / 'data' / 'audit' / 'vald_inventory'
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(LL))

# Reuse the RYA-381 intake scan + verdict logic verbatim (the ticket: "reuse
# intake_solar_nonoptical_rya381.py").
from intake_solar_nonoptical_rya381 import (scan, hfs_verdict,          # noqa: E402
                                            accept_verdict, OPTICAL_REF)
from vald_parse import verify_extraction_threshold                       # noqa: E402

# The delivered 0.001 re-extraction (RyanSchmitt.019704–019708). (file, lo, hi).
DELIVERIES = [
    ('vald_solar_fuv_1150_2000_hfson_raw.txt',     1150.0,  2000.0),
    ('vald_solar_nearuv_2000_3780_hfson_raw.txt',  2000.0,  3780.0),
    ('vald_solar_redopt_6910_9500_hfson_raw.txt',  6910.0,  9500.0),
    ('vald_solar_ir_9500_17000_hfson_raw.txt',     9500.0, 17000.0),
    ('vald_solar_ir_17000_25000_hfson_raw.txt',   17000.0, 25000.0),
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, lo, hi in DELIVERIES:
        p = LL / name
        if not p.exists():
            rows.append({'asset': name, 'verdict': 'MISSING'})
            continue
        r = scan(p, lo, hi)
        r['hfs'] = hfs_verdict(r)
        hard = accept_verdict(r, r['hfs'])                 # RYA-381 hard gates
        thr_verdict, thr_msg, eff = verify_extraction_threshold(p)   # RYA-389 gate
        r['threshold_389'] = f"{thr_verdict} (eff={eff})"
        if hard.startswith('ACCEPT') and thr_verdict == 'ACCEPT':
            r['verdict'] = 'ACCEPT'
        elif not hard.startswith('ACCEPT'):
            r['verdict'] = hard
        else:
            r['verdict'] = f'REJECT — RYA-389 threshold {thr_verdict}: {thr_msg}'
        rows.append(r)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'rya387_intake.csv', index=False)

    pd.set_option('display.width', 240, 'display.max_columns', 40, 'display.max_colwidth', 52)
    print("\n=========== RYA-387 SOLAR NON-OPTICAL RE-INTAKE (0.001) ===========")
    show = ['asset', 'req_lo_A', 'req_hi_A', 'actual_lo_A', 'actual_hi_A',
            'n_data_lines', 'coverage_frac', 'truncated', 'applied_MH',
            'hfs', 'min_depth', 'threshold_389', 'verdict']
    print(df[show].to_string(index=False))

    n_accept = sum(1 for r in rows if str(r.get('verdict', '')) == 'ACCEPT')
    print(f"\n  {n_accept}/{len(DELIVERIES)} wings ACCEPT (RYA-381 hard gates + RYA-389 0.001 gate)")
    print(f"  [out] {OUT / 'rya387_intake.csv'}")
    if n_accept != len(DELIVERIES):
        sys.exit("INTAKE FAILED — not all wings ACCEPT; assembly blocked.")


if __name__ == '__main__':
    main()
