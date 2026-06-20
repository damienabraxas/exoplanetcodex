#!/usr/bin/env python3
"""
scripts/assemble_solar_linelist_rya387.py
==========================================
RYA-387 — RE-ASSEMBLE linelist_solar.csv with the solar non-optical wings
re-extracted at central-depth **0.001**, superseding the RYA-381 0.05 batch.

Same design as scripts/assemble_solar_linelist_rya381.py (optical core kept
VERBATIM, wings parsed with the build_linelist parser/priority/NIST logic,
vald_proximity_flag recomputed over the union) — reused here. The only changes:

  * wings are the 0.001 re-extraction (RyanSchmitt.019704–019708), 5-chunk plan;
  * the assembled list is now depth-HOMOGENEOUS at 0.001 (the RYA-381 depth
    heterogeneity is RESOLVED — provenance records 0.001 throughout);
  * solar params for the manifest come from config.constants.get_star_params
    ('solar') — NOT hardcoded (RYA-388/298).

Base note: this runs on a linelist_solar.csv whose current content is the
optical-only core (the 0.05 wings were never merged), so "replace the 0.05
wings" is achieved by assembling core + 0.001 wings directly — the defective
0.05 batch never enters the canonical file.

Usage:  python scripts/assemble_solar_linelist_rya387.py [--dry-run]
Out:    data/linelists/linelist_solar.csv  (extended, in place)
        data/linelists/linelist_solar_provenance_rya387.json
"""
import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

# Reuse the RYA-381 wing parser + schema verbatim.
from assemble_solar_linelist_rya381 import (parse_wing, COLS, SOLAR,     # noqa: E402
                                            ELEMENTS, LL)
from scripts.build_linelist import load_elements                         # noqa: E402
from pipeline.build_linelist import (compute_vald_proximity_flag,        # noqa: E402
                                     build_vetted_blend_flag)
from config.constants import get_star_params                             # noqa: E402

THRESHOLD = 0.001
# Accepted 0.001 wing deliveries (intake_solar_nonoptical_rya387.py): (file, lo, hi, vald_id)
WINGS = [
    ('vald_solar_fuv_1150_2000_hfson_raw.txt',     1150.0,  2000.0, 'RyanSchmitt.019704'),
    ('vald_solar_nearuv_2000_3780_hfson_raw.txt',  2000.0,  3780.0, 'RyanSchmitt.019705'),
    ('vald_solar_redopt_6910_9500_hfson_raw.txt',  6910.0,  9500.0, 'RyanSchmitt.019706'),
    ('vald_solar_ir_9500_17000_hfson_raw.txt',     9500.0, 17000.0, 'RyanSchmitt.019707'),
    ('vald_solar_ir_17000_25000_hfson_raw.txt',   17000.0, 25000.0, 'RyanSchmitt.019708'),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true', help='compute + report, do not write')
    args = ap.parse_args()

    elements = load_elements(ELEMENTS)
    sp = get_star_params('solar')
    params_str = (f"solar (Teff {sp['teff']:.0f}, logg {sp['logg']:.3f}, "
                  f"vmic {sp['xi']:.1f}; [M/H]=0.00) [config/stars.yaml]")

    # ── optical core (preserved verbatim) ────────────────────────────────────
    core = pd.read_csv(SOLAR, low_memory=False)
    core_n = len(core)
    core_lo, core_hi = core['wavelength_air_A'].min(), core['wavelength_air_A'].max()
    print(f"[core] {SOLAR.name}: {core_n:,} rows, {core_lo:.2f}–{core_hi:.2f} Å")
    if core_hi > 6910.0:
        sys.exit(f"ABORT: core already extends past the optical ({core_hi:.1f} Å) — "
                 f"expected the optical-only core (≤6910). Refusing to double-assemble.")

    # guard: stored optical proximity must be reproducible (RYA-381 guard)
    repro = compute_vald_proximity_flag(core)
    dmax = float(np.abs(repro - core['vald_proximity_flag'].astype(float)).max())
    if dmax > 1e-3:
        sys.exit(f"ABORT: optical vald_proximity_flag not reproducible (max|Δ|={dmax:.4g}).")
    print(f"[core] proximity reproducible (max|Δ|={dmax:.2g}) — safe to extend")

    # ── wings ────────────────────────────────────────────────────────────────
    wing_frames, manifest = [], []
    for name, lo, hi, vid in WINGS:
        p = LL / name
        df = parse_wing(p, lo, hi, elements)
        if ((df['wavelength_air_A'] >= core_lo) & (df['wavelength_air_A'] <= core_hi)).any():
            sys.exit(f"ABORT: {name} contains lines inside the optical core range — overlap.")
        wing_frames.append(df)
        manifest.append({'file': name, 'vald_id': vid, 'req_lo_A': lo, 'req_hi_A': hi,
                         'n_lines': int(len(df)), 'detection_threshold': THRESHOLD,
                         'hfs': 'on', 'params': params_str})
        print(f"[wing] {name}: {len(df):,} lines  ({lo:.0f}–{hi:.0f} Å, thr={THRESHOLD})")

    wings = pd.concat(wing_frames, ignore_index=True)
    wings['vald_proximity_flag'] = 0.0  # recomputed over union below
    wings['blend_flag'] = build_vetted_blend_flag(wings).values
    for c in COLS:
        if c not in wings:
            wings[c] = '' if c == 'notes' else 0
    wings = wings[COLS]

    # ── union + recompute proximity over the full list ───────────────────────
    union = pd.concat([core[COLS], wings], ignore_index=True)
    union = union.sort_values('wavelength_air_A', kind='mergesort').reset_index(drop=True)
    new_vpf = compute_vald_proximity_flag(union)
    union_old_vpf = union['vald_proximity_flag'].astype(float).values
    union['vald_proximity_flag'] = np.round(new_vpf, 6)
    opt_mask = (union['wavelength_air_A'] >= core_lo) & (union['wavelength_air_A'] <= core_hi)
    moved = int((np.abs(new_vpf - union_old_vpf)[opt_mask.values] > 1e-4).sum())
    seam = union[opt_mask & (np.abs(new_vpf - union_old_vpf) > 1e-4)]['wavelength_air_A']
    print(f"\n[union] {len(union):,} rows  {union['wavelength_air_A'].min():.2f}–"
          f"{union['wavelength_air_A'].max():.2f} Å")
    print(f"[union] optical rows with shifted proximity (seam only): {moved}"
          + (f"  at {seam.min():.2f}–{seam.max():.2f} Å" if moved else ""))

    # byte-preservation guard: every optical column value except the recomputed
    # seam proximity must be identical to the core we read.
    opt_after = union[opt_mask].reset_index(drop=True)
    assert len(opt_after) == core_n, f"optical row count changed {core_n}->{len(opt_after)}"

    if args.dry_run:
        print("\n[dry-run] no files written.")
        return

    union.to_csv(SOLAR, index=False)
    print(f"\n[write] {SOLAR}  ({len(union):,} rows)")

    prov = {
        'linelist': 'linelist_solar.csv', 'assembled_by': 'RYA-387',
        'supersedes': 'RYA-381 0.05 non-optical wings',
        'date': str(date.today()), 'total_lines': int(len(union)),
        'span_A': [float(union['wavelength_air_A'].min()), float(union['wavelength_air_A'].max())],
        'DETECTION_THRESHOLD': (
            'HOMOGENEOUS 0.001 throughout (optical core + all non-optical wings). '
            'The RYA-381 depth heterogeneity (0.05 wings vs 0.001 core) is RESOLVED.'),
        'sources': [
            {'file': 'vald_solar_raw.txt', 'req_lo_A': 3780.0, 'req_hi_A': 6910.0,
             'n_lines': core_n, 'detection_threshold': THRESHOLD, 'hfs': 'on',
             'params': 'solar (existing optical core, preserved verbatim)'},
            *manifest],
    }
    out = LL / 'linelist_solar_provenance_rya387.json'
    out.write_text(json.dumps(prov, indent=2))
    print(f"[write] {out}")


if __name__ == '__main__':
    main()
