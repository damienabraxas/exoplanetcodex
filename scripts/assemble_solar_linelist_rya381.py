#!/usr/bin/env python3
"""
scripts/assemble_solar_linelist_rya381.py
==========================================
RYA-381 — ASSEMBLE the accepted non-optical solar VALD deliveries into the
extended solar line list.

Extends data/linelists/linelist_solar.csv from the optical-only 3780–6910 Å pull
to the full delivered 1150–25000 Å span by APPENDING the six RYA-381 wing
deliveries (FUV / near-UV / red-optical / IR×3), all solar-param + HFS-on
(gated by intake_solar_nonoptical_rya381.py).

Design — preserve the curated optical core, add wings homogeneously:
  * The existing optical rows (3780–6910 Å) are kept VERBATIM — they carry the
    curated atomic data (RYA-365 Ni 6300 / RYA-368 gf adjudications, NIST
    injections, vetted blend_flag). This script never edits them.
  * Wing rows are parsed with the SAME parser/priority/NIST logic the optical
    core was built with (scripts/build_linelist) so the schema + derived columns
    are homogeneous.
  * vald_proximity_flag is recomputed over the FULL union (pipeline/build_linelist
    formula). It reproduces the stored optical values exactly except at the
    genuine 3780 Å seam, where near-UV neighbours now (correctly) contribute.
  * blend_flag (vetted, RYA-209/358) is preserved verbatim for optical rows;
    wing rows get the vetted logic (no vetted blend falls in the wings → False).

DETECTION-THRESHOLD HETEROGENEITY (RYA-381 spec §3 — flagged, not silently mixed):
  the wings were extracted at central-depth threshold 0.05, the optical core at
  0.001 (50× shallower). The assembled list is therefore depth-INHOMOGENEOUS:
  dense in the optical, shallow in the wings. This is recorded per-source in the
  provenance sidecar and the file header comment; downstream weak-line / blend
  work in the wings needs a 0.001 re-extraction (owed item, see Linear).

Usage:  python scripts/assemble_solar_linelist_rya381.py [--dry-run]
Out:    data/linelists/linelist_solar.csv  (extended, in place)
        data/linelists/linelist_solar_provenance_rya381.json  (per-source manifest)
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
LL = ROOT / 'data' / 'linelists'

from scripts.build_linelist import (parse_vald, load_elements, filter_and_assign,  # noqa: E402
                                     crosscheck_nist)
from pipeline.build_linelist import (compute_vald_proximity_flag,                  # noqa: E402
                                     build_vetted_blend_flag)

SOLAR = LL / 'linelist_solar.csv'
ELEMENTS = ROOT / 'data' / 'config' / 'elements_master.json'
NIST = LL / 'nist_crosscheck.csv'

# Accepted wing deliveries (intake_solar_nonoptical_rya381.py): (file, lo, hi, vald_id, thr)
WINGS = [
    ('vald_solar_fuv_1150_2000_hfson_raw.txt',    1150.0,  2000.0, 'RyanSchmitt.019691', 0.05),
    ('vald_solar_nearuv_2000_3780_hfson_raw.txt', 2000.0,  3780.0, 'RyanSchmitt.019692', 0.05),
    ('vald_solar_redopt_6910_9500_hfson_raw.txt', 6910.0,  9500.0, 'RyanSchmitt.019693', 0.05),
    ('vald_solar_ir_9500_14000_hfson_raw.txt',    9500.0, 14000.0, 'RyanSchmitt.019694', 0.05),
    ('vald_solar_ir_14000_18000_hfson_raw.txt',  14000.0, 18000.0, 'RyanSchmitt.019695', 0.05),
    ('vald_solar_ir_18000_25000_hfson_raw.txt',  18000.0, 25000.0, 'RyanSchmitt.019696', 0.05),
]
OPTICAL_THRESHOLD = 0.001  # vald_solar_raw.txt (existing core)
COLS = ['element', 'ion', 'wavelength_air_A', 'excitation_potential_eV', 'log_gf',
        'loggf_source', 'nist_grade', 'damping_rad', 'damping_stark', 'damping_vdW',
        'central_depth', 'vald_proximity_flag', 'blend_flag', 'priority', 'notes']


def parse_wing(path: Path, lo: float, hi: float, elements: dict) -> pd.DataFrame:
    """Parse one wing delivery into the linelist schema (priority + NIST grade)."""
    rows = parse_vald(path)                       # VALD3 → base dicts (blend_flag=False, priority=0)
    rows = filter_and_assign(rows, elements)      # priority from elements_master.json
    rows = crosscheck_nist(rows, NIST)            # NIST grades where they apply
    df = pd.DataFrame(rows)
    # keep only this delivery's own λ window (deliveries are non-overlapping by design,
    # but clip defensively so a stray edge line can't double-count across chunks)
    df = df[(df['wavelength_air_A'] >= lo) & (df['wavelength_air_A'] < hi)].copy()
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true', help='compute + report, do not write')
    args = ap.parse_args()

    elements = load_elements(ELEMENTS)

    # ── optical core (preserved verbatim) ────────────────────────────────────
    core = pd.read_csv(SOLAR, low_memory=False)
    core_n = len(core)
    core_lo, core_hi = core['wavelength_air_A'].min(), core['wavelength_air_A'].max()
    print(f"[core] {SOLAR.name}: {core_n:,} rows, {core_lo:.2f}–{core_hi:.2f} Å")

    # guard: the optical core's stored proximity must be reproducible (else we do not
    # understand the file → refuse to rewrite it).
    repro = compute_vald_proximity_flag(core)
    dmax = float(np.abs(repro - core['vald_proximity_flag'].astype(float)).max())
    if dmax > 1e-3:
        sys.exit(f"ABORT: optical vald_proximity_flag not reproducible (max|Δ|={dmax:.4g}); "
                 f"refusing to rewrite a core I cannot regenerate.")
    print(f"[core] proximity reproducible (max|Δ|={dmax:.2g}) — safe to extend")

    # ── wings ────────────────────────────────────────────────────────────────
    wing_frames, manifest = [], []
    for name, lo, hi, vid, thr in WINGS:
        p = LL / name
        df = parse_wing(p, lo, hi, elements)
        # ensure no overlap with the optical core or another wing's range
        if ((df['wavelength_air_A'] >= core_lo) & (df['wavelength_air_A'] <= core_hi)).any():
            sys.exit(f"ABORT: {name} contains lines inside the optical core range — overlap.")
        wing_frames.append(df)
        manifest.append({'file': name, 'vald_id': vid, 'req_lo_A': lo, 'req_hi_A': hi,
                         'n_lines': int(len(df)), 'detection_threshold': thr, 'hfs': 'on',
                         'params': 'solar (Teff 5777, logg 4.44, vmic 1.0; [M/H]=0.00)'})
        print(f"[wing] {name}: {len(df):,} lines  ({lo:.0f}–{hi:.0f} Å, thr={thr})")

    wings = pd.concat(wing_frames, ignore_index=True)
    # wing blend_flag via the vetted logic (RYA-209/358) — none of the vetted blends
    # fall in the wings, so all False; this proves it rather than assuming.
    wings['vald_proximity_flag'] = 0.0  # placeholder; recomputed over union below
    wings['blend_flag'] = build_vetted_blend_flag(wings).values
    for c in COLS:
        if c not in wings:
            wings[c] = '' if c == 'notes' else 0
    wings = wings[COLS]

    # ── union + recompute proximity over the full list ───────────────────────
    union = pd.concat([core[COLS], wings], ignore_index=True)
    union = union.sort_values('wavelength_air_A', kind='mergesort').reset_index(drop=True)
    new_vpf = compute_vald_proximity_flag(union)
    # report how many optical rows shift (only the 3780 seam should move)
    union_old_vpf = union['vald_proximity_flag'].astype(float).values
    union['vald_proximity_flag'] = np.round(new_vpf, 6)
    opt_mask = (union['wavelength_air_A'] >= core_lo) & (union['wavelength_air_A'] <= core_hi)
    moved = int((np.abs(new_vpf - union_old_vpf)[opt_mask.values] > 1e-4).sum())
    seam = union[opt_mask & (np.abs(new_vpf - union_old_vpf) > 1e-4)]['wavelength_air_A']
    print(f"\n[union] {len(union):,} rows  {union['wavelength_air_A'].min():.2f}–"
          f"{union['wavelength_air_A'].max():.2f} Å")
    print(f"[union] optical rows with shifted proximity (seam only): {moved}"
          + (f"  at {seam.min():.2f}–{seam.max():.2f} Å" if moved else ""))

    if args.dry_run:
        print("\n[dry-run] no files written.")
        return

    union.to_csv(SOLAR, index=False)
    print(f"\n[write] {SOLAR}  ({len(union):,} rows)")

    prov = {
        'linelist': 'linelist_solar.csv', 'assembled_by': 'RYA-381',
        'date': str(date.today()), 'total_lines': int(len(union)),
        'span_A': [float(union['wavelength_air_A'].min()), float(union['wavelength_air_A'].max())],
        'DEPTH_HETEROGENEITY': (
            'optical core threshold 0.001; non-optical wings threshold 0.05 (50x '
            'shallower). Assembled list is NOT depth-homogeneous. Weak-line/blend '
            'work beyond the optical needs a 0.001 wing re-extraction.'),
        'sources': [
            {'file': 'vald_solar_raw.txt', 'req_lo_A': 3780.0, 'req_hi_A': 6910.0,
             'n_lines': core_n, 'detection_threshold': OPTICAL_THRESHOLD, 'hfs': 'on',
             'params': 'solar (existing optical core, preserved verbatim)'},
            *manifest],
    }
    out = LL / 'linelist_solar_provenance_rya381.json'
    out.write_text(json.dumps(prov, indent=2))
    print(f"[write] {out}")


if __name__ == '__main__':
    main()
