#!/usr/bin/env python3
"""
scripts/drop_fe_leg_placeholders_rya417.py
==========================================
RYA-417 Part C (data fix) — drop the LIVE-VERIFIED TRUE_PLACEHOLDER lines from the
committed MPIA Bergemann Fe leg (Fe_Bergemann_MPIA.csv) and the stale Si grid
(Si_Bergemann_MPIA.csv). Drops ONLY lines the live MPIA re-query confirmed as
TRUE_PLACEHOLDER (served-but-literal-0 at every node); GENUINE_NEAR_ZERO is never
dropped. Reads the verdict table written by audit_fe_leg_placeholders_rya417.py.

Each dropped line is recorded (with its live-query verdict) in the grid's .prov.json.
Idempotent: re-running on an already-clean grid drops nothing.

    python scripts/drop_fe_leg_placeholders_rya417.py            # write clean grids
    python scripts/drop_fe_leg_placeholders_rya417.py --dry-run  # report only
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
GRIDS = _REPO / 'data' / 'nlte_grids'
AUDIT = _REPO / 'data' / 'audit' / 'fe_leg_placeholders_rya417.csv'
_TOL = 1e-3   # Å — round-match grid wave to the audited wave


def _drop_one(grid_csv: Path, drops: pd.DataFrame, dry: bool) -> int:
    df = pd.read_csv(grid_csv)
    has_ion = 'ion' in df.columns
    mask = pd.Series(False, index=df.index)
    hit_lines = []
    for _, r in drops.iterrows():
        w = float(r['wave_A'])
        m = (df['wave_A'] - w).abs() < _TOL
        if has_ion and not pd.isna(r.get('ion')):
            m &= (df['ion'].astype(str) == str(r['ion']))
        if m.any():
            mask |= m
            hit_lines.append((str(r.get('ion', '')), round(w, 3), int(m.sum())))
    if not hit_lines:
        print(f"  {grid_csv.name}: already clean (no verified placeholders present)")
        return 0
    print(f"  {grid_csv.name}: dropping {len(hit_lines)} verified TRUE_PLACEHOLDER line(s), "
          f"{int(mask.sum())} rows:")
    for ion, w, n in hit_lines:
        print(f"      {ion+' ' if ion else ''}{w}  ({n} nodes)")
    if dry:
        return len(hit_lines)
    df[~mask].to_csv(grid_csv, index=False)

    prov = grid_csv.with_suffix('.prov.json')
    meta = json.loads(prov.read_text()) if prov.exists() else {}
    meta.setdefault('rya417_dropped_placeholders', [])
    meta['rya417_dropped_placeholders'] += [
        {'ion': ion, 'wave_A': w, 'nodes_dropped': n,
         'verdict': 'TRUE_PLACEHOLDER (live MPIA re-query: literal 0 at every served node)',
         'basis': 'RYA-417 audit_fe_leg_placeholders_rya417.py'}
        for ion, w, n in hit_lines]
    prov.write_text(json.dumps(meta, indent=2))
    print(f"      provenance → {prov.name}")
    return len(hit_lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    if not AUDIT.exists():
        sys.exit(f"missing {AUDIT} — run audit_fe_leg_placeholders_rya417.py (live) first")
    audit = pd.read_csv(AUDIT)
    verified = audit[audit['verdict'] == 'TRUE_PLACEHOLDER']
    if len(verified) != len(audit):
        # surface anything NOT confirmed — never dropped, but reported
        other = audit[audit['verdict'] != 'TRUE_PLACEHOLDER']
        print(f"NOTE: {len(other)} candidate(s) NOT TRUE_PLACEHOLDER — KEPT (real physics / not served):")
        for _, r in other.iterrows():
            print(f"  KEEP {r['element']} {r['ion']} {r['wave_A']} → {r['verdict']}")

    fe = verified[verified['element'] == 'Fe']
    si = verified[verified['element'] == 'Si']
    n = 0
    print(f"RYA-417 drop (verified TRUE_PLACEHOLDER only):")
    n += _drop_one(GRIDS / 'Fe_Bergemann_MPIA.csv', fe, args.dry_run)
    n += _drop_one(GRIDS / 'Si_Bergemann_MPIA.csv', si, args.dry_run)
    print(f"{'(dry-run) ' if args.dry_run else ''}total verified placeholder lines dropped: {n}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
