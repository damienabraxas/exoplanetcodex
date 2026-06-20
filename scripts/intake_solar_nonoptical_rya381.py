#!/usr/bin/env python3
"""
scripts/intake_solar_nonoptical_rya381.py
=========================================
RYA-381 — INTAKE GATE for the new non-optical solar VALD extractions.

Ryan extracted the solar line list beyond the optical (3780–6910 Å) per the
RYA-376 revisit chunk plan — FUV 1150–2000, near-UV 2000–3780, red-optical
6910–9500, IR 9500–25000 (split 9500–14000 / 14000–18000 / 18000–25000) — all
at solar parameters with hyperfine structure (HFS) ON. This script gates each
delivery through the codex-vald-extraction intake checks BEFORE assembly:

  * completeness — actual λ-span reaches BOTH ends of the requested window and
    every selected line was delivered (n_data == n_selected_hdr); verified by
    span + count, NOT the 100k line-count heuristic (RYA-376 finding).
  * truncation   — VALD web 100,000-transition cap (RYA-64): explicit warning on
    line 1, or n_data at/above the cap below n_selected_hdr.
  * coverage     — delivered range == requested chunk.
  * HFS-on       — content signal: VALD emits an 'hfs:' source tag in the
    reference field ONLY when hyperfine splitting is on. We require hfs: tags to
    be present when HFS-prone species (Sc/V/Mn/Co/Cu/Ba/Eu/La/Pr/Nd) fall in the
    delivered range; absence-with-no-such-species is reported as INDETERMINATE,
    never silently passed.
  * solar params — applied [M/H] from the composition block must be ~0.00
    (parse_vald_applied_metallicity); the structure grid node is reported but
    NOT gated (VALD substitutes the nearest Castelli node — see vald_parse note).
  * detection threshold homogeneity — min central_depth delivered, compared to
    the existing optical pull (vald_solar_raw.txt) so a threshold mismatch is
    flagged, not silently mixed.

Read-only audit (writes a CSV report only). ACCEPT / REJECT per delivery.

Usage:  python scripts/intake_solar_nonoptical_rya381.py
Out:    data/audit/vald_inventory/rya381_intake.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LL = ROOT / 'data' / 'linelists'
OUT = ROOT / 'data' / 'audit' / 'vald_inventory'
sys.path.insert(0, str(LL))
from vald_parse import (read_vald_header, is_vald_data_line,          # noqa: E402
                        parse_vald_applied_metallicity, parse_vald_structure_grid)

VALD_WEB_CAP = 100000
MIN_COVERAGE_FRAC = 0.95   # delivered span must cover >=95% of the requested window
                           # (a truncated pull stops ~halfway; a sparse high-threshold
                           #  pull reaches ~99.8% — this separates the two)
MAX_INTERIOR_GAP_FRAC = 0.10   # no single inter-line gap may exceed this fraction of the
                               # window width (catches a true mid-range cutoff)
METAL_TOL = 0.05        # applied [M/H] must be within this of solar (0.00)
REF_THRESHOLD = 0.001   # the existing optical pull's detection threshold (vald_solar_raw)

# Species with significant hyperfine structure in VALD (odd nuclear spin). If any
# fall in the delivered range, an HFS-on extraction MUST carry hfs: tags for them.
HFS_PRONE = {'Sc', 'V', 'Mn', 'Co', 'Cu', 'Ba', 'Eu', 'La', 'Pr', 'Nd', 'Nb'}

# The new solar non-optical deliveries and the optical reference (for threshold
# homogeneity). (filename, requested_lo, requested_hi). Optical ref last.
DELIVERIES = [
    ('vald_solar_fuv_1150_2000_hfson_raw.txt',    1150.0,  2000.0),
    ('vald_solar_nearuv_2000_3780_hfson_raw.txt', 2000.0,  3780.0),
    ('vald_solar_redopt_6910_9500_hfson_raw.txt', 6910.0,  9500.0),
    ('vald_solar_ir_9500_14000_hfson_raw.txt',    9500.0, 14000.0),
    ('vald_solar_ir_14000_18000_hfson_raw.txt',  14000.0, 18000.0),
    ('vald_solar_ir_18000_25000_hfson_raw.txt',  18000.0, 25000.0),
]
OPTICAL_REF = ('vald_solar_raw.txt', 3780.0, 6910.0)


def scan(path: Path, req_lo: float, req_hi: float) -> dict:
    """Full intake scan of one VALD long-format delivery."""
    hdr = read_vald_header(path)
    wls, depths = [], []
    hfs_tag_count = 0
    species_present = set()
    for line in open(path, errors='replace'):
        if 'hfs:' in line:
            hfs_tag_count += 1
        hit = is_vald_data_line(line)
        if hit is None:
            continue
        sp, rest = hit
        fields = rest.split(',')
        try:
            wls.append(float(fields[1]))
            depths.append(float(fields[13]))
        except (ValueError, IndexError):
            continue
        species_present.add(sp.split()[0])
    wls = np.array(wls)
    depths = np.array(depths)
    n = len(wls)

    truncated = bool(hdr['truncated']) or (n >= VALD_WEB_CAP and n < hdr['n_selected'])
    actual_lo = float(wls.min()) if n else np.nan
    actual_hi = float(wls.max()) if n else np.nan
    width = req_hi - req_lo
    # completeness by SPAN, not line count (RYA-376): the delivered lines must cover
    # essentially the whole requested window with no large interior gap. A sparse
    # high-threshold pull legitimately reaches ~99.8%; a web-cap truncation stops
    # partway (low coverage_frac); a mid-range cutoff shows as one oversized gap.
    coverage_frac = float((actual_hi - actual_lo) / width) if n else 0.0
    if n > 1:
        gaps = np.diff(np.sort(wls))
        max_gap_frac = float(gaps.max() / width)
    else:
        max_gap_frac = 1.0
    span_ok = coverage_frac >= MIN_COVERAGE_FRAC and max_gap_frac <= MAX_INTERIOR_GAP_FRAC
    all_delivered = (n == hdr['n_selected'])
    hfs_species_in_range = sorted(species_present & HFS_PRONE)

    applied_mh = parse_vald_applied_metallicity(path)
    grid = parse_vald_structure_grid(path)

    return {
        'asset': path.name,
        'req_lo_A': req_lo, 'req_hi_A': req_hi,
        'actual_lo_A': round(actual_lo, 2) if n else np.nan,
        'actual_hi_A': round(actual_hi, 2) if n else np.nan,
        'n_selected_hdr': hdr['n_selected'], 'n_data_lines': n,
        'n_processed': hdr['n_processed'],
        'coverage_frac': round(coverage_frac, 4),
        'max_gap_frac': round(max_gap_frac, 4),
        'span_ok': bool(span_ok),
        'all_delivered': bool(all_delivered), 'truncated': truncated,
        'min_depth': round(float(depths.min()), 4) if n else np.nan,
        'vmicro': hdr['vmicro'],
        'applied_MH': None if applied_mh is None else round(applied_mh, 3),
        'struct_grid_MH': grid,
        'hfs_tag_count': hfs_tag_count,
        'hfs_prone_species': ','.join(hfs_species_in_range),
        'n_elements': len(species_present),
    }


def hfs_verdict(row) -> str:
    if row['hfs_tag_count'] > 0:
        return 'ON'
    if row['hfs_prone_species']:
        return 'OFF?? (HFS-prone species present but NO hfs: tags)'
    return 'INDETERMINATE (no HFS-prone species in range)'


def accept_verdict(row, hfs) -> str:
    """Hard intake gates (completeness/coverage/HFS/solar). The detection-threshold
    HOMOGENEITY mismatch is a separate FLAG (see threshold_flag), not a REJECT — the
    delivery is internally valid, it just isn't depth-homogeneous with the optical
    core (RYA-381 spec: 'flag, don't silently mix')."""
    reasons = []
    if row['truncated']:
        reasons.append('TRUNCATED at web cap')
    if not row['span_ok']:
        reasons.append(f"INCOMPLETE span (coverage={row['coverage_frac']:.3f}, "
                       f"max_gap={row['max_gap_frac']:.3f})")
    if not row['all_delivered']:
        reasons.append(f"n_data {row['n_data_lines']} != n_selected {row['n_selected_hdr']}")
    if row['applied_MH'] is None:
        reasons.append('no composition block (applied [M/H] unreadable)')
    elif abs(row['applied_MH']) > METAL_TOL:
        reasons.append(f"applied [M/H]={row['applied_MH']:+.2f} not solar")
    if hfs.startswith('OFF'):
        reasons.append('HFS appears OFF')
    return 'ACCEPT' if not reasons else 'REJECT — ' + '; '.join(reasons)


def threshold_flag(row) -> str:
    """Detection-threshold homogeneity vs the optical core (REF_THRESHOLD=0.001)."""
    md = row.get('min_depth')
    if md is None or (isinstance(md, float) and np.isnan(md)):
        return ''
    if abs(md - REF_THRESHOLD) <= 1e-4:
        return 'homogeneous (0.001)'
    return f'MISMATCH (delivered {md:g} vs optical {REF_THRESHOLD:g}, {md/REF_THRESHOLD:.0f}x shallower)'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, lo, hi in DELIVERIES + [OPTICAL_REF]:
        p = LL / name
        if not p.exists():
            rows.append({'asset': name, 'verdict': 'MISSING'})
            continue
        r = scan(p, lo, hi)
        r['hfs'] = hfs_verdict(r)
        r['threshold'] = threshold_flag(r)
        r['verdict'] = ('REFERENCE (existing optical)' if name == OPTICAL_REF[0]
                        else accept_verdict(r, r['hfs']))
        rows.append(r)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'rya381_intake.csv', index=False)

    pd.set_option('display.width', 240, 'display.max_columns', 40, 'display.max_colwidth', 48)
    print("\n================ RYA-381 SOLAR NON-OPTICAL INTAKE ================")
    show = ['asset', 'req_lo_A', 'req_hi_A', 'actual_lo_A', 'actual_hi_A',
            'n_data_lines', 'coverage_frac', 'truncated', 'applied_MH',
            'hfs', 'verdict']
    print(df[show].to_string(index=False))
    print("\n---- detection-threshold homogeneity (RYA-381 spec §3) ----")
    print(df[['asset', 'min_depth', 'threshold']].to_string(index=False))
    print(f"\n  [out] {OUT / 'rya381_intake.csv'}")

    n_accept = sum(1 for r in rows if str(r.get('verdict', '')).startswith('ACCEPT'))
    n_mismatch = sum(1 for r in rows if 'MISMATCH' in str(r.get('threshold', '')))
    print(f"\n  {n_accept}/{len(DELIVERIES)} new deliveries ACCEPT (hard gates); "
          f"{n_mismatch} carry a detection-threshold MISMATCH flag.")


if __name__ == '__main__':
    main()
