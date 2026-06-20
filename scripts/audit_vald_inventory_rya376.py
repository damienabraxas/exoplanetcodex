"""
scripts/audit_vald_inventory_rya376.py
======================================
RYA-376 — AUDIT the VALD line-list inventory: coverage matrix + near-IR readiness.

Read-only. Inventories every VALD asset under data/linelists/ — raw VALD3 "Extract
Stellar" deliveries (vald_*.txt) and assembled line lists (linelist_*.csv,
canonical_gf.csv) — and answers: what wavelength / element / band coverage do we
hold, and does any near-IR (0.95–2.5 µm) coverage exist anywhere?

Per raw delivery it runs the codex-vald-extraction intake gate:
  * truncation — the VALD web 100,000-transition cap (RYA-64); a capped delivery
    carries the warning on line 1 (read_vald_header) AND/OR n_selected > 100k with a
    short data block. Flagged either way.
  * coverage   — requested vs actual wavelength span + per-band reach.
  * HFS        — hyperfine splitting on/off (from the filename convention + record count).
  * ACCEPT / REJECT verdict.

Bands (air Å): optical <9500 · Y 9500–11000 · J 11000–14000 · H 14000–18000 ·
K 18000–25000 (0.95 / 1.10 / 1.40 / 1.80 / 2.50 µm).

Writes: data/audit/vald_inventory/{assets.csv, coverage_matrix.csv}.
Usage:  python scripts/audit_vald_inventory_rya376.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LL = ROOT / 'data' / 'linelists'
OUT = ROOT / 'data' / 'audit' / 'vald_inventory'
sys.path.insert(0, str(LL))
from vald_parse import read_vald_header, is_vald_data_line   # noqa: E402

VALD_WEB_CAP = 100000
BANDS = [('optical', 3000, 9500), ('Y', 9500, 11000), ('J', 11000, 14000),
         ('H', 14000, 18000), ('K', 18000, 25000)]


def band_of(wl_A: float) -> str:
    for name, lo, hi in BANDS:
        if lo <= wl_A < hi:
            return name
    return 'other'


# ── Raw VALD deliveries ───────────────────────────────────────────────────────

def scan_raw(path: Path) -> dict:
    """Lightweight scan: header + per-record (species, wavelength) only."""
    hdr = read_vald_header(path)
    wls, elems, ions, species = [], set(), set(), set()
    n = 0
    for line in open(path, errors='replace'):
        hit = is_vald_data_line(line)
        if hit is None:
            continue
        sp, rest = hit
        try:
            wl = float(rest.split(',')[1])
        except (ValueError, IndexError):
            continue
        n += 1
        wls.append(wl)
        parts = sp.split()
        elems.add(parts[0])
        ions.add(sp)
        species.add(parts[0])
    wls = np.array(wls)
    # truncation: explicit warning, OR delivered count at/above the web cap while the
    # header claims more selected than were delivered.
    truncated = bool(hdr['truncated']) or (n >= VALD_WEB_CAP and n < hdr['n_selected'])
    over_cap = hdr['n_selected'] > VALD_WEB_CAP
    band_counts = {b: int(((wls >= lo) & (wls < hi)).sum()) for b, lo, hi in BANDS}
    name = path.name
    hfs = ('off' if 'hfsoff' in name else 'on' if 'hfson' in name else 'unknown')
    return {
        'asset': name, 'type': 'raw', 'star': _star_of(name),
        'req_lo_A': hdr['wl_start'], 'req_hi_A': hdr['wl_end'],
        'actual_lo_A': round(float(wls.min()), 1) if n else np.nan,
        'actual_hi_A': round(float(wls.max()), 1) if n else np.nan,
        'n_selected_hdr': hdr['n_selected'], 'n_data_lines': n,
        'truncated': truncated, 'over_100k_cap': over_cap,
        'hfs': hfs, 'quarantine': 'quarantine' in name,
        'n_elements': len(elems), 'elements': ','.join(sorted(elems)),
        'gf_source': 'VALD3 (raw)',
        **{f'band_{b}': band_counts[b] for b, _, _ in BANDS},
        'verdict': _raw_verdict(truncated, n, hdr, 'quarantine' in name),
    }


def _raw_verdict(truncated, n, hdr, quarantine):
    if quarantine:
        return 'QUARANTINE (excluded by curator — HFS-off / superseded)'
    if truncated:
        return 'REJECT — truncated at web cap (re-extract in sub-ranges)'
    if n == 0:
        return 'REJECT — no parseable transitions'
    return 'ACCEPT'


def _star_of(name: str) -> str:
    for key in ('55cnc', 'alpha_cen_a', 'alpha_cen_b', 'procyon', 'solar'):
        if key in name:
            return key
    return '?'


# ── Assembled line lists ──────────────────────────────────────────────────────

def scan_csv(path: Path) -> dict:
    cols = pd.read_csv(path, nrows=0, comment='#').columns.tolist()
    wcol = 'wavelength_air_A' if 'wavelength_air_A' in cols else None
    if wcol is None:
        return {'asset': path.name, 'type': 'assembled', 'verdict': 'SKIP (no wl column)'}
    usecols = [c for c in (wcol, 'element', 'loggf_source', 'loggf_reference') if c in cols]
    df = pd.read_csv(path, usecols=usecols, comment='#')
    df[wcol] = pd.to_numeric(df[wcol], errors='coerce')
    df = df.dropna(subset=[wcol])
    if df.empty:
        return {'asset': path.name, 'type': 'assembled', 'star': _star_of(path.name),
                'n_data_lines': 0, 'verdict': 'SKIP (no numeric wavelengths — not a line list)',
                **{f'band_{b}': 0 for b, _, _ in BANDS}}
    w = df[wcol].to_numpy(dtype=float)
    gfcol = 'loggf_source' if 'loggf_source' in df else ('loggf_reference' if 'loggf_reference' in df else None)
    gf = (df[gfcol].value_counts().to_dict() if gfcol else {})
    band_counts = {b: int(((w >= lo) & (w < hi)).sum()) for b, lo, hi in BANDS}
    return {
        'asset': path.name, 'type': 'assembled', 'star': _star_of(path.name),
        'req_lo_A': np.nan, 'req_hi_A': np.nan,
        'actual_lo_A': round(float(w.min()), 1), 'actual_hi_A': round(float(w.max()), 1),
        'n_selected_hdr': np.nan, 'n_data_lines': len(df),
        'truncated': np.nan, 'over_100k_cap': np.nan,
        'hfs': 'n/a', 'quarantine': 'quarantine' in path.name,
        'n_elements': int(df['element'].nunique()) if 'element' in df else np.nan,
        'elements': ','.join(sorted(df['element'].dropna().unique())) if 'element' in df else '',
        'gf_source': '; '.join(f'{k}:{v}' for k, v in list(gf.items())[:4]),
        **{f'band_{b}': band_counts[b] for b, _, _ in BANDS},
        'verdict': ('has NIR' if band_counts['Y'] + band_counts['J'] + band_counts['H'] + band_counts['K'] > 0
                    else 'optical-only'),
        '_df': df.rename(columns={wcol: 'wl'}),
    }


# ── Coverage matrix (element × band) over the NIR-bearing assets ──────────────

def coverage_matrix(csv_assets: dict, targets: list) -> pd.DataFrame:
    """For each target element, which bands are present, aggregated over the assembled
    assets that carry NIR (proves the VALD atomic data is held & extractable)."""
    rows = []
    for el in targets:
        row = {'element': el}
        for b, _, _ in BANDS:
            present = set()
            for a in csv_assets.values():
                df = a.get('_df')
                if df is None or 'element' not in df:
                    continue
                lo, hi = next((lo, hi) for bn, lo, hi in BANDS if bn == b)
                sub = df[(df['element'] == el) & (df['wl'] >= lo) & (df['wl'] < hi)]
                if len(sub):
                    present.add(a['asset'].replace('linelist_', '').replace('.csv', ''))
            row[b] = ('✓' + (f"({','.join(sorted(present))})" if b != 'optical' else '')) if present else '·'
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    targets = ['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'K', 'Ca', 'Sc', 'Ti',
               'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Sr', 'Y', 'Zr', 'Ba',
               'Eu', 'Li']

    raw_assets = {}
    for p in sorted(LL.glob('vald_*.txt')):
        try:
            raw_assets[p.name] = scan_raw(p)
        except Exception as exc:                              # surface, never skip
            raw_assets[p.name] = {'asset': p.name, 'type': 'raw', 'verdict': f'ERROR: {exc}'}

    csv_assets = {}
    for p in sorted(LL.glob('*.csv')):
        try:
            csv_assets[p.name] = scan_csv(p)
        except Exception as exc:
            csv_assets[p.name] = {'asset': p.name, 'type': 'assembled', 'verdict': f'ERROR: {exc}'}

    # Asset table (drop heavy internal cols)
    rows = []
    for a in list(raw_assets.values()) + list(csv_assets.values()):
        rows.append({k: v for k, v in a.items() if not k.startswith('_')})
    assets = pd.DataFrame(rows)
    assets.to_csv(OUT / 'assets.csv', index=False)

    matrix = coverage_matrix(csv_assets, targets)
    matrix.to_csv(OUT / 'coverage_matrix.csv', index=False)

    # ── Console report ────────────────────────────────────────────────────────
    pd.set_option('display.width', 200, 'display.max_columns', 40, 'display.max_colwidth', 28)
    print("\n================ RAW VALD DELIVERIES ================")
    rcols = ['asset', 'star', 'req_lo_A', 'req_hi_A', 'actual_hi_A', 'n_data_lines',
             'truncated', 'over_100k_cap', 'hfs', 'band_Y', 'band_J', 'band_H', 'band_K', 'verdict']
    print(assets[assets.type == 'raw'][rcols].to_string(index=False))

    print("\n================ ASSEMBLED LINE LISTS ================")
    acols = ['asset', 'star', 'actual_lo_A', 'actual_hi_A', 'n_data_lines',
             'band_Y', 'band_J', 'band_H', 'band_K', 'gf_source', 'verdict']
    print(assets[assets.type == 'assembled'][acols].to_string(index=False))

    print("\n================ COVERAGE MATRIX (element × band, assembled assets) ================")
    print(matrix.to_string(index=False))

    # Near-IR readiness headline
    nir_csv = [a for a in csv_assets.values()
               if sum(a.get(f'band_{b}', 0) for b in ('Y', 'J', 'H', 'K')) > 0]
    print("\n================ NEAR-IR READINESS ================")
    for a in nir_csv:
        print(f"  {a['asset']}: {a['actual_lo_A']:.0f}-{a['actual_hi_A']:.0f} Å  "
              f"NIR lines Y/J/H/K = {a['band_Y']}/{a['band_J']}/{a['band_H']}/{a['band_K']}  ({a['verdict']})")
    print(f"\n  [out] {OUT/'assets.csv'}  +  {OUT/'coverage_matrix.csv'}")


if __name__ == '__main__':
    main()
