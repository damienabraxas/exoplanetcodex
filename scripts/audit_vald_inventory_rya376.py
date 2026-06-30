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

# ── 27-element IR triage (re-baselined 2026-06-30, per RYA-376 ticket) ─────────
# Tag every target element by the VALUE the IR adds, so the downstream extraction
# scope (RYA-427) is driven by value and the IR-vs-UV/near-blue leg split is explicit.
#   IR-PRIMARY    Group A — IR is a new or best channel (extraction priority)
#   IR-CROSSCHECK Group B — optical owns the value; IR confirms (Fe = the IR sigma ruler)
#   IR-NONE       Group C — IR does not help; element argues for the UV/near-blue leg
# 'prior' is the ticket's expectation; this script CONFIRMS it against actual holdings.
IR_TRIAGE = {
    # Group A — IR-PRIMARY
    'C':  ('A', 'IR-PRIMARY',    'C I near-IR multiplets (flagship)'),
    'N':  ('A', 'IR-PRIMARY',    'N I red + NH (flagship)'),
    'O':  ('A', 'IR-PRIMARY',    'O I 777 triplet + OH (flagship)'),
    'S':  ('A', 'IR-PRIMARY',    'S I 1.045 um >> lone weak optical S I 6757 (RYA-488)'),
    'P':  ('A', 'IR-PRIMARY',    'P I 1.05 um is essentially the ONLY real P channel'),
    'K':  ('A', 'IR-PRIMARY',    'K I 1.17 um sidesteps the O2 A-band telluric on optical 7699'),
    'Na': ('A', 'IR-PRIMARY',    'near-IR subordinate lines UNSATURATED where optical resonance saturates (RYA-465)'),
    'Mg': ('A', 'IR-PRIMARY',    'near-IR subordinate lines UNSATURATED where optical resonance saturates (RYA-465)'),
    'Al': ('A', 'IR-PRIMARY',    'near-IR subordinate lines UNSATURATED where optical resonance saturates (RYA-465)'),
    # Group B — IR-CROSSCHECK (optical owns the value)
    'Si': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms'),
    'Ca': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms'),
    'Ti': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms'),
    'Cr': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms'),
    'Mn': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms (HFS element)'),
    'V':  ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms (HFS element)'),
    'Co': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms (HFS element)'),
    'Ni': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms'),
    'Sc': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms (HFS element)'),
    'Cu': ('B', 'IR-CROSSCHECK', 'optical owns the value; IR confirms (HFS element)'),
    'Fe': ('B', 'IR-CROSSCHECK', 'IR SIGMA RULER — most IR lines by far; characterize IR scatter on Fe, apply to sparse CNOPS'),
    'Zn': ('B', 'IR-CROSSCHECK', '27th element, not in ticket triage; optical Zn I 4722/4810/6362 owns value, IR subordinate'),
    # Group C — IR-NONE (argue for the UV/near-blue leg, not IR)
    'Ba': ('C', 'IR-NONE',       'neutron-capture heavy; ionized lines live in blue/near-UV, IR barren'),
    'Y':  ('C', 'IR-NONE',       'neutron-capture heavy; ionized lines live in blue/near-UV, IR barren'),
    'Zr': ('C', 'IR-NONE',       'neutron-capture heavy; ionized lines live in blue/near-UV, IR barren'),
    'Sr': ('C', 'IR-NONE',       'neutron-capture heavy; blue Sr II doublet (RYA-428/433), IR barren'),
    'Eu': ('C', 'IR-NONE',       'neutron-capture heavy; ionized lines live in blue/near-UV, IR barren'),
    'Li': ('C', 'IR-NONE',       'stays optical (Li I 6707)'),
}

# Named Group-A IR diagnostics to CONFIRM (reach + depth) against held solar data.
# (element, ion, target air-Å, band, channel label). OH/NH/CO are MOLECULAR — not in
# the VALD atomic list — so they are reported separately, not looked up here.
GROUP_A_DIAGNOSTICS = [
    ('S', 'I', 10455.45, 'Y', 'S I 1.045 um triplet'),
    ('S', 'I', 10456.76, 'Y', 'S I 1.045 um triplet'),
    ('S', 'I', 10459.41, 'Y', 'S I 1.045 um triplet'),
    ('P', 'I', 10511.59, 'Y', 'P I 1.05 um quartet'),
    ('P', 'I', 10529.52, 'Y', 'P I 1.05 um quartet'),
    ('P', 'I', 10581.58, 'Y', 'P I 1.05 um quartet'),
    ('P', 'I', 10596.90, 'Y', 'P I 1.05 um quartet'),
    ('K', 'I', 11690.22, 'J', 'K I 1.17 um doublet'),
    ('K', 'I', 11769.64, 'J', 'K I 1.17 um doublet'),
    ('C', 'I', 10691.25, 'Y', 'C I near-IR'),
    ('O', 'I', 7771.94,  'optical', 'O I 777 triplet'),
    ('O', 'I', 11302.38, 'J', 'O I near-IR'),
    ('Mg', 'I', 15740.71, 'H', 'Mg I 1.57 um (unsaturated subordinate)'),
    ('Al', 'I', 16718.96, 'H', 'Al I 1.67 um (unsaturated subordinate)'),
]


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


def ir_triage_table(matrix: pd.DataFrame) -> pd.DataFrame:
    """Tag each target element IR-PRIMARY / IR-CROSSCHECK / IR-NONE and CONFIRM the tag
    against actual NIR holdings (does any Y/J/H/K coverage exist for it?)."""
    held = {r['element']: any(str(r[b]).startswith('✓') for b in ('Y', 'J', 'H', 'K'))
            for _, r in matrix.iterrows()}
    rows = []
    for el, (group, tag, why) in IR_TRIAGE.items():
        nir = held.get(el, False)
        # Confirmation logic: Group A/B want NIR present; Group C expects it barren.
        if group in ('A', 'B'):
            confirm = 'CONFIRMED (NIR held)' if nir else 'GAP (tag says IR but NO NIR held)'
        else:
            confirm = 'CONFIRMED (IR barren → UV/blue leg)' if not nir else 'NIR present but low-value (heavy)'
        rows.append({'element': el, 'group': group, 'tag': tag,
                     'nir_held': nir, 'confirm': confirm, 'rationale': why})
    order = {'A': 0, 'B': 1, 'C': 2}
    return pd.DataFrame(rows).sort_values(['group', 'element'], key=lambda s: s.map(order) if s.name == 'group' else s)


def groupA_depth_confirm() -> pd.DataFrame:
    """Confirm the named Group-A IR diagnostics exist with sane depth in the held solar
    list (linelist_solar.csv). Reach + depth come from the data, never from memory."""
    path = LL / 'linelist_solar.csv'
    if not path.exists():
        return pd.DataFrame([{'note': 'linelist_solar.csv absent'}])
    df = pd.read_csv(path, comment='#', low_memory=False,
                     usecols=['element', 'ion', 'wavelength_air_A', 'log_gf', 'central_depth'])
    df['wavelength_air_A'] = pd.to_numeric(df['wavelength_air_A'], errors='coerce')
    df['central_depth'] = pd.to_numeric(df['central_depth'], errors='coerce')
    rows = []
    for el, ion, wl, band, label in GROUP_A_DIAGNOSTICS:
        sub = df[(df['element'] == el) & (df['ion'] == ion)
                 & (df['wavelength_air_A'].between(wl - 0.6, wl + 0.6))]
        if len(sub):
            # HFS-split features deliver many components at one wl — take the deepest.
            best = sub.loc[sub['central_depth'].idxmax()]
            rows.append({'diagnostic': label, 'species': f'{el} {ion}', 'target_A': wl,
                         'band': band, 'covered': True,
                         'matched_A': round(float(best['wavelength_air_A']), 3),
                         'n_components': len(sub),
                         'max_central_depth': round(float(best['central_depth']), 3)})
        else:
            rows.append({'diagnostic': label, 'species': f'{el} {ion}', 'target_A': wl,
                         'band': band, 'covered': False, 'matched_A': np.nan,
                         'n_components': 0, 'max_central_depth': np.nan})
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

    triage = ir_triage_table(matrix)
    triage.to_csv(OUT / 'ir_triage.csv', index=False)

    depth = groupA_depth_confirm()
    depth.to_csv(OUT / 'groupA_depth_confirm.csv', index=False)

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

    print("\n================ 27-ELEMENT IR TRIAGE (tag × confirmation) ================")
    print(triage[['element', 'group', 'tag', 'nir_held', 'confirm']].to_string(index=False))

    print("\n================ GROUP-A IR DIAGNOSTIC DEPTH CONFIRM (held solar list) ================")
    print(depth.to_string(index=False))

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
