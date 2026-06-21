"""
scripts/build_ba_nlte_cds.py
============================
Build the Ba II NLTE-correction grid (RYA-396) from the published, machine-readable
CDS catalogue **Korotin et al. 2015, A&A 581, A70** (`J/A+A/581/A70`) — verified to
EXIST on VizieR (unlike the Ti/Cr grids, 404 per RYA-241) and to cover solar params.

This is NOT a scraper: VizieR serves the full grid as a table. Columns `d{i}{baFe}` are
the NLTE corrections (A_NLTE − A_LTE) for the four Ba II lines
  line 1 = 4554.029 Å, 2 = 5853.668 Å, 3 = 6141.713 Å, 4 = 6496.897 Å
at [Ba/Fe] = −0.2 / +0.1 / +0.4, over (Teff, logg, [Fe/H], Vt, α-flag).

Reduction to the RYA-235 (Teff,logg,feh) interpolator schema:
  * non-α-enhanced rows (disk stars); Vt = 1 km/s (solar microturbulence; the Vt
    dependence is ≲0.02 dex);
  * the [Ba/Fe] strength axis is interpolated to **[Ba/Fe] = 0** (solar-scaled — we
    report differential [Ba/Fe]); a real interpolation between published nodes, never a
    faked axis (the RYA-242 sin). The strength axis can be re-introduced later for the
    multi-star phase (analogous to the INSPECT A_LTE axis).

Output: data/nlte_grids/Ba_Korotin2015.csv (+ .prov.json), schema
  element, ion, wave_A, teff_K, logg, feh, delta_nlte.

Usage:  python scripts/build_ba_nlte_cds.py
"""
import json
import sys
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'nlte_grids' / 'Ba_Korotin2015.csv'
SRC = 'J/A+A/581/A70'
TSV_URL = ('https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/A%2BA/581/A70'
           '&-out.all&-out.max=unlimited')

LINES = {1: 4554.029, 2: 5853.668, 3: 6141.713, 4: 6496.897}   # Ba II, air Å
BAFE_NODES = (-0.2, 0.1, 0.4)
VT_FIX = 1.0                  # km/s — solar microturbulence node
SOLAR = (5772.0, 4.438, 0.0)


def _download() -> pd.DataFrame:
    req = urllib.request.Request(TSV_URL, headers={'User-Agent': 'Mozilla/5.0'})
    raw = urllib.request.urlopen(req, timeout=60).read().decode('utf-8', 'replace')
    lines = [l for l in raw.splitlines() if l and not l.startswith('#')]
    hdr = lines[0].split('\t')
    rows = [l.split('\t') for l in lines[3:] if l.strip() and not l.startswith('-')]
    df = pd.DataFrame(rows, columns=hdr)
    for c in ('Teff', 'logg', 'Vt', '[Fe/H]'):
        df[c] = df[c].astype(float)
    return df


def _interp_bafe_zero(row, i):
    """Interpolate line i's NLTE correction to [Ba/Fe]=0 between the published nodes."""
    cols = [f'd{i}-0.2', f'd{i}+0.1', f'd{i}+0.4']   # [Ba/Fe] = -0.2, +0.1, +0.4
    y = [float(row[c]) for c in cols]
    return float(np.interp(0.0, BAFE_NODES, y))


def build():
    print(f'=== Ba II NLTE grid (CDS {SRC}, Korotin+2015) ===')
    df = _download()
    print(f'downloaded {len(df)} rows; α-flag values {df["a"].value_counts().to_dict()}')

    # disk stars: non-α-enhanced (blank flag); solar microturbulence node
    sel = df[(df['a'].str.strip() == '') & (df['Vt'] == VT_FIX)].copy()
    print(f'non-α, Vt={VT_FIX}: {len(sel)} rows  '
          f'Teff[{sel.Teff.min():.0f},{sel.Teff.max():.0f}] '
          f'logg[{sel.logg.min():.1f},{sel.logg.max():.1f}] '
          f'feh[{sel["[Fe/H]"].min():.1f},{sel["[Fe/H]"].max():.1f}]')

    out = []
    for _, r in sel.iterrows():
        for i, wave in LINES.items():
            d0 = _interp_bafe_zero(r, i)
            out.append({'element': 'Ba', 'ion': 'II', 'wave_A': wave,
                        'teff_K': int(r['Teff']), 'logg': float(r['logg']),
                        'feh': float(r['[Fe/H]']), 'delta_nlte': round(d0, 4)})
    grid = pd.DataFrame(out)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    grid.to_csv(OUT, index=False)
    print(f'saved {len(grid)} rows → {OUT}')

    # solar anchor (line 2 = 5853, our measured Ba II line)
    g2 = grid[grid['wave_A'] == LINES[2]]
    f = LinearNDInterpolator(g2[['teff_K', 'logg', 'feh']].values, g2['delta_nlte'].values)
    solar_delta = float(f([SOLAR])[0])
    print(f'solar anchor Δ (5853 Å) = {solar_delta:+.4f} dex '
          f'(Teff={SOLAR[0]}, logg={SOLAR[1]}, feh={SOLAR[2]})')
    if np.isnan(solar_delta):
        print('CRITICAL: grid returns NaN at solar params'); sys.exit(1)

    prov = {
        'element': 'Ba', 'ion': 'II', 'reference': 'Korotin et al. 2015, A&A 581, A70',
        'source': f'CDS VizieR {SRC}', 'source_url': TSV_URL,
        'convention': 'delta_nlte = A_NLTE - A_LTE (dex)',
        'lines_air_A': LINES,
        'reduction': {'alpha': 'non-enhanced (disk)', 'Vt_kms': VT_FIX,
                      'BaFe': 'interpolated to [Ba/Fe]=0 between nodes (-0.2,+0.1,+0.4)'},
        'extraction_date': date.today().isoformat(),
        'n_rows': int(len(grid)), 'n_lines': len(LINES),
        'coverage': {'teff_K': [int(grid.teff_K.min()), int(grid.teff_K.max())],
                     'logg': [float(grid.logg.min()), float(grid.logg.max())],
                     'feh': [float(grid.feh.min()), float(grid.feh.max())]},
        'solar_anchor': {'teff_K': SOLAR[0], 'logg': SOLAR[1], 'feh': SOLAR[2],
                         'line_A': LINES[2], 'delta_nlte': round(solar_delta, 4)},
    }
    OUT.with_suffix('.prov.json').write_text(json.dumps(prov, indent=2))
    print(f'saved provenance → {OUT.with_suffix(".prov.json")}')
    return solar_delta


if __name__ == '__main__':
    build()
