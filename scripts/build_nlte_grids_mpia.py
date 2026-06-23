"""
scripts/build_nlte_grids_mpia.py
=================================
Scrape MPIA Spectrum Tools (https://nlte.mpia.de) to build param-driven NLTE
abundance-correction grids. Ca/Ti/Cr (RYA-245) plus Mn/Co/Si (RYA-396) — all hosted
on the MPIA form, confirmed by probing the live `<select>` line dropdowns
(mnlines[]/colines[]/silines[]). Species codes follow the unambiguous SIU `Z.01`
neutral convention (Z = atomic number): Ca 20.01 / Ti 22.01 / Cr 24.01 / Mn 25.01 /
Co 27.01 / Si 14.01 / Mg 12.01.

Usage:
    python scripts/build_nlte_grids_mpia.py --element Mn
    python scripts/build_nlte_grids_mpia.py --element Mn --resume

Form structure (verified by inspection, RYA-245/396):
  POST https://nlte.mpia.de/gui-siuAC_secE.php
  Fields:
    model      : 'mafags-os'  (covers logg 1.0-5; MARCS only goes to 3.5)
    user_input : multiline 'Name Teff logg [Fe/H] vmic_kms'
    linelist   : '0'
    lines_input: multiline 'wavelength_A species_code'
  Result in: <div id='result'><pre><b>StarName   delta</b></pre>...

Line source (RYA-396): the committed, citable measured set
data/measured/sol_ew_results_v1.csv (NOT the gitignored data/processed/solar_ew.csv).
A provenance sidecar `{grid}.prov.json` is written next to each grid.

Science decision (RYA-243/244/245): MAFAGS-OS model grid used throughout.
"""

import argparse
import json
import time
import sys
import re
from datetime import date
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from pathlib import Path
from scipy.interpolate import LinearNDInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config.constants import PATHS

# ── Config ────────────────────────────────────────────────────────────────────
MPIA_URL   = 'https://nlte.mpia.de/gui-siuAC_secE.php'
MODEL      = 'mafags-os'     # logg 1.0-5, covers solar logg=4.44
LINELIST   = '0'
VMIC       = 1.0             # km/s, fixed for the grid
LINE_BATCH = 20              # lines per HTTP request
STAR_BATCH = 12              # stellar param sets per HTTP request
SLEEP_S    = 2.0             # seconds between requests

# Committed, citable measured-line source (RYA-396) — never the gitignored solar_ew.csv
MEASURED_LINES = Path(__file__).resolve().parents[1] / 'data' / 'measured' / 'sol_ew_results_v1.csv'

# SIU Z.01 neutral codes (Z = atomic number); Mn/Co/Si/Mg hosted, probed RYA-396.
SPECIES = {'Ca': '20.01', 'Ti': '22.01', 'Cr': '24.01',
           'Mn': '25.01', 'Co': '27.01', 'Si': '14.01', 'Mg': '12.01'}

# MPIA per-element line-select dropdown names (probed live, RYA-396).
SELECT_NAME = {'Ca': 'calines[]', 'Ti': 'tilines[]', 'Cr': 'crlines[]',
               'Mn': 'mnlines[]', 'Co': 'colines[]', 'Si': 'silines[]', 'Mg': 'mglines[]'}

OUTFILES = {
    'Ca': 'Ca_Mashonkina2017.csv',
    'Ti': 'Ti_Bergemann2011_MPIA.csv',
    'Cr': 'Cr_Bergemann2010_MPIA.csv',
    'Mn': 'Mn_Bergemann_MPIA.csv',
    'Co': 'Co_Bergemann_MPIA.csv',
    'Si': 'Si_Bergemann_MPIA.csv',
    'Mg': 'Mg_Bergemann_MPIA.csv',
}

REFS = {
    'Ca': 'Mashonkina et al. 2017 (A&A 606, A147)',
    'Ti': 'Bergemann 2011 (MNRAS 413, 2184)',
    'Cr': 'Bergemann & Cescutti 2010 (A&A 522, A9)',
    'Mn': 'Bergemann group (MPIA SpectrumTools MAFAGS-OS)',
    'Co': 'Bergemann group (MPIA SpectrumTools MAFAGS-OS)',
    'Si': 'Bergemann group (MPIA SpectrumTools MAFAGS-OS)',
    'Mg': 'Bergemann group (MPIA SpectrumTools MAFAGS-OS)',
}

SOLAR = (5772.0, 4.438, 0.0)   # RYA-298: canonical solar anchor (IAU/GBS)

# Grid points
TEFF_GRID = [5000, 5500, 5750, 5772, 6000, 6500]   # RYA-298: 5777→5772
LOGG_GRID = [3.5, 4.0, 4.4, 4.5]
FEH_GRID  = [-0.5, 0.0, 0.3]

MATCH_TOL_A = 0.15   # Å — EW file to MPIA wavelength matching tolerance
ERROR_CODES = {10, 20, 30}  # lineformation error, not converged, too weak

OUT_DIR = Path(__file__).resolve().parents[1] / 'data' / 'nlte_grids'
RAW_DIR = OUT_DIR / 'raw'


# ── Build stellar parameter grid ──────────────────────────────────────────────
def _build_grid_params():
    rows = []
    for teff in TEFF_GRID:
        for logg in LOGG_GRID:
            for feh in FEH_GRID:
                name = f'T{teff:04d}g{int(logg*10):02d}f{int(abs(feh)*10):02d}{"m" if feh<0 else "p"}'
                rows.append({'name': name, 'teff_K': teff, 'logg': logg, 'feh': feh})
    return rows


GRID_PARAMS = _build_grid_params()


# ── Match solar EW lines to MPIA wavelengths ──────────────────────────────────
def _get_matched_lines(element, mpia_opts):
    if not MEASURED_LINES.exists():
        raise SystemExit(f'CRITICAL: committed measured-line source missing: {MEASURED_LINES}')
    ew = pd.read_csv(MEASURED_LINES)
    our = ew[(ew['element'] == element) & (ew['ion'] == 'I')]['wavelength_air_A'].tolist()
    matched = {}
    for w in our:
        diffs = [(abs(w - m), m) for m in mpia_opts]
        best_diff, best_mpia = min(diffs)
        if best_diff <= MATCH_TOL_A and best_mpia not in matched:
            matched[best_mpia] = w
    return sorted(matched.keys())


# ── GET MPIA form to extract available wavelengths ────────────────────────────
def _get_mpia_options(element):
    sel_name = SELECT_NAME[element]
    r = requests.get(MPIA_URL, timeout=30)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, 'html.parser')
    sel = soup.find('select', {'name': sel_name})
    return [float(o.text.strip()) for o in sel.find_all('option')]


# ── Submit one batch and parse result ─────────────────────────────────────────
def _submit_batch(lines_angstrom, species_code, grid_params):
    """
    lines_angstrom : list of float MPIA wavelengths
    grid_params    : list of dicts {name, teff_K, logg, feh}
    Returns        : dict {star_name: {wave_A: delta_nlte}}
    """
    user_input = '\n'.join(
        f'{p["name"]} {p["teff_K"]} {p["logg"]:.2f} {p["feh"]:+.2f} {VMIC}'
        for p in grid_params
    )
    lines_input = '\n'.join(f'{w:.3f} {species_code}' for w in lines_angstrom)

    data = {
        'model': MODEL,
        'user_input': user_input,
        'linelist': LINELIST,
        'lines_input': lines_input,
    }
    resp = requests.post(MPIA_URL, data=data, timeout=180)
    resp.raise_for_status()

    soup = BeautifulSoup(resp.text, 'html.parser')
    result_div = soup.find('div', {'id': 'result'})
    if not result_div:
        raise RuntimeError('No <div id="result"> in MPIA response — form may have failed')

    bold_tags = [b.text.strip() for b in result_div.find_all('b')]
    if not bold_tags:
        raise RuntimeError('Empty result div')

    # First bold tag is the header: '#lines(A)   wav1   wav2   ...'
    header = bold_tags[0]
    wav_cols = [float(w) for w in header.split()[1:]]

    results = {}
    for tag in bold_tags[1:]:
        parts = tag.split()
        if len(parts) < 2:
            continue
        star_name = parts[0]
        deltas = {}
        for i, val_str in enumerate(parts[1:]):
            if i >= len(wav_cols):
                break
            try:
                val = float(val_str)
                # Error codes: 10, 20, 30 → NaN
                deltas[wav_cols[i]] = np.nan if val in ERROR_CODES else val
            except ValueError:
                deltas[wav_cols[i]] = np.nan
        results[star_name] = deltas

    return results


# ── Build full grid for one element ──────────────────────────────────────────
def build_grid(element, resume=False):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    raw_path  = RAW_DIR / f'MPIA_{element}_raw.csv'
    final_path = OUT_DIR / OUTFILES[element]

    if resume and raw_path.exists():
        print(f'Resuming from {raw_path}')
        raw_df = pd.read_csv(raw_path)
        done_waves = set(raw_df['wave_A'].unique())
    else:
        raw_df = pd.DataFrame()
        done_waves = set()

    print(f'\n=== {element} I NLTE grid (MPIA MAFAGS-OS) ===')

    # Get available MPIA wavelengths for this element
    print('Fetching MPIA wavelength options...')
    mpia_opts = _get_mpia_options(element)
    print(f'MPIA has {len(mpia_opts)} {element} I lines')

    # Match to our solar EW lines
    matched_waves = _get_matched_lines(element, mpia_opts)
    print(f'Matched {len(matched_waves)} lines from solar_ew.csv to MPIA options')
    if not matched_waves:
        print(f'CRITICAL STOP: No {element} lines matched MPIA options.')
        raise SystemExit(1)

    # Filter already-done waves when resuming
    todo_waves = [w for w in matched_waves if w not in done_waves]
    print(f'Lines to scrape: {len(todo_waves)} (skipping {len(done_waves)} done)')

    species_code = SPECIES[element]
    all_rows = list(raw_df.to_dict('records')) if not raw_df.empty else []

    # Process in batches of LINE_BATCH wavelengths × STAR_BATCH stellar params
    line_batches = [todo_waves[i:i+LINE_BATCH]
                    for i in range(0, len(todo_waves), LINE_BATCH)]
    star_batches = [GRID_PARAMS[i:i+STAR_BATCH]
                    for i in range(0, len(GRID_PARAMS), STAR_BATCH)]

    total_batches = len(line_batches) * len(star_batches)
    batch_n = 0
    first = True

    for lb_idx, line_batch in enumerate(line_batches):
        for sb_idx, star_batch in enumerate(star_batches):
            batch_n += 1
            print(f'  Batch {batch_n}/{total_batches}: '
                  f'{len(line_batch)} lines × {len(star_batch)} stars...')

            try:
                results = _submit_batch(line_batch, species_code, star_batch)
            except Exception as e:
                print(f'  ERROR in batch: {e}')
                print('  Saving checkpoint and stopping.')
                if all_rows:
                    pd.DataFrame(all_rows).to_csv(raw_path, index=False)
                    print(f'  Checkpoint saved to {raw_path}')
                raise SystemExit(1)

            for param in star_batch:
                name = param['name']
                if name not in results:
                    continue
                for wav, delta in results[name].items():
                    all_rows.append({
                        'element'   : element,
                        'wave_A'    : wav,
                        'teff_K'    : param['teff_K'],
                        'logg'      : param['logg'],
                        'feh'       : param['feh'],
                        'delta_nlte': delta,
                        'star_name' : name,
                    })

            # Save checkpoint after each batch
            pd.DataFrame(all_rows).to_csv(raw_path, index=False)

            if not first:
                time.sleep(SLEEP_S)
            first = False

    df = pd.DataFrame(all_rows)
    print(f'\n{element} raw grid: {len(df)} rows')

    # ── Placeholder-zero guard (RYA-413) ──────────────────────────────────────
    # MPIA OFFERS some lines in its dropdown but RETURNS a literal 0 for them (it serves
    # but does not NLTE-model the line; e.g. Ca 6166). A line identically 0 across ALL
    # nodes is NOT a physical correction — DROP it loudly, never persist a zero as if it
    # were a correction. (A real correction varies over the teff/logg/feh grid.)
    from pipeline.nlte_corrections import detect_placeholder_zero_lines
    dropped_placeholders = detect_placeholder_zero_lines(df)
    if dropped_placeholders:
        print(f'  [PLACEHOLDER_ZERO] dropping {len(dropped_placeholders)} MPIA-served-zero '
              f'line(s) (0.0 across all nodes, not NLTE-modelled): {dropped_placeholders}')
        df = df[~df['wave_A'].round(3).isin(dropped_placeholders)].copy()

    # ── Solar anchor verification ─────────────────────────────────────────────
    print('\nSolar anchor verification...')
    valid = df[df['delta_nlte'].notna()].copy()
    if valid.empty:
        print(f'CRITICAL STOP: {element} grid has no valid (non-NaN) corrections.')
        raise SystemExit(1)

    # Average delta over all lines at each grid point, then interpolate at solar
    grid_mean = (valid.groupby(['teff_K', 'logg', 'feh'])['delta_nlte']
                      .mean().reset_index())
    pts  = grid_mean[['teff_K', 'logg', 'feh']].values
    vals = grid_mean['delta_nlte'].values

    interp = LinearNDInterpolator(pts, vals)
    solar_delta = float(interp([SOLAR])[0])

    if np.isnan(solar_delta):
        print(f'CRITICAL STOP: {element} grid returns NaN at solar params '
              f'(Teff={SOLAR[0]}, logg={SOLAR[1]}, [Fe/H]={SOLAR[2]}). '
              f'Grid does not cover solar parameter space.')
        raise SystemExit(1)

    print(f'{element} solar anchor delta: {solar_delta:+.4f} dex '
          f'(Teff={SOLAR[0]}, logg={SOLAR[1]}, [Fe/H]={SOLAR[2]})')

    # ── Print coverage summary ────────────────────────────────────────────────
    print(f'\n{element} NLTE grid: {len(df)} rows, '
          f'Teff=[{df["teff_K"].min():.0f}, {df["teff_K"].max():.0f}] K, '
          f'logg=[{df["logg"].min():.1f}, {df["logg"].max():.1f}], '
          f'[Fe/H]=[{df["feh"].min():.1f}, {df["feh"].max():.1f}]')

    # ── Save final grid ───────────────────────────────────────────────────────
    out_cols = ['element', 'wave_A', 'teff_K', 'logg', 'feh', 'delta_nlte']
    df[out_cols].to_csv(final_path, index=False)
    print(f'Saved: {final_path}')

    # ── Provenance sidecar (RYA-396) — loader reads the bare CSV, so provenance
    # lives next to it as JSON (not inline comments). ──────────────────────────
    n_valid = int(df['delta_nlte'].notna().sum())
    prov = {
        'element': element, 'ion': 'I', 'species_code': SPECIES[element],
        'reference': REFS[element], 'model': MODEL,
        'source_url': MPIA_URL, 'mpia_select': SELECT_NAME[element],
        'convention': 'delta_nlte = A_NLTE - A_LTE (dex)',
        'line_source': str(MEASURED_LINES.relative_to(Path(__file__).resolve().parents[1])),
        'extraction_date': date.today().isoformat(),
        'grid_axes': {'teff_K': TEFF_GRID, 'logg': LOGG_GRID, 'feh': FEH_GRID},
        'vmic_kms': VMIC,
        'n_lines': int(df['wave_A'].nunique()),
        'n_rows': int(len(df)), 'n_valid_delta': n_valid,
        'solar_anchor': {'teff_K': SOLAR[0], 'logg': SOLAR[1], 'feh': SOLAR[2],
                         'delta_nlte': round(solar_delta, 4)},
        'dropped_placeholder_zero_lines': dropped_placeholders,   # RYA-413: MPIA-served-zero, not modelled
    }
    prov_path = final_path.with_suffix('.prov.json')
    prov_path.write_text(json.dumps(prov, indent=2))
    print(f'Saved provenance: {prov_path}')

    return solar_delta


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--element', required=True,
                        choices=['Ca', 'Ti', 'Cr', 'Mn', 'Co', 'Si', 'Mg'])
    parser.add_argument('--resume', action='store_true',
                        help='Skip wavelengths already in raw checkpoint')
    args = parser.parse_args()

    delta = build_grid(args.element, resume=args.resume)
    print(f'\nDone. {args.element} solar anchor: {delta:+.4f} dex')
