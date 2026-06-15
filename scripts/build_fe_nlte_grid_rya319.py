#!/usr/bin/env python3
"""
scripts/build_fe_nlte_grid_rya319.py
====================================
RYA-319: build the Fe I + Fe II NLTE correction grid from MPIA Spectrum Tools
(MAFAGS-OS, 1D NLTE) over the full benchmark parameter box, so every Codex
benchmark/target resolves IN-GRID (no Amarsi clamp). Cleared the 2-node smoke
gate (Procyon Fe I Δ ≈ +0.02, not −0.3 → grid-edge confirmed).

Reuses the proven POST/parse (_submit_batch) from build_nlte_grids_mpia.py.
Scrapes both ions (Fe I=26.01 via felines[], Fe II=26.02 via fe2lines[]),
resumable per (ion, wave) checkpoint, and writes
data/nlte_grids/Fe_Bergemann_MPIA.csv with columns:
  element, ion, wave_A, teff_K, logg, feh, delta_nlte

Usage:  python scripts/build_fe_nlte_grid_rya319.py [--resume]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from config.constants import PATHS
from scripts.build_nlte_grids_mpia import _submit_batch, MPIA_URL, LINE_BATCH, STAR_BATCH, SLEEP_S

MATCH_TOL_A = 0.15
ION = {'I': ('felines[]', '26.01'), 'II': ('fe2lines[]', '26.02')}

# Benchmark box (all five targets strictly inside the convex hull):
#   Sun 5772/4.44/0.0 · Procyon 6554/4.00/+0.01 · αCenA 5792/4.30/+0.20
#   αCenB 5231/4.53/+0.20 · 55Cnc 5196/4.45/+0.32
TEFF_GRID = [5000, 5200, 5500, 5772, 6000, 6300, 6500, 6600]
LOGG_GRID = [4.0, 4.3, 4.44, 4.5, 4.6]
FEH_GRID  = [-0.5, -0.2, 0.0, 0.2, 0.35]

OUT_DIR = REPO / 'data' / 'nlte_grids'
RAW = OUT_DIR / 'raw' / 'Fe_MPIA_raw.csv'
FINAL = OUT_DIR / 'Fe_Bergemann_MPIA.csv'


def grid_params():
    rows = []
    for t in TEFF_GRID:
        for g in LOGG_GRID:
            for f in FEH_GRID:
                # MPIA rejects '+'/'-' in the starname ("Bad starname"); use the
                # m/p suffix convention from build_nlte_grids_mpia.py (no sign chars).
                tag = f'T{t:04d}g{int(g*100):03d}f{int(abs(f)*100):03d}{"m" if f < 0 else "p"}'
                rows.append({'name': tag, 'teff_K': t, 'logg': g, 'feh': f})
    return rows


def mpia_options(select_name):
    r = requests.get(MPIA_URL, timeout=60); r.raise_for_status()
    sel = BeautifulSoup(r.text, 'html.parser').find('select', {'name': select_name})
    return [float(o.text.strip()) for o in sel.find_all('option')]


def matched_waves(ion, opts):
    our = pd.read_csv(PATHS['solar_ew'])
    our = our[(our.element == 'Fe') & (our.ion == ion)]['wavelength_air_A'].tolist()
    out, used = [], set()
    for w in sorted(our):
        d, m = min((abs(w - o), o) for o in opts)
        if d <= MATCH_TOL_A and m not in used:
            used.add(m); out.append(m)
    return sorted(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()
    RAW.parent.mkdir(parents=True, exist_ok=True)

    params = grid_params()
    print(f"RYA-319 Fe NLTE grid build — {len(params)} nodes "
          f"(Teff {TEFF_GRID[0]}-{TEFF_GRID[-1]}, logg {LOGG_GRID[0]}-{LOGG_GRID[-1]}, "
          f"[Fe/H] {FEH_GRID[0]}..{FEH_GRID[-1]})", flush=True)

    rows = []
    done = set()
    if args.resume and RAW.exists():
        prev = pd.read_csv(RAW)
        rows = prev.to_dict('records')
        done = {(r['ion'], round(float(r['wave_A']), 3)) for r in rows}
        print(f"  resume: {len(done)} (ion,wave) already done", flush=True)

    star_batches = [params[i:i+STAR_BATCH] for i in range(0, len(params), STAR_BATCH)]
    t0 = time.time()
    for ion, (sel, code) in ION.items():
        opts = mpia_options(sel)
        waves = [w for w in matched_waves(ion, opts) if (ion, round(w, 3)) not in done]
        print(f"\nFe {ion}: {len(waves)} lines to scrape × {len(params)} nodes", flush=True)
        line_batches = [waves[i:i+LINE_BATCH] for i in range(0, len(waves), LINE_BATCH)]
        nb = len(line_batches) * len(star_batches); bi = 0
        for lb in line_batches:
            for sb in star_batches:
                bi += 1
                try:
                    res = _submit_batch(lb, code, sb)
                except Exception as e:
                    print(f"  Fe {ion} batch {bi}/{nb} ERROR {e} — checkpoint+continue", flush=True)
                    pd.DataFrame(rows).to_csv(RAW, index=False)
                    time.sleep(SLEEP_S * 3); continue
                for p in sb:
                    for w, d in res.get(p['name'], {}).items():
                        rows.append({'element': 'Fe', 'ion': ion, 'wave_A': w,
                                     'teff_K': p['teff_K'], 'logg': p['logg'],
                                     'feh': p['feh'], 'delta_nlte': d})
                pd.DataFrame(rows).to_csv(RAW, index=False)   # per-batch checkpoint
                if bi % 10 == 0:
                    print(f"  Fe {ion} {bi}/{nb}  ({time.time()-t0:.0f}s, {len(rows)} rows)", flush=True)
                time.sleep(SLEEP_S)

    df = pd.DataFrame(rows)
    df[['element', 'ion', 'wave_A', 'teff_K', 'logg', 'feh', 'delta_nlte']].to_csv(FINAL, index=False)
    valid = df[df.delta_nlte.notna()]
    print(f"\nDONE: {len(df)} rows ({len(valid)} valid), "
          f"Fe I={int((df.ion=='I').sum())} Fe II={int((df.ion=='II').sum())} "
          f"→ {FINAL.relative_to(REPO)}  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
