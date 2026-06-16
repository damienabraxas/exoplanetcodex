#!/usr/bin/env python3
"""
scripts/smoke_fe_nlte_mpia_rya319.py
====================================
RYA-319 GATE: 2-node MPIA Fe NLTE smoke test BEFORE the full overnight scrape.

Scrapes ONLY solar (5772/4.44/0.00) + Procyon (6554/4.00/+0.01) for Fe I (26.01)
and Fe II (26.02), then reports the go/no-go checks (Ryan's gate):

  1. Procyon Fe I mean delta_nlte — GO if ~-0.05 (grid-edge confirmed → launch
     full scrape); NO-GO if ~-0.3 (grid was never the cause → fix the app path).
  2. Wavelength frame: MPIA λ vs our AIR line list (match offset distribution).
  3. Species 26.01 / 26.02 parse → Fe I / Fe II.
  4. Per-node line counts non-zero, ballpark.
  5. CSV schema sample (element, ion, wave_A, teff_K, logg, feh, delta_nlte).

Does NOT modify pipeline files or launch the full scrape.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from config.constants import PATHS
from scripts.build_nlte_grids_mpia import _submit_batch, MPIA_URL  # reuse proven POST/parse

N_SAMPLE = 40          # lines per ion for the smoke (speed)
MATCH_TOL_A = 0.15
NODES = [
    {'name': 'sun',     'teff_K': 5772, 'logg': 4.44, 'feh': 0.00},
    {'name': 'procyon', 'teff_K': 6554, 'logg': 4.00, 'feh': 0.01},
]
ION = {'I': ('felines[]', '26.01'), 'II': ('fe2lines[]', '26.02')}


def mpia_options(select_name):
    r = requests.get(MPIA_URL, timeout=60); r.raise_for_status()
    sel = BeautifulSoup(r.text, 'html.parser').find('select', {'name': select_name})
    return [float(o.text.strip()) for o in sel.find_all('option')]


def matched(our_waves, opts):
    """Return [(mpia_wave, our_wave, offset)] best-matched within tol."""
    out, used = [], set()
    for w in sorted(our_waves):
        d, m = min((abs(w - o), o) for o in opts)
        if d <= MATCH_TOL_A and m not in used:
            used.add(m); out.append((m, w, m - w))
    return out


def main():
    ew = pd.read_csv(PATHS['solar_ew'])
    print(f"{'='*70}\nRYA-319 GATE — 2-node MPIA Fe NLTE smoke (solar + Procyon)\n{'='*70}")
    rows = []
    for ion, (sel, code) in ION.items():
        opts = mpia_options(sel)
        our = ew[(ew.element == 'Fe') & (ew.ion == ion)]['wavelength_air_A'].tolist()
        m = matched(our, opts)
        offs = np.array([o for _, _, o in m]) if m else np.array([0.])
        print(f"\nFe {ion}: MPIA opts={len(opts)}  our lines={len(our)}  matched={len(m)}"
              f"  | match offset (MPIA−our, Å): median={np.median(offs):+.3f} "
              f"max|·|={np.abs(offs).max():.3f}  → {'AIR-frame (≈0)' if np.abs(offs).max()<0.1 else 'CHECK air/vac (offset large)'}")
        waves = [mw for mw, _, _ in m][:N_SAMPLE]
        res = _submit_batch(waves, code, NODES)   # {node: {wave: delta}}
        for node in NODES:
            d = res.get(node['name'], {})
            vals = np.array([v for v in d.values() if v is not None and np.isfinite(v)])
            print(f"   node {node['name']:8} Fe {ion}: n_delta={len(vals)}  "
                  f"mean Δ_nlte={vals.mean():+.4f}  median={np.median(vals):+.4f}" if len(vals)
                  else f"   node {node['name']:8} Fe {ion}: NO valid deltas")
            for w, v in d.items():
                rows.append({'element': 'Fe', 'ion': ion, 'wave_A': w,
                             'teff_K': node['teff_K'], 'logg': node['logg'],
                             'feh': node['feh'], 'delta_nlte': v, 'node': node['name']})
    df = pd.DataFrame(rows)
    out = REPO / 'data' / 'nlte_grids' / 'raw' / 'Fe_smoke_2node_rya319.csv'
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    # ── GATE verdict ──
    def mean_delta(ion, node):
        s = df[(df.ion == ion) & (df.node == node)]['delta_nlte'].dropna()
        return s.mean() if len(s) else np.nan
    pro_fe1 = mean_delta('I', 'procyon'); pro_fe2 = mean_delta('II', 'procyon')
    sun_fe1 = mean_delta('I', 'sun')
    print(f"\n{'='*70}\nGATE VERDICT")
    print(f"  Procyon Fe I  mean Δ_nlte = {pro_fe1:+.4f}   (GO if ~-0.05, NO-GO if ~-0.3)")
    print(f"  Procyon Fe II mean Δ_nlte = {pro_fe2:+.4f}   (expect ~0)")
    print(f"  Solar   Fe I  mean Δ_nlte = {sun_fe1:+.4f}   (tripwire ~-0.03)")
    go = np.isfinite(pro_fe1) and abs(pro_fe1) < 0.15
    print(f"  → {'GO: launch full scrape' if go else 'NO-GO: app-path RCA (grid not the cause)'}")
    print(f"  schema cols: {list(df.columns)}  → {out.relative_to(REPO)}")
    print('='*70)


if __name__ == '__main__':
    main()
