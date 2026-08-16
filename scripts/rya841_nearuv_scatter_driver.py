#!/usr/bin/env python3
"""
RYA-841 — what actually drives the near-UV line-to-line scatter
===============================================================
RYA-836 measured the split cleanly and named it only loosely:

    RYA-832 pool (40 lines, depth-selected, min_sep 4 A) ... scatter 0.413
    RYA-836 pool (60 lines, whatever the labs measured) ... scatter 0.651
    -> "line selection accounts for 0.238 dex of scatter, the gf for 0.002"

"Line selection" is a description of the difference, not a cause. This names the cause,
and it needs NO new synthesis -- both pools already have per-line abundances on disk.

🔴 THE SELECTOR DOES NOT DO WHAT ITS NAME SUGGESTS. `select_lines(min_sep_A=4.0)` applies
the separation test only BETWEEN LINES IT HAS ALREADY SELECTED (rya759_nearuv_fe_product.py
:127-134). Its own docstring says so -- "keeps the set spread across the band" -- so it is
a SAMPLING rule, not a blend rejection. A selected line may sit directly on top of a strong
neighbour that was never a candidate. In a band whose median line gap is 0.146 A, that is
the normal case, not the exception.

So the physical quantity is not separation between chosen lines. It is how much of the
absorption in a line's own fit window actually belongs to that line:

    blend_depth  = sum of theoretical depth of OTHER lines within the +/-0.40 A fit window
    depth_frac   = own depth / (own depth + blend_depth)

depth_frac = 1 means the window is the line. depth_frac = 0.1 means the fit is mostly
measuring somebody else, and the abundance it returns for Fe is a residual.

WHAT THIS DOES AND DOES NOT LICENCE
  * It reports the correlation for BOTH pools and every line stays in the record (RYA-711).
  * A subset defined by depth_frac has a stated PHYSICAL cause -- "this window is not
    dominated by the target line" -- which is what RYA-161 requires before any line is set
    aside. It is NOT an outlier cut, and it is not chosen by looking at which threshold
    minimises the scatter.
  * It does NOT install a tightened product. Whether the near-UV number should be rebuilt
    on a dominance criterion is a decision, not a measurement, and it is Ryan's.

Usage (Sirius -- needs the near-UV linelist, no synthesis):
    python3 scripts/rya841_nearuv_scatter_driver.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from pipeline.nearuv_synth import build_solar_context, gf_provenance      # noqa: E402
from rya759_nearuv_fe_product import NEARUV_LINELIST, select_lines        # noqa: E402

LO_A, HI_A = 3000.0, 3780.0
HALF_WIDTH_A = 0.40
MIN_SEP_A = 4.0
N_LINES = 40
RESOLVING_POWER = 500000.0

OUT = ROOT / 'data' / 'results' / 'rya841'
POOL_832 = ROOT / 'data' / 'audit' / 'nearuv_fe_product' / 'nearuv_fe_per_line_main.csv'
POOL_836 = ROOT / 'data' / 'results' / 'rya836' / 'rya836_nearuv_lab_gf_per_line.csv'
PER_LINE = OUT / 'rya841_scatter_driver_per_line.csv'
SUMMARY = OUT / 'rya841_scatter_driver.json'

#: Match tolerance when tying a pool line back to its linelist row. Well below the 0.146 A
#: median gap, so it cannot silently pair a line with its neighbour (the RYA-785 failure).
MATCH_TOL_A = 0.005

#: Fixed BEFORE looking at what it selects, and reported as a sweep either way so its
#: leverage is visible rather than asserted. "The target line supplies at least half the
#: window's absorption" is the weakest statement of dominance that means anything.
DOMINANCE_FLOOR = 0.50


def window_blending(linelist, pool_waves: np.ndarray, hw_A: float) -> pd.DataFrame:
    """For each pool line, how much of its fit window's absorption is somebody else's.

    Uses EVERY species in the list, not just Fe I: a window contaminated by a strong
    Ni or Cr line is contaminated regardless of what the blend is made of.
    """
    names = linelist.dtype.names
    w_A = np.asarray(linelist['wave_A'] if 'wave_A' in names
                     else linelist['wave_nm'] * 10.0, dtype=float)
    depth = np.asarray(linelist['theoretical_depth'], dtype=float)
    el = np.asarray([str(x).strip() for x in linelist['element']])
    finite = np.isfinite(w_A) & np.isfinite(depth)
    w_A, depth, el = w_A[finite], depth[finite], el[finite]
    order = np.argsort(w_A)
    w_A, depth, el = w_A[order], depth[order], el[order]

    rows = []
    for wv in pool_waves:
        lo, hi = wv - hw_A, wv + hw_A
        i0, i1 = np.searchsorted(w_A, [lo, hi])
        win_w, win_d, win_e = w_A[i0:i1], depth[i0:i1], el[i0:i1]

        # The target row is the closest line in the window; everything else is blend.
        if len(win_w) == 0:
            rows.append({'wave_A': wv, 'own_depth': np.nan, 'blend_depth': np.nan,
                         'depth_frac': np.nan, 'n_in_window': 0, 'nn_sep_A': np.nan,
                         'nn_depth': np.nan})
            continue
        k = int(np.argmin(np.abs(win_w - wv)))
        own = float(win_d[k])
        others = np.delete(win_d, k)
        others_w = np.delete(win_w, k)
        blend = float(np.nansum(others)) if len(others) else 0.0
        tot = own + blend
        if len(others_w):
            j = int(np.argmin(np.abs(others_w - wv)))
            nn_sep, nn_depth = float(abs(others_w[j] - wv)), float(others[j])
        else:
            nn_sep, nn_depth = np.nan, np.nan

        rows.append({'wave_A': wv, 'own_depth': own, 'blend_depth': blend,
                     'depth_frac': own / tot if tot > 0 else np.nan,
                     'n_in_window': int(len(win_w) - 1),
                     'nn_sep_A': nn_sep, 'nn_depth': nn_depth})
    return pd.DataFrame(rows)


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    """Rank correlation without a scipy dependency; n reported so it can be judged."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return np.nan, int(m.sum())
    # .to_numpy() returns a READ-ONLY view under pandas 3, so the in-place centring below
    # raises rather than silently doing nothing. Copy explicitly.
    rx = np.array(pd.Series(x[m]).rank(), dtype=float)
    ry = np.array(pd.Series(y[m]).rank(), dtype=float)
    rx -= rx.mean()
    ry -= ry.mean()
    den = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return (float((rx * ry).sum() / den) if den > 0 else np.nan), int(m.sum())


def analyse(name: str, df: pd.DataFrame, blend: pd.DataFrame) -> dict:
    """Per-pool: does |A - median| track window dominance?"""
    d = df.merge(blend, on='wave_A', how='left')
    d = d[np.isfinite(d.a_synth)].copy()
    med = float(d.a_synth.median())
    d['resid'] = (d.a_synth - med).abs()

    out = {'pool': name, 'n': int(len(d)), 'median_A': med,
           'scatter': float(d.a_synth.std(ddof=1))}

    print(f'\n=== {name} — n={len(d)}, median {med:.3f}, scatter '
          f'{d.a_synth.std(ddof=1):.3f} ===')
    for col in ('depth_frac', 'blend_depth', 'n_in_window', 'nn_sep_A', 'own_depth'):
        r, n = spearman(d[col].to_numpy(float), d.resid.to_numpy(float))
        out[f'spearman_resid_vs_{col}'] = r
        flag = '  <-- ' if (np.isfinite(r) and abs(r) >= 0.35) else ''
        print(f'  |A-med| vs {col:<13} rho = {r:+.3f}  (n={n}){flag}')

    # Dominance sweep. Reported at several thresholds so no single one is load-bearing.
    print(f'  dominance sweep (depth_frac >= t):')
    sweep = []
    for t in (0.0, 0.25, 0.40, 0.50, 0.60, 0.75):
        sub = d[d.depth_frac >= t]
        if len(sub) < 5:
            sweep.append({'t': t, 'n': int(len(sub)), 'scatter': None, 'median': None})
            print(f'    t={t:.2f}  n={len(sub):3d}  (too few to quote)')
            continue
        sc, mn = float(sub.a_synth.std(ddof=1)), float(sub.a_synth.median())
        sweep.append({'t': t, 'n': int(len(sub)), 'scatter': sc, 'median': mn})
        print(f'    t={t:.2f}  n={len(sub):3d}  scatter={sc:.3f}  median={mn:.3f}')
    out['dominance_sweep'] = sweep
    return out, d


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    prov = gf_provenance(LO_A, HI_A)
    ctx = build_solar_context('Fe', RESOLVING_POWER,
                              linelist_file=str(NEARUV_LINELIST),
                              apply_canonical_gf=prov['apply_canonical_gf'])
    linelist = ctx['linelist']
    print(f'[gf]  {prov["detail"]}')

    # --- the selector's own claim, checked rather than believed ----------------------
    sel = select_lines(linelist, lo_A=LO_A, hi_A=HI_A, n=N_LINES,
                       teff=float(ctx['teff']), min_sep_A=MIN_SEP_A)
    sel_blend = window_blending(linelist, sel.wave_A.to_numpy(float), HALF_WIDTH_A)
    print(f'\n=== is min_sep_A=4.0 a blend rejection? ===')
    print(f'  selected lines: {len(sel)}')
    print(f'  min separation BETWEEN SELECTED lines : '
          f'{np.diff(np.sort(sel.wave_A.to_numpy(float))).min():.3f} A')
    print(f'  median distance to nearest ANY-species neighbour in the FULL list: '
          f'{sel_blend.nn_sep_A.median():.4f} A')
    print(f'  median other-line count inside the +/-0.40 A fit window: '
          f'{sel_blend.n_in_window.median():.0f}')
    print(f'  => the 4.0 A rule spreads the SAMPLE; it does not clear the WINDOW.')

    pools = {}
    frames = []

    if POOL_832.exists():
        p = pd.read_csv(POOL_832)[['wave_A', 'a_synth', 'red_chi2', 'cont_ratio']]
        b = window_blending(linelist, p.wave_A.to_numpy(float), HALF_WIDTH_A)
        res, d = analyse('RYA-832 pool (depth-selected, n=40)', p, b)
        d['pool'] = 'rya832'
        pools['rya832'] = res
        frames.append(d)

    if POOL_836.exists():
        raw = pd.read_csv(POOL_836)
        p = raw.rename(columns={'wavelength_air_A': 'wave_A',
                                'a_control': 'a_synth'})[['wave_A', 'a_synth']]
        b = window_blending(linelist, p.wave_A.to_numpy(float), HALF_WIDTH_A)
        res, d = analyse('RYA-836 pool (lab-covered, n=60, production gf)', p, b)
        d['pool'] = 'rya836'
        pools['rya836'] = res
        frames.append(d)

    # --- RYA-836's owed item 2: lab lines that ALSO pass the 832 selection -----------
    overlap = {}
    if 'rya832' in pools and 'rya836' in pools:
        lab = frames[-1]
        sel_w = sel.wave_A.to_numpy(float)
        keep = [bool(np.any(np.abs(sel_w - w) <= MATCH_TOL_A)) for w in lab.wave_A]
        both = lab[np.array(keep)]
        print(f'\n=== RYA-836 owed item 2: lab lines that ALSO pass 832 selection ===')
        print(f'  {len(both)} of {len(lab)} lab lines are in the 832 selection')
        overlap = {'n': int(len(both))}
        if len(both) >= 5:
            overlap.update({'scatter': float(both.a_synth.std(ddof=1)),
                            'median': float(both.a_synth.median())})
            print(f'  scatter {both.a_synth.std(ddof=1):.3f}, '
                  f'median {both.a_synth.median():.3f}')
        else:
            print(f'  too few to quote a scatter — which is itself the answer 836 '
                  f'predicted: the labs and the clean-line criterion barely overlap.')

    allrows = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    allrows.to_csv(PER_LINE, index=False)

    summary = {
        'ticket': 'RYA-841',
        'measurement': 'what drives the near-UV line-to-line scatter',
        'band_A': [LO_A, HI_A],
        'fit_half_width_A': HALF_WIDTH_A,
        'selector_is_a_sampling_rule_not_a_blend_rejection': {
            'min_sep_between_selected_A': float(
                np.diff(np.sort(sel.wave_A.to_numpy(float))).min()),
            'median_nn_sep_full_list_A': float(sel_blend.nn_sep_A.median()),
            'median_other_lines_in_window': float(sel_blend.n_in_window.median()),
            'note': ('select_lines applies min_sep only between lines it has already '
                     'selected, so it spreads the sample across the band and does not '
                     'clear the fit window of blends.'),
        },
        'pools': pools,
        'lab_and_selected_overlap': overlap,
        'dominance_floor_declared': DOMINANCE_FLOOR,
        'note': ('depth_frac = own theoretical depth / total theoretical depth in the fit '
                 'window. A subset cut on it has a physical cause; it is reported as a '
                 'sweep and NOT installed as a product.'),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2))
    print(f'\n[out] {PER_LINE}\n[out] {SUMMARY}')


if __name__ == '__main__':
    main()
