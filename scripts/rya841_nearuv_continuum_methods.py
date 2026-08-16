#!/usr/bin/env python3
"""
RYA-841 — the continuum-placement systematic, measured by disagreement
======================================================================
The companion script measures the LEVER (dA/d delta). This measures the OTHER half:
how far apart defensible continuum placements actually put the answer.

WHY THIS IS THE MEASUREMENT THAT SETTLES IT. The 0.100 dex term needs a sigma_delta --
how badly the near-UV continuum is placed -- and the repo has none. The one diagnostic
that exists, `cont_ratio` (synth/obs median ratio per window), does NOT measure it:

    across-line slope  dA/d(cont_ratio - 1) = +0.49 +/- 0.38   (n=40, r=+0.21)
    controlled lever   dA/d(delta)          = ~+1.7

Those disagree at about 3 sigma. If cont_ratio tracked a real continuum misplacement the
two would match, because they would be the same physical response. They do not, so
cont_ratio is dominated by MODEL error in the window (missing or wrong blends), not by
normalisation error -- and its 6.8% MAD must not be used as sigma_delta.

So sigma_delta is measured the way it can be: apply several defensible continuum
placements to the SAME line and let them disagree. The spread of that disagreement is
the systematic, and it needs no access to a "true" continuum that is never observed.

THE THREE PLACEMENTS
  atlas    the Kitt Peak atlas normalisation as delivered            (delta = 0)
  local    the window's own upper envelope, the way a human places a local continuum
           (delta = E - 1, E = high percentile of the flux over a wider neighbourhood)
  synth    put the continuum where the model says it is, i.e. force cont_ratio -> 1
           (delta = 1/cont_ratio - 1)

⚠️ `synth` IS PARTLY CIRCULAR AND IS REPORTED ANYWAY. Normalising to the model and then
fitting an abundance to the model rewards agreement. It is included because it is a
method people really use and because its DISTANCE from `atlas` is informative -- but it
is never the recommended placement, and the headline spread is quoted both with and
without it.

Every placement reduces to a per-line delta, so all three go through the SAME tested
`fit_at_delta` as the lever run. One code path, three inputs.

Usage (Sirius):
    python3 scripts/rya841_nearuv_continuum_methods.py
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from pipeline.nearuv_synth import build_solar_context, gf_provenance         # noqa: E402
from rya759_nearuv_fe_product import NEARUV_LINELIST, select_lines           # noqa: E402
from rya759_nearuv_synth import _kp_segments, _load_kp_window                # noqa: E402
from rya841_nearuv_continuum_sensitivity import fit_at_delta                 # noqa: E402

LO_A, HI_A = 3000.0, 3780.0
HALF_WIDTH_A = 0.40
MIN_SEP_A = 4.0
N_LINES = 40
RESOLVING_POWER = 500000.0

OUT = ROOT / 'data' / 'results' / 'rya841'
PER_LINE = OUT / 'rya841_continuum_methods_per_line.csv'
SUMMARY = OUT / 'rya841_continuum_methods.json'
POOL_832 = ROOT / 'data' / 'audit' / 'nearuv_fe_product' / 'nearuv_fe_per_line_main.csv'

#: How wide a neighbourhood the LOCAL placement looks at for its envelope. Wide enough to
#: contain some near-continuum pixels, narrow enough that the continuum slope does not
#: matter. In a band with a 0.146 A median line gap there are no truly line-free pixels,
#: which is exactly why this is a pseudo-continuum and why the method disagrees with the
#: others at all.
LOCAL_NEIGHBOURHOOD_A = 2.0
LOCAL_PERCENTILE = 95.0


def local_envelope_delta(segs, wave_A: float) -> float:
    """delta implied by placing the continuum at the window's own upper envelope."""
    try:
        w, f = _load_kp_window(segs, wave_A, LOCAL_NEIGHBOURHOOD_A)
    except LookupError:
        return float('nan')
    f = np.asarray(f, dtype=float)
    f = f[np.isfinite(f)]
    if f.size < 10:
        return float('nan')
    return float(np.percentile(f, LOCAL_PERCENTILE) - 1.0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--limit', type=int, default=N_LINES)
    ap.add_argument('--from-per-line', type=Path, default=None)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if a.from_per_line:
        report(pd.read_csv(a.from_per_line))
        return

    prov = gf_provenance(LO_A, HI_A)
    ctx = build_solar_context('Fe', RESOLVING_POWER,
                              linelist_file=str(NEARUV_LINELIST),
                              apply_canonical_gf=prov['apply_canonical_gf'])
    segs = _kp_segments()
    cand = select_lines(ctx['linelist'], lo_A=LO_A, hi_A=HI_A, n=a.limit,
                        teff=float(ctx['teff']), min_sep_A=MIN_SEP_A)
    print(f'[sel] {len(cand)} Fe I lines')

    # The atlas placement and its per-window cont_ratio already exist from RYA-759; the
    # `synth` delta is read from them rather than recomputed, which also keeps this run
    # anchored to the published product instead of a private re-derivation.
    ref = pd.read_csv(POOL_832)[['wave_A', 'a_synth', 'cont_ratio']].rename(
        columns={'a_synth': 'a_atlas'})

    tmp_dir = tempfile.mkdtemp(prefix='rya841m_')
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    rows = []
    for i, (_, ln) in enumerate(cand.iterrows(), 1):
        wv = float(ln['wave_A'])
        r = ref[np.isclose(ref.wave_A, wv, atol=0.005)]
        cont = float(r.cont_ratio.iloc[0]) if len(r) else np.nan
        a_atlas = float(r.a_atlas.iloc[0]) if len(r) else np.nan

        deltas = {
            'local': local_envelope_delta(segs, wv),
            'synth': (1.0 / cont - 1.0) if np.isfinite(cont) and cont > 0 else np.nan,
        }
        row = {'wave_A': wv, 'cont_ratio': cont, 'a_atlas': a_atlas,
               'delta_atlas': 0.0}
        for name, dl in deltas.items():
            row[f'delta_{name}'] = dl
            if not np.isfinite(dl):
                row[f'a_{name}'] = np.nan
                continue
            out = fit_at_delta(ctx, segs, wv, HALF_WIDTH_A, dl, tmp_dir)
            row[f'a_{name}'] = out['a_synth']
            # Recorded because a placement can be far outside the regime where a fit
            # means anything, and the abundance alone will not say so.
            row[f'chi2_{name}'] = out['red_chi2']
            row[f'status_{name}'] = out['status']
            print(f'[{i:3d}/{len(cand)}] {wv:9.3f}  {name:5s}  d={dl:+.4f}  '
                  f'A={out["a_synth"]:.3f}  chi2r={out["red_chi2"]:.1f}', flush=True)
        rows.append(row)
        pd.DataFrame(rows).to_csv(PER_LINE, index=False)

    report(pd.DataFrame(rows))


def report(d: pd.DataFrame) -> None:
    cols = {'atlas': 'a_atlas', 'local': 'a_local', 'synth': 'a_synth'}
    print('\n=== per-method medians ===')
    med = {}
    for name, c in cols.items():
        if c not in d.columns:
            continue
        v = d[c].dropna()
        if not len(v):
            continue
        med[name] = float(v.median())
        print(f'  {name:6s} n={len(v):3d}  median={v.median():.3f}  '
              f'scatter={v.std(ddof=1):.3f}  '
              f'median delta={d.get(f"delta_{name}", pd.Series([0.0])).median():+.4f}')

    print('\n=== pairwise disagreement (the systematic) ===')
    pairs = {}
    names = list(med)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            both = d[[cols[n1], cols[n2]]].dropna()
            if not len(both):
                continue
            diff = (both[cols[n1]] - both[cols[n2]])
            pairs[f'{n1}_vs_{n2}'] = {
                'n': int(len(both)),
                'median_diff_dex': float(diff.median()),
                'rms_diff_dex': float(np.sqrt((diff ** 2).mean())),
            }
            print(f'  {n1:6s} - {n2:6s}  n={len(both):3d}  '
                  f'median {diff.median():+.3f}  rms {np.sqrt((diff**2).mean()):.3f}')

    # The headline: how far apart the placements put the ANSWER. Quoted excluding the
    # partly-circular `synth` placement as well, so the number does not depend on it.
    honest = [med[k] for k in ('atlas', 'local') if k in med]
    spread_honest = (max(honest) - min(honest)) if len(honest) > 1 else float('nan')
    spread_all = (max(med.values()) - min(med.values())) if len(med) > 1 else float('nan')
    print(f'\n=== the continuum-placement systematic ===')
    print(f'  spread of the MEDIAN across placements')
    print(f'    atlas vs local (no circularity) : {spread_honest:.3f} dex')
    print(f'    including synth                 : {spread_all:.3f} dex')
    print(f'  assumed by the budget             : 0.100 dex')

    SUMMARY.write_text(json.dumps({
        'ticket': 'RYA-841',
        'measurement': 'continuum-placement systematic by method disagreement',
        'medians': med,
        'pairwise': pairs,
        'spread_atlas_vs_local_dex': spread_honest,
        'spread_all_methods_dex': spread_all,
        'assumed_term_dex': 0.10,
        'local_neighbourhood_A': LOCAL_NEIGHBOURHOOD_A,
        'local_percentile': LOCAL_PERCENTILE,
        'caveat': ('the `synth` placement normalises to the model and then fits the '
                   'model, so it is partly circular; the headline is quoted without it '
                   'as well'),
    }, indent=2))
    print(f'\n[out] {PER_LINE}\n[out] {SUMMARY}')


if __name__ == '__main__':
    main()
