#!/usr/bin/env python3
"""
RYA-841 — what the near-UV pseudo-continuum systematic actually is
=================================================================
RYA-836 made the pseudo-continuum the DOMINANT near-UV systematic (0.147 total, of which
0.100 is this term). This measures it, because it has never been measured.

WHERE THE 0.100 CAME FROM. `pipeline/error_budget.py:170-176` adds it whenever a band's
policy says "pseudo-continuum", with this reasoning:

    "near-UV median flux 0.283-0.805 means the normalisation itself is uncertain at the
     ~10% level in the worst windows"   ->   Term("pseudo-continuum", 0.10, ...)

🔴 THAT STEP EQUATES ~10% IN FLUX WITH 0.10 IN DEX. They are different units, and the
conversion between them -- dA/d(continuum) -- appears NOWHERE in the repo. The budget
therefore assumes, without saying so, that dA/ddelta = 1.0 dex per unit fractional
continuum error.

That assumption is not plausible on its face. A continuum misplaced by delta changes the
depth of EVERY pixel in the fit window, not just the line core: for a weak line,
dEW ~ delta * W_window, so dEW/EW ~ delta * W_window/EW. With W = 0.8 A and EW ~ 50 mA
that is a 16x lever, i.e. dA/ddelta of order several dex per unit delta, not 1.

So this script measures the lever directly, by the only method that answers it -- a
controlled perturbation. For each near-UV line it re-fits the SAME window with the
observed flux renormalised as if the continuum had been placed wrong by a known delta:

    f_used = f_atlas / (1 + delta)          delta = fractional error in adopted continuum

delta > 0 means the continuum was placed too HIGH, which makes lines look deeper and must
drive A UP. The per-line slope dA/ddelta is then read off a symmetric grid, which also
tests linearity (a slope quoted from an asymmetric response would be meaningless).

WHAT THIS DELIBERATELY DOES NOT DO
  * It does not tune the continuum toward the optical or Asplund value. RYA-161. The
    perturbation is symmetric about the delivered atlas normalisation and the zero point
    is never moved.
  * It does not drop lines to shrink anything. RYA-777/RYA-161. Every line is reported.
  * It does not, by itself, license a new budget number. dA/ddelta is only HALF the
    product; the other half is sigma_delta, how badly the continuum is actually placed,
    which is a separate measurement. This script reports the lever and states the
    implication for a range of sigma_delta rather than picking one.

Usage (Sirius -- needs Turbospectrum through iSpec):
    python3 scripts/rya841_nearuv_continuum_sensitivity.py
    python3 scripts/rya841_nearuv_continuum_sensitivity.py --from-per-line <csv>   # re-report
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

from pipeline.nearuv_synth import build_solar_context, gf_provenance      # noqa: E402
from rya759_nearuv_fe_product import NEARUV_LINELIST, select_lines        # noqa: E402
from rya759_nearuv_synth import _kp_segments, _load_kp_window             # noqa: E402

LO_A, HI_A = 3000.0, 3780.0
HALF_WIDTH_A = 0.40
MIN_SEP_A = 4.0
N_LINES = 40
RESOLVING_POWER = 500000.0

OUT = ROOT / 'data' / 'results' / 'rya841'
PER_LINE = OUT / 'rya841_continuum_sensitivity_per_line.csv'
SUMMARY = OUT / 'rya841_summary.json'

#: Symmetric on purpose: the symmetry is what lets non-linearity be SEEN rather than
#: assumed away. +/-4% brackets the "~10% in the worst windows" claim the 0.100 rests on
#: while staying small enough that a linear slope is meaningful.
DELTAS = (-0.04, -0.02, 0.0, +0.02, +0.04)

#: The value under test. NOT used in any fit -- only in the comparison at the end.
ASSUMED_DEX = 0.10


def fit_at_delta(ctx: dict, segs, wave_A: float, hw_A: float, delta: float,
                 tmp_dir: str) -> dict:
    """Flux-fit A(Fe) in one window with the continuum deliberately misplaced by `delta`.

    The ONLY difference from RYA-759's `fit_one` is the renormalisation of the observed
    flux. Everything else -- atmosphere, linelist, window, bounds, broadening -- is held
    fixed, which is what makes the difference attributable to the continuum alone.
    """
    from pipeline.abundances_derive import _fit_synth_flux

    lo_A, hi_A = wave_A - hw_A, wave_A + hw_A
    try:
        ow_A, of = _load_kp_window(segs, wave_A, hw_A + 0.4)
    except LookupError as e:
        return {'status': 'no_atlas', 'reason': str(e), 'a_synth': np.nan,
                'red_chi2': np.nan}

    of = np.asarray(of, dtype=float) / (1.0 + delta)

    a_solar = float(ctx['solar_A'])
    r = _fit_synth_flux(
        ow_A / 10.0, of, ctx['atmosphere'],
        ctx['teff'], ctx['logg'], ctx['feh'], ctx['vturb'],
        ctx['linelist'], ctx['isotopes'], ctx['solar_abund'], 'Fe',
        int(ctx['atom_code']), lo_A / 10.0, hi_A / 10.0,
        max(a_solar - 3.0, 1.0), a_solar + 5.0,
        float(ctx['resolving_power']), float(ctx['macroturbulence']),
        float(ctx['vsini']), tmp_dir=tmp_dir)

    return {'status': r['status'], 'reason': r.get('reason', ''),
            'a_synth': float(r.get('A_X', np.nan)),
            'red_chi2': float(r.get('red_chi2', np.nan)),
            'n_pix': int(r.get('n_pix', 0))}


def slope_for_line(sub: pd.DataFrame) -> dict:
    """dA/ddelta for one line, with the linearity check that makes it quotable."""
    ok = sub[np.isfinite(sub.a_synth)].sort_values('delta')
    if len(ok) < 3:
        return {'slope': np.nan, 'r2': np.nan, 'asym': np.nan, 'n_delta': len(ok)}

    x, y = ok.delta.to_numpy(float), ok.a_synth.to_numpy(float)
    slope, icept = np.polyfit(x, y, 1)
    pred = slope * x + icept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # Asymmetry: if the +delta and -delta arms disagree in magnitude, a single slope is
    # a summary of a curve, and the budget term it implies is direction-dependent.
    asym = np.nan
    zero = ok[ok.delta == 0.0]
    if len(zero) == 1:
        a0 = float(zero.a_synth.iloc[0])
        up = ok[ok.delta > 0]
        dn = ok[ok.delta < 0]
        if len(up) and len(dn):
            asym = abs((float(up.a_synth.iloc[-1]) - a0)
                       + (float(dn.a_synth.iloc[0]) - a0))

    return {'slope': float(slope), 'r2': float(r2), 'asym': float(asym),
            'n_delta': int(len(ok))}


#: RYA-759's unperturbed run, used as a positive control on this harness.
RYA759_PER_LINE = (ROOT / 'data' / 'audit' / 'nearuv_fe_product'
                   / 'nearuv_fe_per_line_main.csv')


def control_against_759(per_line: pd.DataFrame) -> dict:
    """The delta=0 column must reproduce RYA-759, and it tests TWO things at once.

    1. THE HARNESS. If delta=0 does not reproduce 759 line for line, the perturbation is
       not the only thing that changed and no slope from it means anything. A sensitivity
       measured with a drifted baseline would look like a result.
    2. RYA-822's blueward extension. 759 ran with `apply_canonical_gf=False` because
       `canonical_gf.csv` then started at 3780 A. RYA-822 (`7a84c77`) extended it to
       3000.003 A, so THIS run applies canonical gf where 759 did not — 759's recorded
       provenance string is now stale. If the values still agree, 822's adjudication left
       the near-UV gf where VALD had it, which is 822's own "grading buys nothing"
       finding showing up as a null in a different place.

    A disagreement is not a failure to hide -- it is reported, because it would mean the
    near-UV product of record moved when the gf table was extended under it.
    """
    if not RYA759_PER_LINE.exists():
        return {'available': False}

    ref = pd.read_csv(RYA759_PER_LINE)[['wave_A', 'a_synth']].rename(
        columns={'a_synth': 'a_759'})
    mine = per_line[per_line.delta == 0.0][['wave_A', 'a_synth']].rename(
        columns={'a_synth': 'a_841'})
    m = ref.merge(mine, on='wave_A', how='inner')
    m = m[np.isfinite(m.a_759) & np.isfinite(m.a_841)]
    if not len(m):
        return {'available': False}

    d = (m.a_841 - m.a_759).abs()
    out = {'available': True, 'n_matched': int(len(m)),
           'max_abs_diff_dex': float(d.max()),
           'median_abs_diff_dex': float(d.median()),
           'n_differing_gt_0p01': int((d > 0.01).sum())}
    print(f'\n=== control: delta=0 vs RYA-759 (n={len(m)}) ===')
    print(f'  max |diff| = {d.max():.4f} dex, median = {d.median():.4f}, '
          f'{int((d > 0.01).sum())} lines differ by >0.01')
    print('  (759 ran pre-RYA-822 without canonical gf; this run applies it)')
    if d.max() > 0.01:
        print('  ⚠️  BASELINE MOVED — the near-UV product of record shifted when RYA-822 '
              'extended canonical_gf under it. Report this, do not absorb it.')
    return out


def report(per_line: pd.DataFrame) -> dict:
    rows = []
    for wv, sub in per_line.groupby('wave_A'):
        s = slope_for_line(sub)
        s['wave_A'] = float(wv)
        base = sub[sub.delta == 0.0]
        s['a_unperturbed'] = float(base.a_synth.iloc[0]) if len(base) else np.nan
        rows.append(s)
    sl = pd.DataFrame(rows).sort_values('wave_A')

    good = sl[np.isfinite(sl.slope)]
    med = float(np.median(np.abs(good.slope))) if len(good) else np.nan
    p16, p84 = (float(np.percentile(np.abs(good.slope), 16)),
                float(np.percentile(np.abs(good.slope), 84))) if len(good) else (np.nan,) * 2

    print('\n=== dA/ddelta, per line ===')
    for _, r in sl.iterrows():
        lin = '' if (np.isfinite(r.r2) and r.r2 > 0.98) else '   [NON-LINEAR]'
        print(f'  {r.wave_A:9.3f}  A0={r.a_unperturbed:6.3f}  '
              f'dA/dd={r.slope:+7.2f}  r2={r.r2:5.3f}{lin}')

    print(f'\n=== the lever ===')
    print(f'  median |dA/ddelta| = {med:.2f} dex per unit fractional continuum error')
    print(f'  16-84 pct          = {p16:.2f} - {p84:.2f}')
    print(f'  ASSUMED by the budget (implicitly) = 1.00\n')

    print('=== what the term becomes, for a range of sigma_delta ===')
    implied = {}
    for sd in (0.005, 0.01, 0.02, 0.05, 0.10):
        implied[f'{sd:.3f}'] = round(med * sd, 4)
        print(f'  sigma_delta = {sd*100:4.1f}%   ->   {med*sd:.3f} dex'
              f'{"     <-- the assumed 0.100 sits here" if abs(med*sd - ASSUMED_DEX) < 0.02 else ""}')

    nonlin = int((good.r2 <= 0.98).sum())
    print(f'\n  non-linear lines (r2<=0.98): {nonlin} of {len(good)}')

    ctrl = control_against_759(per_line)

    return {
        'control_vs_rya759': ctrl,
        'ticket': 'RYA-841',
        'measurement': 'dA/d(fractional continuum error), controlled perturbation',
        'band_A': [LO_A, HI_A],
        'deltas': list(DELTAS),
        'n_lines': int(len(sl)),
        'n_lines_with_slope': int(len(good)),
        'lever_dex_per_unit_delta': {
            'median_abs': med, 'p16': p16, 'p84': p84,
            'assumed_by_budget': 1.0,
        },
        'implied_term_by_sigma_delta': implied,
        'n_nonlinear': nonlin,
        'assumed_term_dex': ASSUMED_DEX,
        'note': (
            'The 0.100 dex pseudo-continuum term was never measured: error_budget.py '
            'converts "~10% uncertain normalisation" into "0.10 dex" one-to-one, which '
            'silently assumes dA/ddelta = 1. This is that derivative, measured. It does '
            'NOT by itself set a new budget term -- sigma_delta is a separate '
            'measurement and is owed.'),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--limit', type=int, default=N_LINES)
    ap.add_argument('--half-width-A', type=float, default=HALF_WIDTH_A)
    ap.add_argument('--min-sep-A', type=float, default=MIN_SEP_A)
    ap.add_argument('--from-per-line', type=Path, default=None,
                    help='re-report from an existing per-line CSV (no fitting)')
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)

    if a.from_per_line:
        summary = report(pd.read_csv(a.from_per_line))
        SUMMARY.write_text(json.dumps(summary, indent=2))
        print(f'\n[out] {SUMMARY}')
        return

    # gf single-sourcing is NOT available below 3780 A (canonical_gf.csv starts there),
    # so the flag comes from the provenance helper rather than the default (which raises).
    prov = gf_provenance(LO_A, HI_A)
    ctx = build_solar_context('Fe', RESOLVING_POWER,
                              linelist_file=str(NEARUV_LINELIST),
                              apply_canonical_gf=prov['apply_canonical_gf'])
    segs = _kp_segments()
    print(f'[kp]  {len(segs)} atlas segments')
    print(f'[gf]  {prov["detail"]}')

    cand = select_lines(ctx['linelist'], lo_A=LO_A, hi_A=HI_A, n=a.limit,
                        teff=ctx['teff'], min_sep_A=a.min_sep_A)
    print(f'[sel] {len(cand)} Fe I lines, {cand["wave_A"].min():.1f}-'
          f'{cand["wave_A"].max():.1f} A')
    print(f'[pert] deltas = {DELTAS}  ({len(cand)*len(DELTAS)} fits)\n')

    # RYA-824: the fitter writes here and does NOT mkdir. A missing dir kills every fit
    # with "No such file or directory" and an EMPTY table is the only symptom.
    tmp_dir = tempfile.mkdtemp(prefix='rya841_')
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    rows = []
    for i, (_, ln) in enumerate(cand.iterrows(), 1):
        wv = float(ln['wave_A'])
        for delta in DELTAS:
            out = fit_at_delta(ctx, segs, wv, a.half_width_A, delta, tmp_dir)
            rows.append({'wave_A': wv, 'delta': delta,
                         'loggf': float(ln['loggf']), 'ep_eV': float(ln['ep_eV']),
                         'theo_depth': float(ln['theo_depth']), **out})
            print(f'[{i:3d}/{len(cand)}] {wv:9.3f}  d={delta:+.3f}  '
                  f'A={out["a_synth"]:.3f}  chi2r={out["red_chi2"]:.1f}', flush=True)
        # Written every line, not at the end: RYA-836's report crashed after all the
        # expensive fits and only survived because the CSV was already on disk.
        pd.DataFrame(rows).to_csv(PER_LINE, index=False)

    summary = report(pd.DataFrame(rows))
    SUMMARY.write_text(json.dumps(summary, indent=2))
    print(f'\n[out] {PER_LINE}\n[out] {SUMMARY}')


if __name__ == '__main__':
    main()
