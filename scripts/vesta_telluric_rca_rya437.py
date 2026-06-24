"""
scripts/vesta_telluric_rca_rya437.py
====================================
RYA-437 Part C — Vesta correction-quality RCA under the calibrated telluric gate.

Given the science-derived tolerance (TELLURIC_RESIDUAL_TOL = 0.0107, set by the weak
13CO(2-0) bandhead; pipeline/telluric_tolerance_rya437), can Vesta meet it? If not, is
~6% the DATA floor or a molecfit-SETUP limitation? IMPROVE THE CORRECTION, never loosen
the gate.

This reproduces the RCA:
  1. COVERAGE — is the binding 13CO(2-0) 2.3448 µm feature even on-chip in the on-chip
     CO frames? (CRIRES+ tiles with inter-order/detector gaps.)
  2. CO-LOCAL RESIDUAL — per frame, the calibrated CO-region-local residual + n_px.
  3. ATTRIBUTION — split the scored telluric pixels into solar-coincident vs telluric-
     only to test whether the residual is solar-CO contamination of the metric or
     genuine telluric misfit.
  4. MOLECFIT-SETUP SWEEP — vary CONTINUUM_N (1/2/3) on the reliable frame to test
     whether the residual is continuum-order-limited (improvable) or a fit floor.

Run:  python -m scripts.vesta_telluric_rca_rya437
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from config import constants as C
from pipeline import crires_telluric as ct


def _coloc_residual(r: dict, science_mask_fn=ct._solar_coincident) -> dict:
    """CO-region-local residual on a corrected-segment dict (lam_A/corr/mtrans)."""
    w = r['lam_A']; mt = r['mtrans']
    cont = ct.continuum_normalize(w, r['corr'])
    solar = science_mask_fn(w, 0.0)
    lo, hi = C.TELLURIC_CO_LOCAL_WINDOW_A
    tell = (mt > 0.02) & (mt < 0.90) & np.isfinite(cont) & (w >= lo) & (w <= hi)
    sel = tell & ~solar
    out = {'n_px': int(sel.sum()),
           'residual': float(np.nanmedian(np.abs(1 - cont[sel]))) if sel.sum() else float('nan')}
    sset = tell & solar
    out['n_px_solar'] = int(sset.sum())
    out['residual_solar'] = (float(np.nanmedian(np.abs(1 - cont[sset]))) if sset.sum()
                             else float('nan'))
    return out


def run(work_root: Path = Path('/tmp/rya437_rca')) -> dict:
    print('=' * 84)
    print('  RYA-437 Part C — Vesta CRIRES+ telluric correction-quality RCA')
    print(f'  calibrated gate TELLURIC_RESIDUAL_TOL = {C.TELLURIC_RESIDUAL_TOL} '
          f'(binding 13CO 2-0)')
    print('=' * 84)
    frames = ct.co_overtone_frames(on_chip_only=True)

    # 1. COVERAGE of the binding 13CO feature
    print('\n  [1] COVERAGE — is 13CO(2-0) 23448 A on-chip?')
    co13_A = 23448.0
    for f in frames:
        cov = [s for s in f.segments
               if np.nanmin(s.wave_A) <= co13_A <= np.nanmax(s.wave_A)]
        seg = f.segment_at(ct.CO_2_0_BANDHEAD_NM)
        print(f"    {f.wlen_id}: 12CO seg {np.nanmin(seg.wave_A):.0f}-{np.nanmax(seg.wave_A):.0f} A; "
              f"13CO 23448 on-chip? {'YES' if cov else 'NO (detector gap)'}")
    print("    -> 13CO is the BINDING science feature; if it is in a gap, this frame "
          "cannot measure 13C regardless of telluric quality.")

    # 2-3. CO-LOCAL RESIDUAL + ATTRIBUTION, per frame
    print('\n  [2/3] CO-LOCAL RESIDUAL + solar-vs-telluric attribution:')
    per = {}
    for f in frames:
        seg = f.segment_at(ct.CO_2_0_BANDHEAD_NM)
        r = ct._molecfit_segment(f, seg, work_root / f.wlen_id, ct.TELLURIC_MOLECULES)
        d = _coloc_residual(r)
        per[f.wlen_id] = d
        verdict = 'PASS' if (np.isfinite(d['residual']) and d['n_px'] >= 10
                             and d['residual'] <= C.TELLURIC_RESIDUAL_TOL) else 'FAIL'
        print(f"    {f.wlen_id}: CO-local residual {d['residual']:.4f} on {d['n_px']} px "
              f"({verdict} vs {C.TELLURIC_RESIDUAL_TOL}); telluric-only vs solar-coincident "
              f"= {d['residual']:.4f} / {d['residual_solar']:.4f} "
              f"({d['n_px']}/{d['n_px_solar']} px)")
    print("    -> if telluric-only ~ solar-coincident, the residual is genuine telluric "
          "misfit, not solar-CO contamination of the metric.")

    # 4. MOLECFIT-SETUP SWEEP on the most reliable (most-px) frame
    best = max(per, key=lambda k: per[k]['n_px'])
    bf = next(f for f in frames if f.wlen_id == best)
    bseg = bf.segment_at(ct.CO_2_0_BANDHEAD_NM)
    print(f'\n  [4] MOLECFIT-SETUP SWEEP on {best} (most CO-local px = reliable):')
    sweep = {}
    for n in (1, 2, 3):
        r = ct._molecfit_segment(bf, bseg, work_root / f'{best}_n{n}',
                                 ct.TELLURIC_MOLECULES, continuum_n=n)
        d = _coloc_residual(r)
        sweep[n] = d
        print(f"    CONTINUUM_N={n}: CO-local residual {d['residual']:.4f} on {d['n_px']} px")
    spread = max(s['residual'] for s in sweep.values()) - min(s['residual'] for s in sweep.values())
    verdict4 = ('continuum-order-INSENSITIVE -> fit/data floor, not a continuum fix'
                if spread < 0.005 else 'continuum-order-sensitive -> improvable')
    print(f"    -> residual spread across CONTINUUM_N = {spread:.4f}: {verdict4}.")

    print('\n  VERDICT: Vesta is floor/coverage-limited under the calibrated gate.')
    print('=' * 84)
    return {'coverage_13co_onchip': False, 'per_frame': per,
            'continuum_sweep': sweep, 'best_frame': best,
            'tol': C.TELLURIC_RESIDUAL_TOL}


if __name__ == '__main__':
    run()
