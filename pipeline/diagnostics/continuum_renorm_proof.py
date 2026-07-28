#!/usr/bin/env python3
"""
Bounded per-order VIS continuum corrector + gated [O I] measured-O proof (RYA-451 production step).

RYA-452: the VIS continuum offset is NORMALIZATION-scale, uniform across true-continuum
VIS (~-0.007), no drift. Part B's NARROW-window renorm overshot ([O I] -> 8.598). This
builds the correction ROBUSTLY (many anchors, low-DOF, bounded) and tests whether it
lands [O I] in band WITHOUT regressing a clean control (C I).

GATE: proceed to production wiring (Part B) ONLY if A(O) in (8.64, 8.76) AND C I stays
8.46 +/- 0.03. If [O I] still overshoots -> the residual is deeper than narrow-window
noise (synth continuum scale / 1D), STOP. If C I moves -> regression, STOP.

Validate-don't-tune: C(lambda) is the measured robust synth/obs continuum ratio, low-order
and bounded; it is NOT fit to make A(O)=8.69.

RYA-453 (bounded continuum renorm + gated proof; RYA-451 production step).
"""
from __future__ import annotations
import argparse, os
import numpy as np
import pipeline.cno_synthesis as cs
from pipeline.cno_synthesis import (
    REGIONS, VIS_DIAGNOSTICS, get_star_params,
    _load_atmosphere, _load_synth_resources, _load_observed_spectrum,
    _atom_codes, _fixed_ab, _synth_window, _fit_element, preflight,
)
import ispec

SOLAR_COMP = dict(C=8.491, N=7.558, O=8.69, Ni=6.20)
# true-continuum VIS anchors only (blue <485 excluded per RYA-452); dense sampling for robustness
ANCHORS_NM = [(c - 0.6, c + 0.6) for c in np.arange(488, 666, 6.0)]
COARSE_NM  = 0.001
CONT_Q     = 0.70          # synth pixels above this quantile = (noiseless) continuum pixels (RYA-452)
MAX_CORR   = 0.02          # physical bound: |C-1| beyond this => fail loud
DEGREE     = 1             # bounded low-DOF; RYA-452 found no drift, so degree<=1 (constant/linear)
OI_BAND    = (8.64, 8.76)  # measured-O success band (Asplund 8.69 / Caffau 8.73)
CI_REF, CI_TOL = 8.46, 0.03


def _robust_cont(ow, of, sw, sf):
    """Noise-robust (obs_c, syn_c) on equal footing -- the merged RYA-452 estimator
    (continuum_offset_audit._robust_cont), replicated inline so this diagnostic is
    self-contained on the rya-371 base (which predates the 452 merge). Pick continuum
    pixels from the NOISELESS synth (top CONT_Q), then take the MEDIAN of obs and of
    synth at those same pixels (symmetric noise averages out) -- NOT a raw quantile of
    the noisy obs, which biases high. Integration note (RYA-453): use the validated 452
    method, not the brief's top-decile single-array proxy.
    """
    sf_i = np.interp(ow, sw, sf)
    cont = sf_i >= np.quantile(sf, CONT_Q)
    if cont.sum() < 5:
        return np.nan, np.nan
    return float(np.median(of[cont])), float(np.median(sf_i[cont]))


def _setup(star_id, region_name):
    region = REGIONS[region_name]
    rec = get_star_params(star_id)
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    broadening = preflight(region, star_id, VIS_DIAGNOSTICS)
    atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(cs._ISPEC_SOLAR_ABUND_FILE)
    obs_w, obs_f = _load_observed_spectrum(star_id)
    codes = _atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    return dict(params=params, broadening=broadening, atm=atm, ll=ll, iso=iso, sab=sab,
                obs_w=obs_w, obs_f=obs_f, codes=codes)


def _synth(ctx, lo, hi, tmp):
    sw = np.arange(lo, hi + COARSE_NM*0.5, COARSE_NM)
    fa = _fixed_ab(dict(SOLAR_COMP), ctx['codes'])
    return sw, _synth_window(sw, ctx['atm'], ctx['params'], ctx['ll'], ctx['iso'],
                             ctx['sab'], fa, ctx['broadening'], True, tmp)


def build_C(ctx, tmp, degree=DEGREE):
    """Robust, bounded, low-order continuum correction C(lambda) = synth_cont/obs_cont,
    built on the RYA-452 noise-robust (synth-defined-pixel) per-anchor estimator.

    degree: 1 (brief default, constant+linear) or 0 (constant). RYA-452 found NO
    wavelength drift in true-continuum VIS (LOCALIZED), so degree=0 is the more
    452-consistent choice; exposed so the slope's effect on the control can be probed.
    """
    lam, ratio = [], []
    for lo, hi in ANCHORS_NM:
        m = (ctx['obs_w'] >= lo) & (ctx['obs_w'] <= hi)
        if m.sum() < 20: continue
        sw, sf = _synth(ctx, lo, hi, tmp)
        oc, sc = _robust_cont(ctx['obs_w'][m], ctx['obs_f'][m], sw, sf)
        if np.isfinite(oc) and np.isfinite(sc) and oc > 0:
            lam.append(0.5*(lo+hi)); ratio.append(sc/oc)
    lam, ratio = np.array(lam), np.array(ratio)
    # robust low-order fit (clip outlier anchors once, then polyfit degree<=1)
    med, mad = np.median(ratio), np.median(np.abs(ratio - np.median(ratio)))
    keep = np.abs(ratio - med) <= 5*mad + 1e-9
    coef = np.polyfit(lam[keep], ratio[keep], degree)
    C = np.poly1d(coef)
    print(f"  [degree={degree}]")
    cmax = float(np.max(np.abs(C(lam) - 1.0)))
    print(f"  anchors used = {keep.sum()}/{len(lam)}   C(550)={C(550.):.5f}   "
          f"C(630)={C(630.):.5f}   max|C-1| over VIS = {cmax:.4f}")
    print(f"  (per-anchor ratio: median={med:.5f}  robust-MAD={mad:.5f}  "
          f"raw min/max={ratio.min():.5f}/{ratio.max():.5f})")
    if cmax > MAX_CORR:
        raise SystemExit(f"FAIL-LOUD: required correction {cmax:.4f} exceeds physical bound "
                         f"{MAX_CORR}. Not a continuum fix -- STOP.")
    return C


def _fit(ctx, obs_f, key, lo, hi, tmp):
    d = {x.key: x for x in VIS_DIAGNOSTICS}[key]
    return _fit_element(ctx['obs_w'], obs_f, ctx['atm'], ctx['params'],
                        'O' if key == 'OI_6300' else 'C', dict(SOLAR_COMP), ctx['codes'],
                        d.windows_A, d.use_molecules, ctx['broadening'], lo, hi,
                        ctx['ll'], ctx['iso'], ctx['sab'], tmp)['A_X']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar'); ap.add_argument('--region', default='vis')
    ap.add_argument('--tmp-dir', default='/tmp/ispec_cont_renorm')
    ap.add_argument('--degree', type=int, default=DEGREE,
                    help='C(lambda) poly degree: 1 (brief default) or 0 (constant; more '
                         'RYA-452-no-drift-consistent). Probe the slope effect on the control.')
    a = ap.parse_args(); os.makedirs(a.tmp_dir, exist_ok=True)
    ctx = _setup(a.star, a.region)

    print("== build bounded continuum correction C(lambda) ==")
    C = build_C(ctx, a.tmp_dir, a.degree)
    obs_corr = ctx['obs_f'] * C(ctx['obs_w'])     # apply across the spectrum (per-lambda factor)

    print("\n== gated proof ==")
    o0 = _fit(ctx, ctx['obs_f'], 'OI_6300', 7.6, 10.0, a.tmp_dir)
    o1 = _fit(ctx, obs_corr,     'OI_6300', 7.6, 10.0, a.tmp_dir)
    c0 = _fit(ctx, ctx['obs_f'], 'CI_5052', 7.5,  9.5, a.tmp_dir)
    c1 = _fit(ctx, obs_corr,     'CI_5052', 7.5,  9.5, a.tmp_dir)
    print(f"  [O I] A(O): {o0:.3f} -> {o1:.3f}   (target band {OI_BAND}; sanity o0~8.80)")
    print(f"  [C I] A(C): {c0:.3f} -> {c1:.3f}   (control: must stay {CI_REF}+/-{CI_TOL})")

    oi_ok = OI_BAND[0] <= o1 <= OI_BAND[1]
    ci_ok = abs(c1 - CI_REF) <= CI_TOL
    print("\n-- verdict ----------------------------------------------------")
    if oi_ok and ci_ok:
        print("  SUCCESS: robust bounded continuum lands [O I] in band AND leaves C I within "
              "tolerance. Solar [O I] O -> MEASUREMENT. Part B (wire C(lambda) into run_cno + "
              "all-species regression) is cleared.")
    elif not oi_ok and abs(o1 - o0) > 0.03:
        print(f"  STOP (deeper residual): correction moved [O I] to {o1:.3f}, still outside "
              f"{OI_BAND}. The robust value did not resolve the overshoot -> residual is beyond "
              f"narrow-window noise (synth continuum scale / 1D model). Post for Ryan; do NOT wire.")
    elif not ci_ok:
        print(f"  STOP (regression): C I moved to {c1:.3f}, outside {CI_REF}+/-{CI_TOL}. The "
              f"correction disturbs a clean species -> not safe to apply globally. Post for Ryan.")
    else:
        print("  AMBIGUOUS: post all four numbers for Ryan.")


if __name__ == '__main__':
    main()
