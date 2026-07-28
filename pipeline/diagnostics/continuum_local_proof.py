#!/usr/bin/env python3
"""
Per-window robust-local continuum proof (RYA-453 Path A).

RYA-453: a GLOBAL flat/linear C(lambda) can't thread "fix [O I]" vs "don't regress C I"
(degree-1 landed [O I] at 8.662 but over-corrected C I at 505). RYA-452 classified
true-continuum VIS as LOCALIZED. This corrects EACH line by ITS OWN local continuum under
a FIXED, uniform neighborhood rule, and checks [O I] lands in band AND C I 5052/5380 hold.

ANTI-TUNING (structural): W is one uniform constant for all lines; identical estimator for
all; a FIXED W-set {1.5,2.5,4.0} nm with STABILITY required (no selecting a passing W --
W-sensitivity is a STOP); multiple controls. Each correction is the MEASURED local
synth/obs ratio, never chosen to hit 8.69.

RYA-454 (per-window robust-local continuum proof; RYA-453 Path A).
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
COARSE_NM = 0.001
W_SET     = (1.5, 2.5, 4.0)     # FIXED neighborhood half-widths (nm). Do NOT add/select to pass.
CONT_Q    = 0.90                # synth-defined continuum pixels = top decile of synth flux
# (diag key, fit element, fit bounds, ('band',(lo,hi)) | ('ref',(val,tol)))
TESTS = [
    ('OI_6300', 'O', (7.6, 10.0), ('band', (8.64, 8.76))),
    ('CI_5052', 'C', (7.5,  9.5), ('ref',  (8.46, 0.03))),
    ('CI_5380', 'C', (7.5,  9.5), ('ref',  (8.46, 0.03))),
]


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


def _local_corr(ctx, lam0, W, tmp):
    """MEASURED local continuum correction = median(synth)/median(obs) at synth-defined
    continuum pixels in [lam0-W, lam0+W]. Identical rule for every line."""
    lo, hi = lam0 - W, lam0 + W
    m = (ctx['obs_w'] >= lo) & (ctx['obs_w'] <= hi)
    ow, of = ctx['obs_w'][m], ctx['obs_f'][m]
    if len(ow) < 20: return np.nan, 0
    sw, sf = _synth(ctx, lo, hi, tmp)
    sfo = np.interp(ow, sw, sf)
    cont = sfo >= np.quantile(sfo, CONT_Q)          # synth says these are continuum pixels
    if cont.sum() < 10: return np.nan, int(cont.sum())
    return float(np.median(sfo[cont]) / np.median(of[cont])), int(cont.sum())


def _fit(ctx, obs_f, d, el, lo, hi, tmp):
    return _fit_element(ctx['obs_w'], obs_f, ctx['atm'], ctx['params'], el,
                        dict(SOLAR_COMP), ctx['codes'], d.windows_A, d.use_molecules,
                        ctx['broadening'], lo, hi, ctx['ll'], ctx['iso'], ctx['sab'], tmp)['A_X']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar'); ap.add_argument('--region', default='vis')
    ap.add_argument('--tmp-dir', default='/tmp/ispec_cont_local')
    a = ap.parse_args(); os.makedirs(a.tmp_dir, exist_ok=True)
    ctx = _setup(a.star, a.region)
    diags = {x.key: x for x in VIS_DIAGNOSTICS}

    print(f"{'line':>9s} {'W':>5s} {'corr':>8s} {'npx':>5s} {'A_bef':>7s} {'A_aft':>7s} {'pass':>5s}")
    stable, sane = {}, True
    for key, el, (flo, fhi), (kind, ref) in TESTS:
        d = diags[key]
        lam0 = 0.5*(d.windows_A[0][0] + d.windows_A[0][1]) / 10.0      # window center, nm
        a_bef = _fit(ctx, ctx['obs_f'], d, el, flo, fhi, a.tmp_dir)
        if key == 'OI_6300' and abs(a_bef - 8.80) > 0.03:
            sane = False
        passes = []
        for W in W_SET:
            corr, npx = _local_corr(ctx, lam0, W, a.tmp_dir)
            oc = ctx['obs_f'].copy()
            wlo, whi = d.windows_A[0][0]/10.0 - 0.2, d.windows_A[0][1]/10.0 + 0.2
            oc[(ctx['obs_w'] >= wlo) & (ctx['obs_w'] <= whi)] *= corr
            a_aft = _fit(ctx, oc, d, el, flo, fhi, a.tmp_dir)
            ok = (ref[0] <= a_aft <= ref[1]) if kind == 'band' else (abs(a_aft - ref[0]) <= ref[1])
            passes.append(bool(ok))
            print(f"{key:>9s} {W:5.1f} {corr:8.5f} {npx:5d} {a_bef:7.3f} {a_aft:7.3f} {str(ok):>5s}")
        stable[key] = passes

    print("\n-- verdict ----------------------------------------------------")
    if not sane:
        print("  STOP (sanity): A(O) before != ~8.80 -> setup diverges from production. Reconcile.")
        return
    all_pass = all(all(p) for p in stable.values())
    unstable = any(any(p) and not all(p) for p in stable.values())
    if all_pass:
        print("  SUCCESS: every test line lands at its literature value under the fixed uniform "
              "local rule, STABLY across the W-set. LOCALIZED confirmed -- each line corrected by "
              "its OWN local continuum threads [O I] AND the clean controls. This IS the production "
              "continuum approach + measured solar O. Part B (next brief) wires per-window local.")
    elif unstable:
        print("  STOP (W-sensitive): pass/fail flips across the fixed W-set -> local continuum not "
              "robustly defined. Do NOT select a passing W. Post the full table for Ryan.")
    else:
        print("  STOP (stable fail): a line misses its gate under the uniform local rule -> "
              "per-window local does not cleanly thread. If it's [O I] that can't land, it is "
              "continuum-limited; accept Caffau-cited and stop drilling. Post for Ryan.")


if __name__ == '__main__':
    main()
