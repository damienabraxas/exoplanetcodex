"""
scripts/procyon_oi777_continuum_localize_rya478.py
==================================================
RYA-478 Phase 2b — continuum-localized O I 777 re-fit + literature-pinned broadening.

RYA-348 Phase 2 gave A(O) 1D-LTE 9.506, NLTE 8.961, [O/H] +0.23, chi2r 147 — not
science-grade. Root cause (diagnosed): the UVES arm applied BERV but NOT the stellar
systemic RV, so the O I 777 triplet sat blueshifted ~0.12 A (Procyon RV) against the
rest-frame synthesis -> misaligned cores -> inflated A(O) + chi2r 147. Plus the global
continuum is slightly sloped under the strong triplet. This re-fit:
  1. RV-to-rest is now fixed in the arm loader itself (_load_procyon_uves_arm, RYA-478).
  2. CONTINUUM-LOCALIZE: re-normalize a local window around the triplet to a robust
     linear local continuum (not the global level).
  3. Broadening PINNED to Procyon literature (cited, not fit): vsini 2.8 (Bruntt+2010),
     adopted RT vmac 6.0 (Gray scale) — stars.yaml. (This was already F-star, not solar.)
  4. Re-fit A(O) 1D-LTE, re-apply the Amarsi 2019 1D-NLTE delta self-consistently,
     report A(O), chi2r, [O/H] vs our solar O 8.735.

SOLAR CONTROL: the same localized-continuum + pinned-broadening method on a solar
O I 777 (Kitt Peak flux atlas, RYA-459/460) must NOT move the Sun — report solar
A(O) under global vs local continuum (method-stability) and vs the established 8.735.

validate-don't-tune: broadening cited, NLTE vendored (Amarsi). A clean-fit value that
is still offset is a finding, not a number to massage. STOP at verdict, no merge.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pipeline.cno_synthesis as cs                              # noqa: E402
from pipeline import nlte_cno                                    # noqa: E402
from config.constants import get_star_params                     # noqa: E402

SOLAR_OURS_O = 8.735
TMP = '/tmp/ispec_cno_rya478'
RESULTS = ROOT / 'data' / 'results'

# O I 777 geometry (air Angstrom). Fit windows = the OI_777 diagnostic's; the local
# continuum is fit on the shoulders OUTSIDE the triplet span.
OI777_FIT_WINDOWS = ((7771.0, 7772.9), (7773.2, 7775.0), (7774.6, 7776.3))
OI777_TRIPLET_SPAN = (7770.3, 7776.6)        # excluded from the continuum fit
OI777_CONT_WINDOW = (7762.0, 7784.0)         # local continuum region

# Procyon broadening — literature-pinned, CITED (stars.yaml). vsini Bruntt+2010; vmac
# adopted radial-tangential on the Gray scale (FT estimator deferred to RYA-316).
PROCYON_VMAC, PROCYON_VSINI = 6.0, 2.8
PROCYON_BROAD_REF = ('Bruntt et al. 2010 A&A 519 A51 (vsini 2.8 km/s) + adopted RT vmac '
                     '6.0 km/s (Gray scale) — config/stars.yaml; already F-star, not solar')
SOLAR_VMAC, SOLAR_VSINI = 3.8, 1.8            # solar, stars.yaml (RYA-287)
UVES_ANCHOR_R = 107200.0                      # RED760 oi_anchor resolution (registry)
KPNO_R = 350000.0                             # Kurucz 1984 NSO FTS flux atlas (high-res)


def _rule(c='='):
    print(c * 78)


def localize_continuum(wave_A, flux, cont_window=OI777_CONT_WINDOW,
                       exclude=(OI777_TRIPLET_SPAN,), deg=1, n_iter=3):
    """Local continuum: robust low-order polynomial through the line-free shoulders of
    `cont_window` (the triplet span excluded), iteratively rejecting absorption. Returns
    (wave_A_win, flux_localized) over the window. deg=1 removes a sloped pseudo-continuum."""
    lo, hi = cont_window
    m = (wave_A >= lo) & (wave_A <= hi)
    w, f = wave_A[m], flux[m]
    cm = np.ones(w.shape, bool)
    for (a, b) in exclude:
        cm &= ~((w >= a) & (w <= b))
    ww, ff = w[cm], f[cm]
    coef = np.polyfit(ww, ff, deg)
    for _ in range(n_iter):                  # reject absorption dips below the fit
        resid = ff - np.polyval(coef, ww)
        s = np.std(resid) or 1e-6
        keep = resid > -1.0 * s
        if keep.sum() < deg + 2:
            break
        ww, ff = ww[keep], ff[keep]
        coef = np.polyfit(ww, ff, deg)
    return w, f / np.polyval(coef, w)


def _resources():
    ll, iso, chem = cs._load_synth_resources()
    sab = cs.ispec.read_solar_abundances(cs._ISPEC_SOLAR_ABUND_FILE)
    codes = cs._atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    return ll, iso, sab, codes


def fit_oi777(obs_w_nm, obs_f, params, broadening, ll, iso, sab, codes, center=8.8):
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    r = cs._fit_element(obs_w_nm, obs_f, atm, params, 'O',
                        {'C': 8.46, 'N': 7.83, 'O': center, 'Ni': 6.2}, codes,
                        OI777_FIT_WINDOWS, False, broadening,
                        center - 1.2, center + 1.2, ll, iso, sab, TMP)
    return r


def nlte_delta(params, a_lte):
    lab = nlte_cno.resolve_line('OI', 7773.0)
    d = nlte_cno.cno_nlte_delta('OI', lab, params['teff_K'], params['logg'],
                                params['feh'], params['vturb_kms'], a_lte)
    if np.isfinite(d):
        nlte_cno.assert_cno_sign('OI', lab, d)
    return lab, float(d)


def run_procyon(ll, iso, sab, codes):
    _rule(); print("  PROCYON O I 777 — RV-to-rest arm + continuum-localized + pinned broadening"); _rule()
    rec = get_star_params('procyon')
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    arm = cs.star_arm_registry('procyon')['uves']
    w_nm, fl = cs.resolve_arm_spectrum('procyon', arm)        # RV-corrected (RYA-478) + globally normed
    wA = w_nm * 10.0
    # continuum-localize the triplet window, splice back
    wl_win, fl_loc = localize_continuum(wA, fl)
    fl2 = fl.copy()
    idx = np.searchsorted(wA, wl_win)
    inb = (idx >= 0) & (idx < wA.size)
    fl2[idx[inb]] = fl_loc[inb]
    print(f"  broadening (pinned, cited): vmac {PROCYON_VMAC} vsini {PROCYON_VSINI} km/s "
          f"[{PROCYON_BROAD_REF}]")
    r = fit_oi777(w_nm, fl2, params, (UVES_ANCHOR_R, PROCYON_VMAC, PROCYON_VSINI),
                  ll, iso, sab, codes, center=8.8)
    a_lte = r['A_X']
    lab, d = nlte_delta(params, a_lte)
    a_nlte = round(a_lte + d, 3)
    print(f"  A(O) 1D-LTE  = {a_lte:.3f}   chi2r = {r['red_chi2']:.2f}  [{r['status']}]  "
          f"(prior Phase-2: 9.506 / chi2r 147)")
    print(f"  Amarsi 1D-NLTE delta ({lab}) = {d:+.3f}  (self-consistent at A_lte)")
    print(f"  A(O) 1D-NLTE = {a_nlte:.3f}   [O/H] = {a_nlte - SOLAR_OURS_O:+.3f}  vs our Sun {SOLAR_OURS_O}")
    return dict(star='procyon', a_lte=a_lte, chi2r=r['red_chi2'], delta=round(d, 3),
                a_nlte=a_nlte, oh=round(a_nlte - SOLAR_OURS_O, 3), status=r['status'],
                vmac=PROCYON_VMAC, vsini=PROCYON_VSINI, broad_ref=PROCYON_BROAD_REF)


def run_solar_control(ll, iso, sab, codes):
    _rule(); print("  SOLAR CONTROL — Kitt Peak O I 777 (RYA-459/460): method must not move the Sun"); _rule()
    rec = get_star_params('solar')
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    kp = pd.read_csv(ROOT / 'data' / 'solar_reference' / 'kpno_flux_atlas' / 'kpno_OI_777_triplet.csv',
                     comment='#')
    wA = kp['wavelength_air_A'].to_numpy(float)
    fl = kp['residual_flux'].to_numpy(float)            # atlas continuum (already normalized)
    w_nm = wA / 10.0
    broad = (KPNO_R, SOLAR_VMAC, SOLAR_VSINI)
    # (a) atlas/global continuum as-is
    r_glob = fit_oi777(w_nm, fl, params, broad, ll, iso, sab, codes, center=8.9)
    # (b) localized continuum
    wl_win, fl_loc = localize_continuum(wA, fl)
    fl2 = fl.copy()
    idx = np.searchsorted(wA, wl_win); inb = (idx >= 0) & (idx < wA.size)
    fl2[idx[inb]] = fl_loc[inb]
    r_loc = fit_oi777(w_nm, fl2, params, broad, ll, iso, sab, codes, center=8.9)
    lab, d = nlte_delta(params, r_loc['A_X'])
    a_nlte_loc = round(r_loc['A_X'] + d, 3)
    print(f"  broadening (solar, stars.yaml): vmac {SOLAR_VMAC} vsini {SOLAR_VSINI} km/s")
    print(f"  A(O) LTE  global-continuum = {r_glob['A_X']:.3f}  chi2r {r_glob['red_chi2']:.2f}")
    print(f"  A(O) LTE  local-continuum  = {r_loc['A_X']:.3f}  chi2r {r_loc['red_chi2']:.2f}")
    print(f"  method shift (local - global) = {r_loc['A_X'] - r_glob['A_X']:+.3f} dex  "
          f"(must be ~0 -> method does not move the Sun)")
    print(f"  Amarsi NLTE delta = {d:+.3f}  -> A(O) NLTE (local) = {a_nlte_loc:.3f}  "
          f"(established solar O {SOLAR_OURS_O}; KP-ESPRESSO zero-point ~0.04, RYA-460)")
    return dict(star='solar_kpno', a_lte_global=r_glob['A_X'], a_lte_local=r_loc['A_X'],
                method_shift=round(r_loc['A_X'] - r_glob['A_X'], 3),
                chi2r_local=r_loc['red_chi2'], delta=round(d, 3), a_nlte_local=a_nlte_loc)


def main():
    Path(TMP).mkdir(parents=True, exist_ok=True)
    ll, iso, sab, codes = _resources()
    sol = run_solar_control(ll, iso, sab, codes)
    pro = run_procyon(ll, iso, sab, codes)
    _rule(); print("  VERDICT INPUTS (interpretation is Ryan's + Claude's)"); _rule()
    bankable = (pro['chi2r'] < 10 and abs(pro['oh']) < 0.10 and abs(sol['method_shift']) < 0.03)
    print(f"  Procyon O I 777: chi2r {pro['chi2r']:.1f} (prior 147), A(O)_NLTE {pro['a_nlte']}, "
          f"[O/H] {pro['oh']:+.3f}")
    print(f"  Solar control method shift: {sol['method_shift']:+.3f} (Sun stable if ~0)")
    print(f"  -> first-pass read: {'BANKABLE' if bankable else 'still a finding/owed (see chi2r and [O/H])'}")
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'procyon_oi777_relfit_rya478.json').write_text(
        json.dumps({'ticket': 'RYA-478 Phase 2b', 'procyon': pro, 'solar_control': sol},
                   indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'procyon_oi777_relfit_rya478.json'}")
    return pro, sol


if __name__ == '__main__':
    main()
