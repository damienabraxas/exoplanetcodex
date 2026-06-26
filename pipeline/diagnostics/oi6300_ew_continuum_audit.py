#!/usr/bin/env python3
"""
[O I] 6300 1D-LTE seam isolation -- observed EW + continuum audit.

RYA-447: 3D differential ~0. RYA-448: Ni partition ruled out (dA(O)/dA(Ni)=-0.38;
even A(Ni)=6.25 buys only ~0.02 dex; the -2.11 Ni gf is Johansson 2003 lab, not an
A(Ni) basis -- RYA-450). The ~0.11 dex excess in our 1D-LTE
[O I] 6300 O lives in the MEASUREMENT level. This audit bifurcates: is our
observed [O I]+Ni feature too STRONG (continuum/normalization/telluric/blend ->
measurement artifact), or NORMAL while our synthesis over-converts EW->A (1D
model/conversion problem)?

Cited solar anchors (comparison only, never a fit target):
  total [O I]+Ni feature EW, Sun:  5.3-5.5 mA (Nissen 2002 5.5; Bedell/Nissen 5.3;
                                               Schuler 2006 5.5; KH96 synth 5.4)
  sensitivity:  +/-0.5 mA -> +/-0.06 dex  => dA/dEW ~ 0.12 dex/mA  (Schuler 2006)
  a CORRECT 1D MARCS analysis of a 5.5 mA feature gives A(O) ~ 8.74 (Nissen 2002),
  NOT 8.80 -- a normal EW that still yields 8.80 indicts the conversion.

RYA-449 (EW + continuum audit; RYA-448 follow-up).
"""
from __future__ import annotations
import argparse
import os
import numpy as np

import pipeline.cno_synthesis as cs
from pipeline.cno_synthesis import (
    REGIONS, VIS_DIAGNOSTICS, SOLAR_ASPLUND2021, get_star_params,
    _load_atmosphere, _load_synth_resources, _load_observed_spectrum,
    _atom_codes, _fixed_ab, _synth_window, preflight,
)
import ispec

LIT_EW_MA = (5.3, 5.5)     # total [O I]+Ni feature EW, Sun (cited above)
DA_DEW    = 0.12           # dex per mA (Schuler 2006)
LIT_1D_AO = 8.74           # Nissen 2002: correct 1D MARCS, 5.5 mA feature
OUR_AO    = 8.80
ASPLUND_O = 8.69

# Feature = 0.35 A centered 6300.31 (Nissen 2002 window), in nm. Clean flanks for continuum.
FEAT_NM = (630.0135, 630.0485)
# RYA-449: the spec's flanks (629.960-630.000, 630.060-630.100) BOTH catch absorption
# in this HARPS-solar spectrum -- a strong line at ~6299.6 (flux ~0.88) and the
# 6300.65 line (flux ~0.94) -- pulling the "continuum" to 0.9806 and garbling the EW.
# Nudged to the genuinely line-free flanks bracketing the feature, inside the clean
# 6299.6-6301.0 region the spec told us to stay within (note logged on the issue).
CONT_NM = [(630.000, 630.015), (630.040, 630.056)]


def _ew_mA(w_nm, f, cont, lo_nm, hi_nm) -> float:
    m = (w_nm >= lo_nm) & (w_nm <= hi_nm)
    if m.sum() < 3:
        return np.nan
    depth = 1.0 - f[m] / cont
    return float(np.trapz(depth, w_nm[m]) * 1.0e4)     # nm -> mA (1 nm = 1e4 mA)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar')
    ap.add_argument('--region', default='vis')
    ap.add_argument('--a-c', type=float, default=8.491)   # converged solar (RYA-448 run)
    ap.add_argument('--a-n', type=float, default=7.558)   # converged solar (RYA-448 run)
    ap.add_argument('--a-ni', type=float, default=None,
                    help='A(Ni) for the synth cross-check (default SOLAR_ASPLUND2021).')
    ap.add_argument('--tmp-dir', default='/tmp/ispec_oi_ew')
    a = ap.parse_args()
    os.makedirs(a.tmp_dir, exist_ok=True)

    # ── observed feature ─────────────────────────────────────────────────────
    obs_w, obs_f = _load_observed_spectrum(a.star)        # obs_w in nm
    cont_vals = [np.median(obs_f[(obs_w >= lo) & (obs_w <= hi)])
                 for lo, hi in CONT_NM if ((obs_w >= lo) & (obs_w <= hi)).any()]
    cont = float(np.median(cont_vals)) if cont_vals else 1.0
    ew_obs    = _ew_mA(obs_w, obs_f, cont, *FEAT_NM)      # vs measured local continuum
    ew_obs_c1 = _ew_mA(obs_w, obs_f, 1.0,  *FEAT_NM)      # vs continuum=1 (isolates norm)

    print(f"observed [O I] 6300 feature audit ({a.star}):")
    print(f"  local continuum (flanks)    = {cont:.4f}   (ideal 1.0; <1 inflates feature)")
    print(f"  EW (vs local continuum)     = {ew_obs:.2f} mA")
    print(f"  EW (vs continuum=1.0)       = {ew_obs_c1:.2f} mA")
    print(f"  literature solar feature EW = {LIT_EW_MA[0]:.1f}-{LIT_EW_MA[1]:.1f} mA "
          f"(Nissen 2002 / Bedell / Schuler / KH96)")
    d_ew = ew_obs - sum(LIT_EW_MA) / 2.0
    implied_dA = d_ew * DA_DEW
    print(f"  EW excess vs lit mid        = {d_ew:+.2f} mA  ->  implied dA(O) = {implied_dA:+.3f} dex")

    # ── synthesis cross-check: does our model reproduce the observed EW, at what A(O)? ──
    ni = float(a.a_ni) if a.a_ni is not None else float(SOLAR_ASPLUND2021['Ni'])
    rec = get_star_params(a.star)
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    region = REGIONS[a.region]
    broadening = preflight(region, a.star, VIS_DIAGNOSTICS)
    atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(cs._ISPEC_SOLAR_ABUND_FILE)
    codes = _atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    sw = np.arange(FEAT_NM[0] - 0.05, FEAT_NM[1] + 0.05, 0.0002)

    def synth_flux(aO):
        fa = _fixed_ab({'C': a.a_c, 'N': a.a_n, 'O': aO, 'Ni': ni}, codes)
        return _synth_window(sw, atm, params, ll, iso, sab, fa, broadening, True, a.tmp_dir)

    sf_our = synth_flux(OUR_AO)
    ew_syn_our = _ew_mA(sw, sf_our, 1.0, *FEAT_NM)
    ew_syn_asp = _ew_mA(sw, synth_flux(ASPLUND_O), 1.0, *FEAT_NM)
    # Synth continuum in the same clean flanks the observed continuum was read from
    # (sw already spans CONT_NM). The fit renders synth on a TRUE continuum ~1.0, so a
    # gap (obs_cont - synth_cont) below ~0 is the normalization offset the direct flux
    # chi2 fit converts into excess line depth -> inflated A(O). RYA-449.
    synth_cont = float(np.median([np.median(sf_our[(sw >= lo) & (sw <= hi)])
                                  for lo, hi in CONT_NM]))
    cont_offset = cont - synth_cont
    print(f"\nsynthesis cross-check (A(Ni)={ni:.2f}, A(C)={a.a_c}, A(N)={a.a_n}):")
    print(f"  synth EW at our A(O)={OUR_AO}      = {ew_syn_our:.2f} mA  (should ~= observed EW vs cont=1.0)")
    print(f"  synth EW at Asplund A(O)={ASPLUND_O}  = {ew_syn_asp:.2f} mA")
    print(f"  synth local continuum (clean flanks) = {synth_cont:.4f}")
    print(f"  continuum offset (obs - synth)       = {cont_offset:+.4f}   "
          f"(<0 => observed continuum normalized low -> fit deepens line -> inflates A(O))")

    # ── verdict ───────────────────────────────────────────────────────────────
    print("\n-- verdict ------------------------------------------------------")
    if (ew_obs > LIT_EW_MA[1] + 0.5) or (cont < 0.99):
        feat_strong = ew_obs > LIT_EW_MA[1] + 0.5
        cause = ("observed feature too strong (EW {ew:.2f} mA > lit {hi:.1f})".format(
                    ew=ew_obs, hi=LIT_EW_MA[1])
                 if feat_strong else
                 "continuum normalized below the synthesis continuum (obs {c:.4f} vs synth "
                 "{s:.4f}, offset {o:+.4f}); EW vs that local continuum is {ew:.2f} mA (~ the "
                 "solar band), so the excess A(O) is a normalization artifact the direct flux "
                 "fit reads as a too-strong feature, NOT a physically strong line".format(
                    c=cont, s=synth_cont, o=cont_offset, ew=ew_obs))
        print(f"  MEASUREMENT seam: {cause}. Suspects in order: (1) continuum/"
              f"normalization [audited obs {cont:.4f} vs synth {synth_cont:.4f}], (2) telluric "
              f"residual in the 6300 window [check BERV/telluric -- codex-data-audit], "
              f"(3) unmodeled blend (CN, usually minor). Fix the measurement; do NOT tune A(O).")
    else:
        print(f"  MODEL/CONVERSION seam: observed EW ({ew_obs:.2f} mA) is within the solar band, yet "
              f"our synthesis needs A(O)={OUR_AO} to reproduce it while a correct 1D MARCS gives "
              f"~{LIT_1D_AO} from the same EW (Nissen 2002). Seam = EW->abundance conversion (1D "
              f"structure / microturbulence / continuum-opacity scale). Document and scope; do NOT tune.")


if __name__ == '__main__':
    main()
