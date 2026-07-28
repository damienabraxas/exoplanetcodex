#!/usr/bin/env python3
"""
[O I] 6300 1D-LTE seam isolation -- A(Ni) partition sensitivity sweep.

RYA-447 STOP established the ~0.07-0.11 dex excess in our 1D-LTE [O I] 6300 O
(our 8.80 vs Caffau 8.73 vs Asplund 8.69) is NOT a 3D-correction problem
(Amarsi 630nm Delta ~= 0). It lives in the 1D-LTE absolute level. [O I] 6300.30
/ Ni I 6300.34 is an unresolved blend; A(O) is partitioned against the PINNED
A(Ni). We pin A(Ni) = SOLAR_ASPLUND2021['Ni'] (6.20, Asplund 2021). The Ni gf
-2.11 (RYA-365) is the Johansson et al. 2003 LABORATORY value (FT emission,
abundance-independent) -- it is NOT tied to any A(Ni) basis, so there is no
"matched gf/abundance pair" to restore. (Allende Prieto, Lambert & Asplund 2001's
*astrophysical* gf was ~ -2.31, ~0.2 dex below the lab value.) The only lever is
the solar A(Ni) itself, and plausible modern values (6.15 Scott 2009/Bedell ->
6.20-6.22 Asplund) all leave O high -- a LOWER A(Ni) weakens Ni and pushes O even
higher. (RYA-450 provenance reconciliation; corrects the RYA-448 brief's framing.)

This MEASURES dA(O)/dA(Ni) by re-fitting A(O) over the real [O I] window at fixed
A(Ni) values, reusing the production fitter (_fit_element). It does NOT change any
pinned value and does NOT tune A(Ni) to an anchor -- adopting a new A(Ni) is a
separate, independently-justified decision for Ryan.

RYA-448 (isolate the [O I] 6300 1D-LTE seam; RYA-447 STOP follow-up).
"""
from __future__ import annotations
import argparse
import os
import numpy as np

import pipeline.cno_synthesis as cs
from pipeline.cno_synthesis import (
    REGIONS, VIS_DIAGNOSTICS, SOLAR_ASPLUND2021, get_star_params,
    _load_atmosphere, _load_synth_resources, _load_observed_spectrum,
    _atom_codes, preflight, _fit_element,
)
import ispec

ASPLUND_O = 8.69
CAFFAU_O  = 8.73
REF_NI_HI = 6.25   # reference high A(Ni) sweep point. The -2.11 Ni gf is Johansson 2003
                   # (lab, abundance-independent) -- NOT tied to any A(Ni) basis. Allende
                   # Prieto 2001's astrophysical gf was ~ -2.31 (~0.2 dex lower).
REF_NI_LO = 6.15   # reference low A(Ni) (Scott 2009 / Bedell) -- brackets the modern solar range.
DEFENSIBLE_NI_HI  = 6.28   # ~ upper edge of plausible solar A(Ni) (Asplund 6.20 +/- ~0.04)


def _setup(star_id, region_name):
    """Mirror run_cno's setup so the O fit is byte-identical to production."""
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
    diag = {d.key: d for d in VIS_DIAGNOSTICS}['OI_6300']
    return dict(params=params, broadening=broadening, atm=atm, ll=ll, iso=iso,
                sab=sab, obs_w=obs_w, obs_f=obs_f, codes=codes, diag=diag)


def _fit_O_at_ni(a_ni, a_c, a_n, o_center, ctx, tmp_dir) -> float:
    d = ctx['diag']
    state = {'C': a_c, 'N': a_n, 'O': o_center, 'Ni': a_ni}   # Ni FIXED at the trial value
    r = _fit_element(ctx['obs_w'], ctx['obs_f'], ctx['atm'], ctx['params'], 'O',
                     state, ctx['codes'], d.windows_A, d.use_molecules, ctx['broadening'],
                     o_center - 1.2, o_center + 1.2, ctx['ll'], ctx['iso'], ctx['sab'], tmp_dir)
    return float(r['A_X']) if np.isfinite(r['A_X']) else np.nan


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar')
    ap.add_argument('--region', default='vis')
    ap.add_argument('--a-c', type=float, default=8.49,
                    help='Converged 1D-LTE A(C) from CH (latest solar run). Affects the '
                         'molecular EOS in the window -- keep matched to the run.')
    ap.add_argument('--a-n', type=float, default=7.90,
                    help='Converged 1D-LTE A(N) from CN (latest solar run).')
    ap.add_argument('--o-center', type=float, default=8.80,
                    help='Our measured 1D-LTE A(O); fit-bracket center + self-check target.')
    ap.add_argument('--ni-grid', default='6.10,6.15,6.20,6.25,6.30,6.35')
    ap.add_argument('--tmp-dir', default='/tmp/ispec_oi_ni_sweep')
    a = ap.parse_args()
    os.makedirs(a.tmp_dir, exist_ok=True)

    ctx = _setup(a.star, a.region)
    ni_pin = float(SOLAR_ASPLUND2021['Ni'])
    print(f"current pinned A(Ni) = SOLAR_ASPLUND2021['Ni'] = {ni_pin:.2f} (Asplund 2021)")
    print(f"Ni gf -2.11 (RYA-365) = Johansson 2003 lab value (abundance-independent; NOT an "
          f"A(Ni) basis). Allende Prieto 2001 astrophysical gf = ~ -2.31.\n")

    ni_grid = sorted({round(float(x), 3) for x in a.ni_grid.split(',')}
                     | {round(ni_pin, 3), REF_NI_HI, REF_NI_LO})

    print(f"{'A(Ni)':>7s}  {'A(O) fit':>9s}  note")
    res = {}
    for ni in ni_grid:
        aO = _fit_O_at_ni(ni, a.a_c, a.a_n, a.o_center, ctx, a.tmp_dir)
        res[ni] = aO
        note = ('<- current pin' if abs(ni - ni_pin) < 1e-6 else
                '<- ref hi (Asplund/AP01 6.25)' if abs(ni - REF_NI_HI) < 1e-6 else
                '<- ref lo (Scott09/Bedell 6.15)' if abs(ni - REF_NI_LO) < 1e-6 else '')
        print(f"{ni:7.3f}  {aO:9.3f}  {note}")

    # self-check: the current pin MUST reproduce production A(O) ~ o_center
    aO_pin = res.get(round(ni_pin, 3), np.nan)
    print(f"\n[self-check] A(O) at current pin {ni_pin:.2f} = {aO_pin:.3f}  "
          f"(production = {a.o_center:.2f})")
    if not np.isfinite(aO_pin) or abs(aO_pin - a.o_center) > 0.02:
        print("CRITICAL: setup does not reproduce production A(O) at the current pin "
              "(>0.02 dex). Sweep INVALID -- match --a-c/--a-n to the latest run or check "
              "the setup. STOP.")
        return

    # local sensitivity dA(O)/dA(Ni) near the pin (finite difference)
    lo = min(res, key=lambda k: abs(k - (ni_pin - 0.05)))
    hi = min(res, key=lambda k: abs(k - (ni_pin + 0.05)))
    if hi != lo and np.isfinite(res[lo]) and np.isfinite(res[hi]):
        print(f"[sensitivity] dA(O)/dA(Ni) ~ {(res[hi]-res[lo])/(hi-lo):+.3f} (near the pin)")

    # A(Ni) implied by each O target (interp across the monotone sweep)
    xs = np.array([k for k in ni_grid if np.isfinite(res[k])])
    ys = np.array([res[k] for k in xs])
    o = np.argsort(ys)
    ni_for = lambda t: float(np.interp(t, ys[o], xs[o])) if ys.min() <= t <= ys.max() else np.nan
    ni_asp, ni_caf = ni_for(ASPLUND_O), ni_for(CAFFAU_O)
    print(f"\nA(Ni) to yield A(O)={ASPLUND_O} (Asplund): {ni_asp:.3f}")
    print(f"A(Ni) to yield A(O)={CAFFAU_O} (Caffau):  {ni_caf:.3f}")

    print("\n-- verdict ------------------------------------------------------")
    if np.isfinite(ni_caf) and ni_caf <= DEFENSIBLE_NI_HI:
        print(f"  Ni-PARTITION LIKELY: A(Ni) ~ {ni_caf:.2f}-{ni_asp:.2f} closes the gap and is "
              f"within the defensible solar range (<= {DEFENSIBLE_NI_HI:.2f}). Seam = the adopted "
              f"solar A(Ni). RECOMMEND (for Ryan, not auto-applied): re-examine the adopted solar "
              f"A(Ni). The -2.11 gf is the Johansson 2003 lab value (fixed, abundance-independent); "
              f"this would be a CITED A(Ni) choice, NOT tuning.")
    else:
        print(f"  Ni-PARTITION INSUFFICIENT: closing the gap needs A(Ni) ~ {ni_caf:.2f}, beyond "
              f"the defensible range (> {DEFENSIBLE_NI_HI:.2f}). Ni is not the whole story -> move "
              f"to EW / continuum / 1D-model. Do NOT force A(Ni). Post for Ryan.")


if __name__ == '__main__':
    main()
