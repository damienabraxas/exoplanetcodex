#!/usr/bin/env python3
# =============================================================================
# THE EXOPLANET CODEX
# =============================================================================
# File:         vmac_sensitivity_rya309.py
# Module:       scripts (RYA-309 rigor rider)
# Description:  A(Fe) sensitivity to the adopted RT macroturbulence: re-derive
#               A(Fe) on a representative strong Fe I line subset at vmac = 5.0
#               and 6.0 km/s (vsini fixed 2.8, RV-corrected obs, converged
#               params) and report the per-line and MAX ΔA(Fe). Quantifies how
#               much the adopted-vs-alternative vmac moves the headline A(Fe)
#               (expected small — vmac couples weakly to abundance), justifying
#               the literature-adopted vmac while the FT estimator waits on
#               RYA-316.
#
# Usage:  python scripts/vmac_sensitivity_rya309.py --feh <converged> --xi <converged>
# =============================================================================

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
import pipeline.abundances_derive as ad   # inserts ISPEC_DIR on sys.path
import ispec
from config.constants import HARPS_R

TMP = '/tmp/ispec_rya309_rider'
TEFF, LOGG = 6554.0, 4.00
VSINI = 2.8
VMAC_PAIR = (5.0, 6.0)
N_LINES = 16
EW_LO, EW_HI = 60.0, 140.0       # strong but unsaturated — vmac-sensitive wings
WL_LO, WL_HI = 4500.0, 6500.0


def strong_fe1(ew_csv):
    e = pd.read_csv(ew_csv)
    m = ((e['element'] == 'Fe') & (e['ion'] == 'I') &
         (~e['blend_flag'].astype(bool)) &
         e['ew_mA'].between(EW_LO, EW_HI) &
         e['wavelength_air_A'].between(WL_LO, WL_HI))
    sel = e[m].sort_values('chi2').head(N_LINES * 3)
    sel = sel.sort_values('wavelength_air_A').iloc[::max(1, len(sel) // N_LINES)]
    return sel.head(N_LINES).reset_index(drop=True)


def fit_A(line, vmac, ctx):
    ow, of, atm, ll, iso, sab, fe_code, fe_solar, feh, xi = ctx
    wave_nm = float(line['wavelength_air_A']) / 10.0
    wbase, wtop = ad._wingwide_window_nm(wave_nm, float(line['ew_mA']))
    r = ad._fit_synth_flux(ow, of, atm, TEFF, LOGG, feh, xi, ll, iso, sab,
                           'Fe', fe_code, wbase, wtop, max(fe_solar - 3, 1.0),
                           fe_solar + 5, HARPS_R, vmac, VSINI, tmp_dir=TMP)
    return r['A_X'] if r['status'] in ('ok', 'edge_pinned') else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--feh', type=float, required=True, help='converged [Fe/H]')
    ap.add_argument('--xi', type=float, required=True, help='converged microturbulence ξ km/s')
    a = ap.parse_args()
    os.makedirs(TMP, exist_ok=True)

    print(f"{'='*70}\nRYA-309 rigor rider — A(Fe) vmac sensitivity ({VMAC_PAIR[0]} vs {VMAC_PAIR[1]} km/s)\n{'='*70}")
    print(f"  params: Teff {TEFF} logg {LOGG} [Fe/H] {a.feh} ξ {a.xi}  vsini {VSINI} (fixed)")
    lines = strong_fe1(REPO_ROOT / 'data' / 'processed' / 'procyon_ew.csv')
    print(f"  strong Fe I subset: {len(lines)} lines "
          f"({lines['wavelength_air_A'].min():.0f}-{lines['wavelength_air_A'].max():.0f} Å, "
          f"EW {lines['ew_mA'].min():.0f}-{lines['ew_mA'].max():.0f} mÅ)")

    atm = ad._load_atmosphere(TEFF, LOGG, a.feh, a.xi)
    ll, iso, chem = ad._load_synth_resources()
    sab = ispec.read_solar_abundances(ad._ISPEC_SOLAR_ABUND_FILE)
    tmp = ispec.create_free_abundances_structure(['Fe'], chem, sab)
    fe_code, fe_solar = int(tmp['code'][0]), float(tmp['Abund'][0]) + 12.036
    ow, of = ad._load_observed_spectrum('procyon')
    ctx = (ow, of, atm, ll, iso, sab, fe_code, fe_solar, a.feh, a.xi)

    rows = []
    for _, ln in lines.iterrows():
        a5 = fit_A(ln, VMAC_PAIR[0], ctx)
        a6 = fit_A(ln, VMAC_PAIR[1], ctx)
        d = a6 - a5 if np.isfinite(a5) and np.isfinite(a6) else np.nan
        rows.append(dict(wl=ln['wavelength_air_A'], ew=ln['ew_mA'],
                         A_vmac5=a5, A_vmac6=a6, dA=d))
        print(f"  {ln['wavelength_air_A']:8.2f} Å  EW={ln['ew_mA']:5.0f}  "
              f"A(5.0)={a5:6.3f}  A(6.0)={a6:6.3f}  ΔA={d:+.3f}")
    df = pd.DataFrame(rows)
    out = REPO_ROOT / 'data' / 'processed' / 'procyon_vmac_sensitivity_rya309.csv'
    df.to_csv(out, index=False)
    dA = df['dA'].dropna().abs()
    print(f"\n{'='*70}")
    print(f"  n={len(dA)}  median |ΔA(Fe)| = {dA.median():.3f} dex  "
          f"MAX |ΔA(Fe)| = {dA.max():.3f} dex  (vmac 5.0↔6.0)")
    print(f"  → {out.relative_to(REPO_ROOT)}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
