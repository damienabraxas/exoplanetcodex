#!/usr/bin/env python3
# =============================================================================
# THE EXOPLANET CODEX
# =============================================================================
# File:         fit_procyon_vmac_rya309.py
# Module:       scripts (RYA-309 §4 synthesis half — RT macroturbulence fit)
# Description:  Fit Procyon's radial-tangential (RT) macroturbulence vmac in
#               synthesis-v2's model, with vsini held FIXED at the Bruntt+2010
#               value (2.8 km/s) per RYA-309 §3.3. Bruntt's own vmac (4.6) is
#               GAUSSIAN and must not be reused as an RT vmac (solar: RT 3.8 vs
#               Bruntt Gaussian 2.4) — so the RT vmac is fit here.
#
#               Method: for a grid of trial vmac, fit each clean Fe I line's
#               abundance with the existing v2 flux fitter (_fit_synth_flux) at
#               (R=HARPS, vsini=2.8, vmac=trial) and record the profile χ²ᵣ.
#               vmac is shape-sensitive (line width), abundance is depth — so the
#               vmac minimizing the summed χ²ᵣ over clean lines is the RT vmac.
#               A coarse grid is refined once around the minimum.
#
#               vmac is fit at pinned Teff/logg (Heiter+2015 6554/4.00) and
#               nominal [Fe/H]/ξ; broadening is only weakly sensitive to those,
#               and the per-line abundance nuisance absorbs depth differences.
#
# Linear issue: RYA-309 §4 (synthesis half)
# =============================================================================

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
import pipeline.abundances_derive as ad   # inserts ISPEC_DIR on sys.path
import ispec                               # importable only after the line above
from config.constants import HARPS_R, get_star_params

TMP = '/tmp/ispec_rya309_vmac'
os.makedirs(TMP, exist_ok=True)

# Pinned fundamentals (RYA-309 §3.4) + nominal solved params for the vmac fit.
TEFF, LOGG = 6554.0, 4.00
FEH_NOMINAL, XI_NOMINAL = -0.04, 1.8

N_CLEAN = 14            # clean Fe I lines to fit vmac on
EW_LO, EW_HI = 25.0, 85.0   # mÅ — unsaturated, well-formed
WL_LO, WL_HI = 4500.0, 6300.0   # Å — good continuum, away from blue forest / Hα


def select_clean_fe1(ew_csv):
    """Deterministic clean Fe I line set: unblended, unsaturated, good Gaussian
    fit, in a clean-continuum band — the lines whose width best constrains vmac."""
    e = pd.read_csv(ew_csv)
    m = ((e['element'] == 'Fe') & (e['ion'] == 'I') &
         (~e['blend_flag'].astype(bool)) &
         e['ew_mA'].between(EW_LO, EW_HI) &
         e['wavelength_air_A'].between(WL_LO, WL_HI) &
         e['ew_err_mA'].notna())
    sel = e[m].copy()
    # best-formed profiles first (lowest fit χ²), spread across the band
    sel = sel.sort_values('chi2').head(N_CLEAN * 3)
    sel = sel.sort_values('wavelength_air_A').iloc[::max(1, len(sel) // N_CLEAN)]
    return sel.head(N_CLEAN).reset_index(drop=True)


def mean_chi2_at_vmac(vmac, lines, ctx):
    obs_w, obs_f, atm, ll, iso, sab, fe_code, fe_solar = ctx
    chi2s = []
    for _, r in lines.iterrows():
        wave_nm = float(r['wavelength_air_A']) / 10.0
        wbase, wtop = ad._wingwide_window_nm(wave_nm, float(r['ew_mA']))
        res = ad._fit_synth_flux(obs_w, obs_f, atm, TEFF, LOGG, FEH_NOMINAL, XI_NOMINAL,
                                 ll, iso, sab, 'Fe', fe_code,
                                 wbase, wtop, max(fe_solar - 3, 1.0), fe_solar + 5,
                                 HARPS_R, vmac, 2.8, tmp_dir=TMP)
        if res['status'] in ('ok', 'edge_pinned') and np.isfinite(res['red_chi2']):
            chi2s.append(res['red_chi2'])
    return (float(np.mean(chi2s)), len(chi2s)) if chi2s else (np.inf, 0)


def main():
    print(f"{'='*70}\nRYA-309 §4 — Procyon RT vmac fit (vsini fixed 2.8, Bruntt+2010)\n{'='*70}")
    rec = get_star_params('procyon')
    print(f"  vsini={rec['vsini']} (fixed)  vmac policy={rec['vmac']}  "
          f"init={rec.get('vmac_init')}  bounds={rec.get('vmac_fit_bounds')}")
    print(f"  pinned Teff/logg={TEFF}/{LOGG}  nominal [Fe/H]/ξ={FEH_NOMINAL}/{XI_NOMINAL}")

    lines = select_clean_fe1(REPO_ROOT / 'data' / 'processed' / 'procyon_ew.csv')
    print(f"  clean Fe I lines for vmac fit: {len(lines)} "
          f"({lines['wavelength_air_A'].min():.0f}-{lines['wavelength_air_A'].max():.0f} Å, "
          f"EW {lines['ew_mA'].min():.0f}-{lines['ew_mA'].max():.0f} mÅ)")

    print("  loading atmosphere + synth resources …", flush=True)
    atm = ad._load_atmosphere(TEFF, LOGG, FEH_NOMINAL, XI_NOMINAL)
    ll, iso, chem = ad._load_synth_resources()
    sab = ispec.read_solar_abundances(ad._ISPEC_SOLAR_ABUND_FILE)
    tmp = ispec.create_free_abundances_structure(['Fe'], chem, sab)
    fe_code, fe_solar = int(tmp['code'][0]), float(tmp['Abund'][0]) + 12.036
    obs_w, obs_f = ad._load_observed_spectrum('procyon')
    ctx = (obs_w, obs_f, atm, ll, iso, sab, fe_code, fe_solar)

    lo, hi = rec.get('vmac_fit_bounds', (3.0, 8.0))
    coarse = np.round(np.arange(lo, hi + 1e-6, 0.5), 2)
    print(f"\n  Coarse grid {coarse[0]}–{coarse[-1]} step 0.5:")
    t0 = time.time()
    results = {}
    for v in coarse:
        mc, n = mean_chi2_at_vmac(v, lines, ctx)
        results[float(v)] = mc
        print(f"    vmac={v:4.2f}  mean χ²ᵣ={mc:8.4f}  (n={n})  [{time.time()-t0:.0f}s]", flush=True)
    best = min(results, key=results.get)

    fine = np.round(np.arange(max(lo, best - 0.5), min(hi, best + 0.5) + 1e-6, 0.25), 2)
    print(f"\n  Refine around {best} (step 0.25):")
    for v in fine:
        if float(v) in results:
            continue
        mc, n = mean_chi2_at_vmac(v, lines, ctx)
        results[float(v)] = mc
        print(f"    vmac={v:4.2f}  mean χ²ᵣ={mc:8.4f}  (n={n})  [{time.time()-t0:.0f}s]", flush=True)
    best = min(results, key=results.get)

    out = REPO_ROOT / 'data' / 'processed' / 'procyon_vmac_fit_rya309.csv'
    pd.DataFrame(sorted(results.items()), columns=['vmac_kms', 'mean_red_chi2']).to_csv(out, index=False)
    print(f"\n{'='*70}")
    print(f"  BEST RT vmac = {best:.2f} km/s  (mean χ²ᵣ={results[best]:.4f})  "
          f"vsini=2.8 fixed")
    print(f"  expected ~5-6 per RYA-309 §3.3 → {'IN RANGE' if 5.0 <= best <= 6.0 else 'OUT OF EXPECTED RANGE — inspect'}")
    print(f"  grid → {out.relative_to(REPO_ROOT)}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
