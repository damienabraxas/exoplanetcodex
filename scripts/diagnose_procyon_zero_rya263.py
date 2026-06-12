"""
scripts/diagnose_procyon_zero_rya263.py
========================================
Trace a single Fe I line through abundance derivation for Procyon.
Print every intermediate value: EW input, model atmosphere loaded,
SPECTRUM call parameters, ispec.determine_abundances output.

RYA-263: all Procyon per-line a_1dlte = 0.0 — root cause investigation.

Target line: Fe I 5379.574 Å (EW ~51.6 mÅ, clean, unblended).
Procyon params: Teff=6530 K, logg=3.96, [Fe/H]=-0.04, vmic=1.66 km/s
"""

import sys
import os
import warnings
import numpy as np
import pandas as pd

# ── Path setup ───────────────────────────────────────────────────────────────
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from config.constants import (
    ISPEC_DIR, PIPELINE, SOLAR_ASPLUND2021,
)
sys.path.insert(0, str(ISPEC_DIR))
import ispec

# ── Config ───────────────────────────────────────────────────────────────────
TARGET_WAVE   = 5379.574   # Å — Fe I, EW~51.6 mÅ, unblended
TARGET_TOL    = 0.05       # Å matching tolerance

STAR = dict(teff_K=6530.0, logg=3.96, feh=-0.04, vturb_kms=1.66)

ATM_GRID = str(ISPEC_DIR / 'input' / 'atmospheres' / 'ATLAS9.Castelli')
SOLAR_ABUND_FILE = str(ISPEC_DIR / 'input' / 'abundances' / 'Asplund.2009' / 'stdatom.dat')
LINELIST_FILE = str(
    ISPEC_DIR / 'input' / 'regions' / '42000_VALD' /
    'limited_but_with_missing_elements_turbospectrum_synth_good_for_abundances_all_extended.txt'
)
EW_CSV = os.path.join(REPO, 'data', 'processed', 'procyon_ew.csv')

SEP = '─' * 60


def section(title: str) -> None:
    print(f'\n{SEP}')
    print(f'  {title}')
    print(SEP)


# ── Checkpoint 1: EW input ────────────────────────────────────────────────────
section('CHECKPOINT 1: EW input')
ew_df = pd.read_csv(EW_CSV)
n_total  = len(ew_df)
n_blends = int((ew_df['blend_flag'] == True).sum()) if 'blend_flag' in ew_df.columns else 0
ew_clean = ew_df[ew_df['blend_flag'] == False] if 'blend_flag' in ew_df.columns else ew_df
ew_clean = ew_clean[(ew_clean['ew_mA'] >= 5) & (ew_clean['ew_mA'] <= 300)]
print(f'EW CSV: {n_total} rows, {n_blends} blend_flag=True, {len(ew_clean)} clean')

target_rows = ew_clean[
    (ew_clean['wavelength_air_A'] - TARGET_WAVE).abs() < TARGET_TOL
]
if target_rows.empty:
    print(f'ERROR: target line {TARGET_WAVE} Å not found in clean EW set — abort')
    sys.exit(1)

target_row = target_rows.iloc[0]
print(f'Target line: {target_row["wavelength_air_A"]:.4f} Å  '
      f'ion={target_row.get("ion","?")}  '
      f'EW={target_row["ew_mA"]:.2f} mÅ  '
      f'blend_flag={target_row.get("blend_flag", "N/A")}')


# ── Checkpoint 2: Model atmosphere ───────────────────────────────────────────
section('CHECKPOINT 2: Model atmosphere')
print(f'Grid: {ATM_GRID}')
print(f'Requested: Teff={STAR["teff_K"]} K  logg={STAR["logg"]}  '
      f'[Fe/H]={STAR["feh"]}  vmic={STAR["vturb_kms"]} km/s')

from pipeline.abundances_derive import _load_atmosphere

pack = ispec.load_modeled_layers_pack(ATM_GRID)
target = {'teff': STAR['teff_K'], 'logg': STAR['logg'], 'MH': STAR['feh'], 'alpha': 0.0}
print(f'Grid pack loaded: {ATM_GRID}')
print(f'Grid target dict: {target}')
in_grid = ispec.valid_atmosphere_target(pack, target)
print(f'ispec.valid_atmosphere_target → {in_grid}')
if not in_grid:
    print('WARNING: Procyon params are OUTSIDE ATLAS9 grid — atmosphere cannot be interpolated!')
else:
    print('Params are within grid bounds.')

with warnings.catch_warnings(record=True) as caught_warnings:
    warnings.simplefilter('always')
    atmosphere = ispec.interpolate_atmosphere_layers(pack, target, code='spectrum')

if caught_warnings:
    print(f'WARNINGS from atmosphere load ({len(caught_warnings)}):')
    for w in caught_warnings:
        print(f'  [{w.category.__name__}] {w.message}')
else:
    print('No warnings from atmosphere interpolation.')

print(f'Atmosphere array shape: {atmosphere.shape}  dtype: {atmosphere.dtype}')
print(f'First layer: {atmosphere[0]}')
print(f'Last  layer: {atmosphere[-1]}')


# ── Checkpoint 3: Build single-line linemask ─────────────────────────────────
section('CHECKPOINT 3: Single-line linemask for ispec.determine_abundances')
solar_abund = ispec.read_solar_abundances(SOLAR_ABUND_FILE)

# Build a linemask for just this one line, using the full linelist
all_regions = ispec.read_line_regions(LINELIST_FILE)
print(f'Full linelist: {len(all_regions)} entries')

match_mask = np.abs(all_regions['wave_A'] - target_row['wavelength_air_A']) < TARGET_TOL
matches = all_regions[match_mask]
if len(matches) == 0:
    print(f'WARNING: target line {target_row["wavelength_air_A"]:.4f} Å '
          f'not found in linelist — cannot build single-line mask')
    single_mask = None
else:
    print(f'Linelist match: {len(matches)} entries near {target_row["wavelength_air_A"]:.4f} Å')
    for m in matches:
        print(f'  wave={m["wave_A"]:.4f}  note={m["note"]}  '
              f'loggf={m["loggf"]:.3f}  EP={m["lower_state_eV"]:.3f}')
    single_mask = matches[:1].copy()
    # Inject measured EW
    single_mask['ew'] = float(target_row['ew_mA'])
    print(f'\nSingle-line mask EW injected: {single_mask["ew"][0]:.2f} mÅ')


# ── Checkpoint 4: ispec.determine_abundances call ────────────────────────────
section('CHECKPOINT 4: ispec.determine_abundances — verbose=1')
if single_mask is None:
    print('Skipping — no linelist match.')
else:
    print(f'Calling determine_abundances with:')
    print(f'  teff={STAR["teff_K"]}  logg={STAR["logg"]}  feh={STAR["feh"]}  '
          f'vmic={STAR["vturb_kms"]}  code=spectrum  verbose=1')
    print()
    with warnings.catch_warnings(record=True) as abund_warnings:
        warnings.simplefilter('always')
        spec_abund, normal_abund, x_over_h, x_over_fe = ispec.determine_abundances(
            atmosphere,
            STAR['teff_K'], STAR['logg'], STAR['feh'], 0.0,
            single_mask, solar_abund,
            microturbulence_vel=STAR['vturb_kms'],
            verbose=1,
            code='spectrum',
            tmp_dir='/tmp/ispec_codex',
        )
    print()
    if abund_warnings:
        print(f'WARNINGS from determine_abundances ({len(abund_warnings)}):')
        for w in abund_warnings:
            print(f'  [{w.category.__name__}] {w.message}')

    print(f'\nRaw return values:')
    print(f'  spec_abund   = {spec_abund}   (SPECTRUM internal scale)')
    print(f'  normal_abund = {normal_abund}  (A(X), log N/N_H + 12 scale)')
    print(f'  x_over_h     = {x_over_h}     ([X/H])')
    print(f'  x_over_fe    = {x_over_fe}    ([X/Fe])')


# ── Checkpoint 5: What does the full pipeline call return? ───────────────────
section('CHECKPOINT 5: Full pipeline _ew_to_abundance for Fe I subset')
# Import the pipeline function directly
from pipeline.abundances_derive import _build_ispec_line_regions, _ew_to_abundance

# Use a small subset of clean Fe I lines only
fe1_clean = ew_clean[
    (ew_clean.get('ion', pd.Series('', index=ew_clean.index)) == 'I') |
    (ew_clean.get('ion', pd.Series('', index=ew_clean.index)) == '')
].copy()

# If 'ion' column not present, filter by looking for Fe I in the linelist
if 'ion' not in fe1_clean.columns:
    print('No ion column — using all clean lines')
    fe1_clean = ew_clean.copy()

stellar_params = {
    'teff_K'   : STAR['teff_K'],
    'logg'     : STAR['logg'],
    'feh'      : STAR['feh'],
    'vturb_kms': STAR['vturb_kms'],
}
print(f'Calling _ew_to_abundance with {len(ew_clean)} clean lines...')
linemasks_out, normal_abund_out, xh_out, xfe_out = _ew_to_abundance(
    ew_clean, stellar_params, atmosphere
)

print(f'\nResult: {len(linemasks_out)} lines returned')
notes = np.array([str(n) for n in linemasks_out['note']])
fe1_idx = np.where(notes == 'Fe 1')[0]
print(f'Fe I lines: {len(fe1_idx)}')
print(f'\nFirst 10 Fe I normal_abund values (A(Fe) absolute):')
for idx in fe1_idx[:10]:
    print(f'  wave={linemasks_out["wave_A"][idx]:.4f}  '
          f'EW={linemasks_out["ew"][idx]:.2f} mÅ  '
          f'normal_abund={normal_abund_out[idx]:.6f}  '
          f'x_over_h={xh_out[idx]:.6f}')

all_fe1_abund = normal_abund_out[fe1_idx]
valid_fe1 = all_fe1_abund[np.isfinite(all_fe1_abund)]
print(f'\nFe I normal_abund stats:')
print(f'  min={np.nanmin(all_fe1_abund):.4f}  max={np.nanmax(all_fe1_abund):.4f}  '
      f'mean={np.nanmean(all_fe1_abund):.4f}  std={np.nanstd(all_fe1_abund):.4f}')
print(f'  n_zero={int((all_fe1_abund == 0.0).sum())}  '
      f'n_finite={int(np.isfinite(all_fe1_abund).sum())}  '
      f'n_nan={int(np.isnan(all_fe1_abund).sum())}')

# ── Checkpoint 6: Solar comparison ───────────────────────────────────────────
section('CHECKPOINT 6: Solar comparison — same call with solar params')
SOLAR_PARAMS = dict(teff_K=5771.0, logg=4.44, feh=0.0, vturb_kms=0.87)
sol_target = {'teff': SOLAR_PARAMS['teff_K'], 'logg': SOLAR_PARAMS['logg'],
              'MH': SOLAR_PARAMS['feh'], 'alpha': 0.0}
atm_solar = ispec.interpolate_atmosphere_layers(pack, sol_target, code='spectrum')
solar_ew_csv = os.path.join(REPO, 'data', 'processed', 'solar_per_line.csv')
if os.path.exists(solar_ew_csv):
    solar_df = pd.read_csv(solar_ew_csv)
    # Check if it has EW data or just abundances
    if 'ew_mA' in solar_df.columns:
        solar_target = solar_df[
            (solar_df['wavelength_air_A'] - TARGET_WAVE).abs() < TARGET_TOL
        ]
        if not solar_target.empty:
            print(f'Solar EW for {TARGET_WAVE} Å: {solar_target.iloc[0]["ew_mA"]:.2f} mÅ')
        else:
            print(f'Target line {TARGET_WAVE} Å not in solar per-line CSV')
    else:
        print(f'solar_per_line.csv has no ew_mA column — columns: {list(solar_df.columns)[:8]}')
else:
    print(f'solar_per_line.csv not found at {solar_ew_csv}')

solar_ew_full = os.path.join(REPO, 'data', 'processed', 'solar_ew.csv')
if os.path.exists(solar_ew_full):
    sol_df = pd.read_csv(solar_ew_full)
    sol_clean = sol_df.copy()
    if 'blend_flag' in sol_clean.columns:
        sol_clean = sol_clean[sol_clean['blend_flag'] == False]
    sol_clean = sol_clean[(sol_clean['ew_mA'] >= 5) & (sol_clean['ew_mA'] <= 300)]
    print(f'Solar EW CSV: {len(sol_df)} total → {len(sol_clean)} clean Fe lines')
    lm_sol, na_sol, xh_sol, xfe_sol = _ew_to_abundance(
        sol_clean, SOLAR_PARAMS, atm_solar
    )
    notes_sol = np.array([str(n) for n in lm_sol['note']])
    fe1_sol = np.where(notes_sol == 'Fe 1')[0]
    sol_abunds = na_sol[fe1_sol]
    print(f'\nSolar Fe I normal_abund first 5:')
    for idx in fe1_sol[:5]:
        print(f'  wave={lm_sol["wave_A"][idx]:.4f}  normal_abund={na_sol[idx]:.6f}')
    print(f'Solar stats: min={np.nanmin(sol_abunds):.4f}  max={np.nanmax(sol_abunds):.4f}  '
          f'mean={np.nanmean(sol_abunds):.4f}  std={np.nanstd(sol_abunds):.4f}')
    print(f'n_zero={int((sol_abunds==0.0).sum())}  n_finite={int(np.isfinite(sol_abunds).sum())}')
else:
    print(f'solar_ew.csv not at {solar_ew_full} — skipping solar comparison')

section('DIAGNOSTIC COMPLETE')
print('Check checkpoint 4 and 5 — if normal_abund=0.0 for all Procyon Fe I:')
print('  Hypothesis 1 (most likely): SPECTRUM bisection converges to [Fe/H]=0.0')
print('     because feh=-0.04 is being passed as initial abundance guess,')
print('     which rounds to 0.0 in relative units (solar - solar = 0).')
print('  Hypothesis 2: SPECTRUM not found / silent failure at Teff=6530K.')
print('  Hypothesis 3: Atmosphere interpolation gives boundary model.')
