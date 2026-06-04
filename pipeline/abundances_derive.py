"""
pipeline/abundances_derive.py
=============================
Derive stellar elemental abundances from measured EWs using iSpec +
Turbospectrum radiative transfer.

Method: 1D LTE (upgradeable to 1D NLTE via Turbospectrum — RYA-165)
Model atmospheres:
  ATLAS9.Castelli — FGK dwarfs (Sun, 55 Cnc A, Alpha Cen A)
  MARCS.GES       — M dwarfs + evolved stars (LHS 3844, Barnard's Star)
Radiative transfer: Turbospectrum (Alvarez & Plez 1998; Plez 2012; Gerber+2023)

Why Turbospectrum:
  - Correct Mg/Si (SPECTRUM has 0.10-0.15 dex bias on some lines, Gaia FGK 2026)
  - Native molecular lines (critical for C and M dwarf targets)
  - Native NLTE path — same code handles LTE → NLTE → 3D NLTE (RYA-165)

EW injection workflow:
  Our solar_ew.csv EWs are matched to iSpec's pre-built GES line regions
  (which carry full atomic data: loggf, EP, Turbospectrum species codes,
  damping constants). Measured EWs replace the region placeholders and
  determine_abundances iterates the abundance A(X) until EW matches.

Solar abundances reference: Asplund.2009 (closest available in iSpec;
  Asplund.2021 values from config/constants.py SOLAR_ASPLUND2021 are used
  for [X/H] reporting).

Linear issue: RYA-167

Calibration notes (first solar run, 2026-06-04)
------------------------------------------------
- ATLAS9.Castelli + Turbospectrum verified working at solar params
- Asplund.2021 not shipped with iSpec; Asplund.2009 used as iSpec reference
- EP slope and REW slope computed for Fe I as stellar-param diagnostics
"""

import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

from config.constants import (
    PATHS, SOLAR_ASPLUND2021, STAR_SOLAR, STAR_55CNC,
    ISPEC_DIR, RADIATIVE_TRANSFER_CODE,
)

# ── iSpec bootstrap ───────────────────────────────────────────────────────────
_ispec_python = ISPEC_DIR / '.venv' / 'bin' / 'python3.12'
sys.path.insert(0, str(ISPEC_DIR))
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import ispec

# ── Grid / file constants ─────────────────────────────────────────────────────
_GRID_DIRS = {
    'ATLAS9': 'ATLAS9.Castelli',
    'MARCS':  'MARCS.GES',
}

_SOLAR_ABUND_FILE = str(
    ISPEC_DIR / 'input' / 'abundances' / 'Asplund.2009' / 'stdatom.dat'
)

# Pre-built line regions with full atomic data for Turbospectrum synthesis.
# The "limited_but_with_missing_elements" file covers the widest element range.
_LINE_REGIONS_TS_ALL = str(
    ISPEC_DIR / 'input' / 'regions' / '42000_GES' /
    'limited_but_with_missing_elements_turbospectrum_synth_good_for_abundances_all_extended.txt'
)
# Fe-only parameter determination regions (higher quality selection)
_LINE_REGIONS_TS_FE = str(
    ISPEC_DIR / 'input' / 'regions' / '42000_GES' /
    'turbospectrum_synth_good_for_params_all_extended.txt'
)


# ── Atmosphere helpers ────────────────────────────────────────────────────────

def _load_pack_and_solar(model_grid: str = 'ATLAS9'):
    """Load model atmosphere pack + iSpec solar abundances reference."""
    grid_dir = ISPEC_DIR / 'input' / 'atmospheres' / _GRID_DIRS[model_grid]
    if not grid_dir.exists():
        raise FileNotFoundError(
            f"Atmosphere grid not found: {grid_dir}\n"
            "Run: cd ispec && tar -xzf input.tar.gz"
        )
    pack = ispec.load_modeled_layers_pack(str(grid_dir))
    solar_abund = ispec.read_solar_abundances(_SOLAR_ABUND_FILE)
    return pack, solar_abund


def _interpolate_atmosphere(pack, teff, logg, feh, alpha=0.0,
                             code=RADIATIVE_TRANSFER_CODE):
    """Interpolate model atmosphere at stellar parameters."""
    target = {'teff': float(teff), 'logg': float(logg),
              'MH': float(feh), 'alpha': float(alpha)}
    if not ispec.valid_atmosphere_target(pack, target):
        raise ValueError(
            f"Stellar params {target} fall outside the atmosphere grid. "
            "Check Teff, logg, [M/H] are within grid bounds."
        )
    return ispec.interpolate_atmosphere_layers(pack, target, code=code)


# ── EW injection ──────────────────────────────────────────────────────────────

def _ours_to_ispec_note(element: str, ion: str) -> str:
    """Convert our 'Fe','I' notation to iSpec 'Fe 1' notation."""
    digit = 1 if ion.strip().upper() == 'I' else 2
    return f"{element} {digit}"


def _inject_ew_into_regions(ew_df: pd.DataFrame,
                             line_regions_path: str,
                             wave_tol: float = 0.15) -> np.ndarray:
    """
    Match our measured EWs to iSpec's pre-built line regions (which carry
    full atomic data) and inject the EW values.

    Returns a structured numpy array (iSpec linemasks dtype) containing
    only the lines for which we have a measured EW.

    Parameters
    ----------
    ew_df       : our solar_ew.csv DataFrame
    line_regions_path : path to iSpec extended line regions file
    wave_tol    : wavelength match tolerance in Å (default 0.15 Å)
    """
    regions = ispec.read_line_regions(line_regions_path)

    # Pre-build lookup: ispec_note → list of (wavelength_A, ew_A, ew_err_A)
    ew_lookup: dict[str, list] = {}
    for _, row in ew_df.iterrows():
        note = _ours_to_ispec_note(str(row['element']), str(row['ion']))
        wav  = float(row['wavelength_air_A'])
        ew_A = float(row['ew_mA']) / 1000.0
        err_A = float(row.get('ew_err_mA', 0.0)) / 1000.0 if 'ew_err_mA' in row.index else 0.0
        ew_lookup.setdefault(note, []).append((wav, ew_A, err_A))

    # Match regions to our EWs
    matched_indices = []
    ew_values  = []
    err_values = []

    for i, region in enumerate(regions):
        note = str(region['note'])
        if note not in ew_lookup:
            continue
        candidates = ew_lookup[note]
        diffs = [abs(wav - float(region['wave_A'])) for wav, _, _ in candidates]
        min_diff = min(diffs)
        if min_diff > wave_tol:
            continue
        best = candidates[int(np.argmin(diffs))]
        matched_indices.append(i)
        ew_values.append(best[1])
        err_values.append(best[2])

    if not matched_indices:
        raise RuntimeError(
            f"No EW measurements matched iSpec line regions "
            f"(tolerance ±{wave_tol} Å). Check element/ion names and wavelengths."
        )

    # Build matched linemasks with injected EWs
    matched = regions[matched_indices].copy()
    for j, (ew_A, err_A) in enumerate(zip(ew_values, err_values)):
        matched[j]['ew']     = ew_A
        matched[j]['ew_err'] = err_A

    return matched


# ── Abundance diagnostics ─────────────────────────────────────────────────────

def _ep_slope(fe1_mask, linemasks, x_over_h):
    """Fe I abundance vs EP slope — diagnostic for Teff."""
    ep  = linemasks['lower_state_eV'][fe1_mask]
    abund = x_over_h[fe1_mask]
    valid = np.isfinite(abund) & np.isfinite(ep)
    if valid.sum() < 5:
        return np.nan
    return float(np.polyfit(ep[valid], abund[valid], 1)[0])


def _rew_slope(fe1_mask, linemasks, x_over_h):
    """Fe I abundance vs reduced EW slope — diagnostic for vturb."""
    ew   = linemasks['ew'][fe1_mask]
    wav  = linemasks['wave_A'][fe1_mask]
    abund = x_over_h[fe1_mask]
    valid = np.isfinite(abund) & (ew > 0) & (wav > 0)
    if valid.sum() < 5:
        return np.nan
    rew = np.log10(ew[valid] / wav[valid])
    return float(np.polyfit(rew, abund[valid], 1)[0])


# ── Core abundance derivation ─────────────────────────────────────────────────

def derive_abundances_ispec(ew_df: pd.DataFrame,
                             stellar_params: dict,
                             model_grid: str = 'ATLAS9',
                             code: str = RADIATIVE_TRANSFER_CODE,
                             line_regions_path: str = _LINE_REGIONS_TS_ALL,
                             wave_tol: float = 0.15) -> pd.DataFrame:
    """
    Derive per-element abundances from measured EWs via iSpec + Turbospectrum.

    Parameters
    ----------
    ew_df            : DataFrame from solar_ew.csv
    stellar_params   : dict with keys teff_K, logg, feh, vturb_kms[, alpha]
    model_grid       : 'ATLAS9' (FGK) or 'MARCS' (M dwarfs)
    code             : radiative transfer code (default: RADIATIVE_TRANSFER_CODE)
    line_regions_path: path to iSpec extended line regions file
    wave_tol         : EW wavelength match tolerance in Å

    Returns
    -------
    DataFrame with columns: element, ion, A_X, A_X_std, n_lines, XH, XFe
    """
    teff  = float(stellar_params['teff_K'])
    logg  = float(stellar_params['logg'])
    feh   = float(stellar_params['feh'])
    vturb = float(stellar_params['vturb_kms'])
    alpha = float(stellar_params.get('alpha', 0.0))

    pack, solar_abund = _load_pack_and_solar(model_grid)
    atm = _interpolate_atmosphere(pack, teff, logg, feh, alpha, code=code)

    linemasks = _inject_ew_into_regions(ew_df, line_regions_path, wave_tol)
    print(f"  Matched {len(linemasks)} EW lines to iSpec regions")

    # Discard zero/negative EWs (unmatched placeholders)
    linemasks = linemasks[linemasks['ew'] > 0]
    print(f"  After EW > 0 filter: {len(linemasks)} lines")

    # ── Call Turbospectrum ────────────────────────────────────────────────────
    spec_abund, normal_abund, x_over_h, x_over_fe = ispec.determine_abundances(
        atm, teff, logg, feh, alpha, linemasks, solar_abund,
        microturbulence_vel=vturb,
        verbose=0,
        code=code,
        tmp_dir='/tmp/ispec_codex',
    )

    # ── Per-element aggregation ───────────────────────────────────────────────
    notes  = np.array([str(n) for n in linemasks['note']])
    a_fe_solar = SOLAR_ASPLUND2021['Fe']

    rows = []
    for note in np.unique(notes):
        mask = notes == note
        xh   = x_over_h[mask]
        xfe  = x_over_fe[mask]
        ab   = spec_abund[mask]
        valid = np.isfinite(xh)
        if valid.sum() == 0:
            continue

        parts = note.split()
        elem = parts[0]
        ion  = 'I' if parts[1] == '1' else 'II'

        a_x_mean = float(np.nanmedian(ab[valid]))
        xh_mean  = float(np.nanmedian(xh[valid]))
        xfe_mean = float(np.nanmedian(xfe[valid]))
        xh_std   = float(np.nanstd(xh[valid])) if valid.sum() > 1 else np.nan

        rows.append({
            'element'  : elem,
            'ion'      : ion,
            'A_X'      : round(a_x_mean, 3),
            'A_X_std'  : round(xh_std, 3),
            'n_lines'  : int(valid.sum()),
            'XH'       : round(xh_mean, 3),
            'XFe'      : round(xfe_mean, 3),
        })

    results = pd.DataFrame(rows).sort_values(['element', 'ion']).reset_index(drop=True)

    # ── Fe diagnostics ────────────────────────────────────────────────────────
    fe1_mask = notes == 'Fe 1'
    fe2_mask = notes == 'Fe 2'
    ep_sl  = _ep_slope(fe1_mask, linemasks, x_over_h)
    rew_sl = _rew_slope(fe1_mask, linemasks, x_over_h)

    fe1_xh = float(np.nanmedian(x_over_h[fe1_mask & np.isfinite(x_over_h)])) if fe1_mask.any() else np.nan
    fe2_xh = float(np.nanmedian(x_over_h[fe2_mask & np.isfinite(x_over_h)])) if fe2_mask.any() else np.nan
    ion_bal = fe1_xh - fe2_xh if np.isfinite(fe1_xh) and np.isfinite(fe2_xh) else np.nan

    print(f"\n  ── Fe diagnostics ──")
    print(f"  EP slope (Teff):   {ep_sl:+.4f} dex/eV  (target |slope| < 0.02)")
    print(f"  REW slope (vturb): {rew_sl:+.4f}         (target flat)")
    print(f"  Ionisation bal.:   {ion_bal:+.4f} dex    (target |diff| < 0.05)")

    return results


# ── Legacy COG stub (preserved for reference — DO NOT USE for science) ────────

def _derive_abundances_cog_legacy(ew_df, stellar_params):
    """
    DEPRECATED — linear COG calibrated to Fe I only.
    Errors of ±0.5–2 dex for non-Fe elements; Fe I excitation slope ~0.10 dex/eV.
    Retained for diagnostic comparison only. Use derive_abundances_ispec() instead.
    """
    raise DeprecationWarning(
        "_derive_abundances_cog_legacy is retired. Use derive_abundances_ispec()."
    )


# ── Public run() entry point ──────────────────────────────────────────────────

def run(star_id: str = 'solar',
        model_grid: str = 'ATLAS9',
        stellar_params_override: dict = None) -> pd.DataFrame:
    """
    Derive abundances for a star and save results.

    Parameters
    ----------
    star_id  : 'solar' or '55cnc'
    model_grid: 'ATLAS9' (default, FGK) or 'MARCS' (M dwarfs)
    stellar_params_override : override dict for uncertainty sensitivity runs (RYA-158)

    Returns
    -------
    DataFrame with element, ion, A_X, A_X_std, n_lines, XH, XFe
    """
    print(f"\n{'='*60}")
    print(f"  abundances_derive — {star_id}  |  {model_grid}  |  {RADIATIVE_TRANSFER_CODE}")
    print(f"{'='*60}\n")

    # ── Stellar parameters ────────────────────────────────────────────────────
    if stellar_params_override:
        params = stellar_params_override
    elif 'solar' in star_id.lower():
        params = STAR_SOLAR.copy()
    else:
        params = STAR_55CNC.copy()

    print(f"[1/4] Stellar params:  Teff={params['teff_K']} K  "
          f"logg={params['logg']}  vturb={params['vturb_kms']} km/s  "
          f"[Fe/H]={params['feh']}")

    # ── Load EW table ─────────────────────────────────────────────────────────
    if 'solar' in star_id.lower():
        ew_path = PATHS['solar_ew']
    else:
        ew_path = PATHS.get(f'{star_id}_ew',
                             Path(PATHS['solar_ew']).parent / f'{star_id}_ew.csv')

    if not Path(str(ew_path)).exists():
        raise FileNotFoundError(
            f"EW table not found: {ew_path}\n"
            "Run pipeline/lines_fit.py first."
        )

    ew_df = pd.read_csv(ew_path)
    ew_df = ew_df[(ew_df['ew_mA'] > 0) & ew_df['ew_mA'].notna()].copy()
    print(f"[2/4] Loaded {len(ew_df)} EW measurements from {Path(str(ew_path)).name}")

    # ── Derive abundances ─────────────────────────────────────────────────────
    print(f"[3/4] Running iSpec + Turbospectrum ({model_grid})...")
    results = derive_abundances_ispec(ew_df, params, model_grid=model_grid)

    # ── Save ──────────────────────────────────────────────────────────────────
    out_path = Path(str(ew_path)).parent / f'{star_id}_abundances.csv'
    results.to_csv(out_path, index=False)
    print(f"\n[4/4] Saved → {out_path.name}  ({len(results)} elements)")

    # ── Solar calibration gate check ──────────────────────────────────────────
    fe_rows = results[results['element'] == 'Fe']
    if len(fe_rows) > 0:
        a_fe = float(fe_rows['A_X'].mean())
        delta = a_fe - 7.46
        status = '✓ PASS' if abs(delta) <= 0.05 else '✗ FAIL'
        print(f"\n  A(Fe)☉ = {a_fe:.3f}  (target 7.46 ± 0.05)  {status}")

    print(f"\n✓  abundances_derive complete.\n")
    return results


if __name__ == '__main__':
    import sys as _sys
    star = _sys.argv[1] if len(_sys.argv) > 1 else 'solar'
    grid = _sys.argv[2] if len(_sys.argv) > 2 else 'ATLAS9'
    run(star, grid)
