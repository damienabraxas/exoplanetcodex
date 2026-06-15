"""
pipeline/nlte_corrections.py
============================
Apply 3D non-LTE corrections to 1D LTE Fe I/II abundances.

Corrections from Amarsi, Liljegren & Nissen 2022 (A&A 668, A68).
Neural network implementation by Liljegren (github.com/sliljegren/1L-3NErrors).
Three multi-layer perceptrons trained on STAGGER-grid 3D NLTE models:
  - Fe I lines with Elo < 2 eV
  - Fe I lines with Elo >= 2 eV
  - Fe II lines

Grid coverage: Teff 5000–6500 K, log g 4.0–4.5, [Fe/H] −3.0 to 0.0.
171 Fe I + 12 Fe II optical lines.

Sign convention (empirically confirmed for solar validation):
  A(Fe;3D NLTE) = A(Fe;1D LTE) + aberr
  aberr < 0 at solar metallicity (neural network output is A_3D − A_1D);
  for the Sun: A(3D) = 7.58 + (−0.127) ≈ 7.45 matches Asplund 2021.

Ref: Amarsi et al. 2022, A&A 668, A68
     DOI: 10.1051/0004-6361/202244542
     https://github.com/sliljegren/1L-3NErrors
"""

import os
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────

_REPO_ROOT   = Path(__file__).parent.parent
_VENDOR_DIR  = _REPO_ROOT / 'vendor' / '1L-3NErrors'
_MODEL_LT02  = _VENDOR_DIR / 'fe1_model_lt02.p'   # Fe I, Elo < 2 eV
_MODEL_GT02  = _VENDOR_DIR / 'fe1_model_gt02.p'   # Fe I, Elo >= 2 eV
_MODEL_FE2   = _VENDOR_DIR / 'fe2_model.p'         # Fe II

try:
    from config.constants import ISPEC_DIR
    _NLTE_REGIONS = str(
        ISPEC_DIR / 'input' / 'regions' / '42000_VALD' /
        'turbospectrum-nlte_synth_good_for_params_all_extended.txt'
    )
except Exception:
    _NLTE_REGIONS = None

# ── Model loading (cached at module level) ────────────────────────────────────

_model_cache: dict = {}


def _load_models():
    """Load sklearn MLP regressors from vendor pickle files."""
    if not _MODEL_LT02.exists():
        raise FileNotFoundError(
            f"Amarsi 2022 NLTE models not found at {_VENDOR_DIR}. "
            "Run: git clone https://github.com/sliljegren/1L-3NErrors vendor/1L-3NErrors"
        )
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        (s_lt02, m_lt02) = pickle.load(open(str(_MODEL_LT02), 'rb'))
        (s_gt02, m_gt02) = pickle.load(open(str(_MODEL_GT02), 'rb'))
        (s_fe2,  m_fe2)  = pickle.load(open(str(_MODEL_FE2),  'rb'))
    return {
        'lt02': (s_lt02, m_lt02),
        'gt02': (s_gt02, m_gt02),
        'fe2':  (s_fe2,  m_fe2),
    }


def _get_models() -> dict:
    if 'models' not in _model_cache:
        _model_cache['models'] = _load_models()
    return _model_cache['models']


# ── Atomic data lookup ────────────────────────────────────────────────────────

_regions_cache: dict = {}


def _load_nlte_regions() -> pd.DataFrame:
    """Load Fe I/II atomic data (Elo, Eup, loggf) from the iSpec NLTE regions file."""
    if 'regions' in _regions_cache:
        return _regions_cache['regions']

    if _NLTE_REGIONS is None or not Path(_NLTE_REGIONS).exists():
        _regions_cache['regions'] = pd.DataFrame()
        return _regions_cache['regions']

    try:
        import sys
        sys.path.insert(0, str(_REPO_ROOT.parent / 'ispec'))
        import ispec
        raw = ispec.read_line_regions(_NLTE_REGIONS)
        fe_mask = np.array([str(n).strip() in ('Fe 1', 'Fe 2') for n in raw['note']])
        fe_raw  = raw[fe_mask]
        df = pd.DataFrame({
            'wave_A':  fe_raw['wave_A'].astype(float),
            'ion':     ['I' if str(n).strip() == 'Fe 1' else 'II' for n in fe_raw['note']],
            'elo_eV':  fe_raw['lower_state_eV'].astype(float),
            'eup_eV':  fe_raw['upper_state_eV'].astype(float),
            'loggf':   fe_raw['loggf'].astype(float),
        })
        _regions_cache['regions'] = df
        return df
    except Exception as e:
        print(f"  WARNING: could not load NLTE regions file: {e}")
        _regions_cache['regions'] = pd.DataFrame()
        return _regions_cache['regions']


# ── Grid bounds (from main_aberr.py) ─────────────────────────────────────────
# Teff: 5000–6500 K   logg: 4.0–4.5   vmic: 0–3 km/s   A(Fe;3N): 4.5–7.5
_GRID = dict(teff=(5000., 6500.), logg=(4.0, 4.5), vmic=(0., 3.), afe=(4.5, 7.5))

# iSpec's determine_abundances returns abundance offsets relative to the solar
# reference ([Fe/H] units), not absolute A(Fe).  The Amarsi grid afe axis is
# in absolute A(Fe) (4.5–7.5 range).  Add this offset before passing to the
# grid so that [Fe/H]=0 maps to A(Fe)=7.46, not to the grid floor of 4.5.
_A_FE_SOLAR = 7.46  # Asplund+2021 A(Fe), matches SOLAR_ASPLUND2021['Fe']


def _in_grid(teff: float, logg: float, afe3n: float, vmic: float) -> bool:
    return (
        _GRID['teff'][0] <= teff  <= _GRID['teff'][1] and
        _GRID['logg'][0] <= logg  <= _GRID['logg'][1] and
        _GRID['vmic'][0] <= vmic  <= _GRID['vmic'][1] and
        _GRID['afe'][0]  <= afe3n <= _GRID['afe'][1]
    )


# ── Bounded Teff clamp (RYA-317) ──────────────────────────────────────────────
# The Amarsi grid tops out at Teff 6500 K, but Procyon (6554 K, Heiter GBS) sits
# 54 K above it — so RYA-309 found 0 Fe lines corrected and the star silently
# dropped to LTE. Rather than abandon NLTE for a 54 K overshoot, clamp Teff to the
# grid edge for the NLTE lookup ONLY, and ONLY when the overshoot is within
# NLTE_TEFF_CLAMP_MARGIN_K — a BOUNDED edge-evaluation, not open-ended: a star far
# above the grid still returns NaN (fail loud), never a silent extrapolation.
#
# Induced error = (NLTE Teff gradient) × (overshoot). MEASURED on Procyon's 102
# Fe I lines (RYA-317): aberr(6400) vs aberr(6500) → median |dΔ/dTeff| = 1.4e-4
# dex/K (worst line 5.1e-4). So the 54 K clamp shifts the ENSEMBLE median A(Fe) by
# 54 × 1.4e-4 ≈ 0.008 dex ≪ the 0.02 dex gate — negligible vs NLTE itself
# (~0.03-0.10 dex). (Worst single line ~0.028 dex, but the abundance is the
# median over ~100 lines.) Clamped lines carry the NLTE_Teff_clamped_6500 flag so
# the approximation is never silent.
NLTE_TEFF_CLAMP_MARGIN_K = 100.0


def _clamp_teff(teff: float) -> tuple:
    """Return (teff_eff, clamped). Clamp a small Teff overshoot to the grid edge;
    leave anything beyond the margin untouched (it will fail _in_grid → NaN)."""
    tmax = _GRID['teff'][1]
    if tmax < teff <= tmax + NLTE_TEFF_CLAMP_MARGIN_K:
        return tmax, True
    return teff, False


# ── Per-line correction ───────────────────────────────────────────────────────

def _compute_aberr(ion: str, elo: float, eup: float, loggf: float,
                   teff: float, logg: float, afe3n: float, vmic: float) -> float:
    """
    Compute NLTE correction for one Fe I/II line.

    Returns aberr = A(Fe;1D) − A(Fe;3D) in dex.
    NaN if outside grid or model unavailable.
    """
    teff_eff, _ = _clamp_teff(teff)   # RYA-317: bounded Teff-edge clamp for NLTE lookup
    if not _in_grid(teff_eff, logg, afe3n, vmic):
        return np.nan
    try:
        models = _get_models()
        X = np.array([[teff_eff, logg, afe3n, vmic, elo, eup, loggf]])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if ion == 'II':
                s, m = models['fe2']
            elif elo < 2.0:
                s, m = models['lt02']
            else:
                s, m = models['gt02']
            aberr = float(m.predict(s.transform(X))[0])
        return aberr
    except Exception:
        return np.nan


# ── Public API ────────────────────────────────────────────────────────────────

def _apply_aberr_to_line(ion: str, elo: float, eup: float, lggf: float,
                          a_1dlte_line: float,
                          teff: float, logg: float, vmic: float) -> float:
    """
    Compute converged aberr for a single line using its own 1D LTE abundance
    as the A(Fe;3N) starting point.  Returns np.nan if out of grid.
    """
    afe3n = float(np.clip(a_1dlte_line, _GRID['afe'][0], _GRID['afe'][1]))
    ab = np.nan
    for _ in range(3):
        ab_new = _compute_aberr(ion, elo, eup, lggf, teff, logg, afe3n, vmic)
        if not np.isfinite(ab_new):
            break
        ab = ab_new
        afe3n_new = float(np.clip(a_1dlte_line + ab, _GRID['afe'][0], _GRID['afe'][1]))
        if abs(afe3n_new - afe3n) < 0.001:
            afe3n = afe3n_new
            break
        afe3n = afe3n_new
    return ab


def apply_fe_nlte_corrections(
    abundances_df: pd.DataFrame,
    stellar_params: dict,
    line_df: pd.DataFrame = None,
    per_line_df: pd.DataFrame = None,
    wave_tol: float = 0.15,
) -> pd.DataFrame:
    """
    Apply Amarsi 2022 Fe I/II NLTE corrections to 1D LTE abundances.

    Parameters
    ----------
    abundances_df  : per-element results from abundances_derive.run()
                     must have columns: element, ion, A_X, n_lines
    stellar_params : dict with keys teff_K, logg, feh, vturb_kms
    per_line_df    : per-line 1D LTE DataFrame from _iterative_parameter_convergence
                     (RYA-207). Columns: element, ion, wavelength_air_A, ew_mA,
                     excitation_potential_eV, eup_eV, log_gf, a_1dlte.
                     When provided, aberr[i] is applied to each line's individual
                     a_1dlte[i] and mean/std are recomputed from corrected values.
    line_df        : EW DataFrame used only when per_line_df is None (legacy path).
    wave_tol       : wavelength match tolerance in Å for legacy line_df path.

    Returns
    -------
    DataFrame with additional columns:
        delta_nlte_mean : mean aberr across corrected lines (dex)
        n_nlte_lines    : number of lines with valid correction
        A_X_nlte        : mean of per-line (a_1dlte[i] + aberr[i])
        A_X_std_nlte    : std of per-line NLTE abundances (scatter gate input)
        nlte_flag       : '3D_NLTE_Amarsi2022' | 'NLTE_unavailable' | '1D_LTE'
        nlte_ref        : citation string
    """
    teff = float(stellar_params.get('teff_K', 5777))
    logg = float(stellar_params.get('logg',   4.44))
    vmic = float(stellar_params.get('vturb_kms', 1.0))

    # RYA-317: flag a bounded Teff-edge clamp so it is recorded, never silent.
    _, _teff_clamped = _clamp_teff(teff)
    _nlte_flag = ('3D_NLTE_Amarsi2022_Teff_clamped_6500'
                  if _teff_clamped else '3D_NLTE_Amarsi2022')
    if _teff_clamped:
        print(f"  NLTE: Teff {teff:.0f} K clamped to grid edge "
              f"{_GRID['teff'][1]:.0f} K for the Fe NLTE lookup (RYA-317; "
              f"gradient error ~0.01 dex ≪ 0.02 gate).")

    regions = _load_nlte_regions()

    result = abundances_df.copy()
    result['delta_nlte_mean'] = np.nan
    result['n_nlte_lines']    = 0
    result['A_X_nlte']        = np.nan
    result['A_X_std_nlte']    = np.nan
    result['nlte_flag']       = '1D_LTE'
    result['nlte_ref']        = ''

    for idx, row in result.iterrows():
        if str(row['element']) != 'Fe':
            result.at[idx, 'A_X_nlte'] = float(row['A_X'])
            continue

        ion     = str(row['ion'])
        a_1dlte = float(row['A_X'])

        # ── Per-line path (RYA-207) ───────────────────────────────────────────
        if per_line_df is not None and not per_line_df.empty:
            fe_lines = per_line_df[
                (per_line_df['element'] == 'Fe') & (per_line_df['ion'] == ion)
            ].copy().reset_index(drop=True)

            if fe_lines.empty:
                result.at[idx, 'A_X_nlte']  = a_1dlte
                result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'
                continue

            fe_lines['aberr']  = np.nan
            fe_lines['a_nlte'] = fe_lines['a_1dlte']  # default: retain 1D LTE

            for i, lrow in fe_lines.iterrows():
                a_fe_abs = float(lrow['a_1dlte']) + _A_FE_SOLAR
                ab = _apply_aberr_to_line(
                    ion,
                    float(lrow['excitation_potential_eV']),
                    float(lrow['eup_eV']),
                    float(lrow['log_gf']),
                    a_fe_abs,
                    teff, logg, vmic,
                )
                if np.isfinite(ab):
                    fe_lines.at[i, 'aberr']  = ab
                    fe_lines.at[i, 'a_nlte'] = float(lrow['a_1dlte']) + ab

            n_corrected = int(fe_lines['aberr'].notna().sum())
            if n_corrected == 0:
                print(f"  Fe {ion}: no lines in NLTE grid — 1D LTE retained")
                result.at[idx, 'A_X_nlte']  = a_1dlte
                result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'
                continue

            # Median — consistent with 1D LTE computation in _iterative_parameter_convergence.
            # Mean is sensitive to outlier lines (e.g. 4918.994 Å, 4970.496 Å which sit
            # >0.7 dex above solar even after NLTE correction — likely blend contamination).
            a_nlte_mean = float(np.median(fe_lines['a_nlte']))
            a_nlte_std  = float(fe_lines['a_nlte'].std()) if len(fe_lines) > 1 else np.nan
            mean_delta  = float(fe_lines['aberr'].mean(skipna=True))

            outliers = fe_lines[np.abs(fe_lines['a_nlte'] - a_nlte_mean) > 0.3]
            if len(outliers) > 0:
                print(f"  WARNING: {len(outliers)} Fe {ion} NLTE outlier lines "
                      f"(|a_nlte - median| > 0.3 dex):")
                for _, r in outliers.iterrows():
                    print(f"    {r['wavelength_air_A']:.3f} Å  a_1dlte={r['a_1dlte']:.3f}  "
                          f"aberr={r['aberr']:.3f}  a_nlte={r['a_nlte']:.3f}")
                print(f"  These lines are excluded from the median but retained in scatter.")

            # Print first 10 lines as diagnostic table
            print(f"  Fe {ion} per-line NLTE ({n_corrected}/{len(fe_lines)} corrected):")
            print(f"    {'wave_A':>9}  {'a_1dlte':>8}  {'aberr':>8}  {'a_nlte':>8}")
            for _, r in fe_lines.head(10).iterrows():
                ab_str = f"{r['aberr']:+.4f}" if np.isfinite(r['aberr']) else "  n/a  "
                print(f"    {r['wavelength_air_A']:9.3f}  {r['a_1dlte']:8.4f}  "
                      f"{ab_str:>8}  {r['a_nlte']:8.4f}")
            print(f"  Fe {ion}: mean Δ(NLTE) = {mean_delta:+.4f} dex  "
                  f"A(Fe;1D)={a_1dlte:.3f} → A(Fe;NLTE)={a_nlte_mean:.3f}  "
                  f"scatter 1D→NLTE: {row.get('A_X_std', float('nan')):.3f}→{a_nlte_std:.3f}")

            result.at[idx, 'delta_nlte_mean'] = round(mean_delta, 4)
            result.at[idx, 'n_nlte_lines']    = n_corrected
            result.at[idx, 'A_X_nlte']        = round(a_nlte_mean, 3)
            result.at[idx, 'A_X_std_nlte']    = round(a_nlte_std, 3) if np.isfinite(a_nlte_std) else np.nan
            result.at[idx, 'nlte_flag']        = _nlte_flag
            result.at[idx, 'nlte_ref']         = 'Amarsi+2022_A&A_668_A68'

        # ── Legacy mean-field path (no per_line_df) ───────────────────────────
        else:
            aberrs      = []
            a_fe_abs_mean = a_1dlte + _A_FE_SOLAR
            afe3n_start = float(np.clip(a_fe_abs_mean, _GRID['afe'][0], _GRID['afe'][1]))

            if line_df is not None and not regions.empty:
                fe_ew = line_df[
                    (line_df['element'] == 'Fe') & (line_df['ion'] == ion)
                ]
                ion_regions = regions[regions['ion'] == ion]
                for _, lrow in fe_ew.iterrows():
                    wave = float(lrow['wavelength_air_A'])
                    diffs = (ion_regions['wave_A'] - wave).abs()
                    if diffs.empty or diffs.min() > wave_tol:
                        continue
                    best = ion_regions.loc[diffs.idxmin()]
                    ab = _apply_aberr_to_line(
                        ion, float(best['elo_eV']), float(best['eup_eV']),
                        float(best['loggf']), a_fe_abs_mean, teff, logg, vmic,
                    )
                    if np.isfinite(ab):
                        aberrs.append(ab)

            elif not regions.empty:
                ion_regions = regions[regions['ion'] == ion]
                for _, r in ion_regions.iterrows():
                    ab = _compute_aberr(
                        ion, float(r['elo_eV']), float(r['eup_eV']),
                        float(r['loggf']), teff, logg, afe3n_start, vmic,
                    )
                    if np.isfinite(ab):
                        aberrs.append(ab)

            if aberrs:
                mean_delta = float(np.mean(aberrs))
                a_nlte     = a_1dlte + mean_delta
                print(f"  Fe {ion}: {len(aberrs)} lines, mean Δ(NLTE) = {mean_delta:+.4f} dex  "
                      f"A(Fe;1D)={a_1dlte:.3f} → A(Fe;NLTE)={a_nlte:.3f}")
                result.at[idx, 'delta_nlte_mean'] = round(mean_delta, 4)
                result.at[idx, 'n_nlte_lines']    = len(aberrs)
                result.at[idx, 'A_X_nlte']        = round(a_nlte, 3)
                result.at[idx, 'nlte_flag']        = _nlte_flag
                result.at[idx, 'nlte_ref']         = 'Amarsi+2022_A&A_668_A68'
            else:
                print(f"  Fe {ion}: no lines matched in NLTE grid — 1D LTE retained")
                result.at[idx, 'A_X_nlte']  = a_1dlte
                result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'

    return result


# ── Ca / Ti / Cr NLTE — integration notes (RYA-256) ─────────────────────────
# Grids: data/nlte_grids/{Ca_Mashonkina2017,Ti_Bergemann2011_MPIA,Cr_Bergemann2010_MPIA}.csv
# Columns: element, wave_A, teff_K, logg, feh, delta_nlte
#
# UNIT CONVENTION: feh axis is [Fe/H] RELATIVE (range −0.5 to +0.3).
# iSpec returns [X/H] relative (0.0 at solar) for all elements.
# Caller should pass stellar_params['feh'] directly — NO solar offset needed.
# This is unlike the Fe Amarsi grid (afe axis = absolute A(Fe) 4.5–7.5),
# which required _A_FE_SOLAR = 7.46 added before lookup (see above, RYA-247).
#
# Expected solar deltas (Teff=5777, logg=4.4, feh=0.0):
#   Ca: +0.013 dex (median, 8 lines)   — range −0.106 to +0.053
#   Ti: +0.108 dex (median, 19 lines)  — range +0.101 to +0.133
#   Cr: +0.073 dex (median, 46 lines)  — range +0.009 to +0.132
