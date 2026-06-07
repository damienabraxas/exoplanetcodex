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


def _in_grid(teff: float, logg: float, afe3n: float, vmic: float) -> bool:
    return (
        _GRID['teff'][0] <= teff  <= _GRID['teff'][1] and
        _GRID['logg'][0] <= logg  <= _GRID['logg'][1] and
        _GRID['vmic'][0] <= vmic  <= _GRID['vmic'][1] and
        _GRID['afe'][0]  <= afe3n <= _GRID['afe'][1]
    )


# ── Per-line correction ───────────────────────────────────────────────────────

def _compute_aberr(ion: str, elo: float, eup: float, loggf: float,
                   teff: float, logg: float, afe3n: float, vmic: float) -> float:
    """
    Compute NLTE correction for one Fe I/II line.

    Returns aberr = A(Fe;1D) − A(Fe;3D) in dex.
    NaN if outside grid or model unavailable.
    """
    if not _in_grid(teff, logg, afe3n, vmic):
        return np.nan
    try:
        models = _get_models()
        X = np.array([[teff, logg, afe3n, vmic, elo, eup, loggf]])
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

def apply_fe_nlte_corrections(
    abundances_df: pd.DataFrame,
    stellar_params: dict,
    line_df: pd.DataFrame = None,
    wave_tol: float = 0.15,
) -> pd.DataFrame:
    """
    Apply Amarsi 2022 Fe I/II NLTE corrections to 1D LTE abundances.

    Parameters
    ----------
    abundances_df  : per-element results from abundances_derive.run()
                     must have columns: element, ion, A_X, n_lines
    stellar_params : dict with keys teff_K, logg, feh, vturb_kms
    line_df        : EW DataFrame (element, ion, wavelength_air_A)
                     used to look up per-line atomic data; if None,
                     applies mean correction over all grid Fe lines
    wave_tol       : wavelength match tolerance in Å (default 0.15)

    Returns
    -------
    DataFrame with additional columns:
        delta_nlte_mean : mean aberr across matched lines (dex)
        n_nlte_lines    : number of lines with valid correction
        A_X_nlte        : A_X − delta_nlte_mean  (3D NLTE abundance)
        nlte_flag       : '3D_NLTE_Amarsi2022' | 'NLTE_unavailable' | '1D_LTE'
        nlte_ref        : citation string
    """
    teff = float(stellar_params.get('teff_K', 5777))
    logg = float(stellar_params.get('logg',   4.44))
    vmic = float(stellar_params.get('vturb_kms', 1.0))

    regions = _load_nlte_regions()

    result = abundances_df.copy()
    result['delta_nlte_mean'] = np.nan
    result['n_nlte_lines']    = 0
    result['A_X_nlte']        = np.nan
    result['nlte_flag']       = '1D_LTE'
    result['nlte_ref']        = ''

    for idx, row in result.iterrows():
        if str(row['element']) != 'Fe':
            result.at[idx, 'A_X_nlte'] = float(row['A_X'])
            continue

        ion     = str(row['ion'])
        a_1dlte = float(row['A_X'])
        afe3n_guess = a_1dlte  # initial guess for iterative convergence

        aberrs = []

        # Initial A(Fe;3N) guess: clamp 1D LTE value to grid ceiling [4.5, 7.5]
        afe3n_start = float(np.clip(a_1dlte, _GRID['afe'][0], _GRID['afe'][1]))

        if line_df is not None and not regions.empty:
            # Per-line corrections using our measured EW wavelengths
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
                elo   = float(best['elo_eV'])
                eup   = float(best['eup_eV'])
                lggf  = float(best['loggf'])

                # Iterate A(Fe;3N) to convergence (2-3 passes)
                # Convention: A(Fe;NLTE) = A(Fe;1D) + aberr (aberr is A_3D − A_1D)
                a3n = afe3n_start
                ab  = np.nan
                for _ in range(3):
                    ab_new = _compute_aberr(ion, elo, eup, lggf, teff, logg, a3n, vmic)
                    if not np.isfinite(ab_new):
                        break
                    ab = ab_new
                    a3n_new = float(np.clip(a_1dlte + ab, _GRID['afe'][0], _GRID['afe'][1]))
                    if abs(a3n_new - a3n) < 0.001:
                        a3n = a3n_new
                        break
                    a3n = a3n_new

                if np.isfinite(ab):
                    aberrs.append(ab)

        elif not regions.empty:
            # Fallback: mean over all grid Fe lines for this ion
            a3n = afe3n_start
            ion_regions = regions[regions['ion'] == ion]
            for _, r in ion_regions.iterrows():
                ab = _compute_aberr(
                    ion, float(r['elo_eV']), float(r['eup_eV']), float(r['loggf']),
                    teff, logg, a3n, vmic
                )
                if np.isfinite(ab):
                    aberrs.append(ab)

        if aberrs:
            mean_delta = float(np.mean(aberrs))
            n_lines    = len(aberrs)
            a_nlte = a_1dlte + mean_delta
            print(f"  Fe {ion}: {n_lines} lines, mean Δ(NLTE) = {mean_delta:+.4f} dex  "
                  f"A(Fe;1D)={a_1dlte:.3f} → A(Fe;NLTE)={a_nlte:.3f}")
            result.at[idx, 'delta_nlte_mean'] = round(mean_delta, 4)
            result.at[idx, 'n_nlte_lines']    = n_lines
            result.at[idx, 'A_X_nlte']        = round(a_nlte, 3)
            result.at[idx, 'nlte_flag']        = '3D_NLTE_Amarsi2022'
            result.at[idx, 'nlte_ref']         = 'Amarsi+2022_A&A_668_A68'
        else:
            print(f"  Fe {ion}: no lines matched in NLTE grid — 1D LTE retained")
            result.at[idx, 'A_X_nlte']  = a_1dlte
            result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'

    return result
