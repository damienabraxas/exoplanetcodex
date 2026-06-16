"""
pipeline/nlte_corrections.py
============================
Apply non-LTE corrections to 1D LTE Fe I/II abundances.

LIVE Fe NLTE source (RYA-319): MPIA Spectrum Tools / Bergemann mafags-os **1D
NLTE** grid (data/nlte_grids/Fe_Bergemann_MPIA.csv) — covers the full benchmark
box (Procyon Teff, α Cen A A(Fe), α Cen B logg, 55 Cnc), interpolated per Fe
line over (teff,logg,feh) via _mpia_fe_delta / apply_fe_nlte_corrections.
The Amarsi 3D-NLTE MLP below (Teff≤6500 ceiling) is ARCHIVED — kept only for the
solar 3D-vs-1D cross-check (RYA-283). Trade: 3D→1D, ≲0.05 dex at metal-rich
targets (Amarsi 2022 Fig 7), accepted for single-methodology coverage.

Archived Amarsi reference (the MLP functions, no longer the live correction):
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


# ── Per-line correction ───────────────────────────────────────────────────────

def _compute_aberr(ion: str, elo: float, eup: float, loggf: float,
                   teff: float, logg: float, afe3n: float, vmic: float) -> float:
    """
    Compute the NLTE correction-to-ADD for one Fe I/II line: A(Fe;3D-NLTE) −
    A(Fe;1D-LTE), in dex, so A_NLTE = A_1D + correction (matches the module
    convention and the MPIA δ convention).

    RYA-339 sign fix: the trained network returns A(1D) − A(3D-NLTE) (e.g. ≈ −0.20
    for Procyon Fe I). Amarsi et al. 2022 (2209.13449v3, §6): for Procyon, 1D LTE is
    0.20 dex SMALLER than 3D-NLTE → the physical correction is +0.20 (a hot F dwarf
    overionizes Fe I, so NLTE needs MORE iron). The prior code returned the network
    value directly and the consumers added it (a_1dlte + aberr), flipping the sign →
    the spurious −0.25 per-line / −0.2 [Fe/H] (RYA-322/339). Negate the network
    output here — the output-convention adapter at the MLP call.

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
            aberr_raw = float(m.predict(s.transform(X))[0])   # network: A(1D) − A(3D)
        return -aberr_raw   # RYA-339: adapt to correction-to-add A(3D) − A(1D)
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
    # RYA-339 sign guard: in the overionization (hot) regime the Fe I NLTE correction
    # MUST be positive (NLTE needs more Fe than LTE; Amarsi 2022 §6: Procyon Fe I =
    # +0.20). Fail loud if a sign flip ever returns. Cooler stars (≲6000 K) can carry
    # small negative Fe I corrections, so only the hot regime is guarded.
    if ion == 'I' and np.isfinite(ab) and teff >= 6000.0 and ab < -0.02:
        raise ValueError(
            f"Fe I NLTE correction {ab:+.3f} is NEGATIVE at Teff={teff:.0f} K "
            f"(overionization regime) — expected positive (Amarsi 2022 §6: Procyon "
            f"Fe I = +0.20). A sign flip has returned (RYA-339).")
    return ab


# ── MPIA/Bergemann Fe NLTE grid (RYA-319) — the LIVE Fe NLTE source ───────────
# Full switch from the Amarsi MLP (3D NLTE, Teff≤6500 ceiling — degraded at the
# edge, gave Procyon Fe I −0.3) to the MPIA Spectrum Tools / Bergemann mafags-os
# 1D NLTE grid, which covers the full benchmark box (Procyon Teff, α Cen A A(Fe),
# α Cen B logg, 55 Cnc). Grid = data/nlte_grids/Fe_Bergemann_MPIA.csv, keyed
# (ion, wave_A, teff, logg, feh) → delta_nlte. For each Fe line we interpolate
# delta over (teff,logg,feh) at the line's own wavelength; outside the convex
# hull → NaN (in-grid handled naturally, no clamp). The Amarsi MLP functions
# above (_load_models/_compute_aberr/_apply_aberr_to_line/_in_grid) are ARCHIVED
# — retained ONLY for the solar 3D-vs-1D cross-check (RYA-283), not the live
# correction. Trade: 3D→1D NLTE, ≲0.05 dex at our metal-rich targets (Amarsi
# 2022 Fig 7), accepted for single-methodology coverage (RYA-319 decision).
_MPIA_FE_GRID = _REPO_ROOT / 'data' / 'nlte_grids' / 'Fe_Bergemann_MPIA.csv'
_mpia_cache: dict = {}


def _load_mpia_fe_grid() -> dict:
    """Build per-(ion, wave) LinearNDInterpolator over (teff, logg, feh)."""
    if 'interp' in _mpia_cache:
        return _mpia_cache
    from scipy.interpolate import LinearNDInterpolator
    df = pd.read_csv(_MPIA_FE_GRID)
    df = df[df['delta_nlte'].notna()]
    interp, waves = {}, {'I': [], 'II': []}
    for (ion, wave), g in df.groupby(['ion', 'wave_A']):
        pts = g[['teff_K', 'logg', 'feh']].values
        if len(pts) < 4:
            continue
        try:
            f = LinearNDInterpolator(pts, g['delta_nlte'].values)
        except Exception:
            continue
        interp[(str(ion), round(float(wave), 3))] = f
        waves[str(ion)].append(round(float(wave), 3))
    _mpia_cache['interp'] = interp
    _mpia_cache['waves'] = {k: np.array(sorted(v)) for k, v in waves.items()}
    _mpia_cache['bounds'] = {
        'teff': (float(df.teff_K.min()), float(df.teff_K.max())),
        'logg': (float(df.logg.min()), float(df.logg.max())),
        'feh':  (float(df.feh.min()),  float(df.feh.max())),
    }
    return _mpia_cache


def _mpia_fe_delta(ion: str, wave_A: float, teff: float, logg: float,
                   feh: float, tol: float = 0.15) -> float:
    """delta_nlte = A(Fe;1D-NLTE) − A(Fe;1D-LTE) for one Fe line, interpolated
    from the MPIA grid. NaN if the line has no wave match (not in grid) or the
    star is outside the (teff,logg,feh) hull."""
    c = _load_mpia_fe_grid()
    w = c['waves'].get(str(ion), np.array([]))
    if len(w) == 0:
        return np.nan
    j = int(np.abs(w - wave_A).argmin())
    if abs(float(w[j]) - wave_A) > tol:
        return np.nan
    f = c['interp'].get((str(ion), float(w[j])))
    return float(f([[teff, logg, feh]])[0]) if f is not None else np.nan


def fe_grid_in_bounds(teff: float, logg: float, feh: float) -> bool:
    """Preflight (RYA-318): is (teff,logg,feh) inside the MPIA Fe grid box?"""
    b = _load_mpia_fe_grid()['bounds']
    return (b['teff'][0] <= teff <= b['teff'][1] and
            b['logg'][0] <= logg <= b['logg'][1] and
            b['feh'][0]  <= feh  <= b['feh'][1])


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
    feh  = float(stellar_params.get('feh', 0.0))   # RYA-319: MPIA grid axis

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
                # RYA-319: per-line δ from the MPIA grid, interpolated over
                # (teff,logg,feh) at the line's own wavelength. a_1dlte is the
                # absolute 1D-LTE A(Fe) (normal_abund); a_nlte = a_1dlte + δ.
                ab = _mpia_fe_delta(ion, float(lrow['wavelength_air_A']),
                                    teff, logg, feh)
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
            result.at[idx, 'nlte_flag']        = 'NLTE_MPIA_Bergemann_1D'
            result.at[idx, 'nlte_ref']         = 'MPIA SpectrumTools / Bergemann mafags-os 1D NLTE (RYA-319)'

        # ── Legacy mean-field path (no per_line_df) ───────────────────────────
        else:
            # RYA-319: MPIA grid δ, per the star's measured Fe lines if available,
            # else a mean over all grid lines for this ion (generic fallback).
            aberrs = []
            if line_df is not None:
                fe_ew = line_df[
                    (line_df['element'] == 'Fe') & (line_df['ion'] == ion)
                ]
                for _, lrow in fe_ew.iterrows():
                    ab = _mpia_fe_delta(ion, float(lrow['wavelength_air_A']),
                                        teff, logg, feh)
                    if np.isfinite(ab):
                        aberrs.append(ab)
            else:
                for w in _load_mpia_fe_grid()['waves'].get(ion, np.array([])):
                    ab = _mpia_fe_delta(ion, float(w), teff, logg, feh)
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
                result.at[idx, 'nlte_flag']        = '3D_NLTE_Amarsi2022'
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
