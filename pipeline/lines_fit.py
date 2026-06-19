"""
pipeline/lines_fit.py
=====================
Measure equivalent widths (EWs) of absorption lines in the HARPS solar spectrum.

For each line: extract local window → re-normalise continuum → fit Voigt
(Gaussian fallback) → integrate profile → EW in mÅ → propagate uncertainty.

Special handling
----------------
O I  6300.304 : Si I 6299.6 contaminates ±2 Å window → narrow ±0.15 Å window;
                Ni I 6300.336 contribution predicted via COG from clean Ni I lines
                and subtracted (Allende Prieto et al. 2001 method).
Ba II 5853.668 : 23 HFS components unresolved at HARPS R; fit single profile = total EW
Eu II 6645.127 : 30 HFS components span 6645.07–6645.16 Å; fit centred on HFS
                centroid at 6645.127, ±0.30 Å narrow window to capture full spread
Li I  6707.840 : doublet midpoint (components at 6707.76 + 6707.91); ±0.25 Å narrow
                window avoids contamination at 6707.46; flagged CN_blend_possible
Na, O, C, Li   : NLTE flag appended to notes

Fixes: RYA-102 (Eu II centroid), RYA-103 (Li contamination), RYA-104 (O I COG Ni sub)

Scope
-----
Main sample  : priority>0, blend_flag=False, depth >= PIPELINE['min_fit_depth']
Special cases: O I, Ni I (blend pair), Ba II, Eu II, Li I, Na I, P I
               force-measured regardless of blend_flag

Output
------
data/processed/solar_ew.csv          (gitignored — 22 MB class outputs)
results/plots/solar_ew_diagnostic.png (grid of Tier 1 Fe line profiles)

Linear issue: RYA-45
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from scipy.special import voigt_profile

from config.constants import (PIPELINE, PATHS, STAR_SOLAR, STAR_PROCYON,
                               get_star_params,
                               STAR_LINELISTS, EW_FIT_PARAMS,
                               NI6300_COG, LINE_SCORE_PARAMS)


def _ew_scores(ew_mA: float, ew_err: float, chi2: float) -> tuple:
    """Compute ew_snr_score, fit_chi2_score, saturation_score for one line."""
    P = LINE_SCORE_PARAMS
    # SNR score
    if np.isfinite(ew_err) and ew_err > 0:
        snr_score = float(min(ew_mA / ew_err / P['snr_saturation'], 1.0))
    else:
        snr_score = 0.5
    # Chi2 score
    if np.isfinite(chi2):
        if chi2 <= P['chi2_clean_max']:
            chi2_score = 1.0
        else:
            chi2_score = float(max(0.0, 1.0 - (chi2 - P['chi2_clean_max']) /
                                   (P['chi2_floor'] - P['chi2_clean_max'])))
    else:
        chi2_score = 0.5
    # Saturation score
    lo, hi = P['saturation_ew_low_mA'], P['saturation_ew_high_mA']
    sat_score = float(np.clip((ew_mA - lo) / (hi - lo), 0.0, 1.0))
    return round(snr_score, 4), round(chi2_score, 4), round(sat_score, 4)


# ── Special-case line registry ────────────────────────────────────────────────

# Lines force-measured regardless of blend_flag: (element, ion, nominal_wav_A)
SPECIAL_MEASURES = [
    # O I handled separately after main loop (COG Ni subtraction) — NOT in worklist
    ('O',  'I',  6363.776),   # forbidden line; CN blend present (flag only, don't subtract)
    ('Ba', 'II', 5853.668),   # HFS feature; total EW of all components
    ('Eu', 'II', 6645.127),   # HFS centroid (was 6645.064 = lowest component, wrong)
    ('Li', 'I',  6707.840),   # doublet midpoint (was 6707.756 = first component only)
    ('Na', 'I',  5682.633),
    ('Na', 'I',  5688.205),
    ('P',  'I',  6034.040),
    ('P',  'I',  6043.120),
    # Tier 1 Fe anchors — blend_flag=True in linelist but required for A(Fe) validation
    ('Fe', 'I',  5576.089),
    ('Fe', 'I',  6065.482),
    ('Fe', 'I',  6136.615),
    ('Fe', 'II', 6149.258),
    ('Fe', 'II', 6247.557),
]

# Per-line fit window overrides (Å half-width) for crowded / HFS lines.
# Default = PIPELINE['fit_window_A'] (2.0 Å).
LINE_WINDOWS = {
    ('O',  'I',  6300.304): 0.30,  # anchors at ±[0.225–0.30] Å (7–8 px, flux 0.98–0.99);
                                   # Si I 6299.599 at −0.705 Å and Sc II HFS at +0.38 Å
                                   # both safely outside; right edge clears Sc II by 0.08 Å
    ('O',  'I',  6363.776): 0.25,  # deep photospheric lines at −0.89 Å and +0.60 Å;
                                   # ±0.25 Å anchors [6363.526–6363.651] and [6363.901–6364.026]
                                   # land in CN gaps (p80 filter handles remaining CN features)
    ('Eu', 'II', 6645.127): 0.35,  # anchors at ±[0.26–0.35] Å — left in clean continuum,
                                   # right avoids known contamination beyond 6645.48
    ('Li', 'I',  6707.840): 0.25,  # avoid 6707.46 contamination feature
    ('Mg', 'I',  5711.088): 0.20,  # tightened: anchors at ±[0.15–0.20] Å, isolated clean line
    ('Ca', 'I',  6122.217): 2.78,  # wide continuum window: anchors at ±[2.08–2.78] Å have
                                   # p80~0.989 (genuinely clean); right cutoff ≈ 6125 Å
}

# Per-line profile-fit window — can be narrower than LINE_WINDOWS to avoid including
# broad wings or blends in the Voigt fit while still using wide anchors for the continuum.
# Keyed identically to LINE_WINDOWS; missing = use the same window as LINE_WINDOWS.
LINE_FIT_WINDOWS = {
    ('O',  'I',  6300.304): 0.15,  # fit ±0.15 Å only; continuum from wide ±0.30 Å anchors
                                   # (narrow fit avoids Si I 6299.6; wide continuum gives 7–8 px)
    ('O',  'I',  6363.776): 0.08,  # CN features at 6363.750 and 6363.830 flank the O I core;
                                   # ±0.08 Å fit captures core only; continuum from ±0.25 Å
    ('Ca', 'I',  6122.217): 0.10,  # ±0.1 Å fit only; Ca I 6102/6162 contaminate wider windows
                                   # (was 0.35 Å → measured 243 mÅ vs expected 140 mÅ)
}

_SPECIAL_WAV_TOL = 0.20  # Å tolerance for matching to SPECIAL_MEASURES

def _is_special(elem: str, ion: str, line_wav: float) -> bool:
    return any(e == elem and i == ion and abs(w - line_wav) < _SPECIAL_WAV_TOL
               for e, i, w in SPECIAL_MEASURES)

def _get_fit_window(elem: str, ion: str, line_wav: float,
                    depth: float = 0.0, base_window: float = None) -> float:
    """
    Return the fit half-window (Å) for this line.
    Explicit LINE_WINDOWS overrides take precedence; otherwise the window
    scales with line depth so that continuum anchors clear the Voigt wings.
    base_window overrides PIPELINE['fit_window_A'] for per-star tuning.
    """
    for (e, i, w), win in LINE_WINDOWS.items():
        if e == elem and i == ion and abs(w - line_wav) < _SPECIAL_WAV_TOL:
            return win
    base   = base_window if base_window is not None else PIPELINE['fit_window_A']
    breaks = PIPELINE['cont_window_depth_breaks']
    scales = PIPELINE['cont_window_scale_factors']
    for i, brk in enumerate(breaks):
        if depth < brk:
            return base * scales[i]
    return base * scales[-1]

def _get_profile_fit_window(elem: str, ion: str, line_wav: float,
                            cont_window: float) -> float:
    """
    Return the profile-fit half-window (Å) for this line.
    Always ≤ cont_window so the fit never exceeds the continuum extraction region.
    """
    for (e, i, w), fw in LINE_FIT_WINDOWS.items():
        if e == elem and i == ion and abs(w - line_wav) < _SPECIAL_WAV_TOL:
            return min(fw, cont_window)
    return cont_window


# Elements requiring NLTE note in output
NLTE_ELEMENTS = {'O', 'Na', 'C', 'Li'}

# Elements with HFS (single-profile fit = total HFS EW)
HFS_ELEMENTS = {'Ba', 'Eu', 'Li'}

# Tier 1 anchor lines shown in the diagnostic plot
TIER1_FE = [
    ('Fe', 'I',  5576.089, 'Fe I 5576'),
    ('Fe', 'I',  6065.482, 'Fe I 6065'),
    ('Fe', 'I',  6136.615, 'Fe I 6136'),
    ('Fe', 'II', 6149.258, 'Fe II 6149'),
    ('Fe', 'II', 6247.557, 'Fe II 6247'),
    ('[O', 'I',  6300.304, '[O I] 6300 + Ni'),  # special blend entry for plot
]


# ── Profile models ────────────────────────────────────────────────────────────

def _voigt_abs(x, x0, depth, sigma, gamma):
    """Single Voigt absorption: 1 – depth × V_norm(x – x0)."""
    v_peak = voigt_profile(0.0, sigma, gamma)
    if v_peak == 0:
        return np.ones_like(x)
    return 1.0 - depth * voigt_profile(x - x0, sigma, gamma) / v_peak


def _gauss_abs(x, x0, depth, sigma):
    """Gaussian absorption fallback."""
    return 1.0 - depth * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def _two_voigt_abs(x, x01, d1, s1, g1, x02, d2, s2, g2):
    """Two-component Voigt for O I + Ni I blend."""
    vp1 = voigt_profile(0.0, s1, g1) or 1.0
    vp2 = voigt_profile(0.0, s2, g2) or 1.0
    v1 = voigt_profile(x - x01, s1, g1) / vp1
    v2 = voigt_profile(x - x02, s2, g2) / vp2
    return 1.0 - d1 * v1 - d2 * v2


# ── Local continuum re-normalisation ─────────────────────────────────────────

def _local_renorm(wav: np.ndarray, flux: np.ndarray,
                  line_wav: float, window: float = None) -> tuple:
    """
    Estimate a linear continuum from the outer cont_edge_frac of the window on
    each side. Uses a percentile filter to exclude absorption-line pixels from
    the polyfit, so nearby lines in the anchor strips do not depress the
    continuum estimate. Returns (renorm_flux, continuum).

    window: half-width in Å; defaults to PIPELINE['fit_window_A'].
    cont_anchor_percentile (e.g. 80): pixels above this percentile in each
      anchor strip are used for the polyfit — they represent the continuum
      peaks between absorption lines rather than the line cores.
    """
    if window is None:
        window = PIPELINE['fit_window_A']
    edge_A = window * PIPELINE['cont_edge_frac']
    pct    = PIPELINE['cont_anchor_percentile']

    left  = (wav >= line_wav - window) & (wav <= line_wav - window + edge_A)
    right = (wav >= line_wav + window - edge_A) & (wav <= line_wav + window)

    if left.sum() < 3 or right.sum() < 3:
        return flux, np.ones_like(flux)

    # Filter each anchor strip to the top (100-pct)% of flux values — these
    # are the continuum-level pixels, not the absorption-line dips.
    tl = np.percentile(flux[left],  pct)
    tr = np.percentile(flux[right], pct)
    ml = flux[left]  >= tl
    mr = flux[right] >= tr

    x_anc = np.concatenate([wav[left][ml],  wav[right][mr]])
    y_anc = np.concatenate([flux[left][ml], flux[right][mr]])

    # Fall back to all anchor pixels if the filter left too few points
    if len(x_anc) < 4:
        x_anc = np.concatenate([wav[left], wav[right]])
        y_anc = np.concatenate([flux[left], flux[right]])

    try:
        coeffs = np.polyfit(x_anc, y_anc, 1)
    except Exception:
        return flux, np.ones_like(flux)

    cont = np.clip(np.polyval(coeffs, wav), 1e-6, None)
    return flux / cont, cont


# ── Single-line profile fitting ───────────────────────────────────────────────

_S0 = 0.025   # initial Gaussian sigma [Å]  (HARPS FWHM ~ 0.06 Å → σ ~ 0.025)
_G0 = 0.010   # initial Lorentzian gamma [Å]


def _fit_profile(wav: np.ndarray, flux: np.ndarray,
                 line_wav: float) -> tuple:
    """
    Attempt Voigt fit; fall back to Gaussian if Voigt fails or chi²_red > 0.05.
    Returns (popt, pcov, profile_type, chi2_red).
    popt indices: [x0, depth, sigma, (gamma for Voigt)].
    Depth is estimated from the ±0.1 Å core, not the window global min, so that
    strong nearby contaminants do not corrupt the initial guess.
    """
    core = np.abs(wav - line_wav) < 0.10
    flux_for_depth = flux[core] if core.sum() > 0 else flux
    depth0 = float(np.clip(1.0 - np.nanmin(flux_for_depth), 0.005, 0.995))

    # ── Voigt ─────────────────────────────────────────────────────────────────
    try:
        pv, cv = curve_fit(
            _voigt_abs, wav, flux,
            p0=[line_wav, depth0, _S0, _G0],
            bounds=([line_wav - 0.30, 0.001, 0.005, 0.0],
                    [line_wav + 0.30, 1.000, 0.40,  0.40]),
            maxfev=2000, method='trf',
        )
        res = flux - _voigt_abs(wav, *pv)
        chi2v = float(np.sum(res ** 2) / max(len(wav) - 4, 1))
        if chi2v < 0.05 and np.all(np.isfinite(cv)):
            return pv, cv, 'voigt', chi2v
    except Exception:
        pass

    # ── Gaussian fallback ─────────────────────────────────────────────────────
    try:
        pg, cg = curve_fit(
            _gauss_abs, wav, flux,
            p0=[line_wav, depth0, _S0],
            bounds=([line_wav - 0.30, 0.001, 0.005],
                    [line_wav + 0.30, 1.000, 0.40]),
            maxfev=2000, method='trf',
        )
        res = flux - _gauss_abs(wav, *pg)
        chi2g = float(np.sum(res ** 2) / max(len(wav) - 3, 1))
        return pg, cg, 'gaussian', chi2g
    except Exception:
        return None, None, 'failed', np.nan


def _fit_two_component(wav: np.ndarray, flux: np.ndarray) -> tuple:
    """
    Two-component Voigt fit for the O I 6300.304 + Ni I 6300.336 blend.
    Ni I is constrained to depth <= 0.08 (weak in the solar spectrum).
    Returns (popt_8, pcov_8, chi2_red) or (None, None, nan) on failure.
    """
    depth0 = float(np.clip(1.0 - np.nanmin(flux), 0.005, 0.99))
    try:
        popt, pcov = curve_fit(
            _two_voigt_abs, wav, flux,
            p0=[6300.304, depth0 * 0.85, _S0, _G0,
                6300.336, depth0 * 0.05, _S0, _G0],
            bounds=([6300.25, 0.001, 0.005, 0.0, 6300.30, 0.001, 0.005, 0.0],
                    [6300.36, 1.000, 0.20,  0.2, 6300.38, 0.080, 0.20,  0.2]),
            maxfev=3000, method='trf',
        )
        res  = flux - _two_voigt_abs(wav, *popt)
        chi2 = float(np.sum(res ** 2) / max(len(wav) - 8, 1))
        return popt, pcov, chi2
    except Exception:
        return None, None, np.nan


# ── Ni I 6300.336 COG prediction (for O I subtraction) ───────────────────────

def _predict_ni6300_ew(ew_df: pd.DataFrame, lines_df: pd.DataFrame) -> float:
    """
    Predict EW of Ni I 6300.336 via a robust linear COG fit to clean Ni I lines.

    Method (Allende Prieto et al. 2001, ApJ 556, L63):
      log(EW) ≈ log_gf − EP × θ + const   (linear COG regime only)
      θ = 5040 / Teff

    Robustness measures applied in order:
      1. Restrict to EW < 80 mÅ (linear part of COG; avoids saturated lines)
      2. 3-pass sigma-clip: reject |resid| > 0.7 dex
      3. Sanity cap: if prediction > STAR_SOLAR['ni6300_ew_lit_mA'] × 3, use
         the literature solar value (Allende Prieto+2001: ~1.1 mÅ) instead.

    Returns predicted EW in mÅ.
    """
    theta   = 5040.0 / get_star_params('solar')['teff']   # RYA-298: single source
    ni_ref  = STAR_SOLAR['ni6300_ew_lit_mA']   # 1.1 mÅ — solar literature value

    ni_ews = ew_df[
        (ew_df['element'] == 'Ni') &
        (ew_df['ion']     == 'I') &
        (ew_df['blend_flag'] == False) &
        (ew_df['ew_mA'] >= PIPELINE['ew_min_mA']) &
        (ew_df['ew_mA'] <  80.0)   # linear COG regime only
    ].copy()

    if len(ni_ews) < 3:
        return ni_ref   # fall back to literature value

    ni_lines = lines_df[
        (lines_df['element'] == 'Ni') & (lines_df['ion'] == 'I')
    ][['wavelength_air_A', 'log_gf', 'excitation_potential_eV']].copy()

    merged = ni_ews.merge(
        ni_lines, on='wavelength_air_A', how='inner'
    ).dropna(subset=['log_gf', 'excitation_potential_eV', 'ew_mA'])

    if len(merged) < 3:
        return ni_ref

    X = (merged['log_gf'] - merged['excitation_potential_eV'] * theta).to_numpy()
    Y = np.log10(merged['ew_mA'].clip(lower=0.1).to_numpy())

    # Iterative sigma-clipping (3 passes, 0.7 dex threshold)
    mask = np.ones(len(X), dtype=bool)
    for _ in range(3):
        if mask.sum() < 3:
            break
        const = float(np.median(Y[mask] - X[mask]))
        resid = np.abs(Y - X - const)
        mask  = resid < 0.7

    if mask.sum() < 3:
        return ni_ref

    const = float(np.median(Y[mask] - X[mask]))
    ni6300_log_ew = NI6300_COG['log_gf'] - NI6300_COG['excitation_potential_eV'] * theta + const
    pred = float(10 ** ni6300_log_ew)

    # Sanity cap: prediction must not exceed 3× the solar literature value
    if pred > ni_ref * 3.0:
        return ni_ref

    return pred


def _measure_oi6300(solar_wav: np.ndarray, solar_flux: np.ndarray,
                    results: list, lines_df: pd.DataFrame,
                    ni_ew_mA: float) -> None:
    """
    Measure O I 6300.304 after subtracting a COG-predicted Ni I 6300.336 model.

    Method (Allende Prieto et al. 2001, ApJ 556, L63):
      1. Caller supplies ni_ew_mA from _predict_ni6300_ew() (COG of clean Ni I lines)
      2. Build a Gaussian Ni I model at HARPS resolution (σ ≈ 0.023 Å) with that EW
      3. Subtract the Ni model → cleaned O I spectrum
      4. Fit O I alone with _fit_profile on the cleaned spectrum
      5. Flag: Ni_blend_subtracted | forbidden_line | NLTE_flag

    Uses narrow ±0.15 Å window to exclude Si I 6299.599.
    Appends one row to `results` in-place.
    """
    WAV_O  = 6300.304
    WAV_NI = 6300.336
    SIGMA  = 0.023   # HARPS σ at 6300 Å (R~115,000 → FWHM~0.055 Å → σ=0.023 Å)

    cont_win = LINE_WINDOWS[('O', 'I', WAV_O)]       # 0.50 Å — wide continuum (12–13 px anchors)
    fit_win  = LINE_FIT_WINDOWS[('O', 'I', WAV_O)]   # 0.15 Å — narrow fit (avoids Si I 6299.6)

    cont_mask = (solar_wav >= WAV_O - cont_win) & (solar_wav <= WAV_O + cont_win)
    if cont_mask.sum() < 10:
        return

    wav_c  = solar_wav[cont_mask]
    edge_A = cont_win * PIPELINE['cont_edge_frac']

    # Anchor pixel diagnostic
    left_n  = int(((solar_wav >= WAV_O - cont_win) &
                   (solar_wav <= WAV_O - cont_win + edge_A)).sum())
    right_n = int(((solar_wav >= WAV_O + cont_win - edge_A) &
                   (solar_wav <= WAV_O + cont_win)).sum())
    print(f"    [O I 6300 continuum] left_anchor={left_n}px  right_anchor={right_n}px  "
          f"(cont_win=±{cont_win}Å  fit_win=±{fit_win}Å)")

    flux_r, _ = _local_renorm(solar_wav[cont_mask], solar_flux[cont_mask],
                              WAV_O, window=cont_win)

    # Build Ni I Gaussian model and remove it (Allende Prieto et al. 2001)
    ni_depth  = ni_ew_mA / (SIGMA * np.sqrt(2 * np.pi) * 1000.0)
    ni_model  = ni_depth * np.exp(-0.5 * ((wav_c - WAV_NI) / SIGMA) ** 2)
    flux_clean = flux_r + ni_model   # restore continuum under Ni feature

    # Narrow to fit window for profile fitting
    fm = np.abs(wav_c - WAV_O) <= fit_win
    wav_fit   = wav_c[fm] if fm.sum() >= 5 else wav_c
    flux_fit  = flux_clean[fm] if fm.sum() >= 5 else flux_clean

    # Fit O I alone on Ni-subtracted, fit-window-sliced spectrum
    popt, pcov, profile_t, chi2 = _fit_profile(wav_fit, flux_fit, WAV_O)
    ew_mA = _integrate_profile(wav_fit, popt, profile_t) if popt is not None else np.nan

    # Uncertainty from continuum-window edge pixels
    edge_mask = (wav_c <= WAV_O - cont_win + edge_A) | (wav_c >= WAV_O + cont_win - edge_A)
    ew_err    = (_ew_error(flux_r[edge_mask], ew_mA, popt, pcov, profile_type=profile_t)
                 if popt is not None else np.nan)

    _ew_v = float(ew_mA) if not np.isnan(ew_mA) else 0.0
    _err_v = float(ew_err) if not np.isnan(ew_err) else np.nan
    _chi_v = float(chi2) if not np.isnan(chi2) else np.nan
    _snr_s, _chi2_s, _sat_s = _ew_scores(_ew_v, _err_v, _chi_v)
    results.append(dict(
        element='O', ion='I', wavelength_air_A=WAV_O,
        ew_mA=round(float(ew_mA), 2) if not np.isnan(ew_mA) else np.nan,
        ew_err_mA=round(float(ew_err), 2) if not np.isnan(ew_err) else np.nan,
        profile_type=profile_t if popt is not None else 'failed',
        chi2=round(chi2, 4) if not np.isnan(chi2) else np.nan,
        ew_snr_score=_snr_s,
        fit_chi2_score=_chi2_s,
        saturation_score=_sat_s,
        blend_flag=True,
        notes=(f'Ni_blend_subtracted; ni_ew_pred={ni_ew_mA:.3f}_mA; '
               'forbidden_line; NLTE_flag; method=Allende_Prieto+2001'),
    ))


# ── EW calculation ────────────────────────────────────────────────────────────

def _integrate_profile(wav: np.ndarray, popt, profile_type: str) -> float:
    """Numerically integrate (1 – model) over the window → EW in mÅ."""
    x = np.linspace(wav[0], wav[-1], 5000)
    if profile_type == 'voigt':
        y = _voigt_abs(x, *popt)
    elif profile_type == 'gaussian':
        y = _gauss_abs(x, *popt)
    else:
        return np.nan
    return float(np.trapz(np.clip(1.0 - y, 0.0, None), x)) * 1000.0  # → mÅ


def _ew_error(flux_edges: np.ndarray, ew_mA: float,
              popt, pcov, profile_type: str = 'gaussian') -> float:
    """
    EW uncertainty = fitting error ⊕ continuum placement error.
    Continuum error: σ_cont × FWHM_fit (line width smeared across continuum scatter).
    Fitting error:   δ(depth) / depth × EW  (depth is popt[1] in all our models).
    """
    if popt is None or pcov is None or not np.all(np.isfinite(pcov)):
        return np.nan
    sigma_cont = float(np.std(flux_edges)) if len(flux_edges) > 2 else 0.01
    if profile_type == 'voigt' and len(popt) == 4:
        sigma = float(popt[2])
        gamma = float(popt[3])
        fg = 2.355 * sigma
        fl = 2.0 * gamma
        fwhm_A = 0.5346 * fl + np.sqrt(0.2166 * fl**2 + fg**2)
    else:
        fwhm_A = float(popt[2]) * 2.355
    cont_err   = sigma_cont * fwhm_A * 1000.0  # → mÅ

    depth_std  = float(np.sqrt(max(pcov[1, 1], 0.0)))
    depth_val  = float(max(abs(popt[1]), 1e-6))
    fit_err    = (depth_std / depth_val) * abs(ew_mA)

    return float(np.sqrt(fit_err ** 2 + cont_err ** 2))


# ── Worklist construction ─────────────────────────────────────────────────────

def _build_worklist(lines_df: pd.DataFrame) -> pd.DataFrame:
    """
    Return the set of absorption features to measure.

    Includes:
      (a) priority>0, blend_flag=False, central_depth >= min_fit_depth
      (b) SPECIAL_MEASURES entries regardless of blend_flag
    Deduplicates on (element, ion, wavelength rounded to 0.1 Å).
    """
    df = lines_df.copy()
    df['blend_flag'] = (df['blend_flag'].astype(str)
                         .str.strip().str.lower() == 'true')

    min_depth = PIPELINE['min_fit_depth']

    clean = df[
        (df['priority'] > 0) &
        (~df['blend_flag']) &
        (df['central_depth'] >= min_depth)
    ].copy()

    specials = []
    for elem, ion, nom_wav in SPECIAL_MEASURES:
        hit = df[
            (df['element'] == elem) &
            (df['ion']     == ion) &
            (df['wavelength_air_A'].between(nom_wav - 0.15, nom_wav + 0.15))
        ]
        if len(hit) == 0:
            continue
        rep = hit.loc[hit['central_depth'].idxmax()].copy()
        rep['wavelength_air_A'] = nom_wav   # use nominal for consistency
        specials.append(rep)

    if specials:
        combined = pd.concat([clean, pd.DataFrame(specials)], ignore_index=True)
        combined['_wkey'] = (combined['wavelength_air_A'].round(1).astype(str) +
                             combined['element'] + combined['ion'])
        combined = combined.drop_duplicates('_wkey').drop(columns=['_wkey'])
    else:
        combined = clean

    return combined.sort_values('wavelength_air_A').reset_index(drop=True)


# ── Main measurement loop ─────────────────────────────────────────────────────

def _measure_all(solar_wav: np.ndarray, solar_flux: np.ndarray,
                 worklist: pd.DataFrame,
                 lines_df: pd.DataFrame,
                 run_oi_cog: bool = True,
                 fit_window_A: float = None) -> pd.DataFrame:
    """Fit every line in worklist; return results DataFrame.
    run_oi_cog: perform [O I] 6300 COG Ni-subtraction (solar only).
    fit_window_A: per-star fit half-window override (None → PIPELINE default).
    """
    ew_min    = PIPELINE['ew_min_mA']
    ew_max    = PIPELINE['ew_max_mA']
    blue_warn = PIPELINE['blue_edge_warn_A']
    edge_frac = PIPELINE['cont_edge_frac']

    results = []

    n = len(worklist)
    print(f"  Fitting {n} lines …")
    for i, row in worklist.iterrows():
        if i % 100 == 0:
            print(f"    {i}/{n}", end='\r', flush=True)

        elem     = str(row['element'])
        ion      = str(row['ion'])
        line_wav = float(row['wavelength_air_A'])
        blend    = bool(row.get('blend_flag', False))

        # O I and Ni I 6300 are handled after the main loop via COG subtraction
        if elem == 'O' and ion == 'I' and abs(line_wav - 6300.304) < 0.15:
            continue
        if elem == 'Ni' and ion == 'I' and abs(line_wav - 6300.336) < 0.15:
            continue

        # ── Depth peek → adaptive fit window ─────────────────────────────────
        peek = np.abs(solar_wav - line_wav) < 0.15
        depth_est = float(np.clip(1.0 - np.nanmin(solar_flux[peek]), 0.0, 1.0)) \
                    if peek.sum() > 0 else 0.0
        window = _get_fit_window(elem, ion, line_wav, depth_est,
                                 base_window=fit_window_A)
        edge_A = window * edge_frac

        # ── Extract and re-normalise window ───────────────────────────────────
        mask = ((solar_wav >= line_wav - window) &
                (solar_wav <= line_wav + window))
        if mask.sum() < 10:
            continue

        wav_w  = solar_wav[mask]
        flux_w = solar_flux[mask]
        flux_r, _ = _local_renorm(wav_w, flux_w, line_wav, window=window)

        # Edge pixels for noise estimate (always from the full continuum window)
        edge_mask = (wav_w <= line_wav - window + edge_A) | (wav_w >= line_wav + window - edge_A)
        flux_edges = flux_r[edge_mask]

        # Apply narrower profile-fit window if specified (e.g. Ca I 6122 ±0.1 Å)
        fit_win = _get_profile_fit_window(elem, ion, line_wav, window)
        if fit_win < window:
            fm = np.abs(wav_w - line_wav) <= fit_win
            wav_fit  = wav_w[fm] if fm.sum() >= 5 else wav_w
            flux_fit = flux_r[fm] if fm.sum() >= 5 else flux_r
        else:
            wav_fit, flux_fit = wav_w, flux_r

        # ── Fit ───────────────────────────────────────────────────────────────
        popt, pcov, profile_t, chi2 = _fit_profile(wav_fit, flux_fit, line_wav)
        if popt is None:
            continue

        ew_mA  = _integrate_profile(wav_fit, popt, profile_t)
        ew_err = _ew_error(flux_edges, ew_mA, popt, pcov, profile_type=profile_t)

        is_sp = _is_special(elem, ion, line_wav)

        if np.isnan(ew_mA):
            continue
        if not is_sp and (ew_mA < ew_min or ew_mA > ew_max):
            continue

        # ── Build notes ───────────────────────────────────────────────────────
        notes = []
        if ew_mA < ew_min:
            notes.append('ew_below_threshold')
        if ew_mA > ew_max:
            notes.append('ew_above_threshold')
        if blend:
            notes.append('blend_flag')
        if elem in NLTE_ELEMENTS:
            notes.append('NLTE_flag')
            if elem == 'Li':
                notes.append('LTE_lower_bound')
        if elem in HFS_ELEMENTS:
            notes.append('hfs_total_ew')
        if elem == 'Li' and ion == 'I':
            notes.append('CN_blend_possible; upper_limit')
        if elem == 'O' and ion == 'I' and abs(line_wav - 6363.776) < 0.20:
            notes.extend(['CN_blend_possible', 'forbidden_line', 'upper_limit'])
        if line_wav < blue_warn:
            notes.append('low_snr_blue_edge')

        _snr_s, _chi2_s, _sat_s = _ew_scores(ew_mA, ew_err, chi2)
        results.append(dict(
            element=elem, ion=ion,
            wavelength_air_A=round(line_wav, 5),
            ew_mA=round(ew_mA, 2),
            ew_err_mA=round(ew_err, 2) if not np.isnan(ew_err) else np.nan,
            profile_type=profile_t,
            chi2=round(chi2, 4) if not np.isnan(chi2) else np.nan,
            ew_snr_score=_snr_s,
            fit_chi2_score=_chi2_s,
            saturation_score=_sat_s,
            blend_flag=blend,
            notes='; '.join(notes),
        ))

    print(f"    {n}/{n} — done.         ")

    # ── O I 6300 + Ni I COG subtraction (after main loop, solar only) ──────────
    if run_oi_cog:
        print(f"  Measuring [O I] 6300 via COG Ni subtraction …")
        ew_df_temp  = pd.DataFrame(results)
        ni_ew_pred  = _predict_ni6300_ew(ew_df_temp, lines_df)
        print(f"    Ni I 6300.336 predicted EW = {ni_ew_pred:.3f} mÅ  (COG, Allende Prieto+2001)")
        _measure_oi6300(solar_wav, solar_flux, results, lines_df, ni_ew_pred)
    else:
        print(f"  Skipping [O I] 6300 COG subtraction (solar-only procedure)")

    return pd.DataFrame(results)


# ── Diagnostic plot ───────────────────────────────────────────────────────────

def _plot_tier1_fe(solar_wav: np.ndarray, solar_flux: np.ndarray,
                   results: pd.DataFrame, out_path,
                   star_key: str = 'solar') -> None:
    """2 × 3 grid of Tier 1 Fe anchor line fits + [O I] 6300 blend panel."""
    ew_params = EW_FIT_PARAMS.get(star_key, EW_FIT_PARAMS['solar'])
    window = ew_params['fit_window_A']

    fig = plt.figure(figsize=(15, 9))
    if star_key == 'solar':
        _f = get_star_params('solar')   # RYA-298: fundamentals from single source
        title = (f"Tier 1 line profile fits — Solar HARPS  |  "
                 f"Teff={_f['teff']:.0f} K  log g={_f['logg']}  "
                 f"ξ={STAR_SOLAR['vturb_kms']} km/s\n"
                 f"{STAR_SOLAR['citation']}  |  ESO {STAR_SOLAR['program_id']}")
    else:
        sp = STAR_PROCYON; _f = get_star_params('procyon')
        title = (f"Tier 1 line profile fits — {sp['name']} HARPS  |  "
                 f"Teff={_f['teff']:.0f} K  log g={_f['logg']}  "
                 f"ξ={sp['vturb_kms']} km/s  |  {sp['hd']}")
    fig.suptitle(title, fontsize=10)
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.32)

    for pi, (elem, ion, line_wav, label) in enumerate(TIER1_FE):
        ax = fig.add_subplot(gs[pi // 3, pi % 3])

        # Resolve [O alias used in TIER1_FE
        elem_real = elem.lstrip('[')

        mask = ((solar_wav >= line_wav - window) &
                (solar_wav <= line_wav + window))
        if mask.sum() < 5:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                    transform=ax.transAxes)
            ax.set_title(label, fontsize=9)
            continue

        wav_w  = solar_wav[mask]
        flux_r, _ = _local_renorm(solar_wav[mask], solar_flux[mask], line_wav)

        ax.plot(wav_w, flux_r, 'o', color='#4a90d9', ms=1.8, alpha=0.7,
                label='Observed')
        ax.axvline(line_wav, color='gray', lw=0.6, ls='--', alpha=0.5)
        ax.axhline(1.0,      color='gray', lw=0.5, ls=':', alpha=0.5)

        # EW from results
        hit = results[
            (results['element'] == elem_real) &
            (results['ion']     == ion) &
            (results['wavelength_air_A'].between(line_wav - 0.1, line_wav + 0.1))
        ]
        ew_str = f"EW = {hit.iloc[0]['ew_mA']:.1f} mÅ" if len(hit) > 0 else 'not measured'

        # Re-fit for plotting curve
        x_fine = np.linspace(wav_w[0], wav_w[-1], 2000)
        if abs(line_wav - 6300.304) < 0.1:
            p2, _, _ = _fit_two_component(wav_w, flux_r)
            if p2 is not None:
                ax.plot(x_fine, _two_voigt_abs(x_fine, *p2),
                        '-', color='tomato', lw=1.2, label='2-comp Voigt')
                ax.axvline(6300.336, color='darkorange', lw=0.8,
                           ls=':', alpha=0.8, label='Ni I 6300.336')
        else:
            popt, _, pt, _ = _fit_profile(wav_w, flux_r, line_wav)
            if popt is not None:
                fn = _voigt_abs if pt == 'voigt' else _gauss_abs
                ax.plot(x_fine, fn(x_fine, *popt),
                        '-', color='tomato', lw=1.2,
                        label=f"{'Voigt' if pt=='voigt' else 'Gauss'} fit")

        ax.set_title(f'{label}\n{ew_str}', fontsize=9)
        ax.set_xlabel('Wavelength (Å)', fontsize=8)
        ax.set_ylabel('Norm. flux', fontsize=8)
        ax.set_ylim(0.0, 1.15)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc='lower center')

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {out_path.name}")


# ── O I 6300 diagnostic plot ──────────────────────────────────────────────────

def _plot_oi6300_diagnostic(solar_wav: np.ndarray, solar_flux: np.ndarray,
                            results: pd.DataFrame, out_path) -> None:
    """
    Three-panel diagnostic for [O I] 6300.304 showing the Ni subtraction.
    Left:   raw spectrum + Ni I model (shaded).
    Centre: Ni-subtracted spectrum + O I fit.
    Right:  fit residuals.
    """
    import re

    WAV_O  = 6300.304
    WAV_NI = 6300.336
    SIGMA  = 0.023
    win    = LINE_WINDOWS[('O', 'I', WAV_O)]  # 0.15 Å

    mask   = (solar_wav >= WAV_O - win) & (solar_wav <= WAV_O + win)
    wav_w  = solar_wav[mask]
    flux_r, _ = _local_renorm(solar_wav[mask], solar_flux[mask], WAV_O, window=win)

    hit = results[
        (results['element'] == 'O') &
        (results['ion']     == 'I') &
        (results['wavelength_air_A'].between(WAV_O - 0.1, WAV_O + 0.1))
    ]
    if len(hit) == 0:
        return

    row       = hit.iloc[0]
    ew_O      = row['ew_mA']
    notes_str = str(row.get('notes', ''))
    m         = re.search(r'ni_ew_pred=([\d.]+)_mA', notes_str)
    ni_ew     = float(m.group(1)) if m else STAR_SOLAR['ni6300_ew_lit_mA']

    ni_depth  = ni_ew / (SIGMA * np.sqrt(2 * np.pi) * 1000.0)
    ni_model  = ni_depth * np.exp(-0.5 * ((wav_w - WAV_NI) / SIGMA) ** 2)
    flux_clean = flux_r + ni_model

    popt, _, pt, chi2 = _fit_profile(wav_w, flux_clean, WAV_O)
    x_fine = np.linspace(wav_w[0], wav_w[-1], 1000)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle(
        f"[O I] 6300.304 — Ni I 6300.336 subtraction diagnostic\n"
        f"O I EW = {ew_O:.2f} mÅ  (target ~4.5 mÅ)  |  "
        f"Ni I predicted = {ni_ew:.3f} mÅ  (Allende Prieto et al. 2001)",
        fontsize=10,
    )

    # Left: raw blend + Ni model
    axes[0].plot(wav_w, flux_r, 'o', color='#4a90d9', ms=3, alpha=0.8, label='Observed')
    axes[0].fill_between(wav_w, 1.0 - ni_model, 1.0, alpha=0.30,
                         color='orange', label=f'Ni model ({ni_ew:.3f} mÅ)')
    axes[0].axvline(WAV_O,  color='purple', lw=0.9, ls='--', alpha=0.7, label='O I 6300.304')
    axes[0].axvline(WAV_NI, color='orange', lw=0.9, ls='--', alpha=0.7, label='Ni I 6300.336')
    axes[0].axhline(1.0, color='gray', lw=0.4, ls=':')
    axes[0].set_title('Raw + Ni model', fontsize=9)
    axes[0].legend(fontsize=7)
    axes[0].set_xlabel('Wavelength (Å)', fontsize=9)
    axes[0].set_ylabel('Renorm. flux', fontsize=9)
    axes[0].tick_params(labelsize=8)

    # Centre: Ni-subtracted + O I fit
    axes[1].plot(wav_w, flux_clean, 'o', color='#4a90d9', ms=3, alpha=0.8,
                 label='Ni-subtracted')
    if popt is not None:
        fn = _voigt_abs if pt == 'voigt' else _gauss_abs
        axes[1].plot(x_fine, fn(x_fine, *popt), '-', color='tomato', lw=1.5,
                     label=f'O I fit ({ew_O:.2f} mÅ)')
    axes[1].axvline(WAV_O, color='purple', lw=0.9, ls='--', alpha=0.7)
    axes[1].axhline(1.0, color='gray', lw=0.4, ls=':')
    axes[1].set_title(f'Ni-subtracted + O I fit\nEW = {ew_O:.2f} mÅ', fontsize=9)
    axes[1].legend(fontsize=7)
    axes[1].set_xlabel('Wavelength (Å)', fontsize=9)
    axes[1].set_ylabel('Renorm. flux', fontsize=9)
    axes[1].tick_params(labelsize=8)

    # Right: residuals
    if popt is not None:
        fn = _voigt_abs if pt == 'voigt' else _gauss_abs
        resid = flux_clean - fn(wav_w, *popt)
        axes[2].plot(wav_w, resid, 'o', color='steelblue', ms=3, alpha=0.8)
        axes[2].axhline(0.0, color='gray', lw=0.8, ls='--')
        axes[2].set_title(f'Residuals\nchi²_red = {chi2:.5f}', fontsize=9)
    else:
        axes[2].text(0.5, 0.5, 'Fit failed', ha='center', va='center',
                     transform=axes[2].transAxes)
        axes[2].set_title('Residuals', fontsize=9)
    axes[2].set_xlabel('Wavelength (Å)', fontsize=9)
    axes[2].set_ylabel('Residuals', fontsize=9)
    axes[2].tick_params(labelsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {out_path.name}")


# ── Ca I 6122 diagnostic plot ─────────────────────────────────────────────────

def _plot_ca6122_diagnostic(solar_wav: np.ndarray, solar_flux: np.ndarray,
                            results: pd.DataFrame, out_path) -> None:
    """
    Two-panel diagnostic for Ca I 6122.217 showing why ±0.1 Å fit window is needed.
    Left: full ±2.78 Å continuum region with anchor strips highlighted.
    Right: ±0.40 Å zoom showing the fitted profile vs ±0.1 Å fit boundary.
    """
    line_wav  = 6122.217
    cont_win  = LINE_WINDOWS[('Ca', 'I', line_wav)]   # 2.78 Å
    fit_win   = LINE_FIT_WINDOWS[('Ca', 'I', line_wav)]  # 0.10 Å
    edge_A    = cont_win * PIPELINE['cont_edge_frac']

    cont_mask = ((solar_wav >= line_wav - cont_win) &
                 (solar_wav <= line_wav + cont_win))
    wav_c   = solar_wav[cont_mask]
    flux_c  = solar_flux[cont_mask]
    flux_r, _ = _local_renorm(wav_c, flux_c, line_wav, window=cont_win)

    # Profile fit on narrow window
    fm = np.abs(wav_c - line_wav) <= fit_win
    wav_fit  = wav_c[fm] if fm.sum() >= 5 else wav_c
    flux_fit = flux_r[fm] if fm.sum() >= 5 else flux_r
    popt, _, pt, chi2 = _fit_profile(wav_fit, flux_fit, line_wav)

    hit = results[
        (results['element'] == 'Ca') &
        (results['ion']     == 'I') &
        (results['wavelength_air_A'].between(line_wav - 0.2, line_wav + 0.2))
    ]
    ew_str = f"EW = {hit.iloc[0]['ew_mA']:.1f} mÅ" if len(hit) > 0 else 'not measured'

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(
        f"Ca I 6122.217 Å — continuum window ±{cont_win} Å, fit window ±{fit_win} Å\n"
        f"{ew_str}  |  profile={pt}  chi²={chi2:.4f}" if popt is not None
        else f"Ca I 6122.217 Å — continuum window ±{cont_win} Å, fit window ±{fit_win} Å",
        fontsize=10,
    )

    # Left: full continuum region
    left_anc  = wav_c <= line_wav - cont_win + edge_A
    right_anc = wav_c >= line_wav + cont_win - edge_A
    ax1.plot(wav_c, flux_r, 'o', color='#4a90d9', ms=1.5, alpha=0.6, label='Observed')
    ax1.plot(wav_c[left_anc],  flux_r[left_anc],  's', color='green',  ms=3, alpha=0.8,
             label='Continuum anchors')
    ax1.plot(wav_c[right_anc], flux_r[right_anc], 's', color='green',  ms=3, alpha=0.8)
    ax1.axvline(line_wav, color='gray', lw=0.8, ls='--', alpha=0.6, label='Ca I 6122')
    ax1.axvspan(line_wav - fit_win, line_wav + fit_win, alpha=0.12,
                color='tomato', label=f'Fit window ±{fit_win} Å')
    ax1.axhline(1.0, color='gray', lw=0.5, ls=':', alpha=0.5)
    ax1.set_xlabel('Wavelength (Å)', fontsize=9)
    ax1.set_ylabel('Renorm. flux', fontsize=9)
    ax1.set_title('Full continuum region', fontsize=9)
    ax1.legend(fontsize=7, loc='lower right')
    ax1.tick_params(labelsize=8)

    # Right: zoom on fit window
    zoom = 0.40
    zoom_mask = np.abs(wav_c - line_wav) <= zoom
    ax2.plot(wav_c[zoom_mask], flux_r[zoom_mask], 'o', color='#4a90d9',
             ms=2.5, alpha=0.8, label='Observed')
    if popt is not None:
        x_fine = np.linspace(line_wav - fit_win, line_wav + fit_win, 1000)
        fn = _voigt_abs if pt == 'voigt' else _gauss_abs
        ax2.plot(x_fine, fn(x_fine, *popt), '-', color='tomato', lw=1.5,
                 label=f"{'Voigt' if pt=='voigt' else 'Gauss'} fit")
    ax2.axvline(line_wav, color='gray', lw=0.8, ls='--', alpha=0.6)
    ax2.axvspan(line_wav - fit_win, line_wav + fit_win, alpha=0.12,
                color='tomato', label=f'Fit ±{fit_win} Å')
    ax2.axhline(1.0, color='gray', lw=0.5, ls=':', alpha=0.5)
    ax2.set_xlim(line_wav - zoom, line_wav + zoom)
    ax2.set_xlabel('Wavelength (Å)', fontsize=9)
    ax2.set_ylabel('Renorm. flux', fontsize=9)
    ax2.set_title(f'Fit zoom  |  {ew_str}', fontsize=9)
    ax2.legend(fontsize=7, loc='lower right')
    ax2.tick_params(labelsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {out_path.name}")


# ── Entry point ───────────────────────────────────────────────────────────────

def run(star_id: str = 'solar') -> pd.DataFrame:
    """
    Called by run_pipeline.py. Returns the EW results DataFrame.
    star_id: any value containing 'solar' routes to the solar calibration path.
    """
    print(f"\n{'='*60}\n  lines_fit — {star_id}\n{'='*60}")

    # ── Route star ────────────────────────────────────────────────────────────
    star_key = 'solar'
    for key in STAR_LINELISTS:
        if key in star_id.lower():
            star_key = key
            break
    is_solar = (star_key == 'solar')
    ew_params = EW_FIT_PARAMS.get(star_key, EW_FIT_PARAMS['solar'])
    print(f"  Star key: {star_key}  |  fit_window_A={ew_params['fit_window_A']} Å  "
          f"|  voigt_threshold={ew_params['voigt_threshold_mA']} mÅ")

    # ── Load inputs ───────────────────────────────────────────────────────────
    print(f"\n[1/5] Loading inputs")
    norm_path  = PATHS[f'{star_key}_normalized']
    lines_path = STAR_LINELISTS[star_key]

    if not norm_path.exists():
        raise FileNotFoundError(
            f"Normalised spectrum not found: {norm_path}\n"
            "Run spectra_normalize.py first."
        )

    norm       = pd.read_csv(norm_path, comment='#')
    spec_wav   = norm['wavelength_air_A'].to_numpy()
    spec_flux  = norm['flux_normalized'].to_numpy()
    print(f"  Spectrum: {len(spec_wav):,} pixels, "
          f"{spec_wav[0]:.1f}–{spec_wav[-1]:.1f} Å")

    lines = pd.read_csv(lines_path, low_memory=False, comment='#')
    lines = lines[lines['priority'] > 0].copy()
    print(f"  Line list: {len(lines):,} priority>0 lines")

    # ── Build worklist ────────────────────────────────────────────────────────
    print(f"\n[2/5] Building worklist")
    worklist = _build_worklist(lines)
    n_special = (worklist['blend_flag'].astype(str).str.lower() == 'true').sum()
    print(f"  {len(worklist)} features to measure ({n_special} special/blended force-included)")
    fe1 = (worklist['element'] == 'Fe') & (worklist['ion'] == 'I')
    fe2 = (worklist['element'] == 'Fe') & (worklist['ion'] == 'II')
    print(f"  Fe I: {fe1.sum()}  Fe II: {fe2.sum()}")

    # ── Measure EWs ───────────────────────────────────────────────────────────
    print(f"\n[3/5] Measuring EWs")
    results = _measure_all(spec_wav, spec_flux, worklist, lines,
                           run_oi_cog=is_solar,
                           fit_window_A=ew_params['fit_window_A'])
    print(f"  {len(results)} lines measured")

    # ── QA summary ────────────────────────────────────────────────────────────
    print(f"\n[4/5] QA summary")
    print(f"  {'Element':<6} {'N':>5}  {'EW median':>10}  {'EW range':>20}")
    print(f"  {'-'*50}")
    for elem, grp in results.groupby('element'):
        ew = grp['ew_mA'].dropna()
        if len(ew) == 0:
            continue
        print(f"  {elem:<6} {len(ew):>5}  "
              f"{ew.median():>9.1f} mÅ  "
              f"{ew.min():.1f}–{ew.max():.1f} mÅ")

    # ── Save CSV ──────────────────────────────────────────────────────────────
    out_csv  = PATHS[f'{star_key}_ew']
    out_plot = PATHS[f'{star_key}_ew_diagnostic']

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    col_order = ['element', 'ion', 'wavelength_air_A', 'ew_mA', 'ew_err_mA',
                 'profile_type', 'chi2', 'blend_flag', 'notes']
    out_df = results[col_order].copy()
    out_df['notes'] = out_df['notes'].fillna('')
    out_df.to_csv(out_csv, index=False)
    print(f"\n  Saved → {out_csv.name}  ({len(results)} rows)")

    # ── Diagnostic plot ───────────────────────────────────────────────────────
    print(f"\n[5/5] Generating diagnostic plots")
    out_plot.parent.mkdir(parents=True, exist_ok=True)
    _plot_tier1_fe(spec_wav, spec_flux, results, out_plot, star_key=star_key)

    if is_solar:
        out_oi_plot = out_plot.parent / 'solar_oi6300_diagnostic.png'
        _plot_oi6300_diagnostic(spec_wav, spec_flux, results, out_oi_plot)

        out_ca_plot = out_plot.parent / 'solar_ca6122_diagnostic.png'
        _plot_ca6122_diagnostic(spec_wav, spec_flux, results, out_ca_plot)

    print(f"\n✓  lines_fit complete.")
    return results


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar')
    run(star_id=ap.parse_args().star)
