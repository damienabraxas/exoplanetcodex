"""
pipeline/lines_fit.py
=====================
Measure equivalent widths (EWs) of absorption lines in the HARPS solar spectrum.

For each line: extract local window → re-normalise continuum → fit Voigt
(Gaussian fallback) → integrate profile → EW in mÅ → propagate uncertainty.

Special handling
----------------
O I  6300.304 : two-component Voigt fit with Ni I 6300.336 blend; both reported
Ba II 5853.668 : 23 HFS components unresolved at HARPS R; fit single profile = total EW
Eu II 6645.064 : 30 HFS components; same approach
Li I  6707.756 : HFS + NLTE flag; EW is LTE lower bound
Na, O, C, Li   : NLTE flag appended to notes

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

from config.constants import PIPELINE, PATHS, STAR_SOLAR


# ── Special-case line registry ────────────────────────────────────────────────

# Lines force-measured regardless of blend_flag: (element, ion, nominal_wav_A)
SPECIAL_MEASURES = [
    ('O',  'I',  6300.304),   # [O I] forbidden — primary oxygen indicator
    ('Ni', 'I',  6300.336),   # O I blend partner — measured jointly
    ('Ba', 'II', 5853.668),   # HFS feature; total EW of all components
    ('Eu', 'II', 6645.064),   # HFS feature; total EW of all components
    ('Li', 'I',  6707.756),   # HFS + NLTE; EW is LTE lower bound
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

_SPECIAL_WAV_TOL = 0.15  # Å tolerance for matching to SPECIAL_MEASURES

def _is_special(elem: str, ion: str, line_wav: float) -> bool:
    return any(e == elem and i == ion and abs(w - line_wav) < _SPECIAL_WAV_TOL
               for e, i, w in SPECIAL_MEASURES)

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
                  line_wav: float) -> tuple:
    """
    Fit a linear continuum through the outer cont_edge_frac of the window on
    each side. Returns (renormalised_flux, continuum_array).
    Falls back to identity if there are too few edge points.
    """
    window = PIPELINE['fit_window_A']
    edge_A = window * PIPELINE['cont_edge_frac']

    left  = (wav >= line_wav - window) & (wav <= line_wav - window + edge_A)
    right = (wav >= line_wav + window - edge_A) & (wav <= line_wav + window)

    if left.sum() < 3 or right.sum() < 3:
        return flux, np.ones_like(flux)

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
    """
    depth0 = float(np.clip(1.0 - np.nanmin(flux), 0.005, 0.995))

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
              popt, pcov) -> float:
    """
    EW uncertainty = fitting error ⊕ continuum placement error.
    Continuum error: σ_cont × FWHM_fit (line width smeared across continuum scatter).
    Fitting error:   δ(depth) / depth × EW  (depth is popt[1] in all our models).
    """
    if popt is None or pcov is None or not np.all(np.isfinite(pcov)):
        return np.nan
    sigma_cont = float(np.std(flux_edges)) if len(flux_edges) > 2 else 0.01
    fwhm_A     = float(popt[2]) * 2.355        # popt[2] = sigma in all models
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
                 worklist: pd.DataFrame) -> pd.DataFrame:
    """Fit every line in worklist; return results DataFrame."""
    window     = PIPELINE['fit_window_A']
    ew_min     = PIPELINE['ew_min_mA']
    ew_max     = PIPELINE['ew_max_mA']
    blue_warn  = PIPELINE['blue_edge_warn_A']
    edge_frac  = PIPELINE['cont_edge_frac']
    edge_A     = window * edge_frac

    results = []
    ni_done = False   # track whether Ni I 6300.336 was already measured

    n = len(worklist)
    print(f"  Fitting {n} lines …")
    for i, row in worklist.iterrows():
        if i % 100 == 0:
            print(f"    {i}/{n}", end='\r', flush=True)

        elem     = str(row['element'])
        ion      = str(row['ion'])
        line_wav = float(row['wavelength_air_A'])
        blend    = bool(row.get('blend_flag', False))
        priority = int(row.get('priority', 0))

        # ── O I 6300.304 → two-component fit ─────────────────────────────────
        is_oi = (elem == 'O' and ion == 'I' and abs(line_wav - 6300.304) < 0.1)
        if is_oi:
            mask = ((solar_wav >= line_wav - window) &
                    (solar_wav <= line_wav + window))
            if mask.sum() < 10:
                continue
            wav_w = solar_wav[mask]
            flux_renorm, _ = _local_renorm(solar_wav[mask],
                                           solar_flux[mask], 6300.304)
            p2, cov2, chi2_2 = _fit_two_component(wav_w, flux_renorm)

            if p2 is not None:
                x_f = np.linspace(wav_w[0], wav_w[-1], 5000)
                # O I component EW (params 0–3)
                vp1 = voigt_profile(0.0, p2[2], p2[3]) or 1.0
                y_oi = 1.0 - p2[1] * voigt_profile(x_f - p2[0], p2[2], p2[3]) / vp1
                ew_oi = float(np.trapz(np.clip(1.0 - y_oi, 0, None), x_f)) * 1000
                # Ni I component EW (params 4–7)
                vp2 = voigt_profile(0.0, p2[6], p2[7]) or 1.0
                y_ni = 1.0 - p2[5] * voigt_profile(x_f - p2[4], p2[6], p2[7]) / vp2
                ew_ni = float(np.trapz(np.clip(1.0 - y_ni, 0, None), x_f)) * 1000
                profile_t = 'voigt_2comp'
            else:
                # Fallback: single profile for total blend
                flux_renorm_s, _ = _local_renorm(solar_wav[mask],
                                                 solar_flux[mask], 6300.304)
                popt, pcov, profile_t, chi2_2 = _fit_profile(
                    solar_wav[mask], flux_renorm_s, 6300.304)
                ew_oi = _integrate_profile(solar_wav[mask], popt, profile_t) if popt is not None else np.nan
                ew_ni = np.nan
                p2 = popt

            results.append(dict(
                element='O', ion='I', wavelength_air_A=6300.304,
                ew_mA=round(ew_oi, 2) if not np.isnan(ew_oi) else np.nan,
                ew_err_mA=np.nan,
                profile_type=profile_t, chi2=round(chi2_2, 4) if not np.isnan(chi2_2) else np.nan,
                blend_flag=True,
                notes='O_Ni_blend; two_component_fit; NLTE_flag',
            ))
            results.append(dict(
                element='Ni', ion='I', wavelength_air_A=6300.336,
                ew_mA=round(ew_ni, 2) if not np.isnan(ew_ni) else np.nan,
                ew_err_mA=np.nan,
                profile_type=profile_t, chi2=round(chi2_2, 4) if not np.isnan(chi2_2) else np.nan,
                blend_flag=True,
                notes='O_Ni_blend; blend_partner_of_O_I_6300.304',
            ))
            ni_done = True
            continue

        # ── Ni I 6300.336 already measured above ─────────────────────────────
        if (elem == 'Ni' and ion == 'I' and abs(line_wav - 6300.336) < 0.1
                and ni_done):
            continue

        # ── Extract and re-normalise window ───────────────────────────────────
        mask = ((solar_wav >= line_wav - window) &
                (solar_wav <= line_wav + window))
        if mask.sum() < 10:
            continue

        wav_w  = solar_wav[mask]
        flux_w = solar_flux[mask]
        flux_r, _ = _local_renorm(wav_w, flux_w, line_wav)

        # Edge pixels for noise estimate
        edge_mask = (wav_w <= line_wav - window + edge_A) | (wav_w >= line_wav + window - edge_A)
        flux_edges = flux_r[edge_mask]

        # ── Fit ───────────────────────────────────────────────────────────────
        popt, pcov, profile_t, chi2 = _fit_profile(wav_w, flux_r, line_wav)
        if popt is None:
            continue

        ew_mA  = _integrate_profile(wav_w, popt, profile_t)
        ew_err = _ew_error(flux_edges, ew_mA, popt, pcov)

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
        if line_wav < blue_warn:
            notes.append('low_snr_blue_edge')

        results.append(dict(
            element=elem, ion=ion,
            wavelength_air_A=round(line_wav, 5),
            ew_mA=round(ew_mA, 2),
            ew_err_mA=round(ew_err, 2) if not np.isnan(ew_err) else np.nan,
            profile_type=profile_t,
            chi2=round(chi2, 4) if not np.isnan(chi2) else np.nan,
            blend_flag=blend,
            notes='; '.join(notes),
        ))

    print(f"    {n}/{n} — done.         ")
    return pd.DataFrame(results)


# ── Diagnostic plot ───────────────────────────────────────────────────────────

def _plot_tier1_fe(solar_wav: np.ndarray, solar_flux: np.ndarray,
                   results: pd.DataFrame, out_path) -> None:
    """2 × 3 grid of Tier 1 Fe anchor line fits + [O I] 6300 blend panel."""
    window = PIPELINE['fit_window_A']

    fig = plt.figure(figsize=(15, 9))
    fig.suptitle(
        f"Tier 1 line profile fits — Solar HARPS  |  "
        f"Teff={STAR_SOLAR['teff_K']:.0f} K  log g={STAR_SOLAR['logg']}  "
        f"ξ={STAR_SOLAR['vturb_kms']} km/s\n"
        f"{STAR_SOLAR['citation']}  |  ESO {STAR_SOLAR['program_id']}",
        fontsize=10,
    )
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
        if elem == '[O':
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


# ── Entry point ───────────────────────────────────────────────────────────────

def run(star_id: str = 'solar') -> pd.DataFrame:
    """
    Called by run_pipeline.py. Returns the EW results DataFrame.
    star_id: any value containing 'solar' routes to the solar calibration path.
    """
    print(f"\n{'='*60}\n  lines_fit — {star_id}\n{'='*60}")

    if 'solar' not in star_id.lower():
        raise NotImplementedError(
            f"lines_fit not yet implemented for '{star_id}'."
        )

    # ── Load inputs ───────────────────────────────────────────────────────────
    print(f"\n[1/5] Loading inputs")
    norm_path  = PATHS['solar_normalized']
    lines_path = PATHS['linelist_solar']

    if not norm_path.exists():
        raise FileNotFoundError(
            f"Normalised spectrum not found: {norm_path}\n"
            "Run spectra_normalize.py first."
        )

    norm = pd.read_csv(norm_path)
    solar_wav  = norm['wavelength_air_A'].to_numpy()
    solar_flux = norm['flux_normalized'].to_numpy()
    print(f"  Spectrum: {len(solar_wav):,} pixels, "
          f"{solar_wav[0]:.1f}–{solar_wav[-1]:.1f} Å")

    lines = pd.read_csv(lines_path, low_memory=False)
    lines = lines[lines['priority'] > 0].copy()
    print(f"  Line list: {len(lines):,} priority>0 lines")

    # ── Build worklist ────────────────────────────────────────────────────────
    print(f"\n[2/5] Building worklist")
    worklist = _build_worklist(lines)
    print(f"  {len(worklist)} features to measure "
          f"({(worklist['blend_flag'].astype(str).str.lower()=='true').sum()} "
          f"special/blended force-included)")

    # ── Measure EWs ───────────────────────────────────────────────────────────
    print(f"\n[3/5] Measuring EWs")
    results = _measure_all(solar_wav, solar_flux, worklist)
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
    out_csv  = PATHS['solar_ew']
    out_plot = PATHS['solar_ew_diagnostic']

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    col_order = ['element', 'ion', 'wavelength_air_A', 'ew_mA', 'ew_err_mA',
                 'profile_type', 'chi2', 'blend_flag', 'notes']
    out_df = results[col_order].copy()
    out_df['notes'] = out_df['notes'].fillna('')
    out_df.to_csv(out_csv, index=False)
    print(f"\n  Saved → {out_csv.name}  ({len(results)} rows)")

    # ── Diagnostic plot ───────────────────────────────────────────────────────
    print(f"\n[5/5] Generating diagnostic plot")
    out_plot.parent.mkdir(parents=True, exist_ok=True)
    _plot_tier1_fe(solar_wav, solar_flux, results, out_plot)

    print(f"\n✓  lines_fit complete.")
    return results


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    run(star_id='solar')
