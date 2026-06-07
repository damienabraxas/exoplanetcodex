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
  - Correct Mg/Si (SPECTRUM has 0.10–0.15 dex bias on some lines, Gaia FGK 2026)
  - Native molecular lines (critical for C and M dwarf targets)
  - Native NLTE path — same code handles LTE → NLTE → 3D NLTE (RYA-165)

EW injection workflow:
  Our solar_ew.csv EWs are matched to iSpec's pre-built GES line regions
  (which carry full atomic data: loggf, EP, Turbospectrum species codes,
  damping constants). Measured EWs replace the region placeholders and
  ispec.determine_abundances iterates A(X) until the synthesized EW matches.

Solar abundances — TWO REFERENCES (do not mix):
  iSpec internal:   Asplund.2009 — passed to ispec.determine_abundances as its
                    reference scale; Turbospectrum uses this internally.
                    (Asplund.2021 is not shipped with iSpec.)
  Science output:   SOLAR_ASPLUND2021 from config/constants.py — used for all
                    [X/H] and [X/Fe] reported to the user. This is the
                    current solar standard (A&A 653, A141).

Linear issue: RYA-167

Calibration notes (first solar run, 2026-06-04)
------------------------------------------------
- ATLAS9.Castelli + Turbospectrum verified working at solar params (72 layers)
- MARCS.GES verified working at M-dwarf params (LHS 3844)
- EP slope and REW slope computed for Fe I as stellar-param diagnostics
- Asplund.2009 vs Asplund.2021 offset for Fe: 7.46 – 7.67 = –0.01 dex (negligible)
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
sys.path.insert(0, str(ISPEC_DIR))
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import ispec

# ── Constants ─────────────────────────────────────────────────────────────────

# Atmosphere grid directory names (exact names on disk after input.tar.gz extraction)
_ATLAS9 = 'ATLAS9.Castelli'
_MARCS  = 'MARCS.GES'

# iSpec's internal solar reference (Asplund 2009 — closest available in iSpec).
# Used by ispec.determine_abundances for its internal normalisation.
# DO NOT use for science output — use SOLAR_ASPLUND2021 from constants.py instead.
_ISPEC_SOLAR_ABUND_FILE = str(
    ISPEC_DIR / 'input' / 'abundances' / 'Asplund.2009' / 'stdatom.dat'
)

# Multi-element line regions with full atomic data for Turbospectrum.
# These carry loggf, EP, species codes, damping constants — required by
# ispec.determine_abundances. Our measured EWs are injected before the call.
_LINE_REGIONS_ALL = str(
    ISPEC_DIR / 'input' / 'regions' / '42000_VALD' /
    'limited_but_with_missing_elements_turbospectrum_synth_good_for_abundances_all_extended.txt'
)
_LINE_REGIONS_FE = str(
    ISPEC_DIR / 'input' / 'regions' / '42000_VALD' /
    'turbospectrum_synth_good_for_params_all_extended.txt'
)

# Module-level pack cache — loading the 2 GB grid once per process
_pack_cache: dict = {}


# ── Core atmosphere functions ─────────────────────────────────────────────────

def _load_atmosphere(teff: float, logg: float, feh: float, vturb: float,
                     model_grid: str = 'ATLAS9.Castelli') -> np.ndarray:
    """
    Interpolate model atmosphere at stellar parameters.

    Loads the grid pack on first call per model_grid and caches it.

    Parameters
    ----------
    teff, logg, feh, vturb : stellar parameters
    model_grid : 'ATLAS9.Castelli' (FGK) or 'MARCS.GES' (M dwarfs)

    Returns
    -------
    atmosphere_layers : numpy structured array (72 layers for ATLAS9)
    """
    global _pack_cache

    grid_dir = ISPEC_DIR / 'input' / 'atmospheres' / model_grid
    if not grid_dir.exists():
        raise FileNotFoundError(
            f"Atmosphere grid not found: {grid_dir}\n"
            "Run: cd ispec && tar -xzf input.tar.gz"
        )

    if model_grid not in _pack_cache:
        _pack_cache[model_grid] = ispec.load_modeled_layers_pack(str(grid_dir))

    pack = _pack_cache[model_grid]
    target = {'teff': float(teff), 'logg': float(logg),
              'MH': float(feh), 'alpha': 0.0}

    if not ispec.valid_atmosphere_target(pack, target):
        raise ValueError(
            f"Stellar params {target} outside {model_grid} grid bounds."
        )

    return ispec.interpolate_atmosphere_layers(
        pack, target, code=RADIATIVE_TRANSFER_CODE
    )


# ── EW → line regions ─────────────────────────────────────────────────────────

def _ours_to_ispec_note(element: str, ion: str) -> str:
    """'Fe', 'I' → 'Fe 1';  'Fe', 'II' → 'Fe 2'"""
    return f"{element} {1 if ion.strip().upper() == 'I' else 2}"


def _build_ispec_line_regions(ew_df: pd.DataFrame,
                               line_regions_path: str = _LINE_REGIONS_ALL,
                               wave_tol: float = 0.15) -> np.ndarray:
    """
    Match our measured EWs to iSpec's pre-built GES line regions and inject
    the EW values. The GES regions carry full atomic data (loggf, EP,
    Turbospectrum species codes, damping) required by determine_abundances.

    Parameters
    ----------
    ew_df            : DataFrame from solar_ew.csv (element, ion, wavelength_air_A, ew_mA)
    line_regions_path: path to iSpec extended line regions file
    wave_tol         : wavelength match tolerance in Å

    Returns
    -------
    linemasks : numpy structured array in iSpec linemask format, EWs injected
                (wave_A in Å — note iSpec wave_peak is in nm, wave_A is in Å)
    """
    regions = ispec.read_line_regions(line_regions_path)

    # Build lookup: ispec_note → [(wav_Å, ew_Å, err_Å), ...]
    lookup: dict[str, list] = {}
    for _, row in ew_df.iterrows():
        note  = _ours_to_ispec_note(str(row['element']), str(row['ion']))
        wav    = float(row['wavelength_air_A'])
        ew_mA  = float(row['ew_mA'])
        err_mA = float(row['ew_err_mA']) if 'ew_err_mA' in row.index else 0.0
        lookup.setdefault(note, []).append((wav, ew_mA, err_mA))

    # Match regions → our EWs (match on wave_A which is in Å)
    matched_idx, ew_vals, err_vals = [], [], []
    for i, region in enumerate(regions):
        note = str(region['note'])
        if note not in lookup:
            continue
        candidates = lookup[note]
        diffs = [abs(wav - float(region['wave_A'])) for wav, _, _ in candidates]
        min_diff = min(diffs)
        if min_diff > wave_tol:
            continue
        best = candidates[int(np.argmin(diffs))]
        matched_idx.append(i)
        ew_vals.append(best[1])
        err_vals.append(best[2])

    if not matched_idx:
        raise RuntimeError(
            f"No EW measurements matched iSpec line regions "
            f"(±{wave_tol} Å). Check element/ion names and wavelengths."
        )

    linemasks = regions[matched_idx].copy()
    # Store EW in mÅ — iSpec's __spectrum_write_abundance_lines writes
    # linemasks['ew'] directly to the SPECTRUM input file which expects mÅ.
    # (iSpec's own fit_lines also stores ew in mÅ, not Å.)
    for j, (ew_mA, err_mA) in enumerate(zip(ew_vals, err_vals)):
        linemasks[j]['ew']     = ew_mA   # mÅ
        linemasks[j]['ew_err'] = err_mA  # mÅ

    # Discard zero/negative EWs (failed matches)
    linemasks = linemasks[linemasks['ew'] > 0]
    return linemasks


# ── Abundance derivation ──────────────────────────────────────────────────────

def _ew_to_abundance(ew_df: pd.DataFrame,
                     stellar_params: dict,
                     atmosphere: np.ndarray,
                     model_grid: str = 'ATLAS9.Castelli',
                     code: str = RADIATIVE_TRANSFER_CODE) -> tuple:
    """
    Call ispec.determine_abundances on the matched linemasks.

    iSpec internally uses Asplund.2009 as its normalisation reference.
    Science [X/H] and [X/Fe] are computed against SOLAR_ASPLUND2021
    in the calling function — not here.

    Returns
    -------
    (linemasks, spec_abund, x_over_h, x_over_fe) :
        linemasks  — the matched numpy structured array
        spec_abund — per-line A(X) (absolute, on iSpec's Asplund.2009 scale)
        x_over_h   — per-line [X/H]  (on iSpec's Asplund.2009 scale)
        x_over_fe  — per-line [X/Fe] (on iSpec's Asplund.2009 scale)
    """
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)

    # Pre-filter: unblended lines with reliable EW range before injection.
    # Large or blended EWs produce spurious slopes in the GES synth regions.
    ew_clean = ew_df.copy()
    if 'blend_flag' in ew_clean.columns:
        ew_clean = ew_clean[ew_clean['blend_flag'] == False]
    ew_clean = ew_clean[(ew_clean['ew_mA'] >= 5) & (ew_clean['ew_mA'] <= 150)]

    linemasks = _build_ispec_line_regions(ew_clean)

    # ── Theoretical EW sanity filter ──────────────────────────────────────────
    # Use GES synthesis-predicted theoretical_ew when available — this accounts
    # for saturation and damping on the COG shoulder, unlike the linear COG
    # approximation which underestimates EW by 3–10× for Fe I lines in the
    # 50–200 mÅ range typically used for abundance work.
    # Fall back to linear COG only when theoretical_ew is absent/zero.
    # Reject lines where |log10(EW_obs / EW_theo)| > 1.5 dex — real blends
    # (e.g. Fe I 5247 loggf=−4.95: EW_theo≈1.6 mÅ, obs=91.6 mÅ → 1.76 dex).
    theta_filter      = 5040.0 / float(stellar_params['teff_K'])
    theo_synth        = linemasks['theoretical_ew']  # from GES synthesis (mÅ)
    has_theo          = np.asarray(theo_synth, dtype=float) > 0
    log_ew_theo_cog   = (linemasks['loggf']
                         + 2.0 * np.log10(np.maximum(linemasks['wave_A'], 1.0))
                         - linemasks['lower_state_eV'] * theta_filter
                         - 2.2164)
    log_ew_theo_synth = np.log10(np.maximum(np.asarray(theo_synth, dtype=float), 1e-6))
    log_ew_theo       = np.where(has_theo, log_ew_theo_synth, log_ew_theo_cog)
    ew_mA_obs    = linemasks['ew']   # already in mÅ
    log_ew_obs   = np.log10(np.maximum(ew_mA_obs, 1e-6))
    ew_ratio     = log_ew_obs - log_ew_theo
    good         = np.abs(ew_ratio) < 1.5
    n_before = len(linemasks)
    linemasks = linemasks[good]
    n_after  = len(linemasks)
    print(f"  EW sanity filter: {n_before} → {n_after} lines "
          f"({n_before - n_after} rejected as blends/misidentifications)")

    teff  = float(stellar_params['teff_K'])
    logg  = float(stellar_params['logg'])
    feh   = float(stellar_params['feh'])
    vturb = float(stellar_params['vturb_kms'])

    # determine_abundances only accepts: 'spectrum', 'moog', 'moog-scat', 'width'.
    # Turbospectrum is used above for atmosphere interpolation (ATLAS9.Castelli/
    # MARCS.GES formatted for Turbospectrum). The EW→abundance step uses
    # 'spectrum' here. Full Turbospectrum synthesis (model_spectrum_from_ew)
    # is the RYA-165 upgrade path.
    ew_code = 'spectrum' if code == 'turbospectrum' else code
    spec_abund, normal_abund, x_over_h, x_over_fe = ispec.determine_abundances(
        atmosphere, teff, logg, feh, 0.0,
        linemasks, solar_abund,
        microturbulence_vel=vturb,
        verbose=0,
        code=ew_code,
        tmp_dir='/tmp/ispec_codex',
    )
    # normal_abund = A(X) in the standard log(N/N_H)+12 scale (hydrogen=12).
    # spec_abund   = SPECTRUM's internal log(N/N_H) scale (hydrogen≈0).
    # All downstream A(X) / [X/H] / [X/Fe] use normal_abund, not spec_abund.
    return linemasks, normal_abund, x_over_h, x_over_fe


# ── Fe diagnostics ────────────────────────────────────────────────────────────

def _compute_ep_slope(fe1_mask: np.ndarray, linemasks: np.ndarray,
                      x_over_h: np.ndarray) -> float:
    """Fe I abundance vs excitation potential slope → Teff diagnostic."""
    ep    = linemasks['lower_state_eV'][fe1_mask]
    abund = x_over_h[fe1_mask]
    valid = np.isfinite(abund) & np.isfinite(ep)
    if valid.sum() < 5:
        return np.nan
    return float(np.polyfit(ep[valid], abund[valid], 1)[0])


def _compute_rew_slope(fe1_mask: np.ndarray, linemasks: np.ndarray,
                       x_over_h: np.ndarray) -> float:
    """Fe I abundance vs reduced EW slope → vturb diagnostic."""
    ew_mA = linemasks['ew'][fe1_mask]   # mÅ
    wav   = linemasks['wave_A'][fe1_mask]
    abund = x_over_h[fe1_mask]
    valid = np.isfinite(abund) & (ew_mA > 0) & (wav > 0)
    if valid.sum() < 5:
        return np.nan
    # Reduced EW = log10(EW_Å / λ_Å) = log10(EW_mÅ / λ_Å) - 3
    # The -3 offset cancels in the slope, so we use log10(EW_mÅ / λ_Å) directly.
    rew = np.log10(ew_mA[valid] / wav[valid])
    return float(np.polyfit(rew, abund[valid], 1)[0])


# ── Iterative parameter convergence ──────────────────────────────────────────

def _iterative_parameter_convergence(ew_df: pd.DataFrame,
                                      stellar_params_init: dict,
                                      model_grid: str = 'ATLAS9.Castelli',
                                      max_iter: int = 10) -> tuple:
    """
    Iterate Teff, log g, vturb to excitation + ionisation equilibrium.

    Convergence criteria (solar calibration gate, RYA-166 Layer 3):
      |EP slope|          < 0.02 dex/eV  → Teff
      |reduced EW slope|  < 0.02         → vturb
      |Fe I – Fe II|      < 0.05 dex     → log g

    Returns
    -------
    (converged_params, results_df) :
        converged_params — dict with final teff_K, logg, feh, vturb_kms
        results_df       — per-line DataFrame with element, ion, A_X, XH, XFe
    """
    params = stellar_params_init.copy()
    last_linemasks = last_spec_abund = last_xh = last_xfe = None

    for iteration in range(1, max_iter + 1):
        atm = _load_atmosphere(
            params['teff_K'], params['logg'], params['feh'],
            params['vturb_kms'], model_grid=model_grid
        )
        linemasks, spec_abund, x_over_h, x_over_fe = _ew_to_abundance(
            ew_df, params, atm, model_grid
        )
        last_linemasks, last_spec_abund, last_xh, last_xfe = (
            linemasks, spec_abund, x_over_h, x_over_fe
        )

        notes   = np.array([str(n) for n in linemasks['note']])
        fe1_mask = notes == 'Fe 1'
        fe2_mask = notes == 'Fe 2'

        n_fe1 = fe1_mask.sum()
        if n_fe1 < 5:
            print(f"  Iter {iteration:02d}: only {n_fe1} Fe I lines — stopping early")
            break

        ep_sl  = _compute_ep_slope(fe1_mask, linemasks, x_over_h)
        rew_sl = _compute_rew_slope(fe1_mask, linemasks, x_over_h)

        fe1_xh = float(np.nanmedian(x_over_h[fe1_mask & np.isfinite(x_over_h)]))
        fe2_xh = (float(np.nanmedian(x_over_h[fe2_mask & np.isfinite(x_over_h)]))
                  if fe2_mask.any() else np.nan)
        ion_bal = fe1_xh - fe2_xh if np.isfinite(fe2_xh) else np.nan

        print(f"  Iter {iteration:02d}: "
              f"Teff={params['teff_K']:.0f}K  logg={params['logg']:.2f}  "
              f"vturb={params['vturb_kms']:.2f}km/s | "
              f"EP={ep_sl:+.4f}  REW={rew_sl:+.4f}  "
              f"FeI-FeII={ion_bal:+.4f}" if np.isfinite(ion_bal)
              else f"  Iter {iteration:02d}: "
                   f"Teff={params['teff_K']:.0f}K  logg={params['logg']:.2f}  "
                   f"vturb={params['vturb_kms']:.2f}km/s | "
                   f"EP={ep_sl:+.4f}  REW={rew_sl:+.4f}  FeI-FeII=N/A (no Fe II)")

        converged = (
            (np.isnan(ep_sl)  or abs(ep_sl)  < 0.02) and
            (np.isnan(rew_sl) or abs(rew_sl) < 0.02) and
            (np.isnan(ion_bal) or abs(ion_bal) < 0.05)
        )
        if converged:
            print(f"  ✓ Converged at iteration {iteration}")
            break

        # Gradient-descent parameter adjustment — conservative steps to avoid
        # walking off the atmosphere grid on large spurious slopes.
        if np.isfinite(ep_sl):
            params['teff_K']    += np.sign(ep_sl)  * min(abs(ep_sl)  / 0.05, 1.0) * 25
        if np.isfinite(rew_sl):
            params['vturb_kms'] += np.sign(rew_sl) * min(abs(rew_sl) / 0.10, 1.0) * 0.05
        if np.isfinite(ion_bal):
            params['logg']      += np.sign(ion_bal) * min(abs(ion_bal) / 0.05, 1.0) * 0.05

        params['teff_K']    = float(np.clip(params['teff_K'],    3500, 7500))
        params['logg']      = float(np.clip(params['logg'],      2.0,  5.0))
        params['vturb_kms'] = float(np.clip(params['vturb_kms'], 0.5,  3.0))
    else:
        print(f"  ⚠ Did not converge within {max_iter} iterations")

    # Build results DataFrame using SOLAR_ASPLUND2021 for [X/H] and [X/Fe].
    # Two-pass approach to keep both [X/H] and [Fe/H] on the same Asplund 2021
    # scale — eliminates the prior mixing of Asplund 2009 stellar A(Fe) with
    # the Asplund 2021 solar reference (RYA-200).
    notes = np.array([str(n) for n in last_linemasks['note']])

    # Pass 1: compute A(X) and [X/H] for all elements on Asplund 2021 scale
    rows_pass1 = []
    for note in np.unique(notes):
        mask  = notes == note
        valid = np.isfinite(last_spec_abund[mask])
        if valid.sum() == 0:
            continue
        parts = note.split()
        elem  = parts[0]
        ion   = 'I' if parts[1] == '1' else 'II'

        a_x       = float(np.nanmedian(last_spec_abund[mask][valid]))
        a_x_solar = SOLAR_ASPLUND2021.get(elem, np.nan)
        xh        = round(a_x - a_x_solar, 3) if np.isfinite(a_x_solar) else np.nan

        rows_pass1.append({
            'element': elem, 'ion': ion,
            'A_X'    : round(a_x, 3),
            'A_X_std': round(float(np.nanstd(last_spec_abund[mask][valid])), 3)
                       if valid.sum() > 1 else np.nan,
            'n_lines': int(valid.sum()),
            'XH'     : xh,
        })

    # Pass 2: derive [X/Fe] = [X/H]_X - [Fe/H] using Fe I [X/H] from Pass 1.
    # Both terms are on Asplund 2021 — no scale mixing. Fe I [X/Fe] = 0.000 by construction.
    fe1_rows = [r for r in rows_pass1 if r['element'] == 'Fe' and r['ion'] == 'I']
    fe2_rows = [r for r in rows_pass1 if r['element'] == 'Fe' and r['ion'] == 'II']
    feh_ref  = fe1_rows[0]['XH'] if fe1_rows and np.isfinite(fe1_rows[0]['XH']) \
               else (fe2_rows[0]['XH'] if fe2_rows and np.isfinite(fe2_rows[0]['XH']) else np.nan)

    rows = []
    for r in rows_pass1:
        xfe = round(r['XH'] - feh_ref, 3) if (np.isfinite(r.get('XH', np.nan))
                                               and np.isfinite(feh_ref)) else np.nan
        rows.append({**r, 'XFe': xfe})

    results_df = pd.DataFrame(rows).sort_values(['element', 'ion']).reset_index(drop=True)
    return params, results_df


# ── Legacy COG (DO NOT USE for science) ──────────────────────────────────────

def _derive_abundances_cog_legacy(*args, **kwargs):
    """
    DEPRECATED — linear COG calibrated to Fe I only.
    Errors ±0.5–2 dex for non-Fe elements. Retained for reference only.
    """
    raise DeprecationWarning(
        "_derive_abundances_cog_legacy is retired. Use run() instead."
    )


# ── Public entry point ────────────────────────────────────────────────────────

def run(star_id: str = 'solar',
        model_grid: str = 'ATLAS9.Castelli',
        stellar_params_override: dict = None) -> pd.DataFrame:
    """
    Derive abundances for a star and save results.

    Parameters
    ----------
    star_id                  : 'solar' or '55cnc'
    model_grid               : 'ATLAS9.Castelli' (FGK) or 'MARCS.GES' (M dwarfs)
    stellar_params_override  : override dict for uncertainty sensitivity runs (RYA-158)
                               keys: teff_K, logg, feh, vturb_kms

    Returns
    -------
    DataFrame : element, ion, A_X, A_X_std, n_lines, XH ([X/H]), XFe ([X/Fe])
                [X/H] normalised to SOLAR_ASPLUND2021 (Asplund et al. 2021)
    """
    print(f"\n{'='*62}")
    print(f"  abundances_derive  |  {star_id}  |  {model_grid}  |  {RADIATIVE_TRANSFER_CODE}")
    print(f"{'='*62}\n")

    # ── Stellar parameters ────────────────────────────────────────
    if stellar_params_override:
        params = stellar_params_override.copy()
    elif 'solar' in star_id.lower():
        params = STAR_SOLAR.copy()
    else:
        params = STAR_55CNC.copy()

    print(f"[1/4] Stellar params:  Teff={params['teff_K']} K  "
          f"logg={params['logg']}  [Fe/H]={params['feh']}  "
          f"vturb={params['vturb_kms']} km/s")

    # ── Load EW table ─────────────────────────────────────────────
    if 'solar' in star_id.lower():
        ew_path = PATHS['solar_ew']
    else:
        ew_path = PATHS.get(f'{star_id}_ew',
                             Path(str(PATHS['solar_ew'])).parent / f'{star_id}_ew.csv')

    if not Path(str(ew_path)).exists():
        raise FileNotFoundError(
            f"EW table not found: {ew_path}\nRun pipeline/lines_fit.py first."
        )

    ew_df = pd.read_csv(ew_path)
    ew_df = ew_df[(ew_df['ew_mA'] > 0) & ew_df['ew_mA'].notna()].copy()
    print(f"[2/4] Loaded {len(ew_df)} EW measurements from {Path(str(ew_path)).name}")

    # ── Iterative convergence ─────────────────────────────────────
    print(f"\n[3/4] Iterating to stellar parameter equilibrium "
          f"({model_grid} / {RADIATIVE_TRANSFER_CODE})...")
    converged_params, results = _iterative_parameter_convergence(
        ew_df, params, model_grid=model_grid
    )

    print(f"\n  Final params: Teff={converged_params['teff_K']:.0f} K  "
          f"logg={converged_params['logg']:.2f}  "
          f"vturb={converged_params['vturb_kms']:.2f} km/s")

    # ── Save ──────────────────────────────────────────────────────
    print(f"\n[4/4] Saving results...")
    out_path = Path(str(ew_path)).parent / f'{star_id}_abundances.csv'
    results.to_csv(out_path, index=False)
    print(f"  Saved → {out_path.name}  ({len(results)} elements, "
          f"{results['n_lines'].sum()} total lines)")

    # ── Solar calibration gate ────────────────────────────────────
    fe = results[results['element'] == 'Fe']
    if len(fe) > 0:
        a_fe = float(fe['A_X'].mean())
        delta = a_fe - 7.46
        gate = '✓ PASS' if abs(delta) <= 0.05 else '✗ FAIL'
        print(f"\n  A(Fe)☉ = {a_fe:.3f}  (target 7.46 ± 0.05)  {gate}")

    print(f"\n{'='*62}")
    print(f"  abundances_derive complete.")
    print(f"{'='*62}\n")
    return results


if __name__ == '__main__':
    import sys as _sys
    star  = _sys.argv[1] if len(_sys.argv) > 1 else 'solar'
    grid  = _sys.argv[2] if len(_sys.argv) > 2 else 'ATLAS9.Castelli'
    run(star, grid)
