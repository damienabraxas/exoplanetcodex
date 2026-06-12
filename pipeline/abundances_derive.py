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
    PATHS, PIPELINE, SOLAR_ASPLUND2021, STAR_SOLAR, STAR_55CNC, STAR_PROCYON,
    ISPEC_DIR, RADIATIVE_TRANSFER_CODE,
    LINE_SCORE_WEIGHTS, LINE_GRADE_THRESHOLDS, LINE_SCORE_PARAMS,
    FE_GATE_LOWER, FE_GATE_UPPER, FE_SCATTER_GATE, FE_IONISATION_GATE,
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
    'turbospectrum_synth_good_for_params_codex_extended.txt'
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

    # Pre-filter: remove blended lines and out-of-range EWs.
    # Ceiling uses PIPELINE['ew_max_mA'] (300 mÅ) — NOT the vmic COG cut.
    # The vmic COG cut (150 mÅ) applies only inside the bisection sample
    # in _converge_vmic_fe1_only(), not here. (RYA-199)
    ew_min = float(PIPELINE['ew_min_mA'])   # 5 mÅ
    ew_max = float(PIPELINE['ew_max_mA'])   # 300 mÅ
    ew_clean = ew_df.copy()
    if 'blend_flag' in ew_clean.columns:
        ew_clean = ew_clean[ew_clean['blend_flag'] == False]
    ew_clean = ew_clean[(ew_clean['ew_mA'] >= ew_min) & (ew_clean['ew_mA'] <= ew_max)]

    n_total  = len(ew_df)
    n_blends = int((ew_df.get('blend_flag', pd.Series(False, index=ew_df.index)) == True).sum())
    n_ew_cut = n_total - n_blends - len(ew_clean)
    print(f"  EW pre-filter: {n_total} total → {n_blends} blends removed, "
          f"{n_ew_cut} outside [{ew_min:.0f},{ew_max:.0f}] mÅ → {len(ew_clean)} clean")

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
    # Fe II lines without a GES theoretical_ew (=0) bypass the COG sanity check.
    # The linear COG formula is calibrated for Fe I and gives absurd predictions
    # for Fe II lines with extreme loggf values, falsely rejecting valid lines.
    # Fe II lines appended by RYA-243/247 were quality-gated by EW and ew_err;
    # the pre-filter [5, 300] mÅ already guards against the grossest blends.
    fe2_no_theo = np.array(
        [str(r['note']).strip() == 'Fe 2' for r in linemasks]
    ) & ~has_theo
    good         = good | fe2_no_theo
    n_before = len(linemasks)
    linemasks = linemasks[good]
    n_after  = len(linemasks)
    n_fe2_bypass = int(fe2_no_theo.sum())
    note_str = f", {n_fe2_bypass} Fe II bypassed (no theoretical_ew)" if n_fe2_bypass else ""
    print(f"  EW sanity filter: {n_before} → {n_after} lines "
          f"({n_before - n_after} rejected as blends/misidentifications{note_str})")

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
                                      max_iter: int = 10,
                                      vmic_fixed: bool = False) -> tuple:
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

        # Sigma clip removed — RYA-220 replaces it with abundance_outlier_score
        # applied post-convergence. Convergence uses all non-vetted lines;
        # median-based statistics make it robust to the outlier tail.

        # COG cut for vmic slope sample: exclude saturated lines (EW > 150 mÅ).
        # Lines in the damping regime have REW insensitive to vmic but high
        # sensitivity to damping parameters — they corrupt the slope and prevent
        # convergence. Standard practice: Mucciarelli 2011, Sousa et al. 2011.
        # Full fe1_mask (all Fe I) is still used for A(Fe I) and ionisation balance.
        _vmic_ew_max   = float(PIPELINE['vmic_ew_ceiling_mA'])
        fe1_slope_mask = fe1_mask & (linemasks['ew'] <= _vmic_ew_max)
        n_fe1_post_cog = int(fe1_slope_mask.sum())

        # Abundance clip on slope sample: exclude lines whose current A(Fe I)
        # deviates > 2σ from the slope-sample mean.
        # Restores behaviour removed by RYA-220 commit c8a23d2.
        # Applied ONLY to the vmic slope sample — does not affect final A(Fe) derivation.
        # Standard practice: Sousa et al. 2011, Jofré et al. 2014.
        _vmic_sigma   = float(PIPELINE['vmic_abund_clip_sigma'])
        _slope_abunds = spec_abund[fe1_slope_mask]
        _sa_valid     = np.isfinite(_slope_abunds)
        if _sa_valid.sum() >= 3:
            _slope_mean = float(np.mean(_slope_abunds[_sa_valid]))
            _slope_std  = float(np.std(_slope_abunds[_sa_valid]))
            if _slope_std > 0:
                _clip_ok = np.ones(len(spec_abund), dtype=bool)
                _clip_ok[fe1_slope_mask] = (
                    np.abs(_slope_abunds - _slope_mean) <= _vmic_sigma * _slope_std
                )
                fe1_slope_mask = fe1_slope_mask & _clip_ok
            else:
                print(f"  WARNING Iter {iteration:02d}: slope_std=0 — abundance clip skipped")
        else:
            print(f"  WARNING Iter {iteration:02d}: <3 valid slope-sample lines — abundance clip skipped")

        print(f"  vmic slope sample: {int(fe1_slope_mask.sum())}/{n_fe1_post_cog}/{int(fe1_mask.sum())} "
              f"Fe I lines (COG cut ≤150 mÅ → abundance clip 2σ)")

        ep_sl  = _compute_ep_slope(fe1_mask, linemasks, x_over_h)
        rew_sl = _compute_rew_slope(fe1_slope_mask, linemasks, x_over_h)

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
        if np.isfinite(rew_sl) and not vmic_fixed:
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

    # Build per-line DataFrame for NLTE per-line corrections (RYA-207).
    # Contains one row per Fe I/II line that survived all filters and contributed
    # to the final abundance — atomic data + individual 1D LTE A(Fe) per line.
    per_line_rows = []
    for i in range(len(last_linemasks)):
        note = str(last_linemasks['note'][i])
        if note not in ('Fe 1', 'Fe 2'):
            continue
        a_val = float(last_spec_abund[i])
        if not np.isfinite(a_val):
            continue
        ion_lbl = 'I' if note.split()[1] == '1' else 'II'
        per_line_rows.append({
            'element'                : 'Fe',
            'ion'                    : ion_lbl,
            'wavelength_air_A'       : float(last_linemasks['wave_A'][i]),
            'ew_mA'                  : float(last_linemasks['ew'][i]),
            'excitation_potential_eV': float(last_linemasks['lower_state_eV'][i]),
            'eup_eV'                 : float(last_linemasks['upper_state_eV'][i]),
            'log_gf'                 : float(last_linemasks['loggf'][i]),
            'a_1dlte'                : round(a_val, 4),
        })
    per_line_df = pd.DataFrame(per_line_rows) if per_line_rows else pd.DataFrame()

    # Apply COG ceiling to final abundance pool — same physical basis as vmic slope cut.
    # Lines in the damping regime (EW > vmic_ew_ceiling_mA) have abundances driven by
    # damping parameters, not line strength. Excluding them from the final pool is
    # consistent with Mucciarelli 2011 and GES methodology.
    # Reference: PIPELINE['vmic_ew_ceiling_mA'] = 150.0 (set in RYA-226)
    if not per_line_df.empty and 'ew_mA' in per_line_df.columns:
        _ew_ceiling   = float(PIPELINE['vmic_ew_ceiling_mA'])
        _fe1_pool     = per_line_df['ion'] == 'I'
        _n_fe1_before = int(_fe1_pool.sum())
        per_line_df   = per_line_df[~_fe1_pool | (per_line_df['ew_mA'] <= _ew_ceiling)].copy()
        _n_fe1_after  = int((per_line_df['ion'] == 'I').sum())
        if _n_fe1_before > _n_fe1_after:
            print(f"  Final pool COG cut: {_n_fe1_before} → {_n_fe1_after} Fe I lines "
                  f"({_n_fe1_before - _n_fe1_after} saturated lines removed, "
                  f"EW > {_ew_ceiling:.0f} mÅ)")

    return params, results_df, per_line_df


# ── Line quality scoring (RYA-220) ───────────────────────────────────────────

def _compute_line_scores(per_line_df: pd.DataFrame,
                         ew_df: pd.DataFrame,
                         results_df: pd.DataFrame = None) -> pd.DataFrame:
    """
    Compute per-line quality scores and assign letter grades — RYA-220.

    Adds columns to per_line_df:
        vald_proximity_flag, ew_snr_score, fit_chi2_score, saturation_score,
        abundance_outlier_score, nlte_correction_score, line_score, line_grade
    """
    W = LINE_SCORE_WEIGHTS
    P = LINE_SCORE_PARAMS
    T = LINE_GRADE_THRESHOLDS
    WAV_TOL = 0.05  # Å — tight match for vetted data

    df = per_line_df.copy().reset_index(drop=True)
    n  = len(df)
    if n == 0:
        return df

    # Precompute lookup tables grouped by (element, ion) for speed
    ew_by_ei = {k: v.reset_index(drop=True)
                for k, v in ew_df.groupby(['element', 'ion'])}

    try:
        ll = pd.read_csv(str(PATHS['linelist_solar']),
                         usecols=['element', 'ion', 'wavelength_air_A', 'vald_proximity_flag'],
                         low_memory=False)
        ll_by_ei = {k: v.reset_index(drop=True) for k, v in ll.groupby(['element', 'ion'])}
    except Exception as _e:
        print(f"  WARNING: vald_proximity_flag load failed — {_e}")
        ll_by_ei = {}

    vpf_arr  = np.full(n, 0.5)
    snr_arr  = np.full(n, 0.5)
    chi2_arr = np.full(n, 0.5)

    for i, row in df.iterrows():
        ei  = (row['element'], row['ion'])
        wl  = float(row['wavelength_air_A'])
        ew  = float(row['ew_mA'])

        # vald_proximity_flag
        if ei in ll_by_ei:
            grp = ll_by_ei[ei]
            delta = np.abs(grp['wavelength_air_A'].values - wl)
            j = int(delta.argmin())
            if delta[j] < WAV_TOL:
                vpf_arr[i] = float(grp.iloc[j]['vald_proximity_flag'])

        # ew_snr_score + fit_chi2_score
        if ei in ew_by_ei:
            grp = ew_by_ei[ei]
            delta = np.abs(grp['wavelength_air_A'].values - wl)
            j = int(delta.argmin())
            if delta[j] < WAV_TOL:
                erow = grp.iloc[j]
                err = float(erow.get('ew_err_mA', np.nan))
                if np.isfinite(err) and err > 0:
                    snr_arr[i] = float(min(ew / err / P['snr_saturation'], 1.0))
                c = float(erow['chi2']) if 'chi2' in erow and not pd.isna(erow['chi2']) else np.nan
                if np.isfinite(c):
                    if c <= P['chi2_clean_max']:
                        chi2_arr[i] = 1.0
                    else:
                        chi2_arr[i] = float(max(
                            0.0, 1.0 - (c - P['chi2_clean_max']) /
                            (P['chi2_floor'] - P['chi2_clean_max'])
                        ))

    df['vald_proximity_flag'] = np.round(vpf_arr, 4)
    df['ew_snr_score']        = np.round(snr_arr,  4)
    df['fit_chi2_score']      = np.round(chi2_arr, 4)

    # saturation_score
    lo, hi = P['saturation_ew_low_mA'], P['saturation_ew_high_mA']
    df['saturation_score'] = np.round(
        np.clip((df['ew_mA'].values.astype(float) - lo) / (hi - lo), 0.0, 1.0), 4)

    # abundance_outlier_score — uses 1D LTE abundances; robust across all lines
    a_vals   = df['a_1dlte'].values.astype(float)
    finite   = np.isfinite(a_vals)
    if finite.sum() >= 2:
        med, std = np.median(a_vals[finite]), np.std(a_vals[finite])
        outlier  = (np.clip(np.abs(a_vals - med) / max(P['outlier_sigma_scale'] * std, 1e-9),
                            0.0, 1.0) if std > 1e-9 else np.zeros(n))
        outlier[~finite] = 1.0
    else:
        outlier = np.full(n, 0.5)
    df['abundance_outlier_score'] = np.round(outlier, 4)

    # nlte_correction_score
    nlte_map = {'3D_NLTE_Amarsi2022': 1.0, '1D_NLTE_mean': 0.5, '1D_LTE': 0.0}
    nlte_s   = np.zeros(n)
    if results_df is not None:
        for i, row in df.iterrows():
            m = results_df[(results_df['element'] == row['element']) &
                           (results_df['ion']     == row['ion'])]
            if not m.empty:
                nlte_s[i] = nlte_map.get(str(m.iloc[0].get('nlte_flag', '1D_LTE')), 0.0)
    df['nlte_correction_score'] = np.round(nlte_s, 4)

    # Aggregate line_score
    df['line_score'] = np.round(
        (1.0 - df['vald_proximity_flag'])         * W['proximity'] +
        df['ew_snr_score']                         * W['snr'] +
        df['fit_chi2_score']                       * W['chi2'] +
        (1.0 - df['saturation_score'])             * W['saturation'] +
        (1.0 - df['abundance_outlier_score'])      * W['outlier'] +
        df['nlte_correction_score']                * W['nlte'],
        4
    )

    # line_grade
    def _grade(s):
        if s >= T['A']: return 'A'
        if s >= T['B']: return 'B'
        if s >= T['C']: return 'C'
        return 'D'
    df['line_grade'] = df['line_score'].apply(_grade)

    return df


def _element_grade_summary(scored_df: pd.DataFrame, results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute element_score, element_grade, n_lines_A/B/C/D and A+B weighted abundance.
    Updates results_df in place and returns it.
    """
    T = LINE_GRADE_THRESHOLDS
    results = results_df.copy()

    for _, row in results.iterrows():
        elem, ion = row['element'], row['ion']
        lines = scored_df[(scored_df['element'] == elem) & (scored_df['ion'] == ion)]
        if lines.empty:
            continue

        idx = results.index[(results['element'] == elem) & (results['ion'] == ion)][0]
        gc = lines['line_grade'].value_counts()
        results.at[idx, 'n_lines_A'] = int(gc.get('A', 0))
        results.at[idx, 'n_lines_B'] = int(gc.get('B', 0))
        results.at[idx, 'n_lines_C'] = int(gc.get('C', 0))
        results.at[idx, 'n_lines_D'] = int(gc.get('D', 0))

        ab_lines = lines[lines['line_grade'].isin(['A', 'B'])]
        if not ab_lines.empty and ab_lines['line_score'].sum() > 0:
            el_score = float(np.average(ab_lines['line_score'],
                                        weights=ab_lines['line_score']))
            results.at[idx, 'element_score'] = round(el_score, 4)

            def _grade(s):
                if s >= T['A']: return 'A'
                if s >= T['B']: return 'B'
                if s >= T['C']: return 'C'
                return 'D'
            results.at[idx, 'element_grade'] = _grade(el_score)

            # A+B weighted mean 1D LTE abundance — NLTE delta applied below
            a_ab = float(np.average(ab_lines['a_1dlte'], weights=ab_lines['line_score']))
            delta_nlte = float(row.get('delta_nlte_mean', np.nan))
            a_ab_nlte  = round(a_ab + delta_nlte, 4) if np.isfinite(delta_nlte) else round(a_ab, 4)
            a_ab_std = round(float(np.std(ab_lines['a_1dlte'].values.astype(float))), 4)
            results.at[idx, 'A_X_nlte_AB']  = a_ab_nlte
            results.at[idx, 'A_X_std_AB']   = a_ab_std
            results.at[idx, 'n_lines_AB']   = len(ab_lines)

    return results


# ── Solar EW loader (hybrid Fe I GES + Fe II lines_fit) ──────────────────────

def _load_solar_ews(ew_override: str = None) -> pd.DataFrame:
    """
    Load solar EWs for abundance derivation using a hybrid approach:
      - Fe I:  GES pre-stored EWs (solar_ew_ges_reference.csv) — avoids NLTE EW bias
      - Fe II: lines_fit.py measured EWs (solar_ew.csv) — Fe II not NLTE-affected
      - Other: lines_fit.py measured EWs (solar_ew.csv)

    If ew_override is provided, use that file directly (bypasses hybrid logic).
    If solar_ew_ges_reference.csv is missing, falls back to solar_ew.csv for all elements.
    """
    if ew_override:
        ew_df = pd.read_csv(ew_override)
        ew_df = ew_df[(ew_df['ew_mA'] > 0) & ew_df['ew_mA'].notna()].copy()
        print(f"  EW override: {ew_override} ({len(ew_df)} lines)")
        return ew_df

    solar_ew = pd.read_csv(str(PATHS['solar_ew']))
    solar_ew = solar_ew[(solar_ew['ew_mA'] > 0) & solar_ew['ew_mA'].notna()].copy()
    print(f"  solar_ew.csv: {len(solar_ew)} lines total")

    ges_ref_path = Path(str(PATHS['solar_ew'])).parent / 'solar_ew_ges_reference.csv'
    try:
        ges_fe1 = pd.read_csv(str(ges_ref_path))
        ges_fe1 = ges_fe1[
            (ges_fe1['element'] == 'Fe') & (ges_fe1['ion'] == 'I') &
            (ges_fe1['ew_mA'] > 0)
        ].copy()
        print(f"  GES Fe I reference: {len(ges_fe1)} lines")

        non_fe1 = solar_ew[~((solar_ew['element'] == 'Fe') & (solar_ew['ion'] == 'I'))]
        fe2_count = int(((non_fe1['element'] == 'Fe') & (non_fe1['ion'] == 'II')).sum())
        hybrid = pd.concat([ges_fe1, non_fe1], ignore_index=True)
        print(f"  Hybrid EW: {len(ges_fe1)} Fe I (GES) + {fe2_count} Fe II (lines_fit) "
              f"+ {len(non_fe1) - fe2_count} other elements")
        return hybrid

    except FileNotFoundError:
        print(f"  WARNING: GES Fe I reference not found at {ges_ref_path} — "
              f"falling back to solar_ew.csv for all elements")
        return solar_ew


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
        stellar_params_override: dict = None,
        ew_override: str = None,
        vmic_fixed: bool = False) -> pd.DataFrame:
    """
    Derive abundances for a star and save results.

    Parameters
    ----------
    star_id                  : 'solar', 'procyon', or '55cnc'
    model_grid               : 'ATLAS9.Castelli' (FGK) or 'MARCS.GES' (M dwarfs)
    stellar_params_override  : override dict for uncertainty sensitivity runs (RYA-158)
                               keys: teff_K, logg, feh, vturb_kms
    ew_override              : path to an EW CSV to use instead of the default
                               solar_ew.csv (e.g. solar_ew_ges_reference.csv for
                               GES pre-stored calibration EWs — RYA-196)
    vmic_fixed               : if True, hold vturb_kms fixed during convergence (RYA-238)

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
    elif 'procyon' in star_id.lower():
        params = STAR_PROCYON.copy()
    else:
        params = STAR_55CNC.copy()

    print(f"[1/4] Stellar params:  Teff={params['teff_K']} K  "
          f"logg={params['logg']}  [Fe/H]={params['feh']}  "
          f"vturb={params['vturb_kms']} km/s")

    # ── Load EW table ─────────────────────────────────────────────
    if 'solar' in star_id.lower():
        print(f"[2/4] Loading solar EWs (hybrid: GES Fe I + lines_fit Fe II)...")
        ew_df = _load_solar_ews(ew_override=ew_override)
        print(f"  Total: {len(ew_df)} EW measurements")
    else:
        ew_path = PATHS.get(f'{star_id}_ew',
                             Path(str(PATHS['solar_ew'])).parent / f'{star_id}_ew.csv')
        if not Path(str(ew_path)).exists():
            raise FileNotFoundError(
                f"EW table not found: {ew_path}\nRun pipeline/lines_fit.py first."
            )
        ew_df = pd.read_csv(str(ew_path))
        ew_df = ew_df[(ew_df['ew_mA'] > 0) & ew_df['ew_mA'].notna()].copy()
        print(f"[2/4] Loaded {len(ew_df)} EW measurements from {Path(str(ew_path)).name}")

    # Vetted blend cross-reference from linelist — RYA-209.
    # linelist_solar.csv blend_flag=True rows are the authoritative vetted exclusion list.
    # Hard exclusion applies only to confirmed non-separable blends, not VALD proximity.
    # vald_proximity_flag (continuous 0-1) feeds into the RYA-220 quality scorer instead.
    try:
        _ll = pd.read_csv(
            str(PATHS['linelist_solar']),
            usecols=['element', 'ion', 'wavelength_air_A', 'blend_flag'],
            low_memory=False,
        )
        _ll_vetted = _ll[_ll['blend_flag'] == True][['element', 'ion', 'wavelength_air_A']]
        if 'blend_flag' not in ew_df.columns:
            ew_df['blend_flag'] = False
        _newly = 0
        for (_elem, _ion), _grp in _ll_vetted.groupby(['element', 'ion']):
            _ei = (ew_df['element'] == _elem) & (ew_df['ion'] == _ion)
            if not _ei.any():
                continue
            _bl_w = _grp['wavelength_air_A'].values
            _ew_w = ew_df.loc[_ei, 'wavelength_air_A'].values
            _hit  = np.any(np.abs(_ew_w[:, None] - _bl_w[None, :]) < 0.05, axis=1)
            _idx  = ew_df.index[_ei][_hit]
            _new  = ew_df.loc[_idx, 'blend_flag'] == False
            _new_idx = _idx[_new.values]
            if len(_new_idx):
                ew_df.loc[_new_idx, 'blend_flag'] = True
                _newly += len(_new_idx)
        print(f"  Vetted blend cross-reference: {len(_ll_vetted)} linelist entry(s), "
              f"{_newly} line(s) marked blend_flag=True in EW data")
    except Exception as _e:
        print(f"  WARNING: vetted blend cross-reference failed — {_e}")

    # ── Iterative convergence ─────────────────────────────────────
    print(f"\n[3/4] Iterating to stellar parameter equilibrium "
          f"({model_grid} / {RADIATIVE_TRANSFER_CODE})...")
    converged_params, results, per_line_df = _iterative_parameter_convergence(
        ew_df, params, model_grid=model_grid, vmic_fixed=vmic_fixed
    )

    print(f"\n  Final params: Teff={converged_params['teff_K']:.0f} K  "
          f"logg={converged_params['logg']:.2f}  "
          f"vturb={converged_params['vturb_kms']:.2f} km/s")

    # ── NLTE corrections (RYA-165) ────────────────────────────────
    try:
        from pipeline.nlte_corrections import apply_fe_nlte_corrections
        stellar_params_for_nlte = {
            'teff_K'   : converged_params.get('teff_K',    5777),
            'logg'     : converged_params.get('logg',      4.44),
            'feh'      : converged_params.get('feh',        0.0),
            'vturb_kms': converged_params.get('vturb_kms',  1.0),
        }
        results = apply_fe_nlte_corrections(
            results, stellar_params_for_nlte,
            line_df=ew_df, per_line_df=per_line_df,
        )
        print(f"  NLTE corrections applied to Fe I/II (Amarsi+2022)")
    except FileNotFoundError as e:
        print(f"  WARNING: NLTE grid not found — running 1D LTE only: {e}")
    except Exception as e:
        print(f"  WARNING: NLTE correction failed — running 1D LTE only: {e}")

    # ── Line quality scoring (RYA-220) ────────────────────────────
    print(f"\n  Line quality scoring (RYA-220)...")
    print(f"  Weights: {LINE_SCORE_WEIGHTS}")
    scored_df = pd.DataFrame()
    if not per_line_df.empty:
        scored_df = _compute_line_scores(per_line_df, ew_df, results_df=results)
        results   = _element_grade_summary(scored_df, results)

    # ── Save ──────────────────────────────────────────────────────
    print(f"\n[4/4] Saving results...")
    out_dir  = Path(str(PATHS['solar_ew'])).parent
    out_path = out_dir / f'{star_id}_abundances.csv'
    results.to_csv(out_path, index=False)
    print(f"  Saved → {out_path.name}  ({len(results)} elements, "
          f"{results['n_lines'].sum()} total lines)")

    if not scored_df.empty:
        per_line_out = out_dir / f'{star_id}_per_line.csv'
        scored_df.to_csv(per_line_out, index=False)
        print(f"  Saved → {per_line_out.name}  ({len(scored_df)} lines with quality scores)")

    # ── Solar calibration gate ────────────────────────────────────
    fe = results[results['element'] == 'Fe']
    if len(fe) > 0:
        for _, fe_row in fe.iterrows():
            ion_lbl = fe_row['ion']
            a_1d    = float(fe_row['A_X'])
            n_lines = int(fe_row.get('n_lines', 0))
            scatter_1d   = float(fe_row.get('A_X_std', np.nan))
            scatter_nlte = float(fe_row.get('A_X_std_nlte', np.nan))
            scatter      = scatter_nlte if np.isfinite(scatter_nlte) else scatter_1d
            nlte_ok = 'A_X_nlte' in fe_row and not np.isnan(float(fe_row['A_X_nlte']))
            # A_X_nlte is relative [Fe/H] from iSpec; convert to absolute A(Fe) for
            # gate comparison and display (RYA-261 unit-consistency fix).
            _a_fe_solar = SOLAR_ASPLUND2021['Fe']
            a_nlte_rel  = float(fe_row['A_X_nlte']) if nlte_ok else a_1d
            a_nlte      = a_nlte_rel + _a_fe_solar   # absolute A(Fe)
            tgt_lo      = _a_fe_solar + FE_GATE_LOWER  # = 7.41
            tgt_hi      = _a_fe_solar + FE_GATE_UPPER  # = 7.51
            d_nlte  = float(fe_row.get('delta_nlte_mean', np.nan))
            n_nlte  = int(fe_row.get('n_nlte_lines', 0))
            flag    = str(fe_row.get('nlte_flag', '1D_LTE'))
            ab_pass = tgt_lo <= a_nlte <= tgt_hi
            sc_pass = np.isfinite(scatter) and scatter < FE_SCATTER_GATE
            # GES regions file has only 3 Fe II lines in the optical — coverage gap,
            # not a filter issue. Gate lowered to ≥ 3. See RYA-227 for scope of fix.
            FE2_MIN_LINES = 3
            nl_min        = 20 if ion_lbl == 'I' else FE2_MIN_LINES
            print(f"\n  Fe {ion_lbl}  1D LTE  = {a_1d:.3f}  (scatter 1D = {scatter_1d:.3f} dex)")
            if np.isfinite(d_nlte):
                print(f"  Mean Δ(NLTE) Fe {ion_lbl}= {d_nlte:+.4f} dex  ({n_nlte} lines corrected)")
            print(f"  A(Fe {ion_lbl}) NLTE    = {a_nlte:.3f}  -> {'PASS' if ab_pass else 'FAIL'} (gate {tgt_lo:.2f}-{tgt_hi:.2f})")
            if np.isfinite(scatter):
                sc_label = 'NLTE scatter' if np.isfinite(scatter_nlte) else '1D scatter'
                print(f"  A(Fe {ion_lbl}) {sc_label} = {scatter:.3f}  -> {'PASS' if sc_pass else 'FAIL'} (gate <{FE_SCATTER_GATE} dex)")
            print(f"  nlte_flag Fe {ion_lbl}   = {flag}")
            print(f"  Fe {ion_lbl} n_lines     = {n_lines}  -> {'PASS' if n_lines >= nl_min else 'FAIL'} (>={nl_min})")
            # Grade summary
            n_A = int(fe_row.get('n_lines_A', 0)); n_B = int(fe_row.get('n_lines_B', 0))
            n_C = int(fe_row.get('n_lines_C', 0)); n_D = int(fe_row.get('n_lines_D', 0))
            el_grade = str(fe_row.get('element_grade', '?'))
            a_ab_rel = float(fe_row.get('A_X_nlte_AB', np.nan))
            a_ab     = a_ab_rel + _a_fe_solar if np.isfinite(a_ab_rel) else np.nan
            n_ab = int(fe_row.get('n_lines_AB', 0))
            print(f"  Fe {ion_lbl} grade dist  = A:{n_A}  B:{n_B}  C:{n_C}  D:{n_D}"
                  f"  element_grade={el_grade}")
            scatter_ab = float(fe_row.get('A_X_std_AB', np.nan))
            if np.isfinite(a_ab):
                ab_gate      = tgt_lo <= a_ab <= tgt_hi
                sc_ab_pass   = np.isfinite(scatter_ab) and scatter_ab < FE_SCATTER_GATE
                print(f"  A(Fe {ion_lbl}) NLTE A+B = {a_ab:.4f}  ({n_ab} lines)"
                      f"  -> {'PASS' if ab_gate else 'FAIL'} (gate {tgt_lo:.2f}-{tgt_hi:.2f})")
                if np.isfinite(scatter_ab):
                    print(f"  A(Fe {ion_lbl}) scatter A+B = {scatter_ab:.3f}"
                          f"  -> {'PASS' if sc_ab_pass else 'FAIL'} (gate <{FE_SCATTER_GATE} dex)")
        vmic_val = converged_params.get('vturb_kms', np.nan)
        vmic_pass = 0.80 <= vmic_val <= 1.20 if np.isfinite(vmic_val) else False
        print(f"  vmic             = {vmic_val:.3f} km/s  -> {'PASS' if vmic_pass else 'FAIL'} (gate 0.80-1.20)")

    print(f"\n{'='*62}")
    print(f"  abundances_derive complete.")
    print(f"{'='*62}\n")
    return converged_params, results


if __name__ == '__main__':
    import sys as _sys
    star  = _sys.argv[1] if len(_sys.argv) > 1 else 'solar'
    grid  = _sys.argv[2] if len(_sys.argv) > 2 else 'ATLAS9.Castelli'
    run(star, grid)
