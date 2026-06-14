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
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d

from config.constants import (
    PATHS, PIPELINE, SOLAR_ASPLUND2021, STAR_SOLAR, STAR_55CNC, STAR_PROCYON,
    STAR_LINELISTS, STAR_PARAMS, get_star_params, HARPS_R,
    ISPEC_DIR, RADIATIVE_TRANSFER_CODE, EW_BASELINE_CODE,
    LINE_SCORE_WEIGHTS, LINE_GRADE_THRESHOLDS, LINE_SCORE_PARAMS,
    FE_GATE_LOWER, FE_GATE_UPPER, FE_SCATTER_GATE, FE_IONISATION_GATE,
)


def _star_linelist(star_id: str):
    """
    Resolve the per-star linelist path from STAR_LINELISTS (RYA-270).
    Substring match on star_id, same convention as stellar-param routing.
    Falls back to the solar list WITH AN EXPLICIT WARNING — a star silently
    analysed against the solar pool is exactly the failure mode RYA-270 fixed.
    """
    for key, path in STAR_LINELISTS.items():
        if key in star_id.lower():
            return path, key
    print(f"  WARNING: no linelist mapping for star '{star_id}' — "
          f"falling back to linelist_solar (add an entry to "
          f"config.constants.STAR_LINELISTS)")
    return PATHS['linelist_solar'], 'solar (fallback)'

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

# Synthesis-EW resource cache (RYA-285) — loaded once per process
_synth_cache: dict = {}

_SYNTH_LINELIST_FILE = str(
    ISPEC_DIR / 'input' / 'linelists' / 'transitions' /
    'GESv6_atom_hfs_iso.420_920nm' / 'atomic_lines.tsv'
)
_SYNTH_ISOTOPE_FILE = str(ISPEC_DIR / 'input' / 'isotopes' / 'SPECTRUM.lst')
_SYNTH_CHEM_FILE    = str(ISPEC_DIR / 'input' / 'abundances' / 'chemical_elements_symbols.dat')


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

def _ew_engine_available(code: str) -> bool:
    """
    True iff iSpec can actually run the requested EW→abundance engine.

    No silent fallback (RYA-289 / RYA-167 Bug 1): the caller RAISES on False
    rather than quietly downgrading to SPECTRUM. The whole point of RYA-289 is
    that the EW baseline must be the engine that was asked for, not whatever the
    synthesis RT code happens to leave enabled.
    """
    code = (code or '').lower()
    checks = {
        'moog'     : ispec.is_moog_support_enabled,
        'moog-scat': ispec.is_moog_scat_support_enabled,
        'spectrum' : ispec.is_spectrum_support_enabled,
        # 'width' is iSpec-internal (no external binary); availability == SPECTRUM build.
        'width'    : ispec.is_spectrum_support_enabled,
    }
    check = checks.get(code)
    return bool(check()) if check is not None else False


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
    # EW baseline engine is chosen INDEPENDENTLY of the synthesis RT code
    # (`code`, used above only for Turbospectrum atmosphere interpolation).
    # RYA-289 fix: the old `ew_code = 'spectrum' if code=='turbospectrum' else code`
    # silently forced SPECTRUM whenever synthesis was on — so RYA-285/287 compared
    # synthesis-vs-SPECTRUM, never synthesis-vs-MOOG. Now the baseline is explicit
    # (EW_BASELINE_CODE) and we REFUSE to run if it is unavailable rather than
    # downgrading silently (RYA-167 Bug 1: no silent engine fallback).
    ew_code = EW_BASELINE_CODE
    if not _ew_engine_available(ew_code):
        raise RuntimeError(
            f"EW baseline engine '{ew_code}' is not available to iSpec. Refusing to "
            f"downgrade to SPECTRUM silently — that downgrade is exactly the bug RYA-289 "
            f"fixes (RYA-285 ran SPECTRUM, never MOOG). Install/compile MOOGSILENT under "
            f"$ISPEC_DIR/synthesizer/moog/ or set EW_BASELINE_CODE explicitly in "
            f"config/constants.py."
        )
    print(f"  EW baseline engine: {ew_code}  (RYA-289 — decoupled from "
          f"RADIATIVE_TRANSFER_CODE='{code}')")
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


# ── Synthesis-EW bisection (RYA-285) ─────────────────────────────────────────
# Turbospectrum synthesis-EW: synthesize a narrow window at a trial A(X),
# integrate (1-flux) to get a synthetic EW, bisect until |EW_synth - EW_obs|
# < 1% of EW_obs or ΔA < 0.01 dex. Linelist, isotopes, and chemical elements
# are loaded once and cached. No silent fallback — failed lines are logged and
# marked as such in the per-line output (RYA-167 Bug 1 prevention).

def _load_synth_resources() -> tuple:
    """Load GES linelist, isotopes, chemical elements (cached per process)."""
    global _synth_cache
    if not _synth_cache:
        print("  [synth] Loading GES linelist (first call — cached)...")
        _synth_cache['linelist']      = ispec.read_atomic_linelist(_SYNTH_LINELIST_FILE)
        _synth_cache['isotopes']      = ispec.read_isotope_data(_SYNTH_ISOTOPE_FILE)
        _synth_cache['chem_elements'] = ispec.read_chemical_elements(_SYNTH_CHEM_FILE)
        print(f"  [synth] GES linelist: {len(_synth_cache['linelist'])} lines loaded")
    return (
        _synth_cache['linelist'],
        _synth_cache['isotopes'],
        _synth_cache['chem_elements'],
    )


# ── Flux-space fitting broadening (RYA-287; per-star sourcing RYA-288) ────────
# v1 EW matching set macro=vsini=0, R=500000 because EW is broadening-invariant.
# v2 fits the synthetic FLUX directly to observed normalized flux, and flux is
# shape-sensitive — the synthesis MUST be broadened to the observed line profile
# or χ² is dominated by a width mismatch. Broadening = instrumental (HARPS) ⊕
# macroturbulence ⊕ rotation. These are FIXED inputs (not free params — the only
# free parameter is the target abundance, per the RYA-287 spec).
#
# RYA-288: broadening is sourced PER-STAR from STAR_PARAMS (vmac/vsini) via the
# single get_star_params() resolver — no parallel dict, no _DEFAULT fallback.
# Because flux-space χ² is shape-sensitive, wrong broadening biases the fitted
# A(X) itself, so a star with no vmac/vsini must FAIL LOUD rather than inherit the
# Sun's profile (the bug this fixes). Instrumental R = HARPS_R (instrument constant).


def _resolve_broadening(star_id: str) -> tuple:
    """(R, vmac, vsini, fit_vmac) for `star_id`, sourced from STAR_PARAMS via
    get_star_params.

    FAIL LOUD — never default to solar (RYA-288 / RYA-287 review). An unknown star
    raises in get_star_params; a known star missing vmac/vsini raises here with a
    clear message. The data lives only in STAR_PARAMS (no parallel broadening dict).

    vmac may be a number (FIXED input, e.g. solar 3.8) or the string 'fit'
    (RYA-309 §3.3: fit vmac in synth-v2's RT model with vsini held fixed). When
    'fit', the returned vmac is the RT initial guess ('vmac_init') and fit_vmac
    is True — the caller must fit vmac, never silently use the guess as a fixed
    value (no-silent-broadening, RYA-288).
    """
    rec = get_star_params(star_id)   # raises KeyError if the star has no record
    if 'vmac' not in rec or 'vsini' not in rec:
        raise KeyError(
            f"No broadening (vmac/vsini) in STAR_PARAMS for star '{star_id}'. "
            f"Add them (vsini + vmac, or vmac='fit') before running synthesis-v2. "
            f"Refusing to default to solar broadening — flux-space χ² is "
            f"shape-sensitive, so wrong broadening biases A(X) itself "
            f"(RYA-288 / RYA-287 review)."
        )
    vmac_spec = rec['vmac']
    if isinstance(vmac_spec, str) and vmac_spec.strip().lower() == 'fit':
        return (float(HARPS_R), float(rec.get('vmac_init', 5.0)),
                float(rec['vsini']), True)
    return float(HARPS_R), float(vmac_spec), float(rec['vsini']), False


def _build_fixed_ab(element: str, atom_code: int, trial_A: float) -> np.recarray:
    """Single free-element fixed-abundance recarray on the iSpec SPECTRUM scale."""
    fixed_ab               = np.recarray(1, dtype=[('code', int), ('Abund', float),
                                                    ('element', '|U30')])
    fixed_ab['code'][0]    = atom_code
    fixed_ab['Abund'][0]   = trial_A - 12.036   # A_X → log(N/Ntot) SPECTRUM scale
    fixed_ab['element'][0] = element
    return fixed_ab


def _synth_flux_at_abund(waveobs_nm: np.ndarray,
                         atmosphere: np.ndarray,
                         teff: float, logg: float, feh: float, vturb: float,
                         linelist: np.ndarray,
                         isotopes: np.ndarray,
                         solar_abund: np.ndarray,
                         element: str, atom_code: int,
                         trial_A: float,
                         R: float = 500000, macroturbulence: float = 0.0,
                         vsini: float = 0.0,
                         tmp_dir: str = '/tmp/ispec_codex_synth') -> np.ndarray:
    """
    Normalized Turbospectrum synthetic flux over `waveobs_nm` with ONLY the target
    element's abundance set to `trial_A` (all others fixed at converged composition).

    This is the shared spectrum generator (RYA-287 Step 1): v1 EW matching and v2
    flux fitting both call it. Returns the flux array as-is — the v1 EW integration
    is what discards the blend/wing info (RYA-285 FLAG 1/FLAG 2), so v2 keeps the flux.

    R / macroturbulence / vsini default to the EW-mode values (R=500000, 0, 0) so
    v1 behaviour is byte-for-byte unchanged; v2 passes the observed broadening.
    """
    fixed_ab = _build_fixed_ab(element, atom_code, trial_A)
    return ispec.generate_spectrum(
        waveobs_nm, atmosphere,
        teff, logg, feh, 0.0,
        linelist, isotopes, solar_abund, fixed_ab,
        microturbulence_vel=vturb,
        macroturbulence=macroturbulence, vsini=vsini, R=R,
        verbose=0, code='turbospectrum',
        tmp_dir=tmp_dir,
    )


def _synth_ew_at_abund(waveobs_nm: np.ndarray,
                        atmosphere: np.ndarray,
                        teff: float, logg: float, feh: float, vturb: float,
                        linelist: np.ndarray,
                        isotopes: np.ndarray,
                        solar_abund: np.ndarray,
                        element: str, atom_code: int,
                        trial_A: float,
                        tmp_dir: str = '/tmp/ispec_codex_synth') -> float:
    """
    Synthesize a narrow window at trial_A and return integrated EW in mÅ (v1).

    EW = integral(1-flux) dλ_nm × 10000 mÅ/nm. Broadening zeroed (EW-invariant),
    R=500000 — these are now the defaults of the shared _synth_flux_at_abund.
    """
    flux = _synth_flux_at_abund(
        waveobs_nm, atmosphere, teff, logg, feh, vturb,
        linelist, isotopes, solar_abund, element, atom_code, trial_A,
        tmp_dir=tmp_dir,
    )
    ew_nm = float(np.trapz(np.maximum(0.0, 1.0 - flux), waveobs_nm))
    return ew_nm * 10000.0   # nm → mÅ


def _bisect_synth_abundance(waveobs_nm: np.ndarray,
                              ew_obs_mA: float,
                              atmosphere: np.ndarray,
                              teff: float, logg: float, feh: float, vturb: float,
                              linelist: np.ndarray,
                              isotopes: np.ndarray,
                              solar_abund: np.ndarray,
                              element: str, atom_code: int,
                              a_lo: float = 4.0, a_hi: float = 13.0,
                              tol: float = 0.01, max_iter: int = 20,
                              tmp_dir: str = '/tmp/ispec_codex_synth') -> tuple:
    """
    Bisect A(X) until synth EW matches EW_obs within tol (relative).

    Returns (A_X, converged: bool, n_iter: int).
    Returns (nan, False, 0) if the bracket fails or synthesis raises.
    """
    _kw = dict(atmosphere=atmosphere, teff=teff, logg=logg, feh=feh, vturb=vturb,
               linelist=linelist, isotopes=isotopes, solar_abund=solar_abund,
               element=element, atom_code=atom_code, tmp_dir=tmp_dir)
    try:
        ew_lo = _synth_ew_at_abund(waveobs_nm, trial_A=a_lo, **_kw)
        ew_hi = _synth_ew_at_abund(waveobs_nm, trial_A=a_hi, **_kw)
    except Exception as exc:
        print(f"      FAILED bracket eval: {exc}")
        return np.nan, False, 0

    # Differential bisection: the synthesis window integrates ALL species in the
    # window, not just the target element.  EW_obs is the Gaussian-fitted EW of
    # only the target line.  Use ew_floor = ew_lo (blend floor at minimum
    # target-element abundance) to align scales: find A_X where
    # ew_synthesis(A_X) = ew_floor + ew_obs.
    ew_floor  = ew_lo
    ew_target = ew_floor + ew_obs_mA

    if ew_hi <= ew_target:
        return np.nan, False, 0
    # Guard: if blend response (ew_hi - ew_lo) is smaller than 10% of EW_obs,
    # the target line is too weak to measure in this window.
    if (ew_hi - ew_lo) < ew_obs_mA * 0.10:
        return np.nan, False, 0

    a_mid = a_lo
    for n_iter in range(1, max_iter + 1):
        a_mid = (a_lo + a_hi) / 2.0
        try:
            ew_mid = _synth_ew_at_abund(waveobs_nm, trial_A=a_mid, **_kw)
        except Exception as exc:
            print(f"      FAILED iter {n_iter}: {exc}")
            return np.nan, False, n_iter
        if (abs(ew_mid - ew_target) / max(ew_target, 1.0) < tol
                or (a_hi - a_lo) < tol):
            return a_mid, True, n_iter
        if ew_mid < ew_target:
            a_lo = a_mid
        else:
            a_hi = a_mid

    return a_mid, False, max_iter


def _run_synthesis_mode(last_linemasks: np.ndarray,
                         converged_params: dict,
                         atmosphere: np.ndarray,
                         star_id: str,
                         out_dir: 'Path') -> tuple:
    """
    Run synthesis-EW bisection at converged_params for every matched line.

    Applies the same Fe I EW ceiling as the gate pool (RYA-279).
    Failed lines are logged explicitly — no silent substitution.

    Saves {star_id}_abundances_synth.csv and {star_id}_per_line_synth.csv.
    Returns (results_synth_df, per_line_synth_df).
    """
    import time
    from pipeline.progress import ProgressReporter

    t0 = time.time()
    print(f"\n  [synth] Starting synthesis-EW bisection — RYA-285 ({star_id})")

    _synth_tmp = '/tmp/ispec_codex_synth'
    Path(_synth_tmp).mkdir(exist_ok=True, parents=True)

    linelist, isotopes, chem_elements = _load_synth_resources()
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)

    teff   = float(converged_params['teff_K'])
    logg   = float(converged_params['logg'])
    feh    = float(converged_params['feh'])
    vturb  = float(converged_params['vturb_kms'])
    ew_ceil = float(PIPELINE['vmic_ew_ceiling_mA'])

    # Pre-compute atom codes and iSpec solar A(X) per unique element
    unique_elements = sorted({str(last_linemasks['note'][i]).split()[0]
                               for i in range(len(last_linemasks))})
    atom_codes:   dict = {}
    solar_A_ispec: dict = {}
    for elem in unique_elements:
        try:
            tmp = ispec.create_free_abundances_structure([elem], chem_elements, solar_abund)
            atom_codes[elem]     = int(tmp['code'][0])
            solar_A_ispec[elem]  = float(tmp['Abund'][0]) + 12.036
        except Exception as exc:
            print(f"  [synth] WARNING: element '{elem}' not in chemical elements table: {exc}")

    per_line_rows = []
    n_ok = n_fail = n_skip = 0
    n_total = len(last_linemasks)

    # RYA-286: stderr-only progress + ETA so the (hours-long) synthesis run can
    # be left unattended over SSH. Does not touch per_line_rows / the CSV.
    prog = ProgressReporter(total=n_total,
                            label=f"{star_id} synthesis abundances")

    for i in range(n_total):
        note      = str(last_linemasks['note'][i])
        parts     = note.split()
        element   = parts[0]
        ion_lbl   = 'I' if parts[1] == '1' else 'II'
        wave_A    = float(last_linemasks['wave_A'][i])
        wave_nm   = wave_A / 10.0
        ew_obs    = float(last_linemasks['ew'][i])

        if note == 'Fe 1' and ew_obs > ew_ceil:
            n_skip += 1
            prog.update(item_label=f"{element} {ion_lbl} {wave_A:.2f} (skip)")
            continue

        if element not in atom_codes:
            print(f"  [synth] {i+1}/{n_total}: {note} {wave_A:.3f} Å — no atom code, skip")
            n_skip += 1
            per_line_rows.append({
                'element': element, 'ion': ion_lbl,
                'wavelength_air_A': wave_A, 'ew_mA': ew_obs,
                'a_synth': np.nan, 'converged': False, 'n_iter': 0,
            })
            prog.update(item_label=f"{element} {ion_lbl} {wave_A:.2f} (skip)")
            continue

        a_code  = atom_codes[element]
        a_solar = solar_A_ispec[element]
        a_lo    = max(a_solar - 3.0, 1.0)
        a_hi    = a_solar + 5.0
        # Use the GES line-region window (wave_base / wave_top in nm) —
        # these are designed to capture the line without excessive blending.
        # Step 0.0002 nm (0.002 Å) gives ≥10 pts per thermal FWHM (~0.034 Å).
        wbase  = float(last_linemasks['wave_base'][i])
        wtop   = float(last_linemasks['wave_top'][i])
        wstep  = 0.0002  # nm (0.002 Å)
        waveobs_nm = np.arange(wbase, wtop + wstep * 0.5, wstep)

        elapsed = time.time() - t0
        print(f"  [synth] {i+1}/{n_total}: {note} {wave_A:.3f} Å  "
              f"EW={ew_obs:.1f} mÅ  [{a_lo:.2f},{a_hi:.2f}]  "
              f"({elapsed:.0f}s elapsed)", end='  ', flush=True)

        a_x, converged, n_iter = _bisect_synth_abundance(
            waveobs_nm, ew_obs, atmosphere,
            teff, logg, feh, vturb,
            linelist, isotopes, solar_abund,
            element, a_code,
            a_lo=a_lo, a_hi=a_hi,
            tmp_dir=_synth_tmp,
        )

        if converged:
            n_ok += 1
            print(f"→ A={a_x:.4f} ({n_iter} iter)")
        else:
            n_fail += 1
            print(f"→ FAILED (n_iter={n_iter}, bracket or synthesis error)")

        per_line_rows.append({
            'element': element, 'ion': ion_lbl,
            'wavelength_air_A': wave_A, 'ew_mA': ew_obs,
            'a_synth': round(float(a_x), 4) if np.isfinite(a_x) else np.nan,
            'converged': converged, 'n_iter': n_iter,
            'saturated': bool(ew_obs > ew_ceil),
        })

        prog.update(item_label=f"{element} {ion_lbl} {wave_A:.2f}")

    prog.finish()

    wall = time.time() - t0
    print(f"\n  [synth] {star_id}: {n_ok} converged, {n_fail} failed, {n_skip} skipped "
          f"— {wall:.1f} s ({wall/60:.1f} min) wall-clock")

    per_line_synth_df = pd.DataFrame(per_line_rows)

    # Aggregate per-element medians from converged, non-saturated lines.
    # The narrow GES synthesis window (≈0.2-0.3 Å) captures the line core but
    # NOT the damping wings, so saturated lines (EW > vmic_ew_ceiling_mA) over-
    # estimate A(X) by 2-4 dex in this method. They are converged but unphysical
    # — exclude from the element aggregate (the per-line CSV keeps them flagged).
    #
    # A_X is reported DIFFERENTIAL vs the iSpec Asplund-2009 reference, matching
    # the live MOOG/SPECTRUM path (whose A_X column is also differential, ≈0 for
    # the Sun). This keeps both engines on the same scale so rya_engine_diff can
    # join on A_X. A_X_abs carries the absolute A(X) for science use.
    ok_df = per_line_synth_df[
        (per_line_synth_df['converged'] == True) &
        (per_line_synth_df['saturated'] == False)
    ].copy()
    n_sat_excl = int(((per_line_synth_df['converged'] == True) &
                      (per_line_synth_df['saturated'] == True)).sum())
    if n_sat_excl:
        print(f"  [synth] Aggregate excludes {n_sat_excl} converged-but-saturated "
              f"line(s) (EW > {ew_ceil:.0f} mÅ — narrow-window method invalid there)")

    rows = []
    if not ok_df.empty:
        for (elem, ion), grp in ok_df.groupby(['element', 'ion']):
            a_vals = grp['a_synth'].dropna().values.astype(float)
            if len(a_vals) == 0:
                continue
            a_med    = float(np.nanmedian(a_vals))
            a_std    = float(np.nanstd(a_vals)) if len(a_vals) > 1 else np.nan
            a_ref_09 = solar_A_ispec.get(elem, np.nan)   # iSpec Asplund-2009 absolute
            a_diff   = a_med - a_ref_09 if np.isfinite(a_ref_09) else np.nan
            a_sol21  = SOLAR_ASPLUND2021.get(elem, np.nan)
            xh = round(a_med - a_sol21, 3) if np.isfinite(a_sol21) else np.nan
            rows.append({'element': elem, 'ion': ion,
                         'A_X'    : round(a_diff, 3) if np.isfinite(a_diff) else np.nan,
                         'A_X_abs': round(a_med, 3),
                         'A_X_std': round(a_std, 3) if np.isfinite(a_std) else np.nan,
                         'n_lines': len(a_vals), 'XH': xh})

    results_synth = pd.DataFrame(rows).sort_values(['element', 'ion']).reset_index(drop=True)

    # [X/Fe] using Fe I [X/H] as reference
    fe1 = results_synth[(results_synth['element'] == 'Fe') & (results_synth['ion'] == 'I')]
    feh_ref = float(fe1.iloc[0]['XH']) if (not fe1.empty
                                            and np.isfinite(fe1.iloc[0]['XH'])) else np.nan
    results_synth['XFe'] = results_synth['XH'].apply(
        lambda xh: round(xh - feh_ref, 3)
        if (np.isfinite(xh) and np.isfinite(feh_ref)) else np.nan
    )

    # RYA-289 scale label: the cross-engine comparison column is A_X_abs (absolute
    # LTE A(X)); A_X stays the legacy Asplund-2009-differential column. Synthesis
    # is pure 1D LTE — no NLTE layer applied here.
    results_synth['scale']        = 'absolute'   # describes A_X_abs
    results_synth['nlte_applied'] = False

    out_ab  = Path(str(out_dir)) / f'{star_id}_abundances_synth.csv'
    out_pl  = Path(str(out_dir)) / f'{star_id}_per_line_synth.csv'
    results_synth.to_csv(out_ab, index=False)
    per_line_synth_df.to_csv(out_pl, index=False)
    print(f"  [synth] Saved → {out_ab.name} ({len(results_synth)} elements)")
    print(f"  [synth] Saved → {out_pl.name} ({len(per_line_synth_df)} lines)")

    return results_synth, per_line_synth_df


# ── Synthesis v2: flux-space fitting (RYA-287) ───────────────────────────────
# Resolves RYA-285 FLAG 1 (EW double-count) and FLAG 2 (wing clipping) with one
# move: stop matching EW numbers, fit the synthetic spectrum directly to the
# observed normalized flux over a wing-wide window. Blends in-window are
# SYNTHESIZED (never floored) → genuine blend-awareness (FLAG 1); the wing-wide
# window includes the damping wings instead of clipping them (FLAG 2).

def _load_observed_spectrum(star_id: str) -> tuple:
    """
    Load the observed normalized spectrum for flux-space fitting (RYA-287).

    Returns (obs_wave_nm, obs_flux) — continuum-normalized (≈1.0 at continuum).
    The normalized CSV has no per-pixel error column, so χ² weighting uses a
    locally-estimated continuum-scatter σ per window (see _fit_synth_flux).
    """
    key = None
    for k in ('solar', 'procyon', '55cnc'):
        if k in star_id.lower():
            key = k
            break
    if key is None:
        raise ValueError(f"No normalized-spectrum mapping for star '{star_id}'")
    norm_path = PATHS.get(f'{key}_normalized')
    if norm_path is None or not Path(str(norm_path)).exists():
        raise FileNotFoundError(
            f"Normalized spectrum not found for {star_id}: {norm_path}. "
            "Run spectra_normalize.py first.")
    df = pd.read_csv(str(norm_path), comment='#')
    wave_nm = df['wavelength_air_A'].to_numpy(dtype=float) / 10.0
    flux    = df['flux_normalized'].to_numpy(dtype=float)
    print(f"  [synth-v2] Observed spectrum: {len(wave_nm):,} px, "
          f"{wave_nm[0]*10:.1f}–{wave_nm[-1]*10:.1f} Å (continuum=1.0)")
    return wave_nm, flux


def _fit_synth_flux(obs_wave_nm: np.ndarray, obs_flux: np.ndarray,
                    atmosphere: np.ndarray,
                    teff: float, logg: float, feh: float, vturb: float,
                    linelist: np.ndarray, isotopes: np.ndarray,
                    solar_abund: np.ndarray, element: str, atom_code: int,
                    wave_base: float, wave_top: float,
                    a_lo: float, a_hi: float,
                    R: float, vmac: float, vsini: float,
                    tmp_dir: str = '/tmp/ispec_codex_synth') -> dict:
    """
    v2 blend-aware abundance: fit the broadened synthetic spectrum directly to
    observed normalized flux over [wave_base, wave_top] (nm) by minimizing reduced
    χ² over the single free abundance. Blends are synthesized in-window, never
    floored (FLAG 1). Window is wing-wide so damping wings are included (FLAG 2).

    No silent fallback: returns status ∈ {ok, edge_pinned, failed}; failed/edge
    lines carry A_X=nan / the edge value and are marked, never substituted.
    """
    mask = (obs_wave_nm >= wave_base) & (obs_wave_nm <= wave_top)
    ow, of = obs_wave_nm[mask], obs_flux[mask]
    if of.size < 5:
        return {'status': 'failed', 'reason': 'too few observed pixels',
                'A_X': np.nan, 'red_chi2': np.nan, 'n_pix': int(of.size)}

    # σ scale (Step 0: the normalized CSV has NO per-pixel error column). Use the
    # spec's constant 0.01-flux placeholder = model-adequacy floor, NOT the tiny
    # photon-noise continuum scatter (~0.001-0.002 at solar SNR), which would make
    # χ²ᵣ ≫ 1 for ANY 1D-LTE synthesis (model imperfection ≫ noise) and obscure
    # the per-line fit quality. CRUCIAL: σ is constant, so the χ² minimum location
    # — i.e. the fitted A(X) — is σ-scale-invariant; only red_chi2's magnitude
    # changes. This keeps red_chi2 interpretable (order unity = good fit) while
    # leaving every abundance identical. Per-pixel weighting awaits an error column.
    sig_scalar = 0.01

    # Synthesis grid: fine, wing-wide window (0.002 Å). Generated once-per-eval.
    wstep = 0.0002  # nm
    sw = np.arange(wave_base, wave_top + wstep * 0.5, wstep)

    _kw = dict(atmosphere=atmosphere, teff=teff, logg=logg, feh=feh, vturb=vturb,
               linelist=linelist, isotopes=isotopes, solar_abund=solar_abund,
               element=element, atom_code=atom_code, R=R, macroturbulence=vmac,
               vsini=vsini, tmp_dir=tmp_dir)

    n_eval = [0]
    fail   = [None]

    def chi2(a_x):
        n_eval[0] += 1
        try:
            sf = _synth_flux_at_abund(sw, trial_A=float(a_x), **_kw)
        except Exception as exc:
            fail[0] = str(exc)
            return 1e30
        interp = interp1d(sw, sf, bounds_error=False, fill_value=1.0)
        sf_i = interp(ow)
        r = (of - sf_i) / sig_scalar
        return float(np.nansum(r * r))

    res = minimize_scalar(chi2, bounds=(a_lo, a_hi), method='bounded',
                          options={'xatol': 1e-3})
    if fail[0] is not None:
        return {'status': 'failed', 'reason': f'synthesis error: {fail[0]}',
                'A_X': np.nan, 'red_chi2': np.nan, 'n_pix': int(of.size)}

    a_best   = float(res.x)
    dof      = max(of.size - 1, 1)
    red_chi2 = chi2(a_best) / dof
    edge_pinned = min(abs(a_best - a_lo), abs(a_best - a_hi)) < 1e-2
    return {'status': 'edge_pinned' if edge_pinned else 'ok',
            'A_X': round(a_best, 3), 'red_chi2': round(float(red_chi2), 3),
            'n_pix': int(of.size), 'n_eval': int(n_eval[0])}


def _wingwide_window_nm(wave_nm: float, ew_mA: float) -> tuple:
    """
    Wing-wide synthesis/fit window (RYA-287 Step 2). Half-width derived from the
    line's expected wing extent — EW is the observable strength proxy: a strong
    damped line has broader Lorentzian wings than a weak line. Not a fixed ±X Å.

      hw_A = 0.12 + 0.0026 · EW_mA   →  20 mÅ ≈ ±0.17 Å,  150 mÅ ≈ ±0.51 Å
    bounded to [0.15, 0.60] Å so even the weakest line keeps enough continuum
    pixels and the strongest accepted line (≤ ceiling) contains its wings.
    """
    hw_A = float(np.clip(0.12 + 0.0026 * ew_mA, 0.15, 0.60))
    hw_nm = hw_A / 10.0
    return wave_nm - hw_nm, wave_nm + hw_nm


def _run_synthesis_v2_mode(last_linemasks: np.ndarray,
                           converged_params: dict,
                           atmosphere: np.ndarray,
                           star_id: str,
                           out_dir: 'Path') -> tuple:
    """
    Flux-space synthesis fit (RYA-287) at converged_params for every matched line.

    Writes {star_id}_abundances_synth_v2.csv and {star_id}_per_line_synth_v2.csv
    (v1 schema + red_chi2 + status). v1 outputs are NOT touched.
    Returns (results_v2_df, per_line_v2_df).
    """
    import time
    t0 = time.time()
    print(f"\n  [synth-v2] Flux-space fitting — RYA-287 ({star_id})")

    _synth_tmp = '/tmp/ispec_codex_synth_v2'
    Path(_synth_tmp).mkdir(exist_ok=True, parents=True)

    obs_wave_nm, obs_flux = _load_observed_spectrum(star_id)
    linelist, isotopes, chem_elements = _load_synth_resources()
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)

    teff   = float(converged_params['teff_K'])
    logg   = float(converged_params['logg'])
    feh    = float(converged_params['feh'])
    vturb  = float(converged_params['vturb_kms'])
    ew_ceil = float(PIPELINE['vmic_ew_ceiling_mA'])  # ceiling from constants (RYA-277 coord)

    # Broadening (fixed inputs — only abundance is free). RYA-288: sourced per-star
    # from STAR_PARAMS via get_star_params(); fail loud if absent (no solar default).
    R_inst, vmac, vsini, fit_vmac = _resolve_broadening(star_id)
    if fit_vmac:
        # RYA-309 §3.3: vmac is to be FIT in the RT model (vsini fixed). The RT
        # vmac fit is implemented in the RYA-309 synthesis-run step, not the
        # parameter prereq — fail loud rather than silently use the init guess as
        # a fixed vmac (no-silent-broadening, RYA-288).
        raise NotImplementedError(
            f"[{star_id}] vmac='fit' (RYA-309 §3.3): vsini fixed = {vsini} km/s, "
            f"RT init = {vmac} km/s (expect ~5-6). The RT vmac fit is implemented "
            f"in the RYA-309 synthesis run, not the param prereq. Refusing to use "
            f"the init as a fixed vmac."
        )
    print(f"[{star_id}] broadening: R={R_inst:.0f}, vmac={vmac} km/s, vsini={vsini} km/s",
          file=sys.stderr)
    print(f"  [synth-v2] Broadening: R={R_inst:.0f}, vmac={vmac} km/s, vsini={vsini} km/s "
          f"(fixed); EW ceiling={ew_ceil:.0f} mÅ")

    unique_elements = sorted({str(last_linemasks['note'][i]).split()[0]
                               for i in range(len(last_linemasks))})
    atom_codes: dict = {}
    solar_A_ispec: dict = {}
    for elem in unique_elements:
        try:
            tmp = ispec.create_free_abundances_structure([elem], chem_elements, solar_abund)
            atom_codes[elem]    = int(tmp['code'][0])
            solar_A_ispec[elem] = float(tmp['Abund'][0]) + 12.036
        except Exception as exc:
            print(f"  [synth-v2] WARNING: element '{elem}' not in chem table: {exc}")

    per_line_rows = []
    n_ok = n_edge = n_fail = n_sat = n_skip = 0
    n_total = len(last_linemasks)

    for i in range(n_total):
        note    = str(last_linemasks['note'][i])
        parts   = note.split()
        element = parts[0]
        ion_lbl = 'I' if parts[1] == '1' else 'II'
        wave_A  = float(last_linemasks['wave_A'][i])
        wave_nm = wave_A / 10.0
        ew_obs  = float(last_linemasks['ew'][i])

        base_row = {'element': element, 'ion': ion_lbl,
                    'wavelength_air_A': wave_A, 'ew_mA': ew_obs}

        # Saturation policy: ceiling pulled from constants (coordinate RYA-277);
        # above it, flag saturated rather than fit (deferred to RYA-277).
        if ew_obs > ew_ceil:
            n_sat += 1
            per_line_rows.append({**base_row, 'a_synth': np.nan,
                                  'red_chi2': np.nan, 'status': 'saturated',
                                  'saturated': True})
            continue

        if element not in atom_codes:
            n_skip += 1
            per_line_rows.append({**base_row, 'a_synth': np.nan,
                                  'red_chi2': np.nan, 'status': 'failed',
                                  'saturated': False})
            print(f"  [synth-v2] {i+1}/{n_total}: {note} {wave_A:.3f} Å — no atom code")
            continue

        a_code  = atom_codes[element]
        a_solar = solar_A_ispec[element]
        a_lo    = max(a_solar - 3.0, 1.0)
        a_hi    = a_solar + 5.0
        wbase, wtop = _wingwide_window_nm(wave_nm, ew_obs)

        elapsed = time.time() - t0
        print(f"  [synth-v2] {i+1}/{n_total}: {note} {wave_A:.3f} Å  "
              f"EW={ew_obs:.1f} mÅ  win=±{(wtop-wave_nm)*10:.2f}Å  "
              f"({elapsed:.0f}s)", end='  ', flush=True)

        r = _fit_synth_flux(obs_wave_nm, obs_flux, atmosphere,
                            teff, logg, feh, vturb,
                            linelist, isotopes, solar_abund, element, a_code,
                            wbase, wtop, a_lo, a_hi, R_inst, vmac, vsini,
                            tmp_dir=_synth_tmp)

        status = r['status']
        if status == 'ok':
            n_ok += 1
            print(f"→ A={r['A_X']:.3f}  χ²ᵣ={r['red_chi2']:.2f} ({r.get('n_eval','?')} ev)")
        elif status == 'edge_pinned':
            n_edge += 1
            print(f"→ EDGE-PINNED A={r['A_X']:.3f}  χ²ᵣ={r['red_chi2']:.2f}")
        else:
            n_fail += 1
            print(f"→ FAILED ({r.get('reason','?')})")

        per_line_rows.append({**base_row,
                              'a_synth': r['A_X'], 'red_chi2': r['red_chi2'],
                              'status': status, 'saturated': False})

    wall = time.time() - t0
    print(f"\n  [synth-v2] {star_id}: {n_ok} ok, {n_edge} edge_pinned, "
          f"{n_fail} failed, {n_sat} saturated, {n_skip} no-code "
          f"— {wall:.1f} s ({wall/60:.1f} min) wall-clock")

    per_line_v2_df = pd.DataFrame(per_line_rows)

    # Aggregate from status=='ok' lines only (edge_pinned/saturated/failed excluded).
    # A_X differential vs iSpec Asplund-2009 ref (same scale convention as v1).
    ok_df = per_line_v2_df[per_line_v2_df['status'] == 'ok'].copy()
    rows = []
    if not ok_df.empty:
        for (elem, ion), grp in ok_df.groupby(['element', 'ion']):
            a_vals = grp['a_synth'].dropna().values.astype(float)
            if len(a_vals) == 0:
                continue
            a_med    = float(np.nanmedian(a_vals))
            a_std    = float(np.nanstd(a_vals)) if len(a_vals) > 1 else np.nan
            a_ref_09 = solar_A_ispec.get(elem, np.nan)
            a_diff   = a_med - a_ref_09 if np.isfinite(a_ref_09) else np.nan
            a_sol21  = SOLAR_ASPLUND2021.get(elem, np.nan)
            xh = round(a_med - a_sol21, 3) if np.isfinite(a_sol21) else np.nan
            med_chi2 = float(np.nanmedian(grp['red_chi2'].values.astype(float)))
            rows.append({'element': elem, 'ion': ion,
                         'A_X'    : round(a_diff, 3) if np.isfinite(a_diff) else np.nan,
                         'A_X_abs': round(a_med, 3),
                         'A_X_std': round(a_std, 3) if np.isfinite(a_std) else np.nan,
                         'n_lines': len(a_vals),
                         'med_red_chi2': round(med_chi2, 3), 'XH': xh})

    results_v2 = pd.DataFrame(rows).sort_values(['element', 'ion']).reset_index(drop=True)
    fe1 = results_v2[(results_v2['element'] == 'Fe') & (results_v2['ion'] == 'I')]
    feh_ref = float(fe1.iloc[0]['XH']) if (not fe1.empty
                                            and np.isfinite(fe1.iloc[0]['XH'])) else np.nan
    results_v2['XFe'] = results_v2['XH'].apply(
        lambda xh: round(xh - feh_ref, 3)
        if (np.isfinite(xh) and np.isfinite(feh_ref)) else np.nan)

    # RYA-289 scale label: comparison column is A_X_abs (absolute LTE A(X)).
    # A_X stays the legacy Asplund-2009-differential column. v2 flux fit is 1D LTE.
    results_v2['scale']        = 'absolute'   # describes A_X_abs
    results_v2['nlte_applied'] = False

    out_ab = Path(str(out_dir)) / f'{star_id}_abundances_synth_v2.csv'
    out_pl = Path(str(out_dir)) / f'{star_id}_per_line_synth_v2.csv'
    results_v2.to_csv(out_ab, index=False)
    per_line_v2_df.to_csv(out_pl, index=False)
    print(f"  [synth-v2] Saved → {out_ab.name} ({len(results_v2)} elements)")
    print(f"  [synth-v2] Saved → {out_pl.name} ({len(per_line_v2_df)} lines)")

    return results_v2, per_line_v2_df


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
                                      vmic_fixed: bool = False,
                                      skip_convergence: bool = False,
                                      solve_params=None) -> tuple:
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

    # RYA-292 pin/solve policy: only params in `solve` are optimized; pinned params
    # are held fixed at their fundamental (GBS) values. Default = solve all three
    # spectroscopic params (backward-compatible with pre-292 behaviour). 'xi' also
    # honours the legacy vmic_fixed flag. A pinned param's equilibrium criterion is
    # treated as satisfied so convergence is judged only on the solved params.
    solve = set(solve_params) if solve_params is not None else {'teff', 'logg', 'xi'}
    solve_teff = 'teff' in solve
    solve_logg = 'logg' in solve
    solve_xi   = ('xi' in solve) and not vmic_fixed

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
        # NOTE (RYA-279): the same ceiling is applied to the final gate-statistic
        # pool (post-convergence). The two cuts share PIPELINE['vmic_ew_ceiling_mA'].
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

        # RYA-292: judge convergence only on the SOLVED params. A pinned param's
        # equilibrium slope is informational (it reflects line/data systematics,
        # not a free param), so it must not block convergence of the solved ones.
        converged = (
            (not solve_teff or np.isnan(ep_sl)  or abs(ep_sl)  < 0.02) and
            (not solve_xi   or np.isnan(rew_sl) or abs(rew_sl) < 0.02) and
            (not solve_logg or np.isnan(ion_bal) or abs(ion_bal) < 0.05)
        )
        if converged:
            print(f"  ✓ Converged at iteration {iteration}")
            break

        # RYA-291: pinned-params mode. Run exactly one pass at the supplied params
        # (abundances above are at those params), log the equilibrium slopes for
        # the record, then break WITHOUT updating Teff/logg/ξ. The results table
        # below is built from this pinned pass, so the EW baseline lands at exactly
        # the params the synthesis-v2 tables used — a pure engine diff, no walk.
        if skip_convergence:
            print(f"  ⏸ skip_convergence (RYA-291) — pinned at "
                  f"Teff={params['teff_K']:.0f}K logg={params['logg']:.2f} "
                  f"[Fe/H]={params['feh']} vturb={params['vturb_kms']:.2f}km/s "
                  f"(no param update; slopes above are diagnostic only)")
            break

        # Gradient-descent parameter adjustment — conservative steps to avoid
        # walking off the atmosphere grid on large spurious slopes. RYA-292: only
        # SOLVED params are stepped; pinned Teff/logg/ξ are held at their
        # fundamental values (no walk — the RYA-290 directive, now data-driven).
        if solve_teff and np.isfinite(ep_sl):
            params['teff_K']    += np.sign(ep_sl)  * min(abs(ep_sl)  / 0.05, 1.0) * 25
        if solve_xi and np.isfinite(rew_sl):
            params['vturb_kms'] += np.sign(rew_sl) * min(abs(rew_sl) / 0.10, 1.0) * 0.05
        if solve_logg and np.isfinite(ion_bal):
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

    # RYA-279: single ceiling-correct Fe I pool — gate statistics AND per-line CSV
    # both draw from the same set of lines. Previously A_X_std was computed over
    # the full Fe I pool (all EW) while per_line_df had the ceiling applied
    # post-hoc, creating a divergence: the gate σ reflected lines the CSV excluded.
    # Fix: apply ceiling once here; use this mask everywhere downstream.
    # Fe II: no ceiling (Fe II lines are fewer, weaker, and not in the damping regime
    # at the EW values seen in practice). Ceiling is Fe I specific.
    _ew_ceiling      = float(PIPELINE['vmic_ew_ceiling_mA'])
    _lm_ew           = np.array([float(last_linemasks['ew'][i])
                                  for i in range(len(last_linemasks))])
    _fe1_ceiling_mask = (notes == 'Fe 1') & (_lm_ew <= _ew_ceiling)
    _n_fe1_all        = int((notes == 'Fe 1').sum())
    _n_fe1_ceiling    = int(_fe1_ceiling_mask.sum())
    if _n_fe1_all > _n_fe1_ceiling:
        print(f"  Gate pool COG cut (RYA-279): {_n_fe1_all} → {_n_fe1_ceiling} Fe I lines "
              f"({_n_fe1_all - _n_fe1_ceiling} EW > {_ew_ceiling:.0f} mÅ removed; "
              f"gate σ and per-line CSV now share the same pool)")

    # Pass 1: compute A(X) and [X/H] for all elements on Asplund 2021 scale.
    # Fe I uses the ceiling-correct pool; all other species use the full pool.
    rows_pass1 = []
    for note in np.unique(notes):
        if note == 'Fe 1':
            mask = _fe1_ceiling_mask
        else:
            mask = notes == note
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
    # Fe I: ceiling-correct pool only (same pool as gate statistics, per RYA-279).
    # Fe II: all lines with finite abundances.
    per_line_rows = []
    for i in range(len(last_linemasks)):
        note = str(last_linemasks['note'][i])
        if note not in ('Fe 1', 'Fe 2'):
            continue
        if note == 'Fe 1' and _lm_ew[i] > _ew_ceiling:
            continue  # excluded from gate pool — also excluded from CSV (RYA-279)
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

    return params, results_df, per_line_df, last_linemasks


# ── Line quality scoring (RYA-220) ───────────────────────────────────────────

def _compute_line_scores(per_line_df: pd.DataFrame,
                         ew_df: pd.DataFrame,
                         results_df: pd.DataFrame = None,
                         linelist_path=None) -> pd.DataFrame:
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

    # RYA-280: load linelist for proximity flag lookup.
    # comment='#' is required — per-star linelists (e.g. linelist_procyon.csv) carry
    # metadata comment blocks before the header row. Without it, pandas reads the
    # first comment line as the header, all four usecols fail to match, and the
    # old bare-except block defaulted vpf_arr to 0.5 silently — scoring was inert.
    _ll_path = linelist_path if linelist_path is not None else PATHS['linelist_solar']
    _ll_required = ['element', 'ion', 'wavelength_air_A', 'vald_proximity_flag']
    _ll_raw = pd.read_csv(str(_ll_path), comment='#', low_memory=False)
    _missing = [c for c in _ll_required if c not in _ll_raw.columns]
    if _missing:
        raise KeyError(
            f"vald_proximity_flag load failed: columns {_missing} not found in "
            f"{_ll_path}. Expected columns: {_ll_required}. "
            f"Actual columns: {list(_ll_raw.columns[:10])}..."
        )
    ll = _ll_raw[_ll_required]
    ll_by_ei = {k: v.reset_index(drop=True) for k, v in ll.groupby(['element', 'ion'])}

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
        vmic_fixed: bool = False,
        engine: str = 'spectrum',
        skip_convergence: bool = False) -> tuple:
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
    engine                   : 'spectrum' (default MOOG/SPECTRUM path) or 'synthesis'
                               (Turbospectrum synthesis-EW bisection — RYA-285).
                               'synthesis' runs the normal convergence first, then
                               re-derives per-line A(X) via synthesis at the converged
                               params, saving to {star_id}_abundances_synth.csv.

    Returns
    -------
    (converged_params, results_df) :
        converged_params — dict with final teff_K, logg, feh, vturb_kms
        results_df       — element, ion, A_X, A_X_std, n_lines, XH, XFe
                           [X/H] normalised to SOLAR_ASPLUND2021
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

    # ── Pin/solve policy (RYA-292) ────────────────────────────────
    # Resolve the per-star fundamental-param record (fail LOUD if absent), pin the
    # fundamental Teff/logg to their GBS values, and pass the solve set to the
    # convergence so pinned params are held fixed. Benchmark calibrators pin
    # Teff+logg (RYA-290 directive as data); program stars without a fundamental
    # logg keep the full spectroscopic solve. Override runs bypass the pin.
    star_policy = get_star_params(star_id)
    solve_params = list(star_policy.get('solve', ['teff', 'logg', 'feh', 'xi']))
    if not stellar_params_override:
        if 'teff' in star_policy.get('pin', []):
            params['teff_K'] = float(star_policy['teff'])
        if 'logg' in star_policy.get('pin', []):
            params['logg'] = float(star_policy['logg'])
    pin_lbl = ','.join(star_policy.get('pin', [])) or '(none)'
    solve_lbl = ','.join(solve_params) or '(none)'
    print(f"  Param policy (RYA-292): PIN [{pin_lbl}]  SOLVE [{solve_lbl}]  "
          f"— {star_policy.get('source', '?')}")

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

    # Per-star linelist routing (RYA-270) — explicit and logged, so a star
    # silently analysed against the solar pool can never go unnoticed again.
    _star_ll_path, _star_ll_key = _star_linelist(star_id)
    try:
        _ll_full = pd.read_csv(str(_star_ll_path), comment='#', low_memory=False)
        _fe_pool = _ll_full[(_ll_full['element'] == 'Fe')]
        if 'priority' in _ll_full.columns:
            _fe_pool = _fe_pool[_fe_pool['priority'] > 0]
        _n_fe1 = int((_fe_pool['ion'] == 'I').sum())
        _n_fe2 = int((_fe_pool['ion'] == 'II').sum())
        print(f"  Linelist loaded for {star_id}: "
              f"{Path(str(_star_ll_path)).relative_to(Path(str(PATHS['data_root'])).parent)} "
              f"({_n_fe1} Fe I, {_n_fe2} Fe II)")
    except Exception as _e:
        raise FileNotFoundError(
            f"Star linelist not loadable: {_star_ll_path} — {_e}")

    # Vetted blend cross-reference from linelist — RYA-209.
    # The star linelist's blend_flag=True rows are the authoritative vetted exclusion list.
    # Hard exclusion applies only to confirmed non-separable blends, not VALD proximity.
    # vald_proximity_flag (continuous 0-1) feeds into the RYA-220 quality scorer instead.
    try:
        _ll = _ll_full[['element', 'ion', 'wavelength_air_A', 'blend_flag']]
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

    # Fe pool-membership filter (RYA-270): for stars with a bespoke linelist,
    # the linelist's priority>0 Fe rows DEFINE the admissible Fe pool — EW
    # measurements of Fe lines outside it (e.g. solar-pool lines contaminated
    # at F-star temperatures) are excluded from convergence and abundance
    # derivation. The solar path is exempt: its EW pool was built FROM
    # linelist_solar by lines_fit, so membership holds by construction and
    # solar behaviour stays bit-for-bit unchanged. Non-Fe rows pass through
    # (Fe selection is the scope of RYA-270; gates are Fe-only).
    if Path(str(_star_ll_path)) != Path(str(PATHS['linelist_solar'])):
        _fe_ew = (ew_df['element'] == 'Fe')
        _pool_w = {ion: _fe_pool.loc[_fe_pool['ion'] == ion,
                                     'wavelength_air_A'].values
                   for ion in ('I', 'II')}
        _member = np.ones(len(ew_df), dtype=bool)
        for _i in ew_df.index[_fe_ew]:
            _ion = ew_df.at[_i, 'ion']
            _w   = float(ew_df.at[_i, 'wavelength_air_A'])
            _pw  = _pool_w.get(_ion, np.array([]))
            _member[ew_df.index.get_loc(_i)] = (
                len(_pw) > 0 and np.abs(_pw - _w).min() < 0.05)
        _n_dropped = int((~_member).sum())
        ew_df = ew_df[_member].copy()
        print(f"  Fe pool-membership filter ({_star_ll_key} linelist): "
              f"{_n_dropped} Fe EW line(s) outside the F-star pool excluded, "
              f"{int((ew_df['element'] == 'Fe').sum())} Fe line(s) retained")

    # ── Pinned-params synthesis (RYA-285 Step 2) ──────────────────
    # skip_convergence: derive at EXACTLY the supplied params (no Teff/logg/ξ
    # iteration), then run synthesis. Used to compare engines at a fixed
    # benchmark anchor where convergence would otherwise walk the params and
    # confound the engine diff with a stellar-parameter difference.
    if skip_convergence and engine in ('synthesis', 'synthesis-v2'):
        print(f"\n[3/4] PINNED params (skip_convergence) — no iteration:")
        print(f"  Teff={params['teff_K']:.0f} K  logg={params['logg']:.2f}  "
              f"[Fe/H]={params['feh']}  vturb={params['vturb_kms']:.2f} km/s")
        _atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'],
                                params['vturb_kms'], model_grid=model_grid)
        last_linemasks, _, _, _ = _ew_to_abundance(ew_df, params, _atm, model_grid)
        out_dir = Path(str(PATHS['solar_ew'])).parent
        if engine == 'synthesis-v2':
            print(f"\n[4/4] Synthesis-v2 flux-space fit — pinned params...")
            results_synth, _ = _run_synthesis_v2_mode(
                last_linemasks, params, _atm, star_id, out_dir)
        else:
            print(f"\n[4/4] Synthesis-EW mode (Turbospectrum bisection) — pinned params...")
            results_synth, _ = _run_synthesis_mode(
                last_linemasks, params, _atm, star_id, out_dir)
        print(f"\n{'='*62}\n  abundances_derive complete (pinned synthesis).\n{'='*62}\n")
        return params, results_synth

    # ── Iterative convergence (or pinned single pass for the EW baseline) ──────
    # RYA-291: skip_convergence on the EW baseline (engine='spectrum') pins the
    # params to the supplied init (= the synthesis-v2 tables' params) so the MOOG
    # baseline and the synthesis are compared at IDENTICAL params — removing the
    # RYA-289 param confound. Same result-building/NLTE/scoring/save path follows.
    if skip_convergence:
        print(f"\n[3/4] PINNED params (skip_convergence) — EW baseline, no iteration:")
    else:
        print(f"\n[3/4] Iterating to stellar parameter equilibrium "
              f"({model_grid} / {RADIATIVE_TRANSFER_CODE})...")
    converged_params, results, per_line_df, last_linemasks = _iterative_parameter_convergence(
        ew_df, params, model_grid=model_grid, vmic_fixed=vmic_fixed,
        skip_convergence=skip_convergence, solve_params=solve_params,
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
        scored_df = _compute_line_scores(per_line_df, ew_df, results_df=results,
                                         linelist_path=_star_ll_path)
        results   = _element_grade_summary(scored_df, results)

    # ── RYA-289: explicit absolute scale labels (no implicit ≡0 differential) ──
    # MOOG returns absolute A(X) uniformly, so A_X here is 1D-LTE absolute A(X)
    # for every star (this is what dissolves the solar A_X≡0 degeneracy that the
    # SPECTRUM path produced). The NLTE layer (Fe only) is stored as a DISTINCT,
    # explicitly-absolute column (A_X_nlte_absolute) instead of being reconstructed
    # implicitly via +A(Fe)_sun at gate time. The scale-aware diff harness joins
    # on the labelled absolute A_X.
    results['scale']        = 'absolute'   # describes A_X (1D LTE absolute)
    results['nlte_applied'] = False        # A_X is pre-NLTE; NLTE in A_X_nlte_absolute
    if 'A_X_nlte' in results.columns:
        def _nlte_to_absolute(row):
            flag = str(row.get('nlte_flag', '1D_LTE'))
            if flag in ('1D_LTE', 'NLTE_unavailable'):
                return np.nan   # NLTE not actually applied → no absolute NLTE value
            ref = SOLAR_ASPLUND2021.get(str(row['element']), np.nan)
            val = row.get('A_X_nlte', np.nan)
            if pd.notna(val) and np.isfinite(ref):
                return round(float(val) + float(ref), 3)
            return np.nan
        results['A_X_nlte_absolute'] = results.apply(_nlte_to_absolute, axis=1)

    # ── Save ──────────────────────────────────────────────────────
    print(f"\n[4/4] Saving results...")
    out_dir  = Path(str(PATHS['solar_ew'])).parent
    out_path = out_dir / f'{star_id}_abundances.csv'
    results.to_csv(out_path, index=False)
    print(f"  Saved → {out_path.name}  ({len(results)} elements, "
          f"{results['n_lines'].sum()} total lines)")

    if not scored_df.empty:
        scored_df['scale'] = 'absolute'   # RYA-289: per-line A_X is 1D-LTE absolute
        per_line_out = out_dir / f'{star_id}_per_line.csv'
        scored_df.to_csv(per_line_out, index=False)
        print(f"  Saved → {per_line_out.name}  ({len(scored_df)} lines with quality scores)")

    # ── Synthesis mode (RYA-285 v1 EW / RYA-287 v2 flux) ──────────
    if engine in ('synthesis', 'synthesis-v2'):
        _atm_synth = _load_atmosphere(
            converged_params['teff_K'], converged_params['logg'],
            converged_params['feh'],    converged_params['vturb_kms'],
            model_grid=model_grid,
        )
        if engine == 'synthesis-v2':
            print(f"\n[4b/4] Synthesis-v2 flux-space fit (Turbospectrum)...")
            _run_synthesis_v2_mode(last_linemasks, converged_params, _atm_synth,
                                   star_id, out_dir)
        else:
            print(f"\n[4b/4] Synthesis-EW mode (Turbospectrum bisection)...")
            _run_synthesis_mode(last_linemasks, converged_params, _atm_synth,
                                star_id, out_dir)

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
            # nlte_flag distinguishes the applied-NLTE case from unavailable fallback.
            # When NLTE is unavailable, nlte_corrections.py stores the raw A_X value
            # (already absolute) in A_X_nlte — adding _a_fe_solar would double-count
            # the solar offset. Only convert when NLTE was actually applied (RYA-267).
            _flag       = str(fe_row.get('nlte_flag', '1D_LTE'))
            nlte_ok     = _flag not in ('NLTE_unavailable', '1D_LTE')
            _a_fe_solar = SOLAR_ASPLUND2021['Fe']
            a_nlte      = float(fe_row['A_X_nlte']) + _a_fe_solar if nlte_ok else a_1d
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
            # A_X_nlte_AB is relative when NLTE applied, absolute when unavailable — same
            # convention as A_X_nlte. Only add solar offset when NLTE was applied (RYA-267).
            a_ab     = (a_ab_rel + _a_fe_solar) if (np.isfinite(a_ab_rel) and nlte_ok) else \
                       (a_ab_rel if np.isfinite(a_ab_rel) else np.nan)
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
    from pathlib import Path as _Path
    _repo = str(_Path(__file__).resolve().parent.parent)
    if _repo not in _sys.path:
        _sys.path.insert(0, _repo)
    _flags = [a for a in _sys.argv[1:] if a.startswith('--')]
    _pos   = [a for a in _sys.argv[1:] if not a.startswith('--')]
    star   = _pos[0] if len(_pos) > 0 else 'solar'
    grid   = _pos[1] if len(_pos) > 1 else 'ATLAS9.Castelli'
    engine = _pos[2] if len(_pos) > 2 else 'spectrum'
    # RYA-291: --skip-convergence pins the params to the star init (no walk), so the
    # MOOG EW baseline is computed at the synthesis-v2 tables' params for a pure diff.
    skip   = '--skip-convergence' in _flags or '--pin' in _flags
    run(star, grid, engine=engine, skip_convergence=skip)
