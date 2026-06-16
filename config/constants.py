"""
config/constants.py
===================
Single source of truth for all constants, paths, and parameters.
NO hardcoded values anywhere else in the pipeline.

Physical constants: NIST CODATA 2022 recommended values
Stellar parameters: CHARA interferometry (von Braun et al. 2011)
Solar reference:    Asplund et al. 2021 (A&A 653, A141)
"""

__version__ = "2026-05-28"

from pathlib import Path

# ── Repository root ───────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent

# ── Physical constants (NIST CODATA 2022) ────────────────────────────────────
PHYSICS = {
    'c_kms'        : 299792.458,        # Speed of light (km/s) — exact
    'c_ms'         : 299792458.0,       # Speed of light (m/s) — exact
    'h_eV'         : 4.135667696e-15,   # Planck constant (eV·s)
    'k_eV'         : 8.617333262e-5,    # Boltzmann constant (eV/K)
    'k_cgs'        : 1.380649e-16,      # Boltzmann constant (erg/K)
    'sigma_sb'     : 5.670374419e-5,    # Stefan-Boltzmann (erg/cm²/s/K⁴)
    'G_cgs'        : 6.67430e-8,        # Gravitational constant (cm³/g/s²)
    'amu_g'        : 1.66053906660e-24, # Atomic mass unit (g)
}

# ── Astronomical constants ────────────────────────────────────────────────────
ASTRO = {
    'Msun_g'       : 1.989e33,          # Solar mass (g)
    'Rsun_cm'      : 6.957e10,          # Solar radius (cm)
    'Lsun_erg'     : 3.828e33,          # Solar luminosity (erg/s)
    'Teff_sun'     : 5778.0,            # Solar effective temperature (K)
    'logg_sun'     : 4.438,             # Solar surface gravity (cgs)
    'pc_cm'        : 3.085677581e18,    # Parsec (cm)
    'AU_cm'        : 1.495978707e13,    # Astronomical unit (cm)
}

# ── Solar abundances: Asplund et al. 2021 (A&A 653, A141) ────────────────────
# A(X) = log(N_X/N_H) + 12
# This is the reference — every [X/H] in the Codex is defined as A(X)_star - A(X)_sun
SOLAR_ASPLUND2021 = {
    'H'  : 12.00,   # by definition
    'He' : 10.914,
    'Li' :  1.05,
    'C'  :  8.46,
    'N'  :  7.83,
    'O'  :  8.69,
    'Na' :  6.24,
    'Mg' :  7.55,
    'Al' :  6.43,
    'Si' :  7.51,
    'P'  :  5.41,
    'S'  :  7.12,
    'K'  :  5.07,
    'Ca' :  6.30,
    'Sc' :  3.14,
    'Ti' :  4.97,
    'V'  :  3.90,
    'Cr' :  5.62,
    'Mn' :  5.42,
    'Fe' :  7.46,   # NOTE: was 7.50 in Lodders 2003 used in 2010 thesis
    'Co' :  4.94,
    'Ni' :  6.20,
    'Cu' :  4.18,   # Asplund et al. 2021, Table 2 — added RYA-201
    'Ba' :  2.27,
    'Sr' :  2.83,   # Asplund et al. 2021, Table 2 — added RYA-201
    'Y'  :  2.21,
    'Zr' :  2.59,
    'Eu' :  0.52,
}

# ── Abundance scale sanity guard (RYA-334) ───────────────────────────────────
# A(X) = log(N_X/N_H) + 12 with A(H) ≡ 12.00 by definition. Every ABSOLUTE
# abundance sits on this scale: hydrogen is the ceiling, everything else below it
# (Fe ≈ 7.46). A value above ~12.5 is physically impossible and is the fingerprint
# of the RYA-267 relative/absolute double-add — re-adding the 7.46 solar offset to
# an already-absolute A_X / A_X_nlte lands ~15 ("above hydrogen"). This tripwire
# makes that class of mistake fail loud at computation/output instead of printing
# ~15 or hiding behind a gate FAIL. Instances were fixed reactively in RYA-320/330;
# the guard stops the next one.
A_X_SCALE_MIN = -2.0    # below the faintest measured absolute A(X)
A_X_SCALE_MAX = 12.5    # just above A(H)=12; a 7.46 double-add lands ~15


def assert_abundance_on_scale(value, context=""):
    """Fail loud if an ABSOLUTE A(X) is off the A(H)=12 scale (RYA-334 tripwire).

    NaN/None pass through (used for 'not applicable' / NLTE-unavailable). Apply to
    absolute abundances only — A_X, A_X_nlte, A_X_nlte_absolute, A_X_nlte_AB — and
    NOT to relative [X/H], which is legitimately negative for metal-poor stars.
    Returns ``value`` so it can be used inline.
    """
    import math
    try:
        v = float(value)
    except (TypeError, ValueError):
        return value
    if math.isnan(v):
        return value
    if not (A_X_SCALE_MIN < v < A_X_SCALE_MAX):
        raise ValueError(
            f"Absolute abundance A(X)={v:.3f} for {context!r} is off the A(H)=12 "
            f"scale (expected {A_X_SCALE_MIN} < A_X < {A_X_SCALE_MAX}). A value above "
            f"~12.5 ('above hydrogen') is the RYA-267 relative/absolute double-add "
            f"signature: an already-absolute A_X/A_X_nlte had the 7.46 solar offset "
            f"re-added. A_X_nlte is stored absolute (RYA-319) — never re-add the solar "
            f"offset to it (cf. RYA-320/330/334)."
        )
    return value


# ── 55 Cancri A stellar parameters ───────────────────────────────────────────
# Primary source: CHARA interferometry, von Braun et al. 2011 (ApJ 729, 63)
# Secondary: Valenti & Fischer 2005, Teske et al. 2013
STAR_55CNC = {
    'name'         : '55 Cancri A',
    'iau_name'     : 'Copernicus',
    'hd'           : 'HD 75732',
    'hip'          : 'HIP 43587',
    'teff_K'       : 5196.0,    # ± 24 K (CHARA)
    'teff_err'     : 24.0,
    'logg'         : 4.41,      # ± 0.02 (CHARA)
    'logg_err'     : 0.02,
    'feh'          : 0.32,      # ± 0.02 (spectroscopic)
    'feh_err'      : 0.02,
    'vturb_kms'    : 0.9,       # km/s initial estimate — refined from Fe lines
    'mass_msun'    : 0.905,     # ± 0.015 (CHARA)
    'radius_rsun'  : 0.943,     # ± 0.010 (CHARA)
    'age_gyr'      : 10.2,      # ± 2.5
    'distance_pc'  : 12.55,
    'vmag'         : 5.95,
    'vsini_kms'    : 2.3,       # Very slow rotator — no rotational broadening issues
    'rv_kms'       : 27.58,     # Systemic radial velocity — for Doppler correction
}

# ── Canonical 27-element target list (RYA-109) ───────────────────────────────
# 26 unique symbols; Fe covers both Fe I (Teff) and Fe II (log g) — distinguished
# by the 'ion' column in all linelists.  Priorities defined in elements_master.json.
TARGET_ELEMENTS = [
    # Priority 1 — science critical
    'Fe', 'C', 'O', 'Mg', 'Si', 'Ca', 'Ti', 'Ni', 'Na', 'P', 'S',
    # Priority 2 — tracers
    'N', 'Co', 'Cr', 'Al', 'K', 'Ba', 'Y', 'V', 'Cu',
    # Priority 3 — supplementary
    'Mn', 'Sc', 'Li', 'Eu', 'Zr', 'Sr',
]

# ── Solar calibration parameters ─────────────────────────────────────────────
# Reflected sunlight from Ceres observed with HARPS as solar reference spectrum
# Citation: Dumusque et al. 2021, arXiv:2009.01945, ESO program 1102.D-0954(A)
STAR_SOLAR = {
    'name'         : 'Sun',
    'program_id'   : '1102.D-0954(A)',
    'citation'     : 'Dumusque et al. 2021, arXiv:2009.01945',
    'teff_K'       : 5777.0,    # IAU 2015 nominal solar effective temperature
    'logg'         : 4.44,      # cgs
    'feh'          : 0.00,      # by definition — solar zero-point
    'vturb_kms'    : 1.0,       # km/s
    'rv_kms'       : 0.0,       # solar zero-point; BERV correction applied per exposure
    'ni6300_ew_lit_mA' : 1.0,   # Ni I 6300.336 EW in solar photosphere (Allende Prieto+2001, ApJ 556 L63)
                                  # VALD3 log_gf=-2.841 predicts ~3.8 mÅ via COG — 0.575 dex too strong.
                                  # COG sanity cap fires; this literature value is used as the fallback.
                                  # log_gf needs cross-check vs Johansson et al. 2003 (separate ticket).
}

# ── Procyon (α CMi A) stellar parameters ─────────────────────────────────────
# RYA-309 §3.4 / RYA-298: STAR_PARAMS['procyon'] is the AUTHORITATIVE single
# source for Teff/logg/[Fe/H]. This legacy-key dict is the adapter the older
# consumers (lines_fit plot titles, abundances_derive convergence seed) still
# read; its Teff/logg/feh are RECONCILED to STAR_PARAMS (Heiter+2015 GBS
# 6554/4.00, Jofré+2014 [Fe/H]) — NOT the former Allende Prieto 2002 6530/3.96.
# The pinned run already pulls Teff/logg from STAR_PARAMS (abundances_derive
# pin step), so these only seed the SOLVED params (feh, ξ) and display.
# vturb_kms is the ξ-solve seed (Allende Prieto 2002 / RYA-284), not a pin.
# Structural removal of this dict is RYA-298.
STAR_PROCYON = {
    'name'      : 'Procyon',
    'hd'        : 'HD 61421',
    'hip'       : 'HIP 37279',
    'teff_K'    : 6554.0,   # Heiter+2015 (mirror of STAR_PARAMS['procyon']['teff'])
    'logg'      : 4.00,     # Heiter+2015 (mirror of STAR_PARAMS['procyon']['logg'])
    'feh'       : 0.03,     # Jofré+2014 (mirror of STAR_PARAMS['procyon']['feh_ref'])
    'vturb_kms' : 1.66,     # km/s — ξ-solve seed (Allende Prieto 2002), not pinned
    'rv_kms'    : -3.0,
}

# ── Fe abundance gate thresholds (RYA-261) ───────────────────────────────────
# Relative [Fe/H] units — matches iSpec output convention (0.0 at solar).
# Absolute equivalent: A(Fe) = [Fe/H] + SOLAR_ASPLUND2021['Fe'] (= [Fe/H] + 7.46)
# FE_GATE_LOWER = -0.05  →  A(Fe) = 7.41
# FE_GATE_UPPER = +0.05  →  A(Fe) = 7.51
FE_GATE_LOWER      = -0.05   # [Fe/H] solar Fe lower bound
FE_GATE_UPPER      = +0.05   # [Fe/H] solar Fe upper bound
FE_SCATTER_GATE    =  0.10   # dex — maximum acceptable σ(Fe I)
FE_IONISATION_GATE =  0.05   # dex — maximum |ΔFe(I−II)|

# ── Ni I 6300.336 COG constants ───────────────────────────────────────────────
# Used in _predict_ni6300_ew() to model the O I 6300 blend contamination.
# Source: VALD3 log_gf; Allende Prieto et al. 2001 (ApJ 556, L63) for EW.
NI6300_COG = {
    'log_gf'               : -2.841,
    'excitation_potential_eV': 4.266,
    'ew_lit_solar_mA'      : 1.1,   # Allende Prieto+2001
}

# ── Pipeline settings ─────────────────────────────────────────────────────────
PIPELINE = {
    # Wavelength range (full HARPS optical range)
    'wav_min_A'         : 3780.0,
    'wav_max_A'         : 6910.0,

    # S/N requirements
    'snr_min_science'   : 200,    # Minimum S/N for abundance work
    'snr_target'        : 400,    # Target S/N after co-adding

    # Gaussian fitting
    'fit_window_A'      : 2.0,    # ± Angstroms around line center for fitting
    'voigt_threshold_mA': 150,    # EW threshold above which to use Voigt profile

    # Continuum normalization
    'continuum_poly_deg'        : 3,     # Polynomial degree for continuum fit
    'continuum_n_iter'          : 5,     # Sigma-clipping iterations for spline fit
    'continuum_knot_spacing_A'  : 100.0, # Window width (Å) for continuum anchor points
    'continuum_upper_percentile': 95,    # Upper percentile of flux used as anchor
    'sigma_clip'                : 3.0,   # Sigma clipping threshold

    # Line quality thresholds
    'ew_min_mA'         : 5.0,    # Minimum detectable EW (mÅ)
    'ew_max_mA'         : 300.0,  # Maximum EW before saturation concerns
    'vmic_ew_ceiling_mA': 100.0,  # COG cut for the Fe I vmic slope sample AND the
                                   # RYA-279 gate pool — lines above this are in the
                                   # damping regime; their REW is insensitive to vmic
                                   # but sensitive to damping params (RYA-226).
                                   # RYA-330: tightened 150 → 100 mÅ. At 150 the solar
                                   # Fe I anchor (GES reference pool) carried 12
                                   # saturated lines (EW→149 mÅ, A as low as 7.20) that
                                   # drove a spurious −0.21 reduced-EW slope (RYA-328).
                                   # The cut applies to ALL stars (the saturation onset
                                   # is ~star-independent), so Sun and program stars now
                                   # clear the same bar. Sweep (solar, NLTE-confirmed):
                                   # 150 → slope −0.209 σ0.169; 120 → −0.086 σ0.154;
                                   # 100 → −0.053 σ0.151 (flat within error), A(Fe I)
                                   # NLTE stable 7.516 (median-robust → anchor barely
                                   # moves, +0.009). Re-routing Fe I to the measured
                                   # pool was tested and REJECTED (A(Fe I) 7.67, σ0.44,
                                   # slope +0.43 — the measured solar EWs are worse).
    'vmic_abund_clip_sigma': 2.0,  # σ threshold for per-iteration abundance clip on vmic
                                   # slope sample. Restores behaviour removed by RYA-220
                                   # commit c8a23d2. Sousa et al. 2011, Jofré et al. 2014.
    'snr_line_min'      : 5.0,    # Minimum line S/N for reliable EW

    # Fe II EW triage (RYA-305) — kills the has_theo==False fail-open. A Fe II line
    # is judged against its PROPER theoretical EW (SPECTRUM at expected A(Fe)=solar+
    # [Fe/H]), not the Fe I-calibrated linear COG.
    'fe2_ew_sanity_dex'    : 0.5,  # |log10(obs_EW / theo_EW)| below this = clean (keep EW).
                                    # Above = blend-contaminated → excluded from the EW pool
                                    # (synthesis recovers it if the line is real).
    'fe2_synth_ew_floor_mA': 5.0,  # measurability floor on the theoretical Fe II-only EW:
                                    # a non-clean line with theo >= floor is a real line
                                    # blended in EW → RECOVER via synthesis; theo < floor
                                    # (e.g. 5376.5 Å, ~0 mÅ real Fe II under blend) → DROP.

    # NIST grade policy
    'min_nist_grade'    : 'B',    # Minimum acceptable grade for science

    # EW measurement (lines_fit.py)
    'cont_edge_frac'        : 0.25,         # Fraction of fit window used for continuum anchors
    'cont_anchor_percentile': 80,            # Keep pixels above this percentile in anchor strips
    'min_fit_depth'         : 0.008,        # Skip features shallower than this before fitting
    'blue_edge_warn_A'      : 3900.0,       # Flag lines below this wavelength as low-SNR

    # Depth-adaptive continuum window (strong lines have broad Voigt wings).
    # Only lines with depth ≥ 0.50 get wider windows; moderate lines unchanged.
    # Window = base × scale_factor where tier is determined by depth breaks.
    # Tiers: [0, break0) → scale[0], [break0, break1) → scale[1], ...
    'cont_window_depth_breaks'  : (0.50, 0.70),         # depth thresholds
    'cont_window_scale_factors' : (1.0, 1.75, 2.5),     # multipliers per tier
}

# ── File paths ────────────────────────────────────────────────────────────────
PATHS = {
    'data_root'         : ROOT / 'data',
    'raw_spectra'       : ROOT / 'data' / 'raw',
    'processed_spectra' : ROOT / 'data' / 'processed',
    'linelists'         : ROOT / 'data' / 'linelists',
    'model_atmospheres' : ROOT / 'data' / 'model_atmospheres',
    'results'           : ROOT / 'results',
    'plots'             : ROOT / 'results' / 'plots',
    'tables'            : ROOT / 'results' / 'tables',

    # Key data files
    'linelist_master'   : ROOT / 'data' / 'linelists' / 'linelist_master.csv',
    'linelist_solar'    : ROOT / 'data' / 'linelists' / 'linelist_solar.csv',
    'linelist_procyon'  : ROOT / 'data' / 'linelists' / 'linelist_procyon.csv',
    'solar_reference'   : ROOT / 'data' / 'linelists' / 'solar_asplund2021.csv',
    'atlas9_grid'       : ROOT / 'data' / 'model_atmospheres' / 'atlas9',

    # Solar calibration spectra (external data — relative to repo parent)
    'solar_spectra'     : ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Solar Calibration' / 'archive',

    # Procyon HARPS spectra (external data — relative to repo parent)
    'procyon_spectra'   : ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Procyon' / 'Procyon Harps',

    # Processed outputs (CSV gitignored — generated by pipeline; PNG committed to results/)
    'solar_normalized'  : ROOT / 'data' / 'processed' / 'solar_normalized.csv',
    'solar_diagnostic'  : ROOT / 'results' / 'plots' / 'solar_continuum_diagnostic.png',
    'solar_ew'          : ROOT / 'data' / 'processed' / 'solar_ew.csv',
    'solar_ew_diagnostic': ROOT / 'results' / 'plots' / 'solar_ew_diagnostic.png',
    'procyon_normalized' : ROOT / 'data' / 'processed' / 'procyon_normalized.csv',
    'procyon_diagnostic' : ROOT / 'results' / 'plots' / 'procyon_continuum_diagnostic.png',
    'procyon_ew'         : ROOT / 'data' / 'processed' / 'procyon_ew.csv',
    'procyon_ew_diagnostic': ROOT / 'results' / 'plots' / 'procyon_ew_diagnostic.png',
}

# ── Star → linelist mapping (RYA-270) ────────────────────────────────────────
# Explicit per-star linelist routing. Stars not listed here fall back to the
# solar list — the fallback is LOGGED at runtime by abundances_derive so a
# star silently running on the solar pool can never go unnoticed again.
# Keys are matched by substring against star_id (same convention as the
# stellar-parameter routing in abundances_derive.run()).
STAR_LINELISTS = {
    'solar'   : PATHS['linelist_solar'],
    'procyon' : PATHS['linelist_procyon'],
}

# ── Per-star EW fitting parameters (RYA-273) ─────────────────────────────────
# Solar values mirror PIPELINE defaults. Procyon fit_window_A is wider:
#   Teff 6530 K (vs 5777 K solar) → ~6% broader thermal line widths; lower
#   logg 3.96 (vs 4.44) deepens Voigt pressure wings. 2.5 Å keeps continuum
#   anchors clear of those wings while remaining selective in the crowded
#   F-star optical spectrum. voigt_threshold_mA unchanged — profile regime
#   decision depends on EW, which scales with abundance, not with Teff.
EW_FIT_PARAMS = {
    'solar': {
        'fit_window_A'      : 2.0,
        'voigt_threshold_mA': 150,
    },
    'procyon': {
        'fit_window_A'      : 2.5,
        'voigt_threshold_mA': 150,
    },
}

# ── Per-star continuum normalization parameters (RYA-274) ─────────────────────
# Solar values mirror the PIPELINE dict (backward-compat). Procyon values differ:
#   knot_spacing_A=50: finer knots to track curvature near Balmer wings (Hβ/Hγ/Hδ)
#     without over-smoothing the continuum between dense line regions.
#   upper_percentile=97: Procyon has broader, more numerous lines than the Sun;
#     pushing to 97th pct keeps the anchor closer to the true stellar continuum
#     and reduces the risk of the iterative fit latching onto a broad absorption wing.
#   sigma_clip / n_iter: same as solar — the iterative rejection logic is star-agnostic.
CONTINUUM_PARAMS = {
    'solar': {
        'knot_spacing_A'   : 100.0,
        'upper_percentile' : 95,
        'n_iter'           : 5,
        'sigma_clip'       : 3.0,
    },
    'procyon': {
        'knot_spacing_A'   : 50.0,
        'upper_percentile' : 97,
        'n_iter'           : 5,
        'sigma_clip'       : 3.0,
    },
}

# ── Per-star fundamental parameters + pin/solve policy (RYA-292) ──────────────
# Authoritative record of fundamental stellar parameters and an EXPLICIT pin/solve
# policy for parameter convergence. Generalizes the RYA-290 solar directive: any
# star whose logg is fixed by fundamental physics (mass + radius / asteroseismology)
# more precisely than an ionization-balance solve should PIN logg (and Teff),
# because floating a well-determined gravity lets data/setup error hide as a
# parameter walk (the solar logg→3.94 walk, RYA-289/290).
#
# Sources (Gaia FGK Benchmark Stars — same reference lineage as the RYA-203 Heiter
# GES loggf cross-match):
#   Teff, logg : Heiter et al. 2015, A&A 582, A49 — interferometric angular
#                diameter + bolometric flux (Teff); mass + radius / seismology (logg).
#   feh_ref    : Jofré et al. 2014, A&A 564, A133 — homogeneous reference [Fe/H].
#
# Policy keys per star:
#   pin   — params held FIXED at the fundamental value during convergence.
#   solve — params optimized by convergence: 'teff' (excitation/EP slope),
#           'logg' (Fe I/II ionization balance), 'xi' (reduced-EW slope).
#           'feh' is derived from A(Fe), not gradient-stepped.
# feh is SOLVED for every star here, so feh_ref is a reference value only (not a pin).
# Program stars WITHOUT a fundamental logg keep the full spectroscopic solve
# (pin=[]) so the benchmark pin never leaks into the unknown-logg path (RYA-290).
#
# Broadening (RYA-288): synthesis-v2 flux fitting is shape-sensitive (unlike EW),
# so each star's astrophysical broadening — vmac (radial-tangential macroturbulence)
# and vsini (projected rotation), km/s — is sourced HERE alongside the converged
# params and resolved through get_star_params(). Wrong broadening biases the fitted
# A(X) itself, not just its error bar, so there is NO default: a star whose synth-v2
# run needs broadening but has no vmac/vsini entry must FAIL LOUD (resolver raises),
# never silently inherit the Sun's profile. Solar vmac/vsini are lifted verbatim from
# the RYA-287 hardcoded set so the solar synth-v2 run is a no-op. Other stars get
# vmac/vsini only when sourced from their converged record (e.g. Procyon from RYA-284)
# — until then the entry is intentionally absent and the guard fires.
# Instrumental resolution is an INSTRUMENT constant, not per-star:
HARPS_R = 115000   # HARPS resolving power; instrumental convolution for synthesis

STAR_PARAMS = {
    # RYA-325: ξ (xi, microturbulence) is PINNED for the four GBS calibrators —
    # same status as Teff/logg (RYA-292 principle: floating a well-determined
    # param lets a systematic hide as a parameter shift; cf. the RYA-284 Procyon
    # ξ undershoot on then-dirty data). ξ is framework-relative (it absorbs
    # line-list/EW/1D-model choices), so the pinned value is the GBS reference on
    # a consistent scale: Jofré et al. 2014 (A&A 564, A133, GBS Paper III) — ξ via
    # the GES vmic relation at the Heiter+2015 (Paper I) fundamental params;
    # Bruntt et al. 2010 (MNRAS 405, 1907) as the independent cross-check. The
    # convergence still reports the reduced-EW slope at the pinned ξ as a
    # guardrail (flat = consistent; sloped = a finding to surface, never re-float).
    'solar': {
        'teff': 5772.0, 'e_teff': 1.0,  'logg': 4.438, 'e_logg': 0.000,
        'feh_ref': 0.00, 'e_feh': 0.03,
        'xi': 1.0,   # km/s — our scale anchor (Jofré+2014; matches STAR_SOLAR)
        'vmac': 3.8, 'vsini': 1.8,   # km/s — RYA-288: verbatim from the RYA-287
                                      # synth-v2 hardcoded set (solar run = no-op)
        'source': 'GBS: Heiter+2015 (Teff/logg), Jofré+2014 ([Fe/H], ξ)',
        'logg_basis': 'IAU nominal solar mass + radius (logg fixed to ~0.001 dex)',
        'pin': ['teff', 'logg', 'xi'], 'solve': ['feh'],
    },
    'procyon': {
        # RYA-309 §3.4: Teff/logg are the Heiter+2015 GBS values (6554/4.00),
        # confirmed by Ryan. The RYA-309 brief's "6530/3.96 (Heiter 2015)" was a
        # mislabel — those are Allende Prieto 2002 (spectroscopic) values.
        'teff': 6554.0, 'e_teff': 84.0, 'logg': 4.00, 'e_logg': 0.02,
        'feh_ref': 0.03, 'e_feh': 0.04,
        # RYA-325: ξ PINNED = 1.8 (GBS/Jofré scale, Ryan-confirmed 2026-06-15;
        # lit consensus 1.7-1.8, ≈ Bruntt+2010 1.69 — NOT the old ~2.0 solve-guess
        # on dirty data). HARD GATE: the reduced-EW slope on clean (post-309/320)
        # Fe data must be flat at 1.8; if it still insists on ~2.0, that's a finding
        # (F-star line-list/EW systematics, RYA-281/206), surfaced not pinned over.
        # RYA-322 2×2 produces the clean slope → confirms or flags.
        'xi': 1.8,
        # Broadening (RYA-309 §3.3): vsini FIXED at the Bruntt+2010 value; vmac is
        # an ADOPTED literature RT macroturbulence (Gray RT scale — same footing as
        # the solar 3.8 RT value). Bruntt's own vmac (4.6) is GAUSSIAN and must NOT
        # be reused as an RT vmac (cf. solar RT 3.8 vs Bruntt Gaussian 2.4).
        # WHY ADOPTED, not fit: at Procyon's SNR (~2272) the per-line χ²ᵣ is
        # dominated by 1D-LTE model error (~0.4%), not noise, so full-profile χ²
        # monotonically biases vmac high (rails to the search bound, not ~5-6) —
        # full-profile χ² is the wrong vmac estimator here. The proper FT/width-
        # based estimator is deferred to RYA-316 (non-blocking). A(Fe) sensitivity
        # to vmac 5.0↔6.0 is reported with the run (RYA-309 rigor rider) and is
        # small (vmac couples weakly to A(Fe)).
        'vsini': 2.8,   # km/s — FIXED. Bruntt et al. 2010, MNRAS 405, 1907
        'vmac': 6.0,    # km/s — ADOPTED RT (Gray scale), Procyon F5 IV-V ~5.5-6.5
        'source': 'GBS: Heiter+2015 (Teff/logg), Jofré+2014 ([Fe/H], ξ=1.8); '
                  'broadening: Bruntt+2010 vsini=2.8 + adopted RT vmac=6.0 '
                  '(Gray scale; FT estimator deferred to RYA-316)',
        'logg_basis': 'astrometric binary mass (Procyon B WD) + interferometric radius',
        'pin': ['teff', 'logg', 'xi'], 'solve': ['feh'],
    },
    'alpha_cen_a': {
        'teff': 5792.0, 'e_teff': 16.0, 'logg': 4.30, 'e_logg': 0.01,
        'feh_ref': 0.20, 'e_feh': 0.04,
        'xi': 1.1,   # km/s — GBS/Jofré scale (Ryan-confirmed; G2V near-solar,
                     # GES relation @ 5792/4.30/+0.20); slope-validate at α Cen run
        'source': 'GBS: Heiter+2015 (Teff/logg), Jofré+2014 ([Fe/H], ξ)',
        'logg_basis': 'dynamical binary mass + VLTI interferometric radius + asteroseismology',
        'pin': ['teff', 'logg', 'xi'], 'solve': ['feh'],
    },
    'alpha_cen_b': {
        'teff': 5231.0, 'e_teff': 20.0, 'logg': 4.53, 'e_logg': 0.03,
        'feh_ref': 0.20, 'e_feh': 0.04,
        'xi': 1.0,   # km/s — GBS/Jofré scale (Ryan-confirmed; cooler K1V,
                     # GES relation @ 5231/4.53/+0.20); slope-validate at α Cen run
        'source': 'GBS: Heiter+2015 (Teff/logg), Jofré+2014 ([Fe/H], ξ)',
        'logg_basis': 'dynamical binary mass + VLTI interferometric radius + asteroseismology',
        'pin': ['teff', 'logg', 'xi'], 'solve': ['feh'],
    },
    # 55 Cnc A (Copernicus) — flagship science-target primary, host of 55 Cnc e
    # (Janssen). RYA-293: graduates from the solve example to a PINNED target now
    # that von Braun et al. 2011 (ApJ 740, 49) supplies a fundamental logg.
    # Provenance caveat (recorded honestly): the radius is directly interferometric
    # but the mass is isochrone-based (model-dependent) — so this logg pin is
    # slightly softer than the dynamical-mass benchmarks (α Cen) or the IAU-nominal
    # Sun, but still far tighter (±0.01) than any spectroscopic solve, so pinning is
    # correct. Pinning the interferometric Teff (5196 K) may surface an excitation-
    # balance tension with our spectroscopic Teff when spectra run — that visibility
    # is intended (same pattern as the solar logg pin surfacing the Fe imbalance);
    # do NOT silently unpin Teff to hide it. feh_ref is a spectroscopic sanity
    # target (Valenti & Fischer 2005), NOT a pinned input — [Fe/H] stays solved.
    '55cnc_a': {
        'teff': 5196.0, 'e_teff': 24.0, 'logg': 4.45, 'e_logg': 0.01,
        'feh_ref': 0.31, 'e_feh': 0.04,
        # RYA-325: 55 Cnc is a SCIENCE TARGET, not a GBS calibrator — ξ stays
        # SOLVED in our framework (reduced-EW slope on our Fe lines, differential
        # to our solar ξ≈1.0). Literature ξ are heterogeneous & framework-relative
        # (Teske+2013 1.17±0.14 derived vs THEIR ξ_sun=1.38; Ecuvillon+2004
        # 0.98±0.07) — used only as the solve INITIAL GUESS + a cross-check
        # tripwire on a drifted solve (RYA-284 style), never pinned.
        'xi_init': 1.0,   # km/s — solve seed on our ~1.0-solar scale
        'xi_xcheck': (0.85, 1.35),   # scale-adjusted sane band (Teske/Ecuvillon)
        'source': 'von Braun et al. 2011, ApJ 740, 49; ξ x-check Teske+2013 / Ecuvillon+2004',
        'logg_basis': 'CHARA interferometric R=0.943 Rsun (theta_LD 0.711 mas + '
                      'van Leeuwen 2007 parallax) + Yonsei-Yale isochrone M=0.905 '
                      'Msun (model-dependent mass)',
        'pin': ['teff', 'logg'], 'solve': ['feh', 'xi'],
    },
    # Synthetic test fixture (RYA-293) — a star WITHOUT a fundamental logg, kept so
    # the spectroscopic-solve path (pin does not leak) stays covered after 55 Cnc A
    # became pinned. Not a real target; never run on spectra. Dummy params.
    'synthetic_no_logg': {
        'teff': 5800.0, 'e_teff': 100.0, 'logg': 4.40, 'e_logg': 0.20,
        'feh_ref': 0.00, 'e_feh': 0.10,
        'source': 'synthetic test fixture — no fundamental logg available',
        'logg_basis': 'no fundamental logg adopted — solved via Fe ionization balance',
        'pin': [], 'solve': ['teff', 'logg', 'feh', 'xi'],
    },
}


def get_star_params(star_id: str) -> dict:
    """
    Resolve a star_id to its STAR_PARAMS record (fundamental params + pin/solve
    policy). Fail LOUD if there is no record — no silent default (RYA-292; same
    discipline as the RYA-288 broadening guard). Matching mirrors the substring
    convention used elsewhere for star routing.
    """
    s = star_id.strip().lower().replace(' ', '_').replace('-', '_')
    if s in STAR_PARAMS:
        return STAR_PARAMS[s]
    for key, rec in STAR_PARAMS.items():
        if key in s or s in key:
            return rec
    raise KeyError(
        f"No STAR_PARAMS record for star_id='{star_id}'. Add a fundamental-params "
        f"+ pin/solve entry in config/constants.py before running it — refusing to "
        f"guess a pin/solve policy or fall back to a default (RYA-292)."
    )


# ── iSpec / Turbospectrum installation ───────────────────────────────────────
# Path to cloned iSpec repo (machine-specific; override via ISPEC_DIR env var).
import os as _os
ISPEC_DIR = Path(_os.environ.get('ISPEC_DIR',
    str(ROOT.parent / 'ispec')
))
RADIATIVE_TRANSFER_CODE = 'turbospectrum'   # 'turbospectrum' | 'spectrum' | 'moog'

# ── EW→abundance baseline engine (RYA-289) ───────────────────────────────────
# The community-standard EW→abundance engine for the synthesis-vs-EW comparison
# (RYA-285). This MUST be chosen INDEPENDENTLY of RADIATIVE_TRANSFER_CODE:
# setting the RT code to 'turbospectrum' for the synthesis path must never
# silently downgrade this baseline to SPECTRUM. That silent downgrade is exactly
# the bug RYA-289 fixes — RYA-285/287 ran SPECTRUM, never MOOG. MOOG returns
# absolute A(X) uniformly for every star (no NLTE-layer reconstruction), which
# also dissolves the solar A_X≡0 / differential-vs-absolute scale split.
EW_BASELINE_CODE = 'moog'   # 'moog' | 'moog-scat' | 'spectrum' | 'width'

# Auto-create all output directories on import
for _key in ('data_root', 'raw_spectra', 'processed_spectra', 'linelists',
             'model_atmospheres', 'results', 'plots', 'tables'):
    PATHS[_key].mkdir(parents=True, exist_ok=True)

# ── Line quality scoring (RYA-220) ───────────────────────────────────────────
# Weights for the per-line quality score aggregation formula.
# Sum = 1.0. Scientific rationale for each weight:
#   proximity (0.20): VALD density is a known source of systematic EW bias; significant
#                     but correctable via synthesis — moderate weight.
#   snr (0.25):       EW precision scales directly with line SNR; primary measurement
#                     quality driver — highest weight.
#   chi2 (0.20):      Profile fit quality captures blending, continuum errors, and
#                     non-Gaussian wings; high diagnostic value.
#   saturation (0.15):saturated lines have non-linear COG; EW→abundance mapping degrades
#                     above ~150 mÅ — meaningful but secondary to SNR.
#   outlier (0.10):   Abundance outliers signal unresolved systematics; lower weight
#                     because the sigma is itself uncertain with small n.
#   nlte (0.10):      NLTE correction coverage is an analysis quality indicator, not a
#                     measurement quality indicator — lowest weight.
LINE_SCORE_WEIGHTS = {
    'proximity':  0.20,
    'snr':        0.25,
    'chi2':       0.20,
    'saturation': 0.15,
    'outlier':    0.10,
    'nlte':       0.10,
}

# Grade boundaries for line_score (0.0–1.0).
LINE_GRADE_THRESHOLDS = {
    'A': 0.80,   # Clean, NLTE corrected. Enters weighted abundance mean.
    'B': 0.60,   # Minor issues. Enters weighted mean with footnote.
    'C': 0.35,   # Significant issues. Measured and reported with caveats.
    'D': 0.00,   # Unusable. Archived only.
}

# Sub-score normalisation parameters. No hardcoded values elsewhere in scoring code.
LINE_SCORE_PARAMS = {
    'snr_saturation':         30.0,  # EW/ew_err at which ew_snr_score = 1.0
    'chi2_clean_max':          1.5,  # reduced χ² below which fit_chi2_score = 1.0
    'chi2_floor':             10.0,  # reduced χ² above which fit_chi2_score = 0.0
    'saturation_ew_low_mA':  100.0,  # EW below which saturation_score = 0.0
    'saturation_ew_high_mA': 200.0,  # EW above which saturation_score = 1.0
    'outlier_sigma_scale':     2.0,  # |A_i − median| / (scale × σ) → outlier score
}

# ── Post-LTE corrections (apply in abundances_derive.py) ─────────────────────
# All values are additive corrections to log N(X)/N(H) derived from 1D LTE.
CORRECTIONS_3D = {
    # O I forbidden lines: 1D → 3D hydrodynamic correction (Caffau et al. 2008, A&A 483, 591)
    # Both lines share the same -0.07 dex offset from granulation effects.
    'O_6300_3d_dex': -0.07,
    'O_6363_3d_dex': -0.07,
}

# ── Model atmosphere parameters ───────────────────────────────────────────────
MODEL = {
    'type'              : 'ATLAS9',          # Castelli-Kurucz
    'geometry'          : 'plane-parallel',  # Valid for dwarfs (log g > 3.5)
    'radiative_transfer': 'LTE',
    'grid_source'       : 'https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas',
    'nlte_corrections'  : 'Amarsi et al. 2020',  # Applied post-LTE for O, Na, C
}


# ── Validation ────────────────────────────────────────────────────────────────

_REQUIRED_KEYS = {
    'PHYSICS'          : ['c_kms', 'c_ms', 'h_eV', 'k_eV', 'k_cgs', 'sigma_sb', 'G_cgs', 'amu_g'],
    'ASTRO'            : ['Msun_g', 'Rsun_cm', 'Lsun_erg', 'Teff_sun', 'logg_sun', 'pc_cm', 'AU_cm'],
    'SOLAR_ASPLUND2021': ['H', 'C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni', 'Cu', 'Sr'],
    'STAR_55CNC'       : ['teff_K', 'logg', 'feh', 'vturb_kms', 'rv_kms', 'hd', 'hip'],
    'STAR_SOLAR'       : ['teff_K', 'logg', 'feh', 'vturb_kms', 'rv_kms', 'program_id', 'citation'],
    'PIPELINE'         : ['wav_min_A', 'wav_max_A', 'snr_min_science', 'ew_min_mA', 'ew_max_mA',
                          'continuum_n_iter', 'continuum_knot_spacing_A', 'continuum_upper_percentile',
                          'cont_edge_frac', 'min_fit_depth', 'blue_edge_warn_A'],
    'PATHS'            : ['data_root', 'raw_spectra', 'results', 'plots',
                          'solar_spectra', 'solar_normalized', 'solar_diagnostic',
                          'procyon_spectra', 'procyon_normalized', 'procyon_diagnostic'],
}

_DICTS = {
    'PHYSICS'          : PHYSICS,
    'ASTRO'            : ASTRO,
    'SOLAR_ASPLUND2021': SOLAR_ASPLUND2021,
    'STAR_55CNC'       : STAR_55CNC,
    'STAR_SOLAR'       : STAR_SOLAR,
    'PIPELINE'         : PIPELINE,
    'PATHS'            : PATHS,
}


def validate_constants():
    """Check all required keys exist and critical values are sane. Raises on failure."""
    missing = []
    for dict_name, required_keys in _REQUIRED_KEYS.items():
        d = _DICTS[dict_name]
        for key in required_keys:
            if key not in d:
                missing.append(f"{dict_name}['{key}']")
    if missing:
        raise KeyError(f"Missing constants: {', '.join(missing)}")

    assert abs(PHYSICS['c_kms'] - 299792.458) < 0.001, "Speed of light mismatch"
    assert abs(SOLAR_ASPLUND2021['Fe'] - 7.46) < 0.01, "Solar Fe abundance mismatch"
    assert 4000 < STAR_55CNC['teff_K'] < 7000, "Teff out of range"

    print("Constants validated ✓")
