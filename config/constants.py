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
# RYA-298: fundamental teff/logg/feh live ONLY in stars.yaml (get_star_params).
# This legacy dict keeps the 55 Cnc-specific NON-fundamental fields consumers read
# (names, mass/radius, broadening, RV) — no stellar-param literals (was 5196/4.41/0.32;
# 4.41/0.32 had drifted from the canonical 4.45/0.31, now single-sourced).
STAR_55CNC = {
    'name'         : '55 Cancri A',
    'iau_name'     : 'Copernicus',
    'hd'           : 'HD 75732',
    'hip'          : 'HIP 43587',
    'vturb_kms'    : 0.9,       # km/s ξ-solve SEED (55 Cnc ξ is solved, not pinned)
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
# RYA-298: fundamental teff/logg/feh live ONLY in stars.yaml (get_star_params).
# The legacy 5777/4.44 here had drifted from the canonical 5772/4.438 (IAU-nominal /
# GBS Heiter+2015) — RYA-355 flagged it; now single-sourced. Non-fundamental solar
# fields (provenance, RV, Ni 6300 literature EW) stay here.
STAR_SOLAR = {
    'name'         : 'Sun',
    'program_id'   : '1102.D-0954(A)',
    'citation'     : 'Dumusque et al. 2021, arXiv:2009.01945',
    'vturb_kms'    : 1.0,       # km/s ξ-solve seed (solar ξ pinned 1.0 in stars.yaml)
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
    'vturb_kms' : 1.66,     # km/s — ξ-solve seed (Allende Prieto 2002), not pinned
    'rv_kms'    : -3.0,
}

# ── Fe abundance gate thresholds (RYA-261) ───────────────────────────────────
# Relative [Fe/H] units — matches iSpec output convention (0.0 at solar).
# Absolute equivalent: A(Fe) = [Fe/H] + SOLAR_ASPLUND2021['Fe'] (= [Fe/H] + 7.46)
# FE_GATE_LOWER = -0.05  →  A(Fe) = 7.41
# FE_GATE_UPPER = +0.05  →  A(Fe) = 7.51
FE_GATE_LOWER      = -0.05   # [Fe/H] solar Fe lower bound
FE_GATE_UPPER      = +0.05   # [Fe/H] solar Fe upper bound  → absolute [7.41, 7.51]
FE_SCATTER_GATE    =  0.10   # dex — maximum acceptable σ(Fe I)
FE_IONISATION_GATE =  0.05   # dex — maximum |ΔFe(I−II)|

# ── Fe I−Fe II ionization arbiter: SYNTHESIS, not the EW path (RYA-406, DECISION 2) ──
# The Fe ionization-balance VERDICT is scored on flux-space SYNTHESIS Fe II (the ratified
# arbiter, RYA-305/341 DECISION 2), NOT the EW-path Fe II. The EW Fe II pool is blend-
# limited and reads HIGH BY DESIGN (RYA-352: clean EW pool 7.696; the EW path is
# deliberately not forced onto Asplund), so its ΔFe(I−II) is a DIAGNOSTIC, not the verdict
# (RYA-405/406). Values are CITED single-source, never tuned: synth Fe II 7.486 (RYA-305)
# … 7.500 (RYA-341 crowned), balanced ΔFe(I−II) = −0.007 (RYA-305) / −0.015 (RYA-341).
# The synth Fe II's OWN uncertainty (RYA-305 7.486 vs RYA-338 7.571 ≈ 0.085 spread; thin /
# core-broken pool RYA-347) is reported; the balance is scored against the principled
# FE_IONISATION_GATE (0.05) — the ratified −0.015 is well within it, so NO threshold is
# loosened to pass (validate-don't-tune). Keyed by star: only the Sun is a ratified
# calibrator; other stars consume their own synth-v2 product or are flagged (no silent EW
# fallback as a verdict).
FE_IONIZATION_SYNTH_ARBITER = {
    'solar': {
        'dFe'             : -0.015,   # Fe I 7.485 − Fe II 7.500 (RYA-341 crowned; RYA-305 = −0.007)
        'fe2_synth'       :  7.500,   # ratified synth Fe II (RYA-341; RYA-305 = 7.486)
        'fe2_uncertainty' :  0.085,   # RYA-305 7.486 vs RYA-338 7.571 method spread
        'provenance'      : 'RYA-305/341/405 (DECISION 2: synthesis is the Fe ionization arbiter)',
    },
}
# Expected EW-vs-synth Fe II blend-bias band: above this, the spread is FLAGGED as a
# finding rather than waved through as "high by design" (RYA-406 Socratic check). Set from
# the documented EW-high (RYA-352 7.696) vs synth (RYA-341 7.500) ≈ 0.20 + the synth-pool
# uncertainty (0.085) headroom.
FE_EW_SYNTH_SPREAD_BAND = 0.30   # dex — |EW Fe II − synth Fe II| above this → flag, not "by design"

# ── Scale-aware solar Fe gate (RYA-336) ──────────────────────────────────────
# The solar Fe verdict rests on SCALE-ROBUST checks: the reduced-EW slope (line-
# formation shape), the Fe I−Fe II ionization balance (cancels a uniform Fe zero-
# point), and the Fe I scatter. None of these carry the 1D-vs-3D abundance term, and
# they are what the differential science actually rests on.
#
# The ABSOLUTE A(Fe) is a scale-aware DIAGNOSTIC, not the verdict. Our NLTE grids
# output on the 1D-NLTE scale, which runs ~+0.05 dex above the 3D-true solar value:
# 3D granulation strengthens Fe I lines, so a 1D atmosphere needs a higher abundance
# to match. The diagnostic centre is therefore the in-repo 3D-true anchor (Asplund
# 2021, SOLAR_ASPLUND2021['Fe'] = 7.46) PLUS the published solar 1D-3D Fe offset
# below — NOT our own output (that would be circular). Both terms are independent
# published quantities, so the centre is honest. The window is WIDE: a sanity bound
# that catches a gross zero-point error (e.g. a loggf scale slip) — which every
# shift-invariant strong gate misses — while passing the legitimate grid spread
# (Amarsi 7.495 / MPIA 7.516, straddling ~7.51).
#
# Upgrade path: when an APPLIED 3D correction lands (RYA-285 synthesis / Magic 2013),
# the absolute output moves onto the true 7.46 scale and FE_GATE [7.41,7.51] applies
# again — this scaffolding retires. (The RYA-334 A_X<12.5 guard is the separate
# gross double-add tripwire, unrelated to this science sanity bound.)
FE_1D3D_SOLAR_OFFSET  = 0.05   # dex — solar A(Fe I) on a 1D atmosphere runs ~+0.05
                               # above the 3D value (granulation term). Literature:
                               # Magic et al. 2013 (Stagger 3D grid, A&A 557 A26);
                               # cf. Amarsi et al. 2019 (A&A 624 A111) & Asplund et
                               # al. 2021 (A&A 653 A141) §3-4. NOT from our output.
FE_ABS_DIAG_HALFWIDTH = 0.07   # dex — wide half-window around the scale-aware centre
                               # (≈7.51 → [7.44, 7.58]); passes grid spread + loggf
                               # zero-point, fails gross errors.
FE_REW_SLOPE_GATE     = 0.10   # |reduced-EW slope| at pinned ξ — "flat within error"
                               # for the solar guardrail (RYA-330: −0.04..−0.09).

# ── Flux-space synthesis (synth-v2) line-acceptance gate (RYA-342) ────────────
# The flux-space fitter accepts/rejects each line on FIT QUALITY (reduced χ²), NOT
# on the EW saturation ceiling (vmic_ew_ceiling_mA — an EW-path concept). Flux-space
# synthesis exists to measure exactly the strong/blended lines EW saturation kills
# (RYA-338 FLAG 2), so gating it on an EW ceiling defeats its purpose and silently
# suppresses real strong lines for every star/element (RYA-341 found this on solar
# Fe II: 6239.943 EW 112 mÅ was dropped as "saturated" yet flux-fits at χ²ᵣ 1.32).
# The χ²ᵣ here is on iSpec's model-adequacy-floor scale (order-unity = synthetic
# matches observed; ≫1 = unmodeled blend / model error, abundance unreliable). The
# solar Fe II distribution is starkly bimodal — clean lines χ²ᵣ ≤ 3.0, blend-
# dominated lines χ²ᵣ ≥ 128, with a wide empty gap (3.0 → 128) — so 10.0 separates
# them robustly without being delicate. Rejected lines are logged with their χ²ᵣ.
SYNTH_CHI2_GATE       = 10.0   # max reduced-χ² for a synth-v2 line to enter the median

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

    # Fe II EW-quality cull (RYA-352) — re-bases the RYA-239/a88ef0f "fe2_solar_cull_
    # rya347" hack onto its REAL basis: EW measurement quality, NOT RYA-347's synth
    # χ²ᵣ verdict. The two paths are independent; a line failing a synth fit is not a
    # reason to drop a clean EW measurement. Criteria (computed at load, no frozen list):
    #   (1) saturation — EW above vmic_ew_ceiling_mA (the COG-derived ceiling, RYA-330;
    #       lines on the flat curve-of-growth over-read because 1D-LTE can't reproduce
    #       their saturation). Expressed per-line as REW for the report.
    #   (2) measurement quality — ew_err/EW above the threshold below.
    #   (3) blend — blend_flag == True (already excluded in _ew_to_abundance; kept here
    #       for a unified, logged criterion set).
    # INTERIM Fe II-only (RYA-199 lesson: a fixed REW ceiling right for Fe II wrongly
    # kills strong Ba/Na lines — the general element/COG-aware EW gate is routed to the
    # non-Fe curation track, RYA-349). Do NOT blanket this across elements.
    'fe2_ew_err_frac_max'  : 0.5,  # cull Fe II EW line if ew_err_mA/ew_mA exceeds this.

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

    # Method-selection policy (RYA-315 / RYA-306): per (star × species) EW-vs-synthesis matrix
    'method_policy'     : ROOT / 'data' / 'method_policy.yaml',

    # Solar calibration spectra (external data — relative to repo parent)
    'solar_spectra'     : ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Solar Calibration' / 'archive',

    # Procyon HARPS spectra (external data — relative to repo parent)
    'procyon_spectra'   : ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Procyon' / 'Procyon Harps',

    # Processed outputs (CSV gitignored — generated by pipeline; PNG committed to results/)
    'solar_normalized'  : ROOT / 'data' / 'processed' / 'solar_normalized.csv',
    'solar_diagnostic'  : ROOT / 'results' / 'plots' / 'solar_continuum_diagnostic.png',
    # RYA-408: solar_ew.csv is the gitignored lines_fit STAGING output ONLY.
    # It must NOT be a gate/abundance EW input — gates read the committed canonical below.
    'solar_ew'          : ROOT / 'data' / 'processed' / 'solar_ew.csv',
    # Committed CANONICAL solar EW measurements — single source of truth for gates,
    # abundance derivation, and stewardship (RYA-408; promotion via scripts/promote_solar_ew.py).
    'solar_ew_canonical': ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv',
    'solar_ew_diagnostic': ROOT / 'results' / 'plots' / 'solar_ew_diagnostic.png',
    'procyon_normalized' : ROOT / 'data' / 'processed' / 'procyon_normalized.csv',
    'procyon_diagnostic' : ROOT / 'results' / 'plots' / 'procyon_continuum_diagnostic.png',
    'procyon_ew'         : ROOT / 'data' / 'processed' / 'procyon_ew.csv',
    'procyon_ew_diagnostic': ROOT / 'results' / 'plots' / 'procyon_ew_diagnostic.png',
}

# ── Telluric data-input stage (RYA-424) ──────────────────────────────────────
# Telluric correction is a STANDING, instrument-aware data-input stage (RYA-424),
# promoting the RYA-373 molecfit driver + RYA-380 per-night-GDAS recipe from
# per-dataset one-offs into a gate every red-optical / IR dataset passes through
# before it can be marked analysis-ready.
#
# THE RULE IS WAVELENGTH-GATED, NOT INSTRUMENT-GATED (RYA-380, codified here):
#   * λ ≳ TELLURIC_LAMBDA_MIN_A  (red-optical + IR): sharp telluric forest
#     (H2O / O2 A-band / CH4 / CO2) — molecfit + observation-night GDAS MANDATORY.
#     Includes red-optical 6910–9500 Å (O I 7772, K I 7699, N I 8216), not only IR.
#   * λ ≲ TELLURIC_BLUE_MAX_A  (blue/UV) and space-UV (HST/STIS, COS): NO telluric
#     step (space-UV is above the atmosphere; blue = broad ozone/Rayleigh handled
#     by normalization).
#   * mid-optical 3800–6800 Å: largely clean; case-by-case (NOT auto-required).
TELLURIC_LAMBDA_MIN_A = 6800.0   # ≳ this (Å) → telluric forest → correction MANDATORY
TELLURIC_BLUE_MAX_A   = 3800.0   # ≲ this (Å) → blue/UV → no telluric-line step

# Wavelength-regime labels (the gate's verdict on a [lo, hi] Å span).
TELLURIC_REGIME_IR          = 'ir'            # fully ≳ 6800 Å (sharp forest)
TELLURIC_REGIME_RED_OPTICAL = 'red_optical'   # crosses/sits in 6800–~10000 Å
TELLURIC_REGIME_MID_OPTICAL = 'mid_optical'   # 3800–6800 Å, case-by-case
TELLURIC_REGIME_BLUE        = 'blue'          # ≲ 3800 Å ground-based
TELLURIC_REGIME_SPACE_UV    = 'space_uv'      # HST/STIS/COS — above the atmosphere

# Instrument → telluric ENGINE, single source of truth (NO silent default — an
# unknown instrument LOUD-FAILS; engine is chosen explicitly, never guessed).
#   'molecfit'      → ESO molecfit (esorex molecfit_model/calctrans) + per-night GDAS
#                     (RYA-373 driver + RYA-380 recipe).
#   'apero_wapiti'  → APERO extraction + Wapiti telluric (SPIRou permanent rule —
#                     NOT molecfit; SPIRou is reduced by its own DRS).
TELLURIC_ENGINE_MOLECFIT     = 'molecfit'
TELLURIC_ENGINE_APERO_WAPITI = 'apero_wapiti'
TELLURIC_ENGINES = {
    'CRIRES':        TELLURIC_ENGINE_MOLECFIT,
    'CRIRES+':       TELLURIC_ENGINE_MOLECFIT,
    'UVES-red':      TELLURIC_ENGINE_MOLECFIT,
    'ESPRESSO-red':  TELLURIC_ENGINE_MOLECFIT,
    'FEROS-red':     TELLURIC_ENGINE_MOLECFIT,
    'CHIRON':        TELLURIC_ENGINE_MOLECFIT,
    'NIRPS':         TELLURIC_ENGINE_MOLECFIT,
    'SPIRou':        TELLURIC_ENGINE_APERO_WAPITI,
}

# Instrument → observatory site for the per-night GDAS profile (molecfit engine).
# GDAS tarballs are per-site; the molecfit GDAS loc code is "C{lon}{lat}" with EAST-
# positive longitude (Paranal = -70.4, -24.6). SPIRou (CFHT) uses APERO+Wapiti, not
# GDAS, so it has no molecfit GDAS site (None).
SITES = {
    'paranal':  {'name': 'ESO Paranal (VLT)', 'lat_deg': -24.6, 'lon_deg': -70.4,
                 'elev_m': 2635.0, 'horizons_code': '309', 'gdas_loc': 'C-70.4-24.6'},
    'la_silla': {'name': 'ESO La Silla',      'lat_deg': -29.26, 'lon_deg': -70.73,
                 'elev_m': 2400.0, 'horizons_code': '809', 'gdas_loc': 'C-70.7-29.3'},
    'ctio':     {'name': 'CTIO (SMARTS)',     'lat_deg': -30.17, 'lon_deg': -70.81,
                 'elev_m': 2200.0, 'horizons_code': '807', 'gdas_loc': 'C-70.8-30.2'},
    'cfht':     {'name': 'CFHT (Maunakea)',   'lat_deg': 19.83, 'lon_deg': -155.47,
                 'elev_m': 4204.0, 'horizons_code': '568', 'gdas_loc': None},
}
TELLURIC_INSTRUMENT_SITE = {
    'CRIRES':        'paranal',
    'CRIRES+':       'paranal',
    'UVES-red':      'paranal',
    'ESPRESSO-red':  'paranal',
    'FEROS-red':     'la_silla',
    'CHIRON':        'ctio',
    'NIRPS':         'la_silla',
    'SPIRou':        'cfht',
}

# Telluric-residual VERIFICATION tolerance (RYA-424 §3 gate). Score is the median
# |1 − corrected_flux/continuum| at telluric-DOMINATED, science-clean pixels in the
# known telluric windows; a frame whose residual exceeds this is FLAGGED, not silently
# passed. 0.05 (5%) inherits the RYA-373 D1 telluric-specific gate (telluric_residual_gate).
TELLURIC_RESIDUAL_TOL = 0.05

# ── Non-Fe NLTE correction registry (RYA-235) ────────────────────────────────
# Element-level NLTE corrections applied PER LINE after the Fe I/II leg, from the
# real MPIA SpectrumTools MAFAGS-OS grids (scraped RYA-244/245; columns
# element,wave_A,teff_K,logg,feh,delta_nlte). The feh axis is [Fe/H] RELATIVE
# (0 at solar) — pass stellar_params['feh'] directly, NO solar offset (unlike the
# absolute-A(Fe) Amarsi Fe grid). Cr II is deliberately excluded (the RYA-232 −0.777
# dex was two saturated lines in the COG damping wing, a line-matching artifact, not a
# real systematic — RYA-240). Ca/Ti/Cr/Mg/Mn/Si are neutral (ion 1); Ba is
# singly-ionised (Ba II, ion 2 — the dominant Ba stage in FGK photospheres).
# The Asplund-2021 solar anchors (SOLAR_ASPLUND2021) are the validation targets,
# but acceptance is gated on a CURATED line pool (RYA-235 reopen: raw 1D-LTE Cr
# already reads high; NLTE on an un-curated pool overshoots — curate first).
# `flag` overrides the default per-source NLTE provenance tag ('NLTE_MPIA_MAFAGS_1D');
# set it where the grid is NOT an MPIA MAFAGS-OS extraction (Na = INSPECT/Lind,
# Ba = Korotin/VizieR) so a non-MPIA grid is never mislabelled as MPIA (RYA-165).
NLTE_CORRECTION_ELEMENTS = {
    # Ca I — STAYS on Mashonkina2017 MPIA. RYA-411 ADJUDICATED the RYA-410 STOP (+0.064
    # vs +0.012): it is NOT a model-atom difference — Amarsi-2020 (+0.017) and Mashonkina
    # (+0.027) AGREE on the only cleanly-computed shared line (5867). The gap was (a)
    # line-set (RYA-410 used 6166/6169/6455, where Amarsi gives real +0.05..+0.08) and
    # (b) a DATA DEFECT: MPIA's Ca 6166 was a placeholder zero (all 72 nodes = 0.000),
    # dragging the MPIA reference down. RYA-413 FIXED it: confirmed live that MPIA OFFERS
    # 6166 but RETURNS a literal 0 (serves but does not NLTE-model the line) -> dropped the
    # line (never carry a zero as a correction) -> clean MPIA solar Ca median = +0.017
    # (was +0.0118). A registration-time placeholder-zero guard now refuses the class.
    # Adopting Amarsi would close the 55 Cnc clamp AND add the missing NLTE on the
    # subordinate 3d-4p lines (~+0.05 solar Ca shift) — a solar-calibration change to RATIFY
    # against a solar-Ca re-validation (RYA-371), not a silent swap here. Ca keeps MPIA.
    'Ca': {'ion': 1, 'grid': 'Ca_Mashonkina2017.csv',
           'ref': 'Mashonkina et al. 2017 (A&A 606, A147), MPIA MAFAGS-OS 1D NLTE',
           'mpia_species': '20.01'},
    'Ti': {'ion': 1, 'grid': 'Ti_Bergemann2011_MPIA.csv',
           'ref': 'Bergemann 2011 (MNRAS 413, 2184), MPIA MAFAGS-OS 1D NLTE',
           'mpia_species': '22.01'},
    'Cr': {'ion': 1, 'grid': 'Cr_Bergemann2010_MPIA.csv',
           'ref': 'Bergemann & Cescutti 2010 (A&A 522, A9), MPIA MAFAGS-OS 1D NLTE',
           'mpia_species': '24.01'},
    # RYA-165 (grids vendored RYA-396): the rescoped Na/Mg/Ba leg + the Mn/Si gap-fill.
    # Li deferred (RYA-103 Li 6707 CN-blend); C/O owned by nlte_cno (Amarsi 2019).
    # Asplund acceptance still gated on a curated pool (Na +0.27 / Ba +0.19 NLTE-sized;
    # Mg's +1.72 is saturation/curation, NLTE here is a rounding-error rung).
    # RYA-410: re-sourced onto the Amarsi-2020 PySME departure grid ([Fe/H] -> +1, so
    # the 55 Cnc +0.31 clamp is closed — the MPIA/INSPECT grids ceilinged at +0.30).
    # Validate-don't-tune: the PySME-Amarsi SOLAR delta REPRODUCED the prior registered
    # value within tol before the swap (Na -0.129 vs INSPECT -0.107; Mg -0.022 vs MPIA
    # +0.013, both ~0; Si -0.013 vs MPIA -0.004). Now the single source per element.
    'Na': {'ion': 1, 'grid': 'Na_Amarsi2020_PySME.csv',
           'ref': 'Amarsi et al. 2020 (A&A 642, A62) departure grid via PySME; Na I model '
                  'atom Lind et al. 2011 (A&A 528, A103); RYA-410 PySME-derived deltas '
                  '(5682/5688; reproduces the INSPECT Lind-2011 solar -0.107)',
           'flag': 'NLTE_Amarsi2020_PySME_1D'},
    'Mg': {'ion': 1, 'grid': 'Mg_Amarsi2020_PySME.csv',
           'ref': 'Amarsi et al. 2020 (A&A 642, A62) departure grid via PySME; RYA-410 '
                  'PySME-derived deltas (5711/4730; Mg NLTE ~0 in FGK, solar -0.022)',
           'flag': 'NLTE_Amarsi2020_PySME_1D'},
    'Ba': {'ion': 2, 'grid': 'Ba_Korotin2015.csv',
           'ref': 'Korotin et al. 2015 (A&A 581, A70), CDS VizieR J/A+A/581/A70',
           'flag': 'NLTE_Korotin2015_1D'},
    # Mn I — STAYS on MPIA. RYA-411 RESOLVED the RYA-410 STOP (+0.018 vs +0.107). Built an
    # HFS-resolved synthesis (pysme_nlte _synth_ew RYA-411) and ran the probes: HFS-collapse
    # was real but MINOR (the 6013/6016/6021 triplet only moves +0.018 -> +0.024 HFS-
    # resolved); the harness had also picked the wrong (low-EP) lines. The DECISIVE test —
    # same lines as MPIA (4998/6304/6306/6867), HFS-resolved — still gives Amarsi-2020
    # ~+0.05 (clean lines 6304/6306/6867 = +0.048/+0.049/+0.058) vs MPIA/Bergemann +0.107.
    # So the residual is a GENUINE Amarsi-2020-vs-Bergemann Mn model-atom difference (~2x),
    # not HFS/line-set. The dedicated Mn-NLTE literature (Bergemann & Gehren 2008) supports
    # the larger correction -> KEEP MPIA (do not adopt the smaller survey-grid value just to
    # close the clamp = would be tuning). Mn stays clamped-but-loud at 55 Cnc (RYA-409).
    'Mn': {'ion': 1, 'grid': 'Mn_Bergemann_MPIA.csv',
           'ref': 'Bergemann group MPIA SpectrumTools MAFAGS-OS (Mn I; cf. Bergemann & Gehren 2008, A&A 492, 823)',
           'mpia_species': '25.01'},
    # Si I — RYA-410: re-sourced onto the Amarsi-2020 PySME grid (closes the 55 Cnc
    # clamp). Cross-check PASSED (-0.013 vs MPIA -0.004). Si NLTE negligible in FGK.
    'Si': {'ion': 1, 'grid': 'Si_Amarsi2020_PySME.csv',
           'ref': 'Amarsi et al. 2020 (A&A 642, A62) departure grid via PySME; RYA-410 '
                  'PySME-derived deltas (5772/6125; Si NLTE ~0 in FGK, solar -0.013)',
           'flag': 'NLTE_Amarsi2020_PySME_1D'},
    # RYA-402 Family-B (Option 2): derived by NLTE synthesis in PySME from the Amarsi
    # 2020 GALAH departure grid (Al I model atom = Nordlander & Lind 2017). Validated
    # vs the published subordinate-Al correction (solar median -0.022, in the
    # near-zero-to--0.04 band) and confirmed LTE-baseline-portable (eps -0.011 vs the
    # TS/EW baseline). Subordinate doublet 6696/6698 (NOT the +0.2 resonance 3961).
    'Al': {'ion': 1, 'grid': 'Al_Amarsi2020_PySME.csv',
           'ref': 'Amarsi et al. 2020 (A&A 642, A62) departure grid via PySME; Al I model '
                  'atom Nordlander & Lind 2017 (A&A 607, A75); RYA-402 PySME-derived deltas',
           'flag': 'NLTE_Amarsi2020_PySME_1D'},
    # RYA-402 Family-B (Option 2): S I multiplet-8 optical lines 6748/6757 (high-excit,
    # small negative NLTE, solar median -0.016) derived via PySME from the Amarsi 2025
    # (A&A 703, A35) S departure grid. Validated vs the small-negative optical-S band;
    # the NIR 9212/9228/9237 triplet (large NLTE) is a separate indicator (RYA-401).
    'S': {'ion': 1, 'grid': 'S_Amarsi2025_PySME.csv',
          'ref': 'Amarsi et al. 2025 (A&A 703, A35) S departure grid via PySME; '
                 'RYA-402 PySME-derived deltas (optical 6748/6757)',
          'flag': 'NLTE_Amarsi2025_PySME_1D'},
}

# ── 3D abundance corrections — the metal 3D leg (RYA-399) ─────────────────────
# RYA-400 routed Si/Ti/Cr as "3D-owed". This registry holds the 3D DIMENSIONAL
# correction (A(3D) − A(1D), at fixed NLTE) applied ON TOP of the existing 1D-NLTE
# pass (NLTE_CORRECTION_ELEMENTS) — i.e. the increment, so it never double-counts
# the NLTE already applied:  A(3D-NLTE) = A(1D-NLTE) + delta_3d.  Mirrors the C/N/O
# 3D leg (pipeline.nlte_cno) which adds a STAGGER-grid Δ to the 1D abundance.
#
# KEY FINDING (RYA-399, validate-don't-tune): the published solar 3D corrections for
# these metals are SMALL (|Δ| ≲ 0.1 dex) — far below the +0.38/+0.50/+0.40 dex
# graded-pool residuals (RYA-398). For Ti/Cr they are POSITIVE (3D raises the
# abundance), so 3D moves them the WRONG way. 3D is therefore NOT the lever that
# closes the Si/Ti/Cr solar residual; the residual is line-data / gf-zero-point
# (RYA-161 differential territory) and is carried forward, never tuned away.
#
# Solar node only (Teff 5772, logg 4.44, [Fe/H] 0). A full per-(Teff,logg,feh) 3D
# grid exists publicly only for Si (Amarsi & Asplund 2017; Amarsi 2020 GALAH);
# Ti/Cr have NO public off-solar 3D-NLTE grid (Fe-peak excluded from Amarsi 2020) —
# off-solar 3D for those is documented as a gap and carried to the multi-star arc.
THREED_CORRECTION_ELEMENTS = {
    'Si': {'ion': 1, 'grid': 'solar3d_metals_rya399.csv',
           'ref': 'Amarsi & Asplund 2017 (MNRAS 464, 264), 3D non-LTE Si; '
                  'Scott et al. 2015 Paper I (A&A 573, A25), A(Si)=7.51',
           'off_solar_grid': 'Amarsi & Asplund 2017 / Amarsi 2020 GALAH (A&A 642, A62)'},
    'Ti': {'ion': 1, 'grid': 'solar3d_metals_rya399.csv',
           'ref': 'Scott et al. 2015 Paper II (A&A 573, A26), solar 3D Ti, A(Ti)=4.93',
           'off_solar_grid': None},   # no public off-solar 3D-NLTE Ti grid
    'Cr': {'ion': 1, 'grid': 'solar3d_metals_rya399.csv',
           'ref': 'Scott et al. 2015 Paper II (A&A 573, A26), solar 3D Cr, A(Cr)=5.62',
           'off_solar_grid': None},   # no public off-solar 3D-NLTE Cr grid
}

# ── Reflected-solar bodies — SINGLE SOURCE for JPL Horizons ids (RYA-394) ─────
# 4 Vesta / 7 Iris are SMALL BODIES: a bare Horizons id="4" resolves to NAIF 4 =
# MARS BARYCENTER, not asteroid 4 Vesta (NAIF 2000004) — the silent-wrong-value bug
# that mislabelled Vesta's reflected RV as Mars's for four tickets (RYA-372→373→380→391).
# Ganymede is a MAJOR BODY (Jovian satellite, NAIF 503), so a blanket "reject major
# bodies" guard would wrongly reject it — the guard asserts resolved-target == intended
# body via `match`, never by body class. Each entry carries its own id + id_type;
# callers pass a KEY, never a raw Horizons id. id_type strings are astroquery-version-
# dependent — confirmed against astroquery 0.4.11 (RYA-394).
REFLECTED_SOLAR_BODIES = {
    'vesta':    {'id': 'Vesta', 'id_type': 'smallbody', 'match': 'Vesta'},
    'iris':     {'id': 'Iris',  'id_type': 'smallbody', 'match': 'Iris'},
    'ganymede': {'id': '503',   'id_type': None,        'match': 'Ganymede'},  # Jovian satellite
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

import yaml as _yaml

# ── Per-star fundamental parameters — SINGLE SOURCE: config/stars.yaml (RYA-298) ──
# Moved out of this module so per-star data is data, not code (editable, diffable,
# self-evidently authoritative). Loaded once here and exposed as STAR_PARAMS for
# backward-compatible importers; get_star_params() is the resolver. The ξ-pin / GBS-
# scale rationale and per-record provenance now live in stars.yaml.
_STARS_YAML = Path(__file__).resolve().parent / 'stars.yaml'


def _load_stars_yaml() -> dict:
    if not _STARS_YAML.exists():
        raise FileNotFoundError(
            f"stars.yaml not found at {_STARS_YAML} — the single source of per-star "
            f"fundamental params (RYA-298). Refusing to fall back (no silent default).")
    with open(_STARS_YAML) as _fh:
        data = _yaml.safe_load(_fh)
    for _rec in data.values():          # xi_xcheck is a (lo, hi) tuple in code
        if isinstance(_rec.get('xi_xcheck'), list):
            _rec['xi_xcheck'] = tuple(_rec['xi_xcheck'])
    return data


STAR_PARAMS = _load_stars_yaml()


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
    # RYA-298: teff_K/logg/feh removed — fundamental params live in stars.yaml.
    'STAR_55CNC'       : ['vturb_kms', 'rv_kms', 'hd', 'hip'],
    'STAR_SOLAR'       : ['vturb_kms', 'rv_kms', 'program_id', 'citation'],
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
    assert 4000 < get_star_params('55cnc_a')['teff'] < 7000, "Teff out of range"

    print("Constants validated ✓")
