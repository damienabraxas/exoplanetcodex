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
FE_SCATTER_GATE    =  0.10   # dex — SUPERSEDED legacy-universal σ(Fe I) ceiling (RYA-277).
                             # Calibrated on the Sun and silently applied as universal; RYA-407
                             # verdict (c) showed it is mis-set for the current GES-EW 1D-LTE
                             # Fe I pool (raw-pool σ=0.1398; no cited mechanism reaches 0.10).
                             # The G-anchor solar gate now reads its ceiling from
                             # ACCEPTANCE_PROFILES['G'] (RYA-446); this value is retained only as
                             # the legacy default for not-yet-profiled types (F/K/M TODO, RYA-277).
                             # Do NOT use for the Sun.
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

# ── Per-spectral-type acceptance profiles (RYA-277) ──────────────────────────
# THE ARCHITECTURAL FACT: acceptance/quality criteria are a FUNCTION of stellar type,
# not universal. Every legacy threshold (σ-clip, COG/EW ceiling, scatter gate, A(Fe)
# window, n_lines floor, continuum strategy, NLTE applicability) was calibrated on the
# Sun and silently labelled "universal" — an implicit, unlabelled "G-star, 5772 K"
# assumption (the same failure family as the solar-only CODE PATHS in RYA-270/273/274,
# one level up: the solar assumption baked into the scientific CUTS themselves).
#
# This structure makes the type explicit. It is a ROUTING/parameterization layer, NOT a
# science change: the G profile = the current solar values, relabelled — routing solar
# through it reproduces the signed-off solar Fe result BIT-FOR-BIT (A(Fe I) NLTE 7.516).
# That regression (scripts/validate_fe_rya238.py --star sun) is the proof the refactor
# changed nothing for the anchor. Prototype patterns followed: CONTINUUM_PARAMS (RYA-274)
# and EW_FIT_PARAMS (RYA-273) are already per-star — those operative pool-construction
# params stay where they live; this profile CROSS-LINKS to them by key and owns the
# per-type ACCEPTANCE GATE thresholds (the validation criteria) as the single source.
#
# Anchors populated from real data (RYA-277 §3):
#   G = Sun (G2V, 5772 K)     — current solar gate values, regression-proven.
#   F = Procyon (F5IV, 6554 K) — Fe I scatter floor σ=0.222 FINAL from RYA-281 (the 4
#       RYA-458 ABUND_OUTLIER lines tested for bad-gf vs GES Heiter+2021: 0/4 confirmed,
#       EW blends excluded, no gf tuned → 0.222 is the genuine 1D-LTE F-star floor, NOT
#       contamination), nlte_available=False (Fe I NLTE grid runs out > 6500 K, RYA-319 →
#       1D-LTE-only → irreducible NLTE-driven floor), A(Fe) window/dFe/n_Fe2/vmic from the
#       RYA-273 F-anchor run. The "Procyon Fe σ FAILED" verdict was a SOLAR gate applied
#       to an F star; under this profile σ=0.222 is at-floor → PASS for an NLTE-unavailable
#       F dwarf. COG ceiling marked PROVISIONAL (RYA-273 COG_FLAG knee ≈ 120 mÅ) — the
#       operative pool cut stays PIPELINE['vmic_ew_ceiling_mA'] (the σ=0.222 anchor was
#       characterized at that cut; see RYA-277 §4 — F/K/M pool thresholds not reset here).
#   K, M = STUBS with the §2 physics notes + TODO; not populated (no target reached). The
#       accessor fails LOUD on a stub (no silent default), same anti-assumption discipline.
#
# σ-floor trend (RYA-407/281, the point of the per-type gate): the Fe I scatter floor
# RISES from cool G to warm F — G 0.1398 (NLTE-covered) → F 0.222 (NLTE-unavailable).
# A single universal number cannot be right for both; the gate must be a function of type.
ACCEPTANCE_PROFILES = {
    'G': {
        'teff_range'      : (5300, 6000),    # G regime (Sun ≈ 5772 K)
        'method'          : 'EW',            # atomic-line EW measurement is valid
        'nlte_available'  : True,            # Fe I NLTE grid covers G (Teff < 6500 K, RYA-319)
        # Acceptance GATE thresholds (single source; the validate script reads these):
        'a_fe_lo'         : SOLAR_ASPLUND2021['Fe'] + FE_GATE_LOWER,  # 7.41 (RYA-261)
        'a_fe_hi'         : SOLAR_ASPLUND2021['Fe'] + FE_GATE_UPPER,  # 7.51
        'dFe_max'         : FE_IONISATION_GATE,                       # 0.05 (RYA-261)
        # Fe I line-to-line scatter ceiling = the cited RYA-407 honest floor of the
        # current GES-EW 1D-LTE solar Fe I pool: σ(1D) = 0.1398 dex (the gate's NLTE σ
        # ≈ 1D σ; RYA-407 Step 0). Characterized POOL floor, derived BLIND to gate
        # pass/fail — NOT a round number, NOT chosen to make the Sun pass. RYA-407
        # verdict (c): no cited mechanism reaches the legacy 0.10.
        'fe1_scatter_max' : 0.1398,          # RYA-407 σ(1D) floor, GES-EW 1D-LTE Fe I pool
        'n_Fe2_min'       : 8,               # RYA-238 solar Fe II pool floor
        'vmic_lit'        : 1.00,            # km/s, solar microturbulence (pinned)
        # Pool-construction config — DOCUMENTED here, OPERATIVE value stays per-star:
        'cog_ceiling_mA'  : 100.0,           # == PIPELINE['vmic_ew_ceiling_mA'] (RYA-330 saturation cut)
        'continuum_key'   : 'solar',         # → CONTINUUM_PARAMS['solar'] (RYA-274)
        'ew_fit_key'      : 'solar',         # → EW_FIT_PARAMS['solar'] (RYA-273)
        'provenance'      : 'RYA-407 honest floor σ(1D)=0.1398 (GES-EW 1D-LTE solar Fe I pool); G=solar anchor (5772 K)',
    },
    'F': {
        'teff_range'      : (6000, 6800),    # F regime (Procyon ≈ 6554 K converged)
        'method'          : 'EW',
        'nlte_available'  : False,           # Fe I NLTE grid runs out > 6500 K (Procyon 6554 K, RYA-319) → 1D-LTE-only
        'a_fe_lo'         : 7.38, 'a_fe_hi'  : 7.54,   # RYA-273 F-anchor window
        'dFe_max'         : 0.08,            # RYA-273 F-anchor ionization gate
        # FINAL F-star Fe I scatter floor (RYA-281): bad-gf disproved (0/4 GES-confirmed),
        # blends excluded (priority 0, no gf tuned), post-NLTE σ=0.222. Sits just above the
        # ~0.18 NLTE-unavailable expectation → the genuine 1D-LTE F-star floor (Procyon
        # above the Fe NLTE grid ceiling), NOT removable contamination. This RISES above the
        # G floor 0.1398 by design — that rise IS why the gate must be per-type.
        'fe1_scatter_max' : 0.222,           # RYA-281 FINAL F-anchor floor (NLTE-unavailable)
        'n_Fe2_min'       : 3,               # RYA-273 F-anchor Fe II pool floor
        'vmic_lit'        : 1.66,            # km/s, Allende Prieto 2002 (Procyon)
        'cog_ceiling_mA'  : 120.0,           # PROVISIONAL F-anchor: RYA-273 COG_FLAG knee ≈ 120 mÅ
                                             # (operative pool cut stays PIPELINE['vmic_ew_ceiling_mA'];
                                             #  the σ=0.222 anchor was characterized at that cut — RYA-277 §4).
        'continuum_key'   : 'procyon',       # → CONTINUUM_PARAMS['procyon'] (RYA-274)
        'ew_fit_key'      : 'procyon',       # → EW_FIT_PARAMS['procyon'] (RYA-273)
        'provenance'      : 'RYA-281 F-star Fe σ floor 0.222 FINAL (bad-gf disproved, blends excluded, '
                            'NLTE-unavailable > 6500 K → genuine 1D-LTE floor); Procyon 6554 K anchor',
    },
    # ── STUBS — physics notes only, NOT populated (RYA-277 §2/§4) ──────────────
    'K': {
        'status'         : 'TODO',
        'method'         : 'EW',             # likely still EW, but molecular-band-aware
        'nlte_available' : None,             # unknown until a K target is reached
        'physics_note'   : ('Cooler than the Sun: sharper, deeper lines (good) BUT molecular '
                            'bands (CN, and toward late-K the onset of TiO) begin contaminating '
                            'regions that are clean in the Sun → a different BLEND regime. The '
                            'continuum-window strategy must avoid molecular bands. Thresholds '
                            'TODO — populate when a K-dwarf target is reached (RYA-277).'),
        'provenance'     : 'STUB — not populated (RYA-277 §2); no K target reached',
    },
    'M': {
        'status'         : 'TODO',
        'method'         : 'synthesis',      # EW likely INVALID as a method — pseudo-continuum forest
        'nlte_available' : None,
        'physics_note'   : ('M dwarfs (Tier-2 targets): severe molecular blanketing (TiO, VO, H2O) — '
                            'arguably NO true optical continuum, the pseudo-continuum is a forest. '
                            'Atomic-line EW measurement as done for the Sun may not be the valid '
                            'METHOD — full-spectrum synthesis on MARCS.GES atmospheres likely '
                            'mandatory. The M profile may specify METHOD, not just thresholds. '
                            'TODO (RYA-277).'),
        'provenance'     : 'STUB — not populated; method=synthesis flagged (RYA-277 §2)',
    },
}

# Star → spectral type, for acceptance-profile routing (RYA-446 seed; RYA-277 anchors).
STAR_SPECTRAL_TYPE = {'solar': 'G', 'procyon': 'F'}

# Keys a populated (non-stub) acceptance profile must carry — used by the loud accessor
# to distinguish a real anchor from a physics-note stub.
_ACCEPTANCE_GATE_KEYS = ('a_fe_lo', 'a_fe_hi', 'dFe_max', 'fe1_scatter_max',
                         'n_Fe2_min', 'vmic_lit')


def get_acceptance_profile(spectral_type):
    """Active acceptance profile for `spectral_type`, read from ACCEPTANCE_PROFILES
    (single source of truth). Fails LOUD on a missing or stub profile (RYA-277 K/M are
    physics-note stubs) — no silent default, same anti-assumption discipline as the
    RYA-270 linelist-routing log."""
    prof = ACCEPTANCE_PROFILES.get(spectral_type)
    if prof is None or not all(k in prof for k in _ACCEPTANCE_GATE_KEYS):
        _note = (f" ({prof['physics_note']})" if isinstance(prof, dict)
                 and 'physics_note' in prof else '')
        raise KeyError(
            f"No acceptance profile populated for spectral type {spectral_type!r}; it is "
            f"a deferred RYA-277 stub.{_note} Refusing to default silently — populate the "
            f"profile (RYA-277) before gating this type.")
    return prof


def acceptance_profile_for_star(star_id):
    """Resolve a star_id → its spectral type → the populated acceptance profile.
    Fails loud if the star has no type mapping or the type's profile is a stub."""
    s = star_id.strip().lower()
    spec = STAR_SPECTRAL_TYPE.get(s)
    if spec is None:
        spec = next((t for k, t in STAR_SPECTRAL_TYPE.items() if k in s or s in k), None)
    if spec is None:
        raise KeyError(
            f"No spectral-type mapping for star_id={star_id!r}; add it to "
            f"STAR_SPECTRAL_TYPE before acceptance-gating it (RYA-277).")
    return spec, get_acceptance_profile(spec)


def fe1_scatter_threshold(spectral_type):
    """Active Fe I line-to-line scatter ceiling for `spectral_type`, read from
    ACCEPTANCE_PROFILES (single source of truth). Fails loud on a missing/stub profile
    (RYA-277 K/M are stubs) — no silent default. Back-compat alias kept for RYA-446
    callers; new code can use get_acceptance_profile()['fe1_scatter_max']."""
    prof = ACCEPTANCE_PROFILES.get(spectral_type)
    if prof is None or 'fe1_scatter_max' not in prof:
        raise KeyError(
            f"No acceptance profile populated for spectral type {spectral_type!r}; the "
            f"RYA-277 K/M profiles are deferred TODO stubs. Refusing to default "
            f"silently — populate the profile (RYA-277) before gating this type.")
    return prof['fe1_scatter_max']

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

    # RYA-469: per-star NAMESPACED pipeline outputs + the frozen gold solar reference.
    # Working per-star products live under data/outputs/{star}/ (gitignored, regenerable);
    # the star is in the PATH so two stars physically cannot collide on a filename. The
    # gold-standard solar differential denominator is the versioned, write-once, immutable
    # data/reference/solar/solar_abundances_v{N}.csv (committed). See pipeline/data_namespace.py
    # for the accessors and docs/design/adr_data_namespacing_and_gold_reference.md for the law.
    'outputs_root'       : ROOT / 'data' / 'outputs',
    'solar_reference_dir': ROOT / 'data' / 'reference' / 'solar',
}

# ── Solar reference library (RYA-459, under the RYA-162 epic) ─────────────────
# The multi-wavelength solar reference atlases that reach lines HARPS-VIS (380-690
# nm) cannot — the anchors for solar N (N I red + NH 3360), the CNO arms, and the
# P/K/Co/Sc DATA-GAP elements. Each source carries an explicit provenance flag.
#
# CRITICAL honesty constraint (RYA-455 discipline): Hubble cannot observe the Sun
# directly, so there is NO direct solar UV spectrum. The UV composite is
# provenance='cited-composite' and must NEVER be presented as a measurement; the
# audit gate (pipeline/audit_solar_reference) fails loud if a UV/composite source
# is tagged 'measured'. Any downstream UV-derived abundance inherits the cited flag.
#
# provenance in {measured, cited-composite, model}.
SOLAR_REFERENCE_SPECTRA = {
    'kpno_flux_atlas': {
        'path': 'data/solar_reference/kpno_flux_atlas',
        'wavelength_coverage_nm': [296.0, 1300.0],
        'resolution': 'FTS high-res (uneven grid, lambda/dlambda ~ 4e5)',
        'flux_units': 'normalized residual flux + irradiance (uW/cm^2/nm)',
        'provenance': 'measured',
        'citation': 'Kurucz, Furenlid, Brault & Testerman 1984, NSO Atlas No. 1 '
                    '("Solar Flux Atlas from 296 to 1300 nm")',
    },
    'uv_composite': {
        'path': 'data/solar_reference/uv_composite',
        'wavelength_coverage_nm': [119.5, 2695.7],
        'resolution': 'low (dlambda ~ 20 A, R ~ 150-300) — flux composite, not a line atlas',
        'flux_units': 'FLAM (erg/s/cm^2/A, absolute)',
        'provenance': 'cited-composite',
        'citation': 'Colina, Bohlin & Castelli 1996 (AJ 112, 307) / Bohlin+2001 '
                    '(CALSPEC sun_reference_stis_002); UV below 4100A = Woods+1996',
    },
    'ir_atlases': {                     # RYA-390 (K-band CO arm) — already staged
        'path': 'data/solar_reference/ir_atlases',
        'wavelength_coverage_nm': [2289.3, 2349.5],   # CO dv=2 segment extracted
        'resolution': 'FTS (ACE-FTS 0.02 cm^-1)',
        'flux_units': 'normalized intensity',
        'provenance': 'measured',
        'citation': 'ACE-FTS (Hase+2010) / NSO photatl (Livingston & Wallace 1991) '
                    '/ Wallace telluric (NOAO/NSO) — RYA-390',
    },
    'ir_atlas_irtf': {                  # Tier-2 — documented + deferred (RYA-459)
        'path': 'data/solar_reference/ir_atlas',
        'wavelength_coverage_nm': [940.0, 2500.0],
        'resolution': 'R ~ 2000 (IRTF SpeX)',
        'flux_units': 'TBD on stage',
        'provenance': 'measured',
        'citation': 'Rayner et al. 2009 (IRTF) — DEFERRED, see README_SOURCE.md',
        'status': 'deferred',
    },
}

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
    # Sr II — RYA-421: closes the RYA-401 ION-MISMATCH GET-DATA gap (our EW set measured
    # Sr I 6617, but the NLTE grid exists for the Sr II 4077/4215 resonance doublet).
    # Delta-type grid (not an Amarsi-2020 b-factor), source-flag pattern (RYA-165):
    # PRIMARY = Mashonkina 2022 / INASAN (cited; demonstrated to [Fe/H]=+0.25; vendor when
    # machine-readable); WORKING grid = Bergemann 2012a via INSPECT (fetchable, agrees with
    # Mashonkina within 0.15 dex). NEAR-LTE at our metallicities (solar 4077 -0.004 /
    # 4215 -0.006; the correction grows only toward metal-poor). Grid [Fe/H] CEILING = +0.0
    # -> 55 Cnc (+0.32) and tau Boo (+0.25) are FLAGGED extrapolations (out-of-hull ->
    # nlte_flag NLTE_unavailable, loud, RYA-409) into the near-LTE regime; acceptable, never
    # silent (RYA-421 Step 3). 4215 carries an Fe I blend -> cool-star discipline below.
    'Sr': {'ion': 2, 'grid': 'Sr_Bergemann2012_INSPECT.csv',
           'ref': 'WORKING grid ON DISK + APPLIED = Bergemann et al. 2012a (A&A 546 A90; '
                  'arXiv 1207.2451) via INSPECT (the applied delta grid). PRIMARY = Mashonkina '
                  'et al. 2022 / INASAN (Mashonkina+2023; DB cite Mashonkina, Sitnova, Pakhomov '
                  '2016 AstL 42,606) -- NOT YET VENDORED: RYA-433 found the INASAN nLTE.cgi '
                  'endpoint WAF-BLOCKED (403) to all programmatic clients -> MANUAL PULL OWED '
                  '(see data/nlte_grids/Sr_Mashonkina2022_INASAN.PENDING.md). Both are metal-poor '
                  'MARCS grids (INASAN [Fe/H] -5..-2; Bergemann -3..0); published agreement '
                  '~0.15 dex (our cross-check pending the pull). Sr II 4077/4215 doublet (RYA-421/433)',
           'flag': 'NLTE_Bergemann2012_INSPECT_1D'},
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
    # RYA-462: K I resonance doublet 7665/7699 (ground-state 4s 2S -> 4p 2P*), the last
    # PRESENT-but-UNWIRED grid the RYA-460 solar polish surfaced. The grid + the full
    # PySME machinery already existed (pipeline.pysme_nlte, RYA-402); the CSV departure
    # grid K_Amarsi2020_PySME.csv was vendored but K was never added to this registry, so
    # the default run left K at 1D-LTE. K I resonance NLTE is LARGE and NEGATIVE (the line
    # is much stronger in NLTE): solar 7665 -0.27 / 7699 -0.31, applied as-is from the
    # vendored grid (validate-don't-tune — the same interpolation path as Ca/Ti/Cr/Na).
    # 7665 sits in the telluric O2 A-band; the clean diagnostic is 7699. Solar leg:
    # Kitt Peak K I 7699 = 5.411 (1D-LTE), -0.312 -> 5.099 (+0.03 vs Asplund 5.07).
    'K': {'ion': 1, 'grid': 'K_Amarsi2020_PySME.csv',
          'ref': 'Amarsi et al. 2020 (A&A 642, A62) departure grid via PySME; K I resonance '
                 '7665/7699; solar anchor Reggiani et al. 2019 (A&A 627, A177) / Andrievsky '
                 'et al. 2006 (severe negative K I resonance NLTE); RYA-402/462 PySME-derived deltas',
          'flag': 'NLTE_Amarsi2020_PySME_1D'},
    # N I red-optical triplet 7468/8216/8683 (RYA-369 strategy; wired RYA-526). High-excitation
    # (~10.3 eV), atomic-primary, DECOUPLED from the C/O CNO synthesis leg. Derived on Sirius via
    # PySME from the Amarsi-2020 GALAH N grid (nlte_N_scatt_pysme.grd, Zenodo 3982506) over the
    # standard 11-node box; the SOLAR node reproduces the RYA-369 load-test EXACTLY (7468 -0.0115,
    # 8216 -0.0145, 8683 -0.0154 — SMALL negative, N I red is NEAR-LTE at the Sun, sign-discipline
    # OK RYA-339). N I red is a WARM-star indicator: 4 weak-line COG rails at cool (Teff<=5172)
    # metal-rich nodes (7468 x3, 8216 x1) were EXCLUDED from the grid, so those fall out-of-hull ->
    # loud-flagged (RYA-409), never a spurious value. Same PySME-derived path + flag as Na/Mg/Si/Al/S/K.
    'N': {'ion': 1, 'grid': 'N_Amarsi2020_PySME.csv',
          'ref': 'Amarsi et al. 2020 (A&A 642, A62) N I departure grid via PySME (Zenodo 3982506); '
                 'N I red triplet 7468/8216/8683 high-excitation atomic; RYA-369 strategy / RYA-526 '
                 'PySME-derived deltas (solar reproduces the RYA-369 load-test)',
          'flag': 'NLTE_Amarsi2020_PySME_1D'},
}

# ── Sr II resonance-doublet line discipline (RYA-421 Step 4) ──────────────────
# Sr II 4215 is blended by an Fe I line (Bergemann 2012a) — the blend grows toward
# cool stars. Per-line role by stellar temperature, enforced at measurement/loader
# time so a known-blended line is NEVER silently averaged into a cool-star Sr abundance:
#   - Sun / warm (tau Boo): BOTH 4077 + 4215 usable.
#   - COOL (55 Cnc, G8V, Teff<~5300): 4077 PRIMARY; 4215 cross-check-only or DROPPED.
SR2_LINES = {
    'primary':    [4077.709],            # clean resonance line, all targets
    'crosscheck': [4215.519],            # Fe I-blended -> cool-star caution
}
SR2_COOLSTAR_TEFF_K = 5300.0             # below this, 4215 is cross-check-only/dropped
SR2_4215_BLEND = 'Fe I (Bergemann 2012a); avoid in cool stars'

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
    # RYA-459: the Tier-1 solar reference sources must stay registered (audit reads these).
    'SOLAR_REFERENCE_SPECTRA': ['kpno_flux_atlas', 'uv_composite'],
}

_DICTS = {
    'PHYSICS'          : PHYSICS,
    'ASTRO'            : ASTRO,
    'SOLAR_ASPLUND2021': SOLAR_ASPLUND2021,
    'STAR_55CNC'       : STAR_55CNC,
    'STAR_SOLAR'       : STAR_SOLAR,
    'PIPELINE'         : PIPELINE,
    'PATHS'            : PATHS,
    'SOLAR_REFERENCE_SPECTRA': SOLAR_REFERENCE_SPECTRA,
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
