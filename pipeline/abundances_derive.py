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

import os
import sys
import warnings
from pipeline import _runtime as _rt   # RYA-514: force-fork + single-thread BLAS pins; import BEFORE numpy/scipy so the pins land before the BLAS backend spins up. Supersedes the RYA-506 inline fork block.
import numpy as np
from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz (RYA-313)
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
    FE_1D3D_SOLAR_OFFSET, FE_ABS_DIAG_HALFWIDTH, FE_REW_SLOPE_GATE,
    CORRECTIONS_3D,                      # RYA-553 tabulated solar Fe 1D→3D correction
    FE_IONIZATION_SYNTH_ARBITER, FE_EW_SYNTH_SPREAD_BAND,
    SYNTH_CHI2_GATE,
    assert_abundance_on_scale,
    ACCEPTANCE_PROFILES, STAR_SPECTRAL_TYPE, fe1_scatter_threshold,
)
from pipeline.species import species_key, species_note
from pipeline.gf_resolver import apply_to_synth_array, apply_to_regions  # RYA-353 single-source gf
from pipeline import data_namespace as ns                               # RYA-469 per-star namespacing


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
    """'Fe', 'I' → 'Fe 1';  'Fe', 'II' → 'Fe 2'.

    RYA-345: routed through the canonical normalizer (pipeline.species), which
    fixes the old `1 if ion=='I' else 2` collapse (every non-neutral ion mapped
    to 2) and raises loudly on unparseable element/ion instead of mis-keying.
    """
    return species_note(element, ion)


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
    # RYA-353 single-source gf: this is the live EW→abundance gf route (store #3).
    # Overwrite each region's loggf from canonical so EW and synth read one value.
    regions = apply_to_regions(regions)

    # Build lookup keyed by the CANONICAL species key (RYA-345) — (Z, ion) —
    # so the GES region 'note' ('Fe 2') and our (element, ion) strings ('Fe','II')
    # resolve to the same key regardless of encoding. A 0-match per species used
    # to be silent (RYA-339/343); see the guard below.
    lookup: dict[tuple, list] = {}
    for _, row in ew_df.iterrows():
        key    = species_key(str(row['element']), str(row['ion']))
        wav    = float(row['wavelength_air_A'])
        ew_mA  = float(row['ew_mA'])
        err_mA = float(row['ew_err_mA']) if 'ew_err_mA' in row.index else 0.0
        lookup.setdefault(key, []).append((wav, ew_mA, err_mA))

    # Match regions → our EWs (match on wave_A which is in Å)
    matched_idx, ew_vals, err_vals = [], [], []
    matched_species: set = set()
    for i, region in enumerate(regions):
        try:
            key = species_key(str(region['note']))
        except ValueError:
            continue                      # molecular / non-atomic region note
        if key not in lookup:
            continue
        candidates = lookup[key]
        diffs = [abs(wav - float(region['wave_A'])) for wav, _, _ in candidates]
        min_diff = min(diffs)
        if min_diff > wave_tol:
            continue
        best = candidates[int(np.argmin(diffs))]
        matched_idx.append(i)
        ew_vals.append(best[1])
        err_vals.append(best[2])
        matched_species.add(key)

    if not matched_idx:
        raise RuntimeError(
            f"No EW measurements matched iSpec line regions "
            f"(±{wave_tol} Å). Check element/ion names and wavelengths."
        )

    # RYA-345 regression guard: a species that HAS measured lines but matched
    # 0 regions is a silent-encoding failure (the exact 339/343 footgun). Warn
    # loud per species so a 0-match never hides as "no such line".
    for key, cands in lookup.items():
        if key not in matched_species:
            warnings.warn(
                f"RYA-345 0-match guard: species {key} has {len(cands)} measured "
                f"line(s) but matched 0 iSpec regions (±{wave_tol} Å). Check the "
                f"region file coverage / wavelength tolerance — NOT an encoding "
                f"mismatch (all sides routed through species_key).",
                RuntimeWarning, stacklevel=2)

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


# RYA-305: proper per-line Fe II theoretical EW (SPECTRUM) for the EW triage.
# The Fe I-calibrated linear COG mis-predicts Fe II by ~2 dex, which is why the
# old code fail-opened Fe II past the sanity check. calculate_theoretical_ew_and_depth
# gives the real EW each Fe II line would have at the expected composition — the
# physical discriminant for clean / blend-contaminated. Cached per (line set, params)
# since it is intrinsic to the lines and does not change across convergence iterations.
_fe2_theo_cache: dict = {}
_fe2_quarantine_state: dict = {'written': False}
_fe1_quarantine_state: dict = {'written': False}   # RYA-309: Fe I triage parity

# RYA-506: the theoretical-EW step is the universal quarantine gate. Its silent-failure
# mode was iSpec returning a zero-initialised array when its child synthesis process died
# (root cause: the macOS 'spawn' default — fixed at module import above). This guard is the
# loud backstop: an all-zero BATCH is the clean signature of a failed synthesis (physically
# impossible for a real Fe pool), so we RAISE and do NOT cache rather than let a failed
# synthesis silently quarantine a real pool. It does NOT loosen the theo<5 mÅ filter: a
# single weak line legitimately near zero still sits inside a batch with non-zero
# neighbours, so `np.any(theo > 0)` separates "synthesis failed" from "one line quarantined".
_TMP_DIR = '/tmp/ispec_codex'


def ensure_ispec_tmp_dir(tmp_dir: str = _TMP_DIR) -> str:
    """Create the iSpec/MOOG scratch dir and return it. Never silently skipped.

    RYA-506 added the makedirs guard to `_fe2_theoretical_ew` only, but three other
    call sites hand the SAME path to `ispec.determine_abundances` without creating it
    (two here, one in curate_nonfe_pools) — MOOG then dies on
    `FileNotFoundError: /tmp/ispec_codex/tmpXXXX` the first time the dir is absent
    (e.g. a fresh boot, or /tmp reaped). Every caller routes through this helper so the
    precondition holds once, in one place, instead of three literals drifting apart.
    """
    os.makedirs(tmp_dir, exist_ok=True)
    return tmp_dir


def _fe2_theoretical_ew(fe2_linemasks, stellar_params, atmosphere) -> np.ndarray:
    """Theoretical Fe II-only EW (mÅ) per line at expected A(Fe)=solar+[Fe/H], via SPECTRUM.

    RYA-506: RAISES (does not cache) if the synthesis returns an all-zero batch — the
    signature of a silently-failed iSpec child — instead of quarantining a real pool on it."""
    teff = float(stellar_params['teff_K']); logg = float(stellar_params['logg'])
    feh  = float(stellar_params['feh']);    vturb = float(stellar_params['vturb_kms'])
    key = (tuple(np.round(np.asarray(fe2_linemasks['wave_A'], dtype=float), 3)),
           round(teff), round(logg, 2), round(feh, 2))
    if key in _fe2_theo_cache:
        return _fe2_theo_cache[key]
    _, isotopes, _ = _load_synth_resources()
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    os.makedirs(_TMP_DIR, exist_ok=True)   # RYA-506: never let a missing tmp dir crash write_atmosphere
    n_lines = int(len(fe2_linemasks))
    out = ispec.calculate_theoretical_ew_and_depth(
        atmosphere, teff, logg, feh, 0.0, fe2_linemasks.copy(), isotopes, solar_abund,
        microturbulence_vel=vturb, verbose=0, tmp_dir=_TMP_DIR)
    theo = np.asarray(out['theoretical_ew'], dtype=float)
    if not np.any(theo > 0):               # all-zero batch → failed synthesis, not physics
        raise RuntimeError(
            f"[RYA-506] theoretical-EW synthesis returned an all-zero batch for {n_lines} "
            f"Fe lines at teff={teff} logg={logg} feh={feh}. An all-zero theoretical EW is "
            f"physically impossible for a real Fe pool — the iSpec child synthesis died "
            f"without enqueueing a result (silent zero-init; root cause = a non-'fork' "
            f"multiprocessing start method — see the module-level fix + docs/OPEN_QUESTIONS.md). "
            f"REFUSING to cache it or quarantine a real pool on a failed synthesis.")
    _fe2_theo_cache[key] = theo
    return theo


# RYA-102/103 fold-in (RYA-371 Phase 1): the special weak/blended lines that are
# their element's ONLY solar indicator. lines_fit force-MEASURES these; the EW→A(X)
# pre-filter must FORCE-INCLUDE them or they vanish into a DATA-GAP. Eu II 6645
# (HFS-total EW, RYA-102) = a real measurement; Li I 6707 (CN-blend, RYA-103) = an
# UPPER LIMIT — both already carry their flags in 'notes' for the per-element verdict.
EW_FORCE_INCLUDE = (('Eu', 'II', 6645.127), ('Li', 'I', 6707.840))   # RYA-102 / RYA-103


def _ew_prefilter(ew_df: pd.DataFrame, ew_min: float, ew_max: float):
    """Pre-filter the EW table for abundance derivation: drop blend-flagged and
    out-of-range lines, EXCEPT the RYA-102/103 force-include set (kept with flags so
    Eu/Li are not silently lost). Returns (ew_clean, n_force_included). Additive —
    non-force-include rows are filtered exactly as before (Fe/metal pools unchanged)."""
    df = ew_df.copy()
    if not len(df):
        return df, 0
    wcol = 'wavelength_air_A' if 'wavelength_air_A' in df.columns else 'wavelength'
    force = pd.Series(False, index=df.index)
    for el, ion, w in EW_FORCE_INCLUDE:
        force |= ((df['element'].astype(str) == el) & (df['ion'].astype(str) == ion)
                  & (df[wcol].astype(float).sub(w).abs() < 0.1))
    keep = pd.Series(True, index=df.index)
    if 'blend_flag' in df.columns:
        keep &= (df['blend_flag'] == False)
    keep &= (df['ew_mA'] >= ew_min) & (df['ew_mA'] <= ew_max)
    n_forced = int((force & ~keep).sum())
    return df[keep | force], n_forced


def _ew_to_abundance(ew_df: pd.DataFrame,
                     stellar_params: dict,
                     atmosphere: np.ndarray,
                     model_grid: str = 'ATLAS9.Castelli',
                     code: str = RADIATIVE_TRANSFER_CODE,
                     star_id: str = 'solar') -> tuple:
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
    ew_clean, n_forced = _ew_prefilter(ew_df, ew_min, ew_max)

    n_total  = len(ew_df)
    n_blends = int((ew_df.get('blend_flag', pd.Series(False, index=ew_df.index)) == True).sum())
    print(f"  EW pre-filter: {n_total} total → {n_blends} blend-flagged, "
          f"cuts [{ew_min:.0f},{ew_max:.0f}] mÅ → {len(ew_clean)} kept "
          f"(incl. {n_forced} force-included RYA-102/103 Eu II 6645 / Li I 6707)")

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

    # ── Fe II triage (RYA-305) — the has_theo==False fail-open DIES ─────────────
    # Old behaviour fail-opened every Fe II line without a GES theoretical_ew, so
    # physically-impossible blends (e.g. 5376.5 Å, loggf −7.9, ~0 mÅ real Fe II
    # under ~80 mÅ of blend → A 10.3) were averaged into the Fe II mean (RYA-290).
    # The Fe I-calibrated linear COG can't vet Fe II (it mis-predicts by ~2 dex),
    # so judge each Fe II line against its PROPER theoretical EW (SPECTRUM, at the
    # expected composition). Keep only lines whose observed EW is consistent with
    # that theoretical EW (clean → EW abundance). Blend-contaminated lines are
    # EXCLUDED from the EW pool (no silent fail-open): if the line is intrinsically
    # measurable (theo >= floor) the blend-aware synthesis recovers it; if it is
    # negligible (theo < floor) it is dropped/quarantined. Anchor = expected A(Fe)
    # (solar+[Fe/H]), keeping the Fe I−Fe II balance an independent check (RYA-305
    # DECISION). The Fe II ionization-balance ARBITER is synthesis (DECISION 2).
    fe2_mask = np.array([str(r['note']).strip() == 'Fe 2' for r in linemasks])
    if fe2_mask.any():
        fe2_idx  = np.where(fe2_mask)[0]
        theo_fe2 = _fe2_theoretical_ew(linemasks[fe2_idx], stellar_params, atmosphere)
        obs_fe2  = np.array([float(linemasks['ew'][i]) for i in fe2_idx])
        logr_fe2 = np.log10(np.maximum(obs_fe2, 1e-6) / np.maximum(theo_fe2, 1e-6))
        clean_fe2 = np.abs(logr_fe2) < float(PIPELINE['fe2_ew_sanity_dex'])
        for j, i in enumerate(fe2_idx):
            good[i] = bool(clean_fe2[j])   # Fe II verdict is the triage, not the linear COG
        _floor = float(PIPELINE['fe2_synth_ew_floor_mA'])
        n_clean = int(clean_fe2.sum())
        n_recov = int(((~clean_fe2) & (theo_fe2 >= _floor)).sum())
        n_drop  = int(((~clean_fe2) & (theo_fe2 <  _floor)).sum())
        print(f"  Fe II triage (RYA-305): {n_clean} clean (EW kept), "
              f"{n_recov} blended→synthesis-recover (EW-excluded), "
              f"{n_drop} dropped/quarantined (theo<{_floor:.0f} mÅ)")
        # Quarantine record (never delete; provenance preserved — RYA-305 DECISION).
        if not _fe2_quarantine_state['written'] and (~clean_fe2).any():
            try:
                q = pd.DataFrame({
                    'wave_A': np.round(obs_fe2 * 0 + np.asarray(
                        [float(linemasks['wave_A'][i]) for i in fe2_idx]), 3),
                    'obs_ew_mA': np.round(obs_fe2, 2),
                    'theo_ew_mA': np.round(theo_fe2, 3),
                    'log_obs_over_theo': np.round(logr_fe2, 2),
                    'verdict': np.where(clean_fe2, 'clean',
                               np.where(theo_fe2 >= _floor, 'recover_synthesis', 'drop')),
                    'reason': np.where(clean_fe2, 'obs~theo (kept in EW)',
                               np.where(theo_fe2 >= _floor,
                                        'real line, EW blend-contaminated → synthesis',
                                        'theo<floor: negligible Fe II under blend → drop')),
                })
                qpath = ns.diagnostics_dir(star_id) / 'fe2_triage_quarantine.csv'  # RYA-469
                q[q.verdict != 'clean'].to_csv(qpath, index=False)
                _fe2_quarantine_state['written'] = True
                print(f"  [Fe II triage] quarantine → {qpath.name}")
            except Exception as _qe:
                print(f"  [Fe II triage] quarantine write skipped: {_qe}")

    # ── Fe I triage (RYA-309 — extend the RYA-305 vetting to Fe I) ─────────────
    # Procyon's Fe I scatter (0.68 dex) and the implausible NLTE correction
    # (RYA-317) are driven by blend-inflated Fe I lines that the linear-COG `good`
    # filter (1.5 dex) lets through — the linear COG under-predicts EW 3-10×, so a
    # heavily blended line can still land inside 1.5 dex. Vet each Fe I line
    # against its PROPER Fe-I-only theoretical EW (SPECTRUM, expected composition)
    # — identical clean/recover/drop logic to the Fe II triage (_fe2_theoretical_ew
    # is ion-general). Blend-contaminated-but-real lines are EXCLUDED from the EW
    # pool and recovered by synthesis (RYA-305 DECISION 2: synthesis is the Fe
    # arbiter); negligible lines (theo<floor) are dropped/quarantined. Same 0.5 dex
    # sanity threshold + floor as Fe II. Clean solar GES Fe I lines pass (obs~theo).
    fe1_mask = np.array([str(r['note']).strip() == 'Fe 1' for r in linemasks])
    if fe1_mask.any():
        fe1_idx  = np.where(fe1_mask)[0]
        theo_fe1 = _fe2_theoretical_ew(linemasks[fe1_idx], stellar_params, atmosphere)
        obs_fe1  = np.array([float(linemasks['ew'][i]) for i in fe1_idx])
        logr_fe1 = np.log10(np.maximum(obs_fe1, 1e-6) / np.maximum(theo_fe1, 1e-6))
        clean_fe1 = np.abs(logr_fe1) < float(PIPELINE['fe2_ew_sanity_dex'])
        for j, i in enumerate(fe1_idx):
            good[i] = bool(clean_fe1[j])
        _floor = float(PIPELINE['fe2_synth_ew_floor_mA'])
        n_clean = int(clean_fe1.sum())
        n_recov = int(((~clean_fe1) & (theo_fe1 >= _floor)).sum())
        n_drop  = int(((~clean_fe1) & (theo_fe1 <  _floor)).sum())
        print(f"  Fe I triage (RYA-309/305): {n_clean} clean (EW kept), "
              f"{n_recov} blended→synthesis-recover (EW-excluded), "
              f"{n_drop} dropped/quarantined (theo<{_floor:.0f} mÅ)")
        if not _fe1_quarantine_state['written'] and (~clean_fe1).any():
            try:
                q = pd.DataFrame({
                    'wave_A': np.round(np.asarray(
                        [float(linemasks['wave_A'][i]) for i in fe1_idx]), 3),
                    'obs_ew_mA': np.round(obs_fe1, 2),
                    'theo_ew_mA': np.round(theo_fe1, 3),
                    'log_obs_over_theo': np.round(logr_fe1, 2),
                    'verdict': np.where(clean_fe1, 'clean',
                               np.where(theo_fe1 >= _floor, 'recover_synthesis', 'drop')),
                    'reason': np.where(clean_fe1, 'obs~theo (kept in EW)',
                               np.where(theo_fe1 >= _floor,
                                        'real line, EW blend-contaminated → synthesis',
                                        'theo<floor: negligible Fe I under blend → drop')),
                })
                qpath = ns.diagnostics_dir(star_id) / 'fe1_triage_quarantine.csv'  # RYA-469
                q[q.verdict != 'clean'].to_csv(qpath, index=False)
                _fe1_quarantine_state['written'] = True
                print(f"  [Fe I triage] quarantine → {qpath.name}")
            except Exception as _qe:
                print(f"  [Fe I triage] quarantine write skipped: {_qe}")

    n_before = len(linemasks)
    linemasks = linemasks[good]
    n_after  = len(linemasks)
    print(f"  EW sanity filter: {n_before} → {n_after} lines "
          f"({n_before - n_after} rejected as blends/misidentifications; "
          f"Fe II via RYA-305 triage)")

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
        tmp_dir=ensure_ispec_tmp_dir(),
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
        # RYA-353 single-source gf: overwrite GES loggf from canonical (the shared
        # ispec/ tsv is read-only, so we rescale after read — branching preserved).
        _synth_cache['linelist']      = apply_to_synth_array(_synth_cache['linelist'])
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
    ew_nm = float(_trapezoid(np.maximum(0.0, 1.0 - flux), waveobs_nm))
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

# Strong, fairly isolated Fe I lines (air, rest Å) for a robust RV centroid.
# Used to bring the observed spectrum to REST before flux-space synthesis fitting
# (RYA-309): the synthesis is rest-frame, so an uncorrected stellar RV misaligns
# the line cores, inflating χ² and biasing the broadening/abundance fit. Solar is
# already ~rest (BERV applied per-exposure) so this is a near-no-op there; Procyon
# carries an uncorrected ≈ -2.3 km/s shift.
_RV_FE1_REF = np.array([
    5497.516, 5501.465, 5569.618, 5615.644, 5624.542, 5934.655, 6024.058,
    6065.482, 6137.692, 6191.558, 6230.723, 6252.555, 6411.649, 6430.846,
])
_C_KMS = 299792.458


def _measure_rv_kms(wave_nm: np.ndarray, flux: np.ndarray) -> float:
    """Median velocity shift (km/s) of the observed spectrum vs rest, from
    parabolic centroids of strong Fe I reference lines. Positive = redshifted."""
    wl = wave_nm * 10.0
    vs = []
    for r0 in _RV_FE1_REF:
        m = (wl > r0 - 0.25) & (wl < r0 + 0.25)
        if m.sum() < 5:
            continue
        w, f = wl[m], flux[m]
        i = int(np.argmin(f))
        if f[i] > 0.85 or i < 1 or i > len(f) - 2:   # need a real, resolved core
            continue
        x, y = w[i - 1:i + 2], f[i - 1:i + 2]
        denom = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
        if denom == 0:
            continue
        A = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / denom
        B = (x[2]**2 * (y[0] - y[1]) + x[1]**2 * (y[2] - y[0]) + x[0]**2 * (y[1] - y[2])) / denom
        lam = -B / (2 * A) if A != 0 else w[i]
        vs.append(_C_KMS * (lam - r0) / r0)
    return float(np.median(vs)) if vs else 0.0


def _load_observed_spectrum(star_id: str) -> tuple:
    """
    Load the observed normalized spectrum for flux-space fitting (RYA-287).

    Returns (obs_wave_nm, obs_flux) — continuum-normalized (≈1.0 at continuum),
    shifted to REST frame via a measured Fe I RV (RYA-309). The normalized CSV
    has no per-pixel error column, so χ² weighting uses a locally-estimated
    continuum-scatter σ per window (see _fit_synth_flux).
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

    # Bring to REST: the synthesis is rest-frame, so an uncorrected stellar RV
    # misaligns line cores and biases the flux-space fit (RYA-309). Measure from
    # Fe I centroids and shift; no-op (<0.5 km/s) for already-rest spectra (solar).
    rv = _measure_rv_kms(wave_nm, flux)
    if abs(rv) > 0.5:
        wave_nm = wave_nm / (1.0 + rv / _C_KMS)
        print(f"  [synth-v2] RV-corrected observed spectrum to rest: {rv:+.2f} km/s "
              f"(Fe I centroids)")
    else:
        print(f"  [synth-v2] RV {rv:+.2f} km/s — within 0.5 km/s, no shift applied")
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
    # RYA-342: the EW saturation ceiling does NOT gate the flux-space path; acceptance
    # is on fit quality (SYNTH_CHI2_GATE) in the aggregation below.

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
          f"(fixed); line acceptance = fit quality χ²ᵣ < {SYNTH_CHI2_GATE} (RYA-342, "
          f"NOT the EW ceiling)")

    # RYA-338: optional element focus for the decisive test (e.g. the solar
    # Fe I−Fe II ionization read under clean flux-space synthesis). Opt-in via env
    # var (comma-separated element symbols); default = all matched lines, no
    # behaviour change. Restricts only which lines get the per-line flux fit, so the
    # full 27-element run and this Fe-focused run share identical machinery/output.
    import os as _os
    _focus = _os.environ.get('SYNTH_V2_ONLY_ELEMENTS', '').strip()
    if _focus:
        _want = {e.strip() for e in _focus.split(',') if e.strip()}
        _keep = np.array([str(last_linemasks['note'][i]).split()[0] in _want
                          for i in range(len(last_linemasks))])
        last_linemasks = last_linemasks[_keep]
        print(f"  [synth-v2] RYA-338 element focus = {sorted(_want)} → "
              f"{int(_keep.sum())} lines (of {len(_keep)} matched)")

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

        # RYA-342: the flux-space path does NOT gate on the EW saturation ceiling
        # (vmic_ew_ceiling_mA) — that is an EW-path concept. Flux-space synthesis is
        # built (RYA-338 FLAG 2, wing-wide windows) to measure exactly the strong/
        # blended lines EW saturation kills; the gate is FIT QUALITY (reduced χ²),
        # applied in the aggregation below. Every line is fit; acceptance is on merit.
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

    # RYA-342: empty-set guard — accessing ['status'] on an empty DataFrame raises
    # KeyError (root cause of the ceiling-150 crash: an upstream filter left 0 lines).
    # Fail soft + loud rather than KeyError; never a bare except.
    if per_line_v2_df.empty or 'status' not in per_line_v2_df.columns:
        print("  [synth-v2] WARNING: 0 synth lines produced (check line set / element "
              "focus) — returning empty tables, no aggregation.")
        empty = pd.DataFrame(columns=['element', 'ion', 'A_X', 'A_X_abs', 'A_X_std',
                                      'n_lines', 'med_red_chi2', 'XH', 'XFe',
                                      'scale', 'nlte_applied'])
        per_line_v2_df.to_csv(Path(str(out_dir)) / f'{star_id}_per_line_synth_v2.csv',
                              index=False)
        return empty, per_line_v2_df

    # ── Line acceptance = FIT QUALITY (RYA-342), not the EW saturation ceiling ──
    # Accept status=='ok' AND reduced-χ² < SYNTH_CHI2_GATE. Lines fit but poorly
    # (blend/model error, χ²ᵣ ≫ gate) are rejected ON MERIT and LOGGED with their χ²ᵣ —
    # never silently dropped. (edge_pinned/failed still excluded by status.)
    _ok    = per_line_v2_df[per_line_v2_df['status'] == 'ok'].copy()
    ok_df  = _ok[_ok['red_chi2'] < SYNTH_CHI2_GATE].copy()
    _rej   = _ok[_ok['red_chi2'] >= SYNTH_CHI2_GATE]
    for _, r in _rej.sort_values('red_chi2', ascending=False).iterrows():
        print(f"  [synth-v2] REJECT(fit) {r['element']} {r['ion']} "
              f"{r['wavelength_air_A']:.3f} Å — χ²ᵣ={r['red_chi2']:.1f} ≥ "
              f"{SYNTH_CHI2_GATE} (a_synth={r['a_synth']:.3f} not used)")
    if len(_rej):
        print(f"  [synth-v2] fit-quality gate: {len(ok_df)} accepted / "
              f"{len(_rej)} rejected (χ²ᵣ ≥ {SYNTH_CHI2_GATE}) of {len(_ok)} fitted")
    # Regression guard (RYA-342): an EW-domain cut must NEVER drop a good-fit synth
    # line. 'saturated' True with a good χ²ᵣ means an EW ceiling has re-gated the synth
    # path — fail loud (the whole point of this fix).
    _bad = per_line_v2_df[(per_line_v2_df['saturated'] == True) &
                          (per_line_v2_df['red_chi2'] < SYNTH_CHI2_GATE)]
    if len(_bad):
        raise AssertionError(
            f"RYA-342: {len(_bad)} synth line(s) with good χ²ᵣ were flagged 'saturated' "
            f"by an EW cut — flux-space synthesis must gate on fit quality only "
            f"(lines: {_bad['wavelength_air_A'].round(3).tolist()})")
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

    # RYA-682: never write a per-line table in which NO row is usable. The RYA-342
    # guard above only catches an EMPTY frame; the failure that actually happened
    # produces a full-length frame where every row is status='failed' (iSpec's
    # atom-code call breaks on numpy >= 2.3, so every element is skipped). That
    # file parses, looks like a successful run, and silently reduces the two-engine
    # floor to one engine downstream. Fail here, where the cause is still visible.
    if n_ok == 0:
        import numpy as _np
        raise SystemExit(
            f"RYA-682 REFUSING TO WRITE AN UNUSABLE ARTIFACT: {star_id} synthesis-v2 "
            f"produced {len(per_line_rows)} lines and ZERO usable ones "
            f"({n_skip} no-atom-code, {n_fail} failed, {n_edge} edge-pinned).\n"
            f"If no-atom-code dominates, this is the numpy break: ispec/abundances.py:132 "
            f"assigns a size-1 array into a scalar slot, which numpy >= 2.3 rejects "
            f"(running numpy {_np.__version__}). docs/SCIENCE_STANDARDS.md (RYA-517) "
            f"ratifies python 3.12 + numpy 2.2.x as the reference stack — regenerate "
            f"there. Writing this file would hand the two-engine driver an Engine-B "
            f"table with nothing in it.")

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


def _fe1_ceiling_pool_mask(notes, ew_mA, ceiling_mA: float) -> np.ndarray:
    """RYA-279: THE single ceiling-correct Fe I pool. A Fe I line enters the pool iff
    its EW ≤ ceiling. This one predicate feeds BOTH the gate scatter statistic
    (A_X_std) and the per-line CSV, so the two can never diverge again — the original
    bug computed σ over the full Fe I pool while the CSV had the ceiling applied
    post-hoc, so the gate reflected lines the artifact excluded. Fe II carries no
    ceiling (fewer, weaker lines, not in the damping regime at the EW seen in
    practice). Returns a boolean mask over the input lines."""
    notes = np.asarray([str(n) for n in notes])
    ew = np.asarray(ew_mA, dtype=float)
    return (notes == 'Fe 1') & (ew <= float(ceiling_mA))


# ── Iterative parameter convergence ──────────────────────────────────────────

def _iterative_parameter_convergence(ew_df: pd.DataFrame,
                                      stellar_params_init: dict,
                                      model_grid: str = 'ATLAS9.Castelli',
                                      max_iter: int = 10,
                                      skip_convergence: bool = False,
                                      solve_params=None,
                                      star_id: str = 'solar') -> tuple:
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
    # spectroscopic params (backward-compatible with pre-292 behaviour). RYA-325:
    # solve_xi is derived PURELY from the pin/solve policy now — the legacy
    # vmic_fixed flag is retired (a pinned ξ / xi_override is applied in run()
    # before this call, which simply drops 'xi' from solve). A pinned param's
    # equilibrium criterion is treated as satisfied so convergence is judged only
    # on the solved params; the reduced-EW slope is still computed + reported at a
    # pinned ξ as the guardrail (flat = consistent, sloped = a finding).
    solve = set(solve_params) if solve_params is not None else {'teff', 'logg', 'xi'}
    solve_teff = 'teff' in solve
    solve_logg = 'logg' in solve
    solve_xi   = 'xi' in solve

    for iteration in range(1, max_iter + 1):
        atm = _load_atmosphere(
            params['teff_K'], params['logg'], params['feh'],
            params['vturb_kms'], model_grid=model_grid
        )
        linemasks, spec_abund, x_over_h, x_over_fe = _ew_to_abundance(
            ew_df, params, atm, model_grid, star_id=star_id
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
    # RYA-279: ONE predicate for the gate pool AND the per-line CSV (below), so they
    # cannot diverge — both index this exact mask array, not a re-derived threshold.
    _fe1_ceiling_mask = _fe1_ceiling_pool_mask(notes, _lm_ew, _ew_ceiling)
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
        if note == 'Fe 1' and not _fe1_ceiling_mask[i]:
            continue  # same mask as the gate pool — never a re-derived threshold (RYA-279)
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


# ── RYA-329: line_score-weighted aggregation + the discrimination GATE ─────────
# Decision (Ryan): the primary A_X/A_X_nlte should be a line_score-WEIGHTED MEDIAN —
# but ONLY if line_score genuinely discriminates good vs bad lines. The primary must
# never be weighted on an unvalidated signal, so the weighting is gated per pool by
# `_line_score_discriminates`. Where the gate fails, the primary stays the plain
# (grade-blind but outlier-robust) median and the scoring repair routes to RYA-220.
# All five estimators are emitted as loud diagnostics regardless of the gate outcome.
_LS_SUBSCORES = ['ew_snr_score', 'fit_chi2_score', 'saturation_score',
                 'abundance_outlier_score', 'nlte_correction_score']
_LS_INERT_STD = 0.02          # a sub-score with std below this is inert (never populated)
_LS_MIN_SPREAD = 0.05         # line_score needs at least this std to weight anything
_LS_MIN_RHO = 0.30            # non-circular |Spearman| vs |a-median| to count as discriminating


def _weighted_median(values, weights):
    """Weighted median: the value where the cumulative weight crosses half the total."""
    v = np.asarray(values, float); w = np.asarray(weights, float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[m], w[m]
    if v.size == 0:
        return np.nan
    order = np.argsort(v)
    v, w = v[order], w[order]
    cw = np.cumsum(w)
    return float(v[np.searchsorted(cw, 0.5 * cw[-1])])


def _line_score_discriminates(lines: pd.DataFrame) -> dict:
    """RYA-329 gate: does line_score separate good vs bad lines on a NON-circular basis?

    abundance_outlier_score is derived from |a-median|, so it must be excluded before
    testing line_score against abundance deviation (else the test is self-fulfilling).
    Verdict DISCRIMINATES only if the non-circular part of the signal has real spread AND
    correlates with deviation. Returns the gate report (never raises)."""
    from scipy import stats
    n = len(lines)
    a = lines['a_1dlte'].astype(float).values
    ls = lines['line_score'].astype(float).values
    dev = np.abs(a - np.median(a))
    inert = [s for s in _LS_SUBSCORES
             if s in lines and float(lines[s].astype(float).std(ddof=0)) < _LS_INERT_STD]
    # non-circular signal = mean of the sub-scores EXCEPT abundance_outlier_score
    noncirc_cols = [s for s in _LS_SUBSCORES
                    if s != 'abundance_outlier_score' and s in lines]
    ls_noncirc = lines[noncirc_cols].astype(float).mean(axis=1).values if noncirc_cols else ls
    rho_circ = rho_nonc = float('nan')
    if n >= 5 and np.ptp(dev) > 0:
        if np.ptp(ls) > 0:
            rho_circ = float(stats.spearmanr(ls, dev).statistic)
        if np.ptp(ls_noncirc) > 0:
            rho_nonc = float(stats.spearmanr(ls_noncirc, dev).statistic)
    spread = float(np.std(ls, ddof=0))
    discriminates = (spread >= _LS_MIN_SPREAD and np.isfinite(rho_nonc)
                     and rho_nonc <= -_LS_MIN_RHO)
    return {'n': n, 'ls_spread': round(spread, 3), 'inert_subscores': inert,
            'rho_circular': round(rho_circ, 3) if np.isfinite(rho_circ) else None,
            'rho_noncircular': round(rho_nonc, 3) if np.isfinite(rho_nonc) else None,
            'discriminates': bool(discriminates)}


def _aggregation_diagnostics(lines: pd.DataFrame, delta_nlte: float) -> dict:
    """The five estimators + n-per-grade for one species pool (RYA-329 loud diagnostics).
    a_1dlte is 1D-LTE; the species delta_nlte is added to give the NLTE-scale twins."""
    a = lines['a_1dlte'].astype(float).values
    w = lines['line_score'].astype(float).values
    d = float(delta_nlte) if np.isfinite(delta_nlte) else 0.0
    gc = lines['line_grade'].value_counts()
    umean = float(np.mean(a)); pmed = float(np.median(a))
    wmean = float(np.average(a, weights=w)) if w.sum() > 0 else float('nan')
    wmed = _weighted_median(a, w)
    ab = lines[lines['line_grade'].isin(['A', 'B'])]
    abcut = (float(np.average(ab['a_1dlte'], weights=ab['line_score']))
             if len(ab) and ab['line_score'].sum() > 0 else float('nan'))
    return {
        'unweighted_mean': round(umean, 4), 'plain_median': round(pmed, 4),
        'weighted_mean': round(wmean, 4), 'weighted_median': round(wmed, 4),
        'ab_cut': round(abcut, 4) if np.isfinite(abcut) else None,
        'wmed_minus_median': round(wmed - pmed, 4),
        'n_A': int(gc.get('A', 0)), 'n_B': int(gc.get('B', 0)),
        'n_C': int(gc.get('C', 0)), 'n_D': int(gc.get('D', 0)),
        'plain_median_nlte': round(pmed + d, 4), 'weighted_median_nlte': round(wmed + d, 4),
    }


def _element_grade_summary(scored_df: pd.DataFrame, results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute element_score, element_grade, n_lines_A/B/C/D and A+B weighted abundance.
    Updates results_df in place and returns it.
    """
    T = LINE_GRADE_THRESHOLDS
    results = results_df.copy()
    _gate_reports = []

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

        # ── RYA-329: aggregation diagnostics + the discrimination gate ─────────
        d_nlte = float(row.get('delta_nlte_mean', np.nan))
        diag = _aggregation_diagnostics(lines, d_nlte)
        gate = _line_score_discriminates(lines)
        for k, v in diag.items():
            results.at[idx, f'agg_{k}'] = v
        results.at[idx, 'ls_discriminates']   = gate['discriminates']
        results.at[idx, 'ls_spread']          = gate['ls_spread']
        results.at[idx, 'ls_rho_noncircular'] = gate['rho_noncircular']
        results.at[idx, 'ls_inert_subscores'] = ';'.join(gate['inert_subscores'])
        # Primary weighting is GATED: only weight on a validated signal (Ryan, RYA-329).
        if gate['discriminates']:
            a_wm = diag['weighted_median']
            results.at[idx, 'A_X']      = round(a_wm, 3)
            results.at[idx, 'A_X_nlte'] = round(a_wm + (d_nlte if np.isfinite(d_nlte) else 0.0), 3)
            results.at[idx, 'primary_estimator'] = 'line_score_weighted_median (RYA-329)'
        else:
            results.at[idx, 'primary_estimator'] = 'plain_median (RYA-329 gate FAIL: line_score not discriminating -> RYA-220)'
        _gate_reports.append((f'{elem} {ion}', diag, gate))

    # ── Loud aggregation-diagnostics summary (RYA-329) ─────────────────────────
    if _gate_reports:
        print("\n  RYA-329 aggregation diagnostics (line_score-weighted median GATED on discrimination):")
        print("    species  n(A/B/C/D)  u.mean  p.median  w.mean  w.median  A+B   wmed-med  ls_spread rho_nc  discriminates")
        n_pass = 0
        for name, dg, gt in _gate_reports:
            n_pass += int(gt['discriminates'])
            print(f"    {name:7s}  {dg['n_A']}/{dg['n_B']}/{dg['n_C']}/{dg['n_D']:<3}"
                  f"  {dg['unweighted_mean']:.3f}  {dg['plain_median']:.3f}  {dg['weighted_mean']:.3f}"
                  f"  {dg['weighted_median']:.3f}  {str(dg['ab_cut'])[:5]:5s}  {dg['wmed_minus_median']:+.3f}"
                  f"    {gt['ls_spread']:.3f}   {str(gt['rho_noncircular']):>6s}   {gt['discriminates']}")
            if gt['inert_subscores']:
                print(f"             ^ INERT sub-scores (never populated -> RYA-220): {gt['inert_subscores']}")
        if n_pass == 0:
            print("    GATE VERDICT: line_score does NOT discriminate on any pool -> PRIMARY STAYS PLAIN MEDIAN "
                  "(unchanged, no re-bank). The apparent line_score signal is circular (abundance_outlier_score "
                  "= f(|a-median|)) and mostly inert -> scoring repair owed to RYA-220 before the primary can be weighted.")

    return results


# ── Solar EW loader (hybrid Fe I GES + Fe II lines_fit) ──────────────────────

def _apply_fe2_ew_quality_cull(ew_df: pd.DataFrame) -> pd.DataFrame:
    """RYA-352 — cull the solar Fe II EW pool on its REAL basis: EW measurement
    quality, computed from the table (no frozen list, no RYA-347 synth verdict).

    A Fe II EW line is culled if it fails ANY of:
      • SATURATION — EW > PIPELINE['vmic_ew_ceiling_mA'] (the COG-derived ceiling,
        RYA-330). Saturated lines on the flat curve of growth over-read because the
        1D-LTE model can't reproduce their saturation.
      • MEASUREMENT — ew_err/EW > PIPELINE['fe2_ew_err_frac_max'].
      • BLEND — blend_flag == True.

    Every drop is logged with its reason (no silent cull). Only Fe II rows are
    touched — this is the INTERIM Fe II-only gate (RYA-199: a fixed REW ceiling
    wrongly kills strong Ba/Na lines; the element-aware general gate is RYA-349).
    Supersedes the misnamed `fe2_solar_cull_rya347.csv` + `_apply_fe2_cull`
    (a88ef0f), which dropped clean EW lines (5991/6084/6456) for failing a *synth*
    fit — an independent path.
    """
    ceiling = PIPELINE['vmic_ew_ceiling_mA']
    err_max = PIPELINE['fe2_ew_err_frac_max']
    is_fe2 = (ew_df['element'] == 'Fe') & (ew_df['ion'] == 'II')
    if not is_fe2.any():
        return ew_df

    keep = pd.Series(True, index=ew_df.index)
    dropped = []
    for i in ew_df.index[is_fe2]:
        ew = float(ew_df.at[i, 'ew_mA'])
        wl = float(ew_df.at[i, 'wavelength_air_A'])
        err = float(ew_df.at[i, 'ew_err_mA']) if 'ew_err_mA' in ew_df.columns else np.nan
        blend = bool(ew_df.at[i, 'blend_flag']) if 'blend_flag' in ew_df.columns else False
        rew = np.log10(ew / 1000.0 / wl) if ew > 0 else np.nan
        reasons = []
        if ew > ceiling:
            reasons.append(f"SAT(EW={ew:.1f}>{ceiling:.0f}mÅ, REW={rew:.2f})")
        if np.isfinite(err) and ew > 0 and err / ew > err_max:
            reasons.append(f"HIERR(err/EW={err/ew:.2f}>{err_max})")
        if blend:
            reasons.append("BLEND")
        if reasons:
            keep[i] = False
            dropped.append((wl, rew, ';'.join(reasons)))

    if dropped:
        print(f"  [Fe II EW-quality cull RYA-352] dropped {len(dropped)} of "
              f"{int(is_fe2.sum())} Fe II lines (ceiling {ceiling:.0f}mÅ, "
              f"err/EW>{err_max}):")
        for wl, rew, why in sorted(dropped):
            print(f"      {wl:9.3f}  {why}")
    return ew_df[keep].reset_index(drop=True)


def _load_solar_ews(ew_override: str = None) -> pd.DataFrame:
    """
    Load solar EWs for abundance derivation using a hybrid approach:
      - Fe I:  GES pre-stored EWs (solar_ew_ges_reference.csv) — avoids NLTE EW bias
      - Fe II: committed canonical EWs (sol_ew_results_v1.csv) — Fe II not NLTE-affected
      - Other: committed canonical EWs (sol_ew_results_v1.csv)

    RYA-408: the non-Fe-I pool is read from the COMMITTED canonical
    PATHS['solar_ew_canonical'] (data/measured/sol_ew_results_v1.csv), NOT the
    gitignored lines_fit staging output data/processed/solar_ew.csv. No gate or
    abundance derivation may take a gitignored, regenerable file as its EW input
    (root cause of the RYA-406 incident). The reviewed staging→canonical promotion
    that preserves vetted blend_flags lives in scripts/promote_solar_ew.py.

    RYA-330: the GES Fe I routing is RETAINED by evidence, not inertia. RYA-328
    flagged that the GES reference pool is unvetted relative to program-star pools,
    and proposed re-routing Fe I through the measured solar_ew.csv "for consistency".
    That was tested and REJECTED: the measured solar Fe I EWs give A(Fe I) NLTE 7.67
    (gate fail), scatter 0.44, and a +0.43 reduced-EW slope — strictly worse than the
    GES reference (7.516 / 0.151 / −0.05). The measured solar EWs carry inflated
    strong-line EWs that the GES curation does not. The actual fix for the −0.21 slope
    RYA-328 found is the saturation cut, not the pool: vmic_ew_ceiling_mA was tightened
    150 → 100 mÅ (config.constants.PIPELINE), which removes the 12 saturated GES lines
    that drove the slope and applies the SAME bar to all stars. Consistency is achieved
    by the shared cut, not by forcing the Sun onto a worse pool.

    If ew_override is provided, use that file directly (bypasses hybrid logic).
    If solar_ew_ges_reference.csv is missing, falls back to the canonical for all elements.
    """
    if ew_override:
        ew_df = pd.read_csv(ew_override)
        ew_df = ew_df[(ew_df['ew_mA'] > 0) & ew_df['ew_mA'].notna()].copy()
        print(f"  EW override: {ew_override} ({len(ew_df)} lines)")
        return _apply_fe2_ew_quality_cull(ew_df)

    # RYA-408: read the COMMITTED canonical, never the gitignored staging file.
    canon_path = Path(str(PATHS['solar_ew_canonical']))
    if not canon_path.exists():
        raise FileNotFoundError(
            f"Canonical solar EWs not found at {canon_path} (RYA-408). The gate/abundance "
            f"path requires the committed canonical sol_ew_results_v1.csv — it must not fall "
            f"back to the gitignored staging file data/processed/solar_ew.csv. Restore the "
            f"canonical or run scripts/promote_solar_ew.py to promote a reviewed staging set."
        )
    solar_ew = pd.read_csv(str(canon_path))
    solar_ew = solar_ew[(solar_ew['ew_mA'] > 0) & solar_ew['ew_mA'].notna()].copy()
    print(f"  canonical sol_ew_results_v1.csv: {len(solar_ew)} lines total")

    # GES Fe I reference is committed in data/processed (force-added past the gitignore).
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
        print(f"  Hybrid EW: {len(ges_fe1)} Fe I (GES) + {fe2_count} Fe II (canonical) "
              f"+ {len(non_fe1) - fe2_count} other elements (canonical)")
        return _apply_fe2_ew_quality_cull(hybrid)

    except FileNotFoundError:
        print(f"  WARNING: GES Fe I reference not found at {ges_ref_path} — "
              f"falling back to canonical sol_ew_results_v1.csv for all elements")
        return _apply_fe2_ew_quality_cull(solar_ew)


# ── RYA-456: curated non-Fe pool → production A(X) ───────────────────────────
# RYA-371 Phase C found the non-Fe EW elements drop at the GES synthesis-region
# match (_build_ispec_line_regions reads the SPARSE pre-built regions file; for
# Eu/Ba/Sr it has 0 regions), and the few that survive (Cr/Ti/Ni) read wildly
# high on the un-curated pool. RYA-395/398 already built the abundance-BLIND,
# graded-gf curation (curate_nonfe_pools) but it was never wired into the run.
#
# The fix routes the non-Fe metals through the curation's OWN line-building path
# (compute_lte_abundances matches the FULL GES synthesis atomic linelist, note +
# wave ≤0.05 Å — the path that actually carries Eu/Ba/Sr), so the production
# baseline uses the SAME kept lines and SAME A(X) the standalone curation kept.
# This is a routing change only: the cull logic stays in curate_nonfe_pools
# (blind); we never re-derive a threshold or pull a value toward Asplund.
#
# Fe (its own curated pool, RYA-279/347), C/N/O (Phase A synthesis), Li (RYA-103
# force-include upper limit) and P (ground-unreachable DATA-GAP) are EXCLUDED —
# everything else in the solar EW pool is curated here.
_CURATED_NONFE_EXCLUDE = frozenset({'Fe', 'C', 'N', 'O', 'Li', 'P'})


def _curate_nonfe_pool(elements):
    """Reuse curate_nonfe_pools end-to-end (NOT reimplemented): the blind RYA-395
    quality cull + the RYA-398 graded-gf firewall, then the post-cull 1D-LTE A(X).
    Returns the per-line pool DataFrame (carries `kept`, `cull_reason`, `A_lte`)."""
    import pipeline.curate_nonfe_pools as cur
    pool = cur.load_pool(elements)
    atomic = ispec.read_atomic_linelist(_SYNTH_LINELIST_FILE)
    pool = cur.add_blend_ratio(pool, atomic)
    pool = cur.apply_cull(pool, grade_restrict=True)        # RYA-398 independent-gf firewall
    cur.assert_no_astrophysical_gf(pool, "RYA-456 non-Fe wiring")  # blindness firewall stays live
    pool = cur.compute_lte_abundances(pool)                 # A(X) read AFTER the blind cull
    return pool


def _curated_nonfe_rows(elements):
    """Per-element curated A(X) for the non-Fe metals, in production results-row
    format. The median is taken over the curation's KEPT lines of the dominant ion
    (the abundance workhorse: neutral for the light metals — matching
    curate_nonfe_pools.element_diagnostic — singly-ionized for the s-process / Eu),
    NLTE-applied with the curation's solar Δ. Carries the curation's BLIND verdict
    (VALIDATED / RESIDUAL / LOW_CONFIDENCE) for the Phase C classifier. Returns
    (rows_df, pool); rows_df has one row per element that produced ≥1 curated line."""
    import pipeline.curate_nonfe_pools as cur
    pool = _curate_nonfe_pool(elements)
    rows = []
    for el in sorted(set(elements)):
        df_el = pool[pool['element'] == el]
        kept = df_el[df_el['kept'] & df_el['A_lte'].notna()]
        if kept.empty:
            continue
        ion = kept['ion'].value_counts().idxmax()           # dominant (most-kept) ion
        sub = kept[kept['ion'] == ion]
        a_lte = float(np.median(sub['A_lte']))
        scatter = float(np.std(sub['A_lte'])) if len(sub) > 1 else np.nan
        n_clean = int(len(sub))
        dnlte = cur.solar_nlte_delta(el)
        dnlte0 = 0.0 if not np.isfinite(dnlte) else float(dnlte)
        a_nlte = a_lte + dnlte0
        asp = SOLAR_ASPLUND2021.get(el, np.nan)
        resid = (a_nlte - asp) if np.isfinite(asp) else np.nan
        tol = cur.validation_tol(el)
        # Verdict vocabulary = curate_nonfe_pools.element_diagnostic(graded=True),
        # read AFTER the blind cull (validation, never a cull input).
        if n_clean < cur.N_CLEAN_FLOOR:
            cverdict = 'LOW_CONFIDENCE'
        elif np.isfinite(resid) and abs(resid) <= tol:
            cverdict = 'VALIDATED'
        else:
            cverdict = 'RESIDUAL'
        applied = np.isfinite(dnlte) and abs(dnlte0) > 0.0
        rows.append({
            'element': el, 'ion': ion,
            'A_X': round(a_lte, 3),
            'A_X_std': round(scatter, 3) if np.isfinite(scatter) else np.nan,
            'n_lines': n_clean,
            'XH': round(a_lte - asp, 3) if np.isfinite(asp) else np.nan,
            'XFe': np.nan,
            'scale': 'absolute',
            'nlte_applied': bool(applied),
            'nlte_flag': 'NLTE' if applied else '1D_LTE',
            'A_X_nlte': round(a_nlte, 3) if np.isfinite(a_nlte) else np.nan,
            'A_X_nlte_absolute': round(a_nlte, 3) if applied else np.nan,
            'A_X_std_nlte': round(scatter, 3) if np.isfinite(scatter) else np.nan,
            'curation_verdict': cverdict,
            'curation_residual_nlte': round(resid, 3) if np.isfinite(resid) else np.nan,
            'nlte_delta_solar': round(dnlte, 4) if np.isfinite(dnlte) else np.nan,
            'curation_source': 'RYA-395/398 graded blind cull (curate_nonfe_pools)',
        })
    return pd.DataFrame(rows), pool


def _wire_curated_nonfe_results(results, ew_df):
    """RYA-456 wiring: replace the un-curated / dropped non-Fe metal rows in the
    production `results` with the curated (RYA-395/398) A(X). Fe, C/N/O, Li and P
    rows are left byte-stable. Returns (results_out, curated_rows, pool)."""
    routed = sorted({str(e) for e in ew_df['element'].unique()} - _CURATED_NONFE_EXCLUDE)
    curated, pool = _curated_nonfe_rows(routed)
    produced = set(curated['element']) if not curated.empty else set()
    kept = results[~results['element'].isin(routed)].copy()
    out = pd.concat([kept, curated], ignore_index=True) if not curated.empty else kept
    out = out.sort_values(['element', 'ion']).reset_index(drop=True)
    thin = sorted(set(routed) - produced)
    print(f"  RYA-456 curated non-Fe wiring: routed {len(routed)} EW-pool metal(s); "
          f"{len(produced)} produced curated A(X), {len(thin)} no curated line "
          f"(thin/over-culled: {thin}).")
    for _, r in curated.sort_values('element').iterrows():
        a_nlte = r['A_X_nlte']
        a_nlte_s = f"{a_nlte:.3f}" if pd.notna(a_nlte) else "  -  "
        print(f"    {r['element']:>3s} {r['ion']:<3s} A(LTE)={r['A_X']:.3f} "
              f"A(NLTE)={a_nlte_s} n={int(r['n_lines'])}  "
              f"Δvs_Asplund={r['curation_residual_nlte']}  [{r['curation_verdict']}]")
    # ── RYA-519 Part B/C: a routed metal that produced no value carries a CITED
    #    registry disposition (problem_children, RYA-463) instead of silently
    #    vanishing from the verdict; an expected non-emitter with NO registry
    #    disposition is a LOUD failure (Part C), never a silent skip. ──
    import pipeline.problem_children as _pc
    disp_rows, unregistered = [], []
    for el in thin:
        d = _pc.disposition_for(el)
        if d is None:
            unregistered.append(el)
            continue
        disp_rows.append({
            'element': el, 'ion': '', 'A_X': np.nan, 'A_X_nlte': np.nan,
            'n_lines': 0, 'scale': 'absolute', 'nlte_flag': 'not_emitted',
            'curation_verdict': d['problem_class'],
            'disposition_class': d['problem_class'],
            'disposition_treatment': d['required_treatment'],
            'disposition_status': d['status'],
            'disposition_tickets': d['governing_tickets'],
            'disposition_source': 'RYA-463 problem_children registry',
        })
        print(f"    {el:>3s}  no EW-pool value -> registry disposition: "
              f"{d['problem_class']} / {d['required_treatment']} [{d['status']}] "
              f"(RYA-{d['governing_tickets']})")
    if unregistered:
        raise RuntimeError(
            f"RYA-519 Part C LOUD-FAIL: expected metal(s) {unregistered} produced no "
            f"value AND carry no problem_children (RYA-463) disposition. A silently-"
            f"omitted element is a defect — add a cited disposition in "
            f"pipeline/problem_children.py or fix the wiring before the verdict counts.")
    if disp_rows:
        out = pd.concat([out, pd.DataFrame(disp_rows)], ignore_index=True)
        out = out.sort_values(['element', 'ion']).reset_index(drop=True)
    return out, curated, pool


# ── RYA-458: EW-verification layer emission ──────────────────────────────────

def _measured_ew_path(star_id: str):
    """Resolve the MEASURED per-line EW pool for the EW-integrity pass (RYA-458/273).

    Solar uses its promoted canonical pool (sol_ew_results_v1.csv, RYA-408); any other
    star uses its lines_fit EW staging output (PATHS['<star>_ew'], e.g. procyon_ew.csv)
    — the file where the per-line chi2 / profile / blend_flag / notes actually live.
    Falls back to the solar canonical only if the star has no resolvable EW product."""
    s = star_id.lower()
    if 'solar' in s or 'sun' in s:
        return PATHS['solar_ew_canonical']
    key = next((k for k in STAR_LINELISTS if k != 'solar' and k in s), None)
    cand = PATHS.get(f'{key}_ew') if key else None
    return cand if (cand is not None and Path(cand).exists()) else PATHS['solar_ew_canonical']


def _emit_ew_integrity(out_dir, star_id, scored_df, curated_pool):
    """Run the RYA-458 EW-integrity QA pass and save {star}_ew_integrity.csv.

    Builds the per-line frame from the star's MEASURED EW pool (where chi2 / notes /
    blend_flag live — solar canonical, or the per-star lines_fit output for Procyon /
    others, RYA-273), attaches the implied A(X) (Fe from the per-line scoring, non-Fe
    from the curated RYA-456 pool when present) for the ABUND_OUTLIER check, flags
    integrity, and asserts no EW was mutated. Flags only — never edits an EW (RYA-451/458).

    For F-stars (Procyon) the value is the BAD_FIT / ABUND_OUTLIER / COG_FLAG catch on
    the broader, saturation-/blend-prone profiles; the solar charter cases (C I 5380 /
    Li 6707 / Eu 6645) are simply ABSENT and reported as such."""
    import pipeline.ew_integrity as ei
    ew_src = _measured_ew_path(star_id)
    measured = pd.read_csv(str(ew_src), comment='#')
    measured = measured[(measured['ew_mA'] > 0) & measured['ew_mA'].notna()].reset_index(drop=True)

    # implied A(X) per line, keyed (element, ion, wave@2dp): Fe from scoring, else curated.
    a_map: dict = {}
    if scored_df is not None and not scored_df.empty and 'a_1dlte' in scored_df.columns:
        for _, r in scored_df.iterrows():
            a_map[(r['element'], r['ion'], round(float(r['wavelength_air_A']), 2))] = float(r['a_1dlte'])
    if curated_pool is not None and 'A_lte' in getattr(curated_pool, 'columns', []):
        for _, r in curated_pool[curated_pool['A_lte'].notna()].iterrows():
            a_map.setdefault((r['element'], r['ion'], round(float(r['wavelength_air_A']), 2)),
                             float(r['A_lte']))
    measured['a_lte'] = [a_map.get((e, i, round(float(w), 2)), np.nan)
                         for e, i, w in zip(measured['element'], measured['ion'],
                                            measured['wavelength_air_A'])]

    flagged = ei.flag_ew_integrity(measured)
    ei.assert_no_ew_mutation(measured, flagged)          # cardinal guard (also asserted inside)
    out_path = out_dir / f'{star_id}_ew_integrity.csv'
    flagged.to_csv(out_path, index=False)

    n_flag = int((flagged['ew_integrity'].astype(str).str.len() > 0).sum())
    n_excl = int(flagged['ew_excluded'].sum())
    print(f"  RYA-458 EW-verify: {len(flagged)} lines scanned, {n_flag} flagged, "
          f"{n_excl} excluded (no EW mutated). Saved → {out_path.name}")
    cs = ei.charter_summary(flagged)
    for key, c in cs.items():
        if not c.get('present'):
            print(f"    charter {key}: ABSENT from the measured pool")
            continue
        print(f"    charter {key}: EW {c['ew_mA']:.2f} mA  flags=[{c['ew_integrity']}]  "
              f"disposition={c['ew_disposition'] or '-'}  excluded={c['ew_excluded']}")
    return flagged


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
        xi_override: float = None,
        engine: str = 'spectrum',
        skip_convergence: bool = False,
        ew_verify: bool = False) -> tuple:
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
    xi_override              : explicit ξ (km/s) for experimental runs (e.g. the
                               RYA-322 ξ sweep) — pins vturb_kms to this value and
                               drops ξ from the solve, logged. Replaces the retired
                               silent vmic_fixed flag (RYA-325). None = use policy.
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
    else:
        # RYA-298: fundamental teff/logg/feh come from the single source (stars.yaml
        # via get_star_params); the legacy dict supplies only NON-fundamental fields
        # (names, RV, ξ-solve seed, Ni 6300 literature EW, broadening).
        if 'solar' in star_id.lower():
            params = STAR_SOLAR.copy()
        elif 'procyon' in star_id.lower():
            params = STAR_PROCYON.copy()
        else:
            params = STAR_55CNC.copy()
        _fund = get_star_params(star_id)
        params['teff_K'] = float(_fund['teff'])
        params['logg'] = float(_fund['logg'])
        params['feh'] = float(_fund['feh_ref'])
        params.setdefault('vturb_kms',
                          float(_fund.get('xi', _fund.get('xi_init', 1.0))))

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
        # RYA-325: ξ pin/solve. xi_override (the ONE explicit experimental knob —
        # e.g. the RYA-322 ξ sweep) supersedes the pin and is logged; it replaces
        # the retired silent vmic_fixed default. Else a pinned ξ is set from
        # STAR_PARAMS (fail loud if pinned without a value); a star whose ξ is
        # neither pinned-with-value nor solved raises (no silent guess).
        if xi_override is not None:
            params['vturb_kms'] = float(xi_override)
            solve_params = [p for p in solve_params if p != 'xi']
            print(f"  ξ OVERRIDE (explicit, RYA-322): vturb pinned to "
                  f"{float(xi_override):.2f} km/s (no re-solve)")
        elif 'xi' in star_policy.get('pin', []):
            if 'xi' not in star_policy:
                raise KeyError(f"{star_id}: ξ is pinned but no 'xi' value in "
                               f"STAR_PARAMS (RYA-325) — refusing to guess.")
            params['vturb_kms'] = float(star_policy['xi'])
        elif 'xi' not in solve_params:
            raise KeyError(f"{star_id}: ξ is neither pinned-with-value nor in "
                           f"'solve' (RYA-325) — refusing to guess.")
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
        last_linemasks, _, _, _ = _ew_to_abundance(ew_df, params, _atm, model_grid, star_id=star_id)
        out_dir = ns.outputs_dir(star_id)                                  # RYA-469 namespaced
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
        ew_df, params, model_grid=model_grid,
        skip_convergence=skip_convergence, solve_params=solve_params, star_id=star_id,
    )

    print(f"\n  Final params: Teff={converged_params['teff_K']:.0f} K  "
          f"logg={converged_params['logg']:.2f}  "
          f"vturb={converged_params['vturb_kms']:.2f} km/s")

    # ── NLTE corrections (RYA-165) — LOUD-FAIL, no silent 1D-LTE fallback (RYA-525 §3) ──
    # nlte_corrections is loud by design (FileNotFoundError / ValueError[PLACEHOLDER_ZERO]
    # / KeyError). The former `except Exception -> print("...1D LTE only") -> continue`
    # swallows here UNDID that one layer up and let an element ride into a PASS on a
    # silently-downgraded LTE value — the exact cop-out the two-engine floor outlaws
    # (RYA-536 Review 6a). Removed: a missing wired grid / any NLTE failure now RAISES.
    # An element that is LTE-BY-DESIGN carries its cited disposition in the RYA-526
    # two-engine coverage ledger (LTE-only-by-design) — it is never a swallowed exception.
    from pipeline.nlte_corrections import (
        apply_fe_nlte_corrections, apply_element_nlte_corrections)
    _sun = get_star_params('solar')  # RYA-298: canonical solar fallback (5772/4.438)
    stellar_params_for_nlte = {
        'teff_K'   : converged_params.get('teff_K',    _sun['teff']),
        'logg'     : converged_params.get('logg',      _sun['logg']),
        'feh'      : converged_params.get('feh',        0.0),
        'vturb_kms': converged_params.get('vturb_kms',  1.0),
    }
    results = apply_fe_nlte_corrections(
        results, stellar_params_for_nlte, line_df=ew_df, per_line_df=per_line_df,
    )
    print(f"  NLTE corrections applied to Fe I/II (Amarsi+2022)")

    # ── Non-Fe element NLTE (RYA-235): Ca I / Ti I / Cr I per-line, MPIA MAFAGS-OS.
    # Composes with the Fe leg (touches only registry rows). NOTE: this is the
    # application layer — Asplund acceptance is gated on a CURATED pool (raw 1D-LTE
    # Cr reads high; NLTE on an un-curated pool overshoots), tracked under RYA-235.
    # LOUD-FAIL as above (RYA-525 §3) — no silent stay-1D-LTE on failure.
    results = apply_element_nlte_corrections(
        results, stellar_params_for_nlte, per_line_df=per_line_df, line_df=ew_df,
    )
    print(f"  NLTE corrections applied to Ca I / Ti I / Cr I (MPIA MAFAGS-OS)")

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
            # RYA-330: A_X_nlte is stored ABSOLUTE by nlte_corrections.py in every
            # path — A_X_nlte = median(a_1dlte[i] + aberr[i]) and a_1dlte = A_X, the
            # absolute MOOG abundance (see nlte_corrections.py:218,271). The earlier
            # `+ ref` here (and in the solar gate below) was a stale RYA-267 relative-
            # scale assumption that re-added the 7.46 solar offset, inflating the
            # absolute NLTE A(Fe) to ~15 and making the solar gate print FAIL on a
            # passing anchor. A_X_nlte is already absolute → return it directly.
            flag = str(row.get('nlte_flag', '1D_LTE'))
            if flag in ('1D_LTE', 'NLTE_unavailable'):
                return np.nan   # NLTE not actually applied → no absolute NLTE value
            val = row.get('A_X_nlte', np.nan)
            return round(float(val), 3) if pd.notna(val) else np.nan
        results['A_X_nlte_absolute'] = results.apply(_nlte_to_absolute, axis=1)

    # ── RYA-456: wire the curated (RYA-395/398) non-Fe pool into the baseline ──
    # The non-Fe EW metals drop / read wild through the GES region match; route
    # them through the curation's own (blind) line-building so the production
    # baseline carries the SAME kept lines + A(X) the standalone curation kept.
    # Fe / C / N / O / Li / P rows are untouched (asserted byte-stable downstream).
    _curated_pool = None
    if 'solar' in star_id.lower():
        results, _curated_nonfe, _curated_pool = _wire_curated_nonfe_results(results, ew_df)

    # ── RYA-334 range-sanity tripwire (output chokepoint) ─────────
    # Every absolute-scale column must sit on the A(H)=12 scale. A double-add of
    # the 7.46 solar offset (RYA-267 class) lands ~15 = above hydrogen → raise loud
    # here instead of silently writing it / printing it / hiding behind a gate.
    # Relative [X/H] (XH) is intentionally NOT checked (legitimately negative).
    for _abs_col in ('A_X', 'A_X_nlte', 'A_X_nlte_absolute', 'A_X_nlte_AB'):
        if _abs_col in results.columns:
            for _, _r in results.iterrows():
                assert_abundance_on_scale(
                    _r[_abs_col], f"{star_id} {_r.get('element','?')} {_r.get('ion','?')} {_abs_col}")

    # ── RYA-469 Deliverable D: pin the differential denominator version ──
    # A target's [X/H] is differential vs OUR measured Sun. Record which FROZEN gold
    # solar version is the denominator in the target's own output, so re-baselining the
    # Sun later (a new solar_abundances_v{N}) never SILENTLY changes this target's
    # already-derived numbers — you re-run the target against the new version on purpose.
    if 'solar' not in star_id.lower():
        try:
            results = ns.stamp_solar_ref_version(results)
            print(f"  [RYA-469] differential denominator pinned: solar_ref_version="
                  f"{results['solar_ref_version'].iloc[0]}")
        except ns.ImmutableReferenceError as _e:
            print(f"  [RYA-469] solar_ref_version unstamped (no gold reference yet): {_e}")

    # ── RYA-520: synthesis-required elements must NOT emit an authoritative
    #    EW-derived abundance. EW is the WRONG method for molecular-coupled C/N/O
    #    (RYA-237) and the HFS heavies — a saturated line on the flat curve-of-
    #    growth turns a tiny EW error into ~1.8 dex (solar C I 5380 -> C=10.26 vs
    #    the synthesis verdict 8.491). The raw EW value is preserved DIAGNOSTIC-
    #    ONLY; the authoritative abundance columns are nulled + labelled SUPERSEDED
    #    _BY_SYNTHESIS so no committed artifact carries a wrong flagship number.
    #    LOUD-FAIL if one slips through as authoritative. ──
    if 'solar' in star_id.lower():
        from config.constants import TARGET_ELEMENTS as _TGT520
        import pipeline.problem_children as _pc520
        _synth_req = {'C', 'N', 'O'} | {
            _e for _e in _TGT520
            if (_d := _pc520.disposition_for(_e)) and _d['required_treatment'] in ('synthesis', 'HFS_sum')}
        _abs_cols = [c for c in ('A_X', 'A_X_nlte', 'A_X_nlte_absolute', 'A_X_nlte_AB')
                     if c in results.columns]
        for _c, _dft in (('ew_diagnostic_A_X', np.nan), ('authoritative', True),
                         ('method', ''), ('superseded_by', '')):
            if _c not in results.columns:
                results[_c] = _dft
        for _i in results.index:
            _el = str(results.at[_i, 'element'])
            if _el not in _synth_req:
                continue
            # RYA-551: Sr's synthesis requirement is ION-SPECIFIC (Sr II resonance
            # doublet + subordinate lines route to synthesis). The clean Sr I 6617
            # EW leg stays authoritative — it is the RYA-422 ionization-balance
            # cross-check. Never suppress Sr I.
            if _el == 'Sr' and str(results.at[_i, 'ion']).strip() != 'II':
                continue
            _ax = results.at[_i, 'A_X'] if 'A_X' in results.columns else np.nan
            if pd.notna(_ax):
                results.at[_i, 'ew_diagnostic_A_X'] = _ax
                for _col in _abs_cols:
                    results.at[_i, _col] = np.nan
                results.at[_i, 'authoritative'] = False
                results.at[_i, 'method'] = 'EW_SUPERSEDED_BY_SYNTHESIS'
                results.at[_i, 'superseded_by'] = 'synthesis (RYA-237 C/N/O; HFS-resolved metals)'
                print(f"  [RYA-520] {_el} raw-EW A(X)={_ax:.3f} SUPPRESSED -> diagnostic-only "
                      f"(synthesis-required; authoritative value is the synthesis verdict).")
        _leak = sorted({str(results.at[_i, 'element']) for _i in results.index
                        if str(results.at[_i, 'element']) in _synth_req
                        and not (str(results.at[_i, 'element']) == 'Sr'
                                 and str(results.at[_i, 'ion']).strip() != 'II')
                        and results.at[_i, 'authoritative']
                        and pd.notna(results.at[_i, 'A_X'] if 'A_X' in results.columns else np.nan)})
        if _leak:
            raise RuntimeError(
                f"RYA-520 LOUD-FAIL: synthesis-required element(s) {_leak} still carry an "
                f"authoritative EW-derived A(X). EW is the wrong method for these "
                f"(molecular-coupled / HFS) — the value must come from synthesis.")

    # ── RYA-519 completeness: NO target element silently absent from the verdict.
    #    Every metal in TARGET_ELEMENTS either emits a value (n_lines>0) or carries
    #    a cited problem_children (RYA-463) disposition; neither = LOUD failure.
    #    Catches elements outside the EW pool (K/Co/Sc) and the CNO/P paths that
    #    this run did not surface — the amendment's "no silent absence" guarantee. ──
    if 'solar' in star_id.lower():
        from config.constants import TARGET_ELEMENTS as _TGT
        import pipeline.problem_children as _pc519
        _present = set(results['element'].astype(str))
        _comp_rows, _comp_unreg = [], []
        for _el in _TGT:
            if _el in _present:
                continue
            _d = _pc519.disposition_for(_el)
            if _d is None:
                _comp_unreg.append(_el)
                continue
            _comp_rows.append({
                'element': _el, 'ion': '', 'A_X': np.nan, 'A_X_nlte': np.nan,
                'n_lines': 0, 'scale': 'absolute', 'nlte_flag': 'not_emitted',
                'curation_verdict': _d['problem_class'],
                'disposition_class': _d['problem_class'],
                'disposition_treatment': _d['required_treatment'],
                'disposition_status': _d['status'],
                'disposition_tickets': _d['governing_tickets'],
                'disposition_source': 'RYA-463 registry (RYA-519 completeness pass)',
            })
            print(f"    {_el:>3s}  not surfaced by this run -> registry disposition: "
                  f"{_d['problem_class']} / {_d['required_treatment']} [{_d['status']}]")
        if _comp_unreg:
            raise RuntimeError(
                f"RYA-519 completeness LOUD-FAIL: target element(s) {_comp_unreg} are "
                f"absent from the solar verdict AND carry no problem_children (RYA-463) "
                f"disposition — a silently-omitted element is a defect.")
        if _comp_rows:
            results = pd.concat([results, pd.DataFrame(_comp_rows)], ignore_index=True)
            results = results.sort_values(['element', 'ion']).reset_index(drop=True)

    # ── Save ──────────────────────────────────────────────────────
    print(f"\n[4/4] Saving results...")
    out_dir  = ns.outputs_dir(star_id)                       # RYA-469 namespaced per-star
    out_path = out_dir / f'{star_id}_abundances.csv'
    # RYA-521: this raw EW file is DIAGNOSTIC-ONLY — the authoritative abundance is
    # the phase_c verdict. Stamp an in-file label so no consumer mistakes it for the
    # reported number (comment-aware readers strip it; see authoritative_channel.py).
    from pipeline.authoritative_channel import DIAGNOSTIC_ONLY_HEADER as _DIAG_HDR
    out_path.write_text(_DIAG_HDR + results.to_csv(index=False))
    print(f"  Saved → {out_path.name}  ({len(results)} elements, "
          f"{results['n_lines'].sum()} total lines)  [RYA-521 DIAGNOSTIC-ONLY]")

    if not scored_df.empty:
        scored_df['scale'] = 'absolute'   # RYA-289: per-line A_X is 1D-LTE absolute
        per_line_out = out_dir / f'{star_id}_per_line.csv'
        scored_df.to_csv(per_line_out, index=False)
        print(f"  Saved → {per_line_out.name}  ({len(scored_df)} lines with quality scores)")

    # ── RYA-458: EW-verification layer (opt-in --ew-verify) ───────
    # A per-line EW-INTEGRITY QA pass: flag BAD_FIT / ABUND_OUTLIER / COG_FLAG /
    # LIT_DEVIATION and assign the charter-case dispositions. It FLAGS only — a hard
    # assert proves no measured EW is mutated. Runs on the star's measured EW pool (where
    # the per-line fit quality + notes live), with the implied A(X) attached from the
    # Fe scoring + the curated non-Fe pool (RYA-456, solar) for the ABUND_OUTLIER check.
    # RYA-273: generalized beyond solar — Procyon's F-star pool runs through the same
    # layer (the curated non-Fe pool is solar-only, so curated_pool is None here and the
    # ABUND_OUTLIER check stands on the Fe per-line scoring, which is the F-star target).
    if ew_verify:
        _emit_ew_integrity(out_dir, star_id, scored_df, _curated_pool)

    # ── Synthesis mode (RYA-285 v1 EW / RYA-287 v2 flux) ──────────
    if engine in ('synthesis', 'synthesis-v2'):
        _atm_synth = _load_atmosphere(
            converged_params['teff_K'], converged_params['logg'],
            converged_params['feh'],    converged_params['vturb_kms'],
            model_grid=model_grid,
        )
        if engine == 'synthesis-v2':
            print(f"\n[4b/4] Synthesis-v2 flux-space fit (Turbospectrum)...")
            # RYA-342: feed the synth the MEASURED line set, NOT the EW-converged
            # last_linemasks. The flux-space synth must not inherit ANY EW-domain line
            # selection (saturation ceiling, blend_flag, or the theoretical-EW sanity
            # filter that drops recover-tier lines like 6239.943 — χ²ᵣ 1.32 — that
            # synthesis can measure). Its only acceptance gate is fit quality (χ²ᵣ).
            # Build from the same valid EW range the convergence uses, but BEFORE the
            # blend/sanity filter. (EW path keeps its converged set — correct there.)
            _emin = float(PIPELINE['ew_min_mA']); _emax = float(PIPELINE['ew_max_mA'])
            _ewm  = ew_df[(ew_df['ew_mA'] >= _emin) & (ew_df['ew_mA'] <= _emax)].copy()
            _synth_lm = _build_ispec_line_regions(_ewm)
            # Regression guard (RYA-342): the measured set must be a SUPERSET of the
            # EW-converged set — no EW-domain filter (ceiling / blend / convergence)
            # may shrink the synth's input. Fail loud if a converged line is missing.
            _cw = {round(float(w), 2) for w in last_linemasks['wave_A']}
            _mw = {round(float(w), 2) for w in _synth_lm['wave_A']}
            _miss = _cw - _mw
            assert not _miss, (
                f"RYA-342: {len(_miss)} EW-converged line(s) absent from the synth "
                f"measured set — an EW-domain filter shrank the synth input "
                f"({sorted(_miss)[:5]})")
            print(f"  [synth-v2] line set: {len(_synth_lm)} measured "
                  f"(vs {len(last_linemasks)} EW-converged) — synth gates on χ²ᵣ only")
            _run_synthesis_v2_mode(_synth_lm, converged_params, _atm_synth,
                                   star_id, out_dir)
        else:
            print(f"\n[4b/4] Synthesis-EW mode (Turbospectrum bisection)...")
            _run_synthesis_mode(last_linemasks, converged_params, _atm_synth,
                                star_id, out_dir)

    # ── Solar calibration gate (RYA-336: scale-aware) ─────────────
    # The PASS/FAIL verdict rests on the SCALE-ROBUST checks — reduced-EW slope,
    # Fe I−Fe II ionization balance, Fe I scatter — which carry no 1D-vs-3D abundance
    # term and are what the differential science rests on. The ABSOLUTE A(Fe) is a
    # scale-aware DIAGNOSTIC against a wide window centred on the 3D-true anchor +
    # the published 1D-3D offset (NOT our own output), there only to catch a gross
    # zero-point error (e.g. a loggf scale slip) that the shift-invariant strong
    # gates all miss. (Distinct from the RYA-334 A_X<12.5 double-add tripwire.)
    fe = results[results['element'] == 'Fe']
    if len(fe) > 0:
        _fe1_df = fe[fe['ion'] == 'I']
        _fe2_df = fe[fe['ion'] == 'II']

        def _abs_nlte_of(row_df):
            if row_df is None or len(row_df) == 0:
                return np.nan
            r = row_df.iloc[0]
            _f = str(r.get('nlte_flag', '1D_LTE'))
            _v = float(r.get('A_X_nlte', np.nan))
            return _v if (_f not in ('NLTE_unavailable', '1D_LTE') and np.isfinite(_v)) \
                   else float(r['A_X'])

        # (1) Fe I reduced-EW slope (line-formation shape; vturb guardrail)
        rew_slope = np.nan
        if not scored_df.empty:
            _f1 = scored_df[(scored_df['element'] == 'Fe') & (scored_df['ion'] == 'I')].copy()
            _f1 = _f1[(_f1['ew_mA'] > 0) & _f1['a_1dlte'].notna()]
            if len(_f1) >= 5:
                _rew = np.log10(_f1['ew_mA'].values / _f1['wavelength_air_A'].values)
                rew_slope = float(np.polyfit(_rew, _f1['a_1dlte'].values, 1)[0])
        slope_pass = np.isfinite(rew_slope) and abs(rew_slope) < FE_REW_SLOPE_GATE
        # (2) Fe I−Fe II ionization balance — scored on the SYNTHESIS arbiter (RYA-406).
        # The EW-path ΔFe is a DIAGNOSTIC: the EW Fe II pool is blend-limited and reads
        # HIGH BY DESIGN (RYA-352), so its imbalance is spurious. The ratified ionization
        # arbiter is flux-space synthesis (DECISION 2, RYA-305/341). Source: a single-
        # source cited constant for the calibrated Sun; otherwise the star's synth-v2
        # product; otherwise loud-fallback to the EW path (never a silent verdict swap).
        a_fe1 = _abs_nlte_of(_fe1_df); a_fe2 = _abs_nlte_of(_fe2_df)
        dfe_ew = (a_fe1 - a_fe2) if (np.isfinite(a_fe1) and np.isfinite(a_fe2)) else np.nan
        _arb = FE_IONIZATION_SYNTH_ARBITER.get(star_id)
        if _arb is not None:                          # ratified synth arbiter (cited)
            dfe = float(_arb['dFe'])
            fe2_synth = float(_arb['fe2_synth'])
            ion_src = f"synth arbiter ({_arb['provenance']})"
            ion_pass = abs(dfe) < FE_IONISATION_GATE
        else:                                          # no ratified arbiter for this star
            # (a synth-v2 product consumer would slot here; absent → EW path, flagged loud)
            dfe = dfe_ew
            fe2_synth = np.nan
            ion_src = "EW path (NO ratified synth arbiter for this star — RYA-406; diagnostic-grade verdict)"
            ion_pass = np.isfinite(dfe) and abs(dfe) < FE_IONISATION_GATE
        # EW-vs-synth Fe II spread — reported diagnostic, flagged if beyond the blend-bias
        # band (Socratic check: confirm it is the EW blend-bias, not an unexplained term).
        ew_synth_spread = (a_fe2 - fe2_synth) if np.isfinite(fe2_synth) else np.nan
        spread_flag = np.isfinite(ew_synth_spread) and abs(ew_synth_spread) > FE_EW_SYNTH_SPREAD_BAND
        # (3) Fe I scatter
        sc1 = np.nan
        if len(_fe1_df):
            _r1 = _fe1_df.iloc[0]
            sc1 = float(_r1.get('A_X_std_nlte', np.nan))
            if not np.isfinite(sc1):
                sc1 = float(_r1.get('A_X_std', np.nan))
        # Fe I scatter ceiling from the active acceptance profile (RYA-446/277): the Sun
        # (G anchor) reads the cited RYA-407 honest floor; types without a populated
        # profile keep the legacy-universal value until RYA-277 builds them (F/K/M TODO).
        # Threshold SOURCE only — the raw σ (sc1) computation is untouched.
        _spec = STAR_SPECTRAL_TYPE.get(star_id)
        if _spec is not None and ACCEPTANCE_PROFILES.get(_spec) is not None:
            _scatter_max = fe1_scatter_threshold(_spec)
            _scatter_src = f"{_spec}-anchor: {ACCEPTANCE_PROFILES[_spec]['provenance']}"
        else:
            _scatter_max = FE_SCATTER_GATE   # legacy-universal (RYA-277 profile pending)
            _scatter_src = f"legacy-universal {FE_SCATTER_GATE} (RYA-277 profile pending)"
        scat_pass = np.isfinite(sc1) and sc1 < _scatter_max
        primary_pass = slope_pass and ion_pass and scat_pass

        print(f"\n  ── Solar Fe gate — PRIMARY (scale-robust, RYA-336/406) ──")
        print(f"  Fe I reduced-EW slope = {rew_slope:+.3f}  -> {'PASS' if slope_pass else 'FAIL'} (|slope| < {FE_REW_SLOPE_GATE})")
        print(f"  Fe I-Fe II ionization = {dfe:+.3f}  -> {'PASS' if ion_pass else 'FAIL'} (|ΔFe| < {FE_IONISATION_GATE} dex)  [{ion_src}]")
        if np.isfinite(ew_synth_spread):
            _band = (f"FLAG > {FE_EW_SYNTH_SPREAD_BAND} (unexplained?)" if spread_flag
                     else f"within blend-bias band <= {FE_EW_SYNTH_SPREAD_BAND}")
            print(f"    · DIAGNOSTIC EW-vs-synth Fe II spread = {ew_synth_spread:+.3f} "
                  f"(EW {a_fe2:.3f} blend-limited HIGH, RYA-352; synth {fe2_synth:.3f}, RYA-341) -> {_band}")
            print(f"    · (EW-path ΔFe(I−II) = {dfe_ew:+.3f} is NOT the verdict — EW Fe II is blend-biased; RYA-405/406)")
        print(f"  Fe I scatter          = {sc1:.3f}  -> {'PASS' if scat_pass else 'FAIL'} (< {_scatter_max} dex)  [{_scatter_src}]")

        # Absolute A(Fe): scale-aware diagnostic (centre = 3D-true 7.46 + published 1D-3D offset)
        _diag_c  = SOLAR_ASPLUND2021['Fe'] + FE_1D3D_SOLAR_OFFSET
        _diag_lo = _diag_c - FE_ABS_DIAG_HALFWIDTH
        _diag_hi = _diag_c + FE_ABS_DIAG_HALFWIDTH
        print(f"\n  ── Absolute A(Fe): scale-aware DIAGNOSTIC ──")
        print(f"  (centre {_diag_c:.2f} = 7.46 + {FE_1D3D_SOLAR_OFFSET} 1D-3D offset; window [{_diag_lo:.2f},{_diag_hi:.2f}]; not the verdict)")
        for _, fe_row in fe.iterrows():
            ion_lbl = fe_row['ion']
            a_1d    = float(fe_row['A_X'])
            n_lines = int(fe_row.get('n_lines', 0))
            scatter_1d   = float(fe_row.get('A_X_std', np.nan))
            flag        = str(fe_row.get('nlte_flag', '1D_LTE'))
            nlte_ok     = flag not in ('NLTE_unavailable', '1D_LTE')
            _a_nlte_raw = float(fe_row.get('A_X_nlte', np.nan))
            a_nlte      = _a_nlte_raw if (nlte_ok and np.isfinite(_a_nlte_raw)) else a_1d
            # RYA-334: guard the gate value (a re-introduced double-add lands ~15)
            assert_abundance_on_scale(a_nlte, f"{star_id} Fe {ion_lbl} gate a_nlte")
            d_nlte  = float(fe_row.get('delta_nlte_mean', np.nan))
            n_nlte  = int(fe_row.get('n_nlte_lines', 0))
            nl_min  = 20 if ion_lbl == 'I' else 3   # Fe II optical coverage gap (RYA-227)
            diag_pass = _diag_lo <= a_nlte <= _diag_hi
            print(f"\n  Fe {ion_lbl}  1D LTE  = {a_1d:.3f}  (scatter 1D = {scatter_1d:.3f} dex)")
            if np.isfinite(d_nlte):
                print(f"  Mean Δ(NLTE) Fe {ion_lbl}= {d_nlte:+.4f} dex  ({n_nlte} lines corrected)")
            print(f"  A(Fe {ion_lbl}) NLTE abs = {a_nlte:.3f}  -> {'PASS' if diag_pass else 'FAIL'} (diagnostic [{_diag_lo:.2f},{_diag_hi:.2f}])")
            # RYA-553: for the SUN the tabulated Magic-2013 1D→3D correction is APPLIED to
            # the reported Fe I anchor (verdict/gold/gate layer). Show the 3D-corrected
            # value and gate it against the real FE_GATE [7.41,7.51]; the scale-aware
            # window above describes the un-corrected 1D-NLTE MEASUREMENT. SOLAR ONLY,
            # Fe I only (Fe II is NOT shifted); off-solar per-star 3D is owed (RYA-550).
            if 'solar' in star_id.lower() and ion_lbl == 'I':
                _fe_3d   = a_nlte + CORRECTIONS_3D['Fe_1D3D_solar_dex']
                _gate_lo = SOLAR_ASPLUND2021['Fe'] + FE_GATE_LOWER
                _gate_hi = SOLAR_ASPLUND2021['Fe'] + FE_GATE_UPPER
                _gate_ok = _gate_lo <= _fe_3d <= _gate_hi
                assert_abundance_on_scale(_fe_3d, f"{star_id} Fe I 3D-corrected reported anchor")
                print(f"  A(Fe I) 3D-NLTE     = {_fe_3d:.3f}  (1D-NLTE {a_nlte:.3f} "
                      f"{CORRECTIONS_3D['Fe_1D3D_solar_dex']:+.3f}, Magic 2013) -> "
                      f"{'PASS' if _gate_ok else 'FAIL'} (FE_GATE [{_gate_lo:.2f},{_gate_hi:.2f}], RYA-553)")
            print(f"  nlte_flag Fe {ion_lbl}   = {flag}")
            print(f"  Fe {ion_lbl} n_lines     = {n_lines}  -> {'PASS' if n_lines >= nl_min else 'FAIL'} (>={nl_min})")
            n_A = int(fe_row.get('n_lines_A', 0)); n_B = int(fe_row.get('n_lines_B', 0))
            n_C = int(fe_row.get('n_lines_C', 0)); n_D = int(fe_row.get('n_lines_D', 0))
            el_grade = str(fe_row.get('element_grade', '?'))
            a_ab_abs = float(fe_row.get('A_X_nlte_AB', np.nan))
            print(f"  Fe {ion_lbl} grade dist  = A:{n_A}  B:{n_B}  C:{n_C}  D:{n_D}  element_grade={el_grade}")
            if np.isfinite(a_ab_abs):
                ab_diag = _diag_lo <= a_ab_abs <= _diag_hi
                print(f"  A(Fe {ion_lbl}) NLTE A+B = {a_ab_abs:.4f}  -> {'PASS' if ab_diag else 'FAIL'} (diagnostic [{_diag_lo:.2f},{_diag_hi:.2f}])")

        vmic_val = converged_params.get('vturb_kms', np.nan)
        vmic_pass = 0.80 <= vmic_val <= 1.20 if np.isfinite(vmic_val) else False
        print(f"\n  vmic = {vmic_val:.3f} km/s  -> {'PASS' if vmic_pass else 'FAIL'} (gate 0.80-1.20)")
        print(f"\n  ►► Solar Fe VERDICT (scale-robust primary): {'PASS' if primary_pass else 'FAIL'}"
              f"  [slope {'✓' if slope_pass else '✗'} · ionization {'✓' if ion_pass else '✗'} · scatter {'✓' if scat_pass else '✗'}]")
        print(f"     1D-NLTE A(Fe) is the measurement; the reported solar anchor is the "
              f"3D-corrected A(Fe I) gated on FE_GATE (Magic 2013 applied, RYA-553).")

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
    # RYA-458: --ew-verify runs the per-line EW-integrity QA pass + charter cases.
    ew_verify = '--ew-verify' in _flags
    run(star, grid, engine=engine, skip_convergence=skip, ew_verify=ew_verify)
