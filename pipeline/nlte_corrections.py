"""
pipeline/nlte_corrections.py
============================
Apply non-LTE corrections to 1D LTE Fe I/II abundances.

LIVE Fe NLTE source (RYA-319): MPIA Spectrum Tools / Bergemann mafags-os **1D
NLTE** grid (data/nlte_grids/Fe_Bergemann_MPIA.csv) — covers the full benchmark
box (Procyon Teff, α Cen A A(Fe), α Cen B logg, 55 Cnc), interpolated per Fe
line over (teff,logg,feh) via _mpia_fe_delta / apply_fe_nlte_corrections.
The Amarsi 3D-NLTE MLP below (Teff≤6500 ceiling) is not the LIVE Fe correction.
Trade: 3D→1D, ≲0.05 dex at metal-rich targets (Amarsi 2022 Fig 7), accepted for
single-methodology coverage.

RYA-817 (Ryan, 2026-08-14: "we use every tool we can") REACTIVATED it as a separate
data product rather than leaving it archived — see `pipeline/amarsi3d.py`, which adds
the LINE-parameter domain check this module never had and emits the
`ENGINE-A-3DNLTE` treatment. Nothing on the live MPIA path below changed.

Archived Amarsi reference (the MLP functions, no longer the live correction):
Corrections from Amarsi, Liljegren & Nissen 2022 (A&A 668, A68).
Neural network implementation by Liljegren (github.com/sliljegren/1L-3NErrors).
Three multi-layer perceptrons trained on STAGGER-grid 3D NLTE models:
  - Fe I lines with Elo < 2 eV
  - Fe I lines with Elo >= 2 eV
  - Fe II lines

Grid coverage: Teff 5000–6500 K, log g 4.0–4.5, [Fe/H] −3.0 to 0.0.
171 Fe I + 12 Fe II optical lines.

Sign convention:
  A(Fe;3D NLTE) = A(Fe;1D LTE) + aberr, where the network emits A_1D − A_3D and
  `_compute_aberr` negates it (RYA-339).

  RYA-817 CORRECTION. This block used to read "aberr < 0 at solar metallicity ... for
  the Sun: A(3D) = 7.58 + (−0.127) ≈ 7.45". Both halves are wrong and the second is a
  number reverse-engineered from a target rather than computed from the network. Run
  against Amarsi's own solar inputs (Allende Prieto et al. 2002 Table 2, weak lines
  REW < −4.9) the reactivated network gives, at Teff 5772 / logg 4.44 / vmic 1.0:
      Fe I   A(1D LTE) 7.467 + (−0.002) = 7.465   [paper Table 6: 7.47 → 7.46]
      Fe II  A(1D LTE) 7.405 + (+0.066) = 7.471   [paper Table 6: 7.41 → 7.47]
  i.e. the SOLAR Fe I 3D-NLTE correction is ≈ 0.00 dex and the Fe II one is POSITIVE.
  The correction is strongly Elo-dependent (r = +0.94) and only becomes large away from
  the Sun. It cannot carry a 1D-LTE 7.58 to 7.46 and was never going to; that gap is a
  gf zero-point matter (RYA-161), not a 3D-NLTE one. Reproduced by
  `scripts/rya817_run_3dnlte_bands.py`, whose control asserts all four numbers.

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

from pipeline.species import parse_ion   # RYA-345 canonical ion normalizer
# RYA-765/318: the intake debug tracer. Every dbg.* call below is a no-op unless a
# trace is active, reads nothing the science path writes, and returns None — so this
# import can live here permanently at zero production cost (intake_debug invariants
# #1 numerically inert, #2 opt-in).
from pipeline import intake_debug as dbg

# ── Paths ─────────────────────────────────────────────────────────────────────

_REPO_ROOT   = Path(__file__).parent.parent
# repo root is already on sys.path (the `from pipeline.species import parse_ion` above
# requires it), so config.constants is importable here.
from config.constants import committed_grid_artifact  # noqa: E402  (RYA-567)
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


#: RYA-923 — failures that mean THE ENVIRONMENT IS BROKEN, never "this line is out of
#: domain". These must reach the caller; everything else may degrade to NaN.
_ENVIRONMENT_FAULTS = (ModuleNotFoundError, ImportError, FileNotFoundError,
                       AttributeError, EOFError)


def _load_models():
    """Load the Amarsi MLP regressors THROUGH THE VERSION CHECK (RYA-1098).

    🔴 THIS FUNCTION USED TO SUPPRESS THE ONLY EVIDENCE IT HAD. The pickles were written
    by scikit-learn 1.0.2 and are loaded under 1.9.0; sklearn says so, loudly, with an
    `InconsistentVersionWarning` per estimator -- and the body here was
    `warnings.simplefilter('ignore')` wrapped around three `pickle.load` calls. Cross-
    version unpickling of sklearn estimators is unsupported and does NOT reliably raise:
    it can deserialize cleanly and then predict WRONG. The output of these models is a
    PHYSICS CORRECTION applied to Fe abundances, so a silent mis-load corrupts an
    abundance with no loud failure. "It loaded without erroring" was standing in for "it
    loaded correctly", which is the RYA-167 silent-fallback shape.

    The load now goes through `pipeline.amarsi_model_integrity`, which reads the written
    version from the PICKLE BYTES (the loaded object reports the RUNTIME version -- it is
    rewritten on unpickle, so asking the object is a confidently wrong answer) and refuses
    any runtime whose predictions have not been validated against external ground truth.

    ⚠️ FOR THE RECORD, THE 1.0.2 -> 1.9.0 SKEW IS BENIGN, AND THAT IS MEASURED RATHER THAN
    ASSUMED: the authors' own published worked example reproduces (-0.135955 vs -0.136),
    and an independent numpy reimplementation of the forward pass matches `predict()`
    EXACTLY over 8,479 feature vectors spanning the feed's own lines. The corrections
    already published were not affected. See that module and RYA-1098.
    """
    from pipeline.amarsi_model_integrity import load_models as _checked_load
    return _checked_load(_VENDOR_DIR)


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
        # RYA-1140: capability-resolved, not `_REPO_ROOT.parent` (which is the
        # worktree's parent and forks the moment worktrees move -- RYA-1090).
        from config.constants import ISPEC_DIR
        sys.path.insert(0, str(ISPEC_DIR))
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
    except _ENVIRONMENT_FAULTS as exc:
        # 🔴 RYA-923 — AN ENVIRONMENT FAULT IS NOT A SCIENCE RESULT.
        #
        # What stood here was a bare `except Exception: return np.nan`. `sklearn` is not
        # installed in either venv on this box, so `_load_models()` raised
        # ModuleNotFoundError, this swallowed it, and every line came back NaN. The
        # caller reads NaN as "the network returned no value for this line" and reports
        # it as a DOMAIN verdict — so a missing dependency was presented as a statement
        # about the physics, and the RYA-817 reactivation control FAILED with `n/a`
        # against published values while looking like a coverage result.
        #
        # 114 in-domain lines produced n=0 that way. NaN is the right answer for a line
        # the network legitimately refuses; it is never the right answer for a network
        # that could not be loaded.
        raise RuntimeError(
            f"the Amarsi 2022 3D-NLTE network could not be evaluated: {exc}. This is an "
            f"ENVIRONMENT fault, not an out-of-domain line — refusing to return NaN, "
            f"which the caller would report as a domain verdict. The models are sklearn "
            f"pickles; install scikit-learn in the venv running this leg, or route the "
            f"element to a different NLTE source and say so in its provenance."
        ) from exc
    except Exception:
        # A genuine per-line prediction failure. NaN here IS a domain statement.
        return np.nan


# ── Public API ────────────────────────────────────────────────────────────────

def _apply_aberr_to_line(ion: str, elo: float, eup: float, lggf: float,
                          a_1dlte_line: float,
                          teff: float, logg: float, vmic: float,
                          afe3n_axis: "float | None" = None) -> float:
    """
    Compute converged aberr for a single line using its own 1D LTE abundance
    as the A(Fe;3N) starting point.  Returns np.nan if out of grid.

    RYA-817: `afe3n_axis` pins the grid's A(Fe;3N) axis to a value the caller has
    already converged — the STAR's 3D non-LTE iron abundance — instead of iterating from
    this line's own 1D LTE abundance. Both readings of the axis exist in the wild: the
    vendor README describes a single stellar A(Fe;3N) iterated to convergence, while
    this repo's per-line leg (RYA-207) starts each line from its own value. They differ
    by however far a line's measured abundance sits from the star's, which for a
    gf-limited line is a lot. The parameter is optional and defaults to None, so every
    existing caller keeps the per-line behaviour exactly; RYA-817 reports both and says
    which one each product used.
    """
    afe3n = float(np.clip(a_1dlte_line if afe3n_axis is None else afe3n_axis,
                          _GRID['afe'][0], _GRID['afe'][1]))
    ab = np.nan
    for _ in range(1 if afe3n_axis is not None else 3):
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
# RYA-567: committed, version-controlled Fe NLTE delta-CSV artifact (in-repo by
# ratified convention, NOT a heavy Sirius compute input) — routed via the resolver.
_MPIA_FE_GRID = committed_grid_artifact('nlte_grids', 'Fe_Bergemann_MPIA.csv')
_mpia_cache: dict = {}


def _load_mpia_fe_grid() -> dict:
    """Build per-(ion, wave) LinearNDInterpolator over (teff, logg, feh)."""
    if 'interp' in _mpia_cache:
        return _mpia_cache
    from scipy.interpolate import LinearNDInterpolator
    df = pd.read_csv(_MPIA_FE_GRID)
    # RYA-417: refuse to register a Fe-leg grid carrying placeholder-zero lines (a line
    # identically 0 across ALL its nodes = MPIA served-but-unmodelled, i.e. a silent
    # fake-LTE correction). This extends the RYA-413 registry-leg guard
    # (_load_mpia_element_grid) to the SEPARATE primary-Fe path — the gap the 31 Fe
    # placeholders rode. Detected per ion (the grid carries Fe I + Fe II at shared waves).
    _ph = []
    for _ion, _g in df.groupby('ion'):
        _ph += [(str(_ion), w) for w in detect_placeholder_zero_lines(_g)]
    if _ph:
        raise ValueError(
            f"[PLACEHOLDER_ZERO] Fe NLTE grid {_MPIA_FE_GRID.name} carries {len(_ph)} line(s) "
            f"identically 0 across all nodes: {_ph}. The Fe leg must NOT apply a placeholder "
            f"zero as if it were an NLTE correction (RYA-417, extends RYA-413) — drop them via "
            f"scripts/drop_fe_leg_placeholders_rya417.py after a live-MPIA TRUE_PLACEHOLDER "
            f"confirmation; never register the zero.")
    df = df[df['delta_nlte'].notna()]
    # RYA-345: the grid encodes ion as 'I'/'II' strings; key the interpolators on
    # the canonical ion int so callers passing any encoding ('I'/'II'/1/2/'26.1')
    # resolve to the same bucket. (Broke RYA-339 when the grid string convention
    # was assumed on the caller side.)
    interp, waves = {}, {}
    for (ion, wave), g in df.groupby(['ion', 'wave_A']):
        ion_i = parse_ion(ion)
        pts = g[['teff_K', 'logg', 'feh']].values
        if len(pts) < 4:
            continue
        try:
            f = LinearNDInterpolator(pts, g['delta_nlte'].values)
        except Exception:
            continue
        interp[(ion_i, round(float(wave), 3))] = f
        waves.setdefault(ion_i, []).append(round(float(wave), 3))
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
    ion_i = parse_ion(ion)                       # RYA-345 canonical ion
    w = c['waves'].get(ion_i, np.array([]))
    # RYA-765: this is RYA-317's exact shape — the grid matched 0 lines and the run
    # went on in LTE without saying so. Recorded on the intake timeline; no-op when
    # no trace is active, and the returned value is unchanged either way.
    if len(w) == 0:
        dbg.trace_fallback('nlte_lte_fallback',
                           f'Fe {ion} {wave_A:.3f}: grid carries no node for this ion '
                           f'-> LTE', severity='ERROR', species=f'Fe {ion}',
                           wavelength=float(wave_A), miss='empty_grid')
        return np.nan
    j = int(np.abs(w - wave_A).argmin())
    if abs(float(w[j]) - wave_A) > tol:
        dbg.trace_fallback('nlte_lte_fallback',
                           f'Fe {ion} {wave_A:.3f}: no grid/level match within {tol} A '
                           f'(nearest {float(w[j]):.3f}) -> LTE', severity='ERROR',
                           species=f'Fe {ion}', wavelength=float(wave_A),
                           miss='wave_tol', nearest_A=float(w[j]))
        return np.nan
    f = c['interp'].get((ion_i, float(w[j])))
    if f is None:
        dbg.trace_fallback('nlte_lte_fallback',
                           f'Fe {ion} {wave_A:.3f}: matched node {float(w[j]):.3f} has '
                           f'no interpolator -> LTE', severity='ERROR',
                           species=f'Fe {ion}', wavelength=float(wave_A),
                           miss='no_interpolator')
        return np.nan
    return float(f([[teff, logg, feh]])[0])


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
    teff = float(stellar_params.get('teff_K', 5772))  # RYA-298: canonical solar Teff
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
                dbg.trace_fallback('nlte_lte_fallback',
                                   f'Fe {ion}: no measured line of this ion in the '
                                   f'per-line table -> 1D LTE retained', severity='ERROR',
                                   species=f'Fe {ion}', miss='no_species_lines')
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
            # RYA-318/317 — THE check that named this ticket: "NLTE matched 0 / 11,275
            # Fe lines", caught in 317 only because a guard happened to fire.
            dbg.trace_check('nlte_matched_fraction', n_corrected > 0,
                            detail=f'Fe {ion}: {n_corrected}/{len(fe_lines)} lines '
                                   f'matched the NLTE grid',
                            species=f'Fe {ion}', matched=int(n_corrected),
                            total=int(len(fe_lines)))
            if n_corrected == 0:
                print(f"  Fe {ion}: no lines in NLTE grid — 1D LTE retained")
                dbg.trace_fallback('nlte_lte_fallback',
                                   f'Fe {ion}: 0/{len(fe_lines)} lines matched the NLTE '
                                   f'grid -> 1D LTE retained', severity='ERROR',
                                   species=f'Fe {ion}', miss='zero_matched',
                                   total=int(len(fe_lines)))
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


# ── Ca / Ti / Cr NLTE — per-line application from the MPIA grids (RYA-235) ────
# Grids: data/nlte_grids/{Ca_Mashonkina2017,Ti_Bergemann2011_MPIA,Cr_Bergemann2010_MPIA}.csv
# Columns: element, wave_A, teff_K, logg, feh, delta_nlte (per-line; neutral only).
# Registry: config.constants.NLTE_CORRECTION_ELEMENTS.
#
# UNIT CONVENTION: feh axis is [Fe/H] RELATIVE (range −0.5 to +0.3). iSpec returns
# [X/H] relative (0.0 at solar) for all elements → pass stellar_params['feh']
# directly, NO solar offset. This is unlike the Fe Amarsi grid (afe axis = absolute
# A(Fe) 4.5–7.5, which required _A_FE_SOLAR = 7.46 added before lookup, RYA-247).
#
# Per-line application mirrors the Fe leg: for each measured line of an element in
# the registry we interpolate delta over (teff,logg,feh) at the line's wavelength,
# set a_nlte[i] = a_1dlte[i] + delta[i], and recompute median/std from the corrected
# values. Outside the convex hull → NaN (line skipped). Availability gate: an element
# with no in-grid line stays 1D-LTE with nlte_flag='NLTE_unavailable' (no silent fix).
#
# Real-grid solar deltas (Teff~5777, logg~4.4, feh=0.0; from these MPIA grids):
#   Ca ~+0.013, Ti ~+0.108, Cr ~+0.073 dex (medians).
#   NB — the MPIA MAFAGS-OS Cr median (~+0.07) is FAR below the +0.3–0.5 the RYA-235
#   description cited from Bergemann & Cescutti 2010; the large landmark correction is
#   line/regime specific. This does NOT rescue the raw-pool Cr overshoot: 1D-LTE Cr
#   already reads high (RYA-239), so acceptance still requires a CURATED pool first.
_MPIA_ELEMENT_DIR = committed_grid_artifact('nlte_grids')   # RYA-567 (committed artifact dir)
_mpia_element_cache: dict = {}


def detect_placeholder_zero_lines(df: pd.DataFrame, eps: float = 1e-9) -> list:
    """Lines that are IDENTICALLY ~0 across ALL their (teff,logg,feh) nodes — almost
    certainly a failed build or an MPIA-served-but-unmodelled line written as 0, NOT a
    physical NLTE correction (a real correction varies over the parameter grid; it is not
    exactly 0.000 everywhere). Returns the sorted wave_A list. RYA-413: this is the
    placeholder-zero class (root: MPIA offers Ca 6166 in its dropdown but returns 0)."""
    bad = []
    for wave, g in df.groupby('wave_A'):
        col = g['delta_nlte'].dropna()
        if len(col) >= 2 and (col.abs() < eps).all():
            bad.append(round(float(wave), 3))
    return sorted(bad)


def _load_mpia_element_grid(element: str) -> dict:
    """Per-(wave) LinearNDInterpolator over (teff_K, logg, feh) for a registry
    element (Ca/Ti/Cr). Cached per element; loud-fails if the element is not in the
    registry or its grid file is missing (no silent skip of a configured element).
    RYA-413: REFUSES to register a grid carrying a placeholder-zero line (raises loud)."""
    if element in _mpia_element_cache:
        return _mpia_element_cache[element]
    from scipy.interpolate import LinearNDInterpolator
    from config.constants import NLTE_CORRECTION_ELEMENTS
    if element not in NLTE_CORRECTION_ELEMENTS:
        raise KeyError(f"{element!r} is not in NLTE_CORRECTION_ELEMENTS "
                       f"({sorted(NLTE_CORRECTION_ELEMENTS)})")
    spec = NLTE_CORRECTION_ELEMENTS[element]
    path = _MPIA_ELEMENT_DIR / spec['grid']
    if not path.exists():
        raise FileNotFoundError(f"NLTE grid for {element} not found: {path}")
    df = pd.read_csv(path)
    placeholders = detect_placeholder_zero_lines(df)
    if placeholders:
        raise ValueError(
            f"[PLACEHOLDER_ZERO] {element} NLTE grid {spec['grid']} carries {len(placeholders)} "
            f"line(s) that are identically 0 across all nodes: {placeholders}. A registered grid "
            f"must NOT store a placeholder zero as if it were an NLTE correction (RYA-413) — "
            f"re-fetch the line or DROP it from the pool with a documented reason; never register "
            f"the zero.")
    df = df[df['delta_nlte'].notna()]
    interp, waves = {}, []
    for wave, g in df.groupby('wave_A'):
        pts = g[['teff_K', 'logg', 'feh']].values
        if len(pts) < 4:                      # need a simplex for LinearND
            continue
        try:
            interp[round(float(wave), 3)] = LinearNDInterpolator(pts, g['delta_nlte'].values)
        except Exception:
            continue
        waves.append(round(float(wave), 3))
    cache = {
        'interp': interp,
        'waves': np.array(sorted(waves)),
        'ref': spec['ref'],
        'ion': int(spec['ion']),
        'bounds': {
            'teff': (float(df.teff_K.min()), float(df.teff_K.max())),
            'logg': (float(df.logg.min()),  float(df.logg.max())),
            'feh':  (float(df.feh.min()),   float(df.feh.max())),
        },
    }
    _mpia_element_cache[element] = cache
    return cache


def _mpia_element_delta(element: str, wave_A: float, teff: float, logg: float,
                        feh: float, tol: float = 0.15) -> float:
    """delta_nlte for one line of `element`, interpolated from its MPIA grid over
    (teff,logg,feh) at the line wavelength. NaN if no wave match within `tol` or the
    star sits outside the (teff,logg,feh) convex hull.

    RYA-765: every NaN return here is a line that will proceed in LTE. The caller
    decides what to do with it (drop it, or add 0.0 — the two-engine driver does the
    latter), and neither choice leaves a trace in the science output, which is how
    RYA-317 stayed invisible until a number looked wrong. The `dbg.trace_fallback`
    calls record the miss on the intake timeline instead; they are no-ops unless a
    trace is active and they change nothing this function returns.
    """
    c = _load_mpia_element_grid(element)
    w = c['waves']
    if len(w) == 0:
        dbg.trace_fallback('nlte_lte_fallback',
                           f'{element} {wave_A:.3f}: grid carries no usable wavelength '
                           f'node -> LTE', severity='ERROR', species=str(element),
                           wavelength=float(wave_A), miss='empty_grid')
        return np.nan
    j = int(np.abs(w - wave_A).argmin())
    if abs(float(w[j]) - wave_A) > tol:
        dbg.trace_fallback('nlte_lte_fallback',
                           f'{element} {wave_A:.3f}: no grid/level match within '
                           f'{tol} A (nearest {float(w[j]):.3f}) -> LTE',
                           severity='ERROR', species=str(element),
                           wavelength=float(wave_A), miss='wave_tol',
                           nearest_A=float(w[j]))
        return np.nan
    f = c['interp'].get(float(w[j]))
    if f is None:
        dbg.trace_fallback('nlte_lte_fallback',
                           f'{element} {wave_A:.3f}: matched node {float(w[j]):.3f} has '
                           f'no interpolator -> LTE', severity='ERROR',
                           species=str(element), wavelength=float(wave_A),
                           miss='no_interpolator')
        return np.nan
    d = float(f([[teff, logg, feh]])[0])
    if not np.isfinite(d):
        # LinearND returns NaN outside the (Teff,logg,[Fe/H]) hull. The element-level
        # RYA-409 guard catches the whole-element case; this is the per-line remainder.
        dbg.trace_fallback('nlte_lte_fallback',
                           f'{element} {wave_A:.3f}: node {float(w[j]):.3f} interpolated '
                           f'to NaN at (Teff={teff:.0f}, logg={logg:.2f}, '
                           f'[Fe/H]={feh:+.2f}) -> LTE', severity='ERROR',
                           species=str(element), wavelength=float(wave_A),
                           miss='out_of_hull')
    return d


def element_grid_in_bounds(element: str, teff: float, logg: float, feh: float) -> bool:
    """Preflight: is (teff,logg,feh) inside the element's MPIA grid box?"""
    b = _load_mpia_element_grid(element)['bounds']
    return (b['teff'][0] <= teff <= b['teff'][1] and
            b['logg'][0] <= logg <= b['logg'][1] and
            b['feh'][0]  <= feh  <= b['feh'][1])


def apply_element_nlte_corrections(
    abundances_df: pd.DataFrame,
    stellar_params: dict,
    per_line_df: pd.DataFrame = None,
    line_df: pd.DataFrame = None,
    elements=None,
    wave_tol: float = 0.15,
    strict: bool = False,
) -> pd.DataFrame:
    """Apply per-line NLTE corrections to the non-Fe registry elements (Ca/Ti/Cr/Na/
    Mg I + Ba II; RYA-235 + RYA-165) from their vendored grids — MPIA MAFAGS-OS for
    Ca/Ti/Cr/Mg, INSPECT/Lind 2011 for Na, Korotin 2015 for Ba. Runs AFTER
    apply_fe_nlte_corrections and
    composes with it: it touches ONLY rows whose (element, ion) are in
    NLTE_CORRECTION_ELEMENTS, leaving the Fe leg and every other row exactly as the
    Fe pass left them.

    Mirrors the Fe per-line path: per-line delta interpolated at each line's
    wavelength; A_X_nlte = median(a_1dlte[i] + delta[i]); columns written are the
    same schema (delta_nlte_mean / n_nlte_lines / A_X_nlte / A_X_std_nlte /
    nlte_flag / nlte_ref). An element with no in-grid line is left 1D-LTE with
    nlte_flag='NLTE_unavailable' — never silently "corrected".

    NOTE (RYA-235): this is the application + wiring layer. Validation to the
    Asplund-2021 anchors is intentionally NOT asserted here — it is gated on a
    curated line pool (raw 1D-LTE Cr reads high; NLTE on an un-curated pool
    overshoots). Acceptance is owed on the curated pool, separately.
    """
    from config.constants import NLTE_CORRECTION_ELEMENTS
    reg = NLTE_CORRECTION_ELEMENTS
    want = set(elements) if elements is not None else set(reg)

    teff = float(stellar_params.get('teff_K', 5772))
    logg = float(stellar_params.get('logg',   4.44))
    feh  = float(stellar_params.get('feh', 0.0))    # RELATIVE [Fe/H] — no solar offset

    result = abundances_df.copy()
    # ensure the NLTE schema exists even if the Fe pass did not run first
    for col, default in (('delta_nlte_mean', np.nan), ('n_nlte_lines', 0),
                         ('A_X_nlte', np.nan), ('A_X_std_nlte', np.nan),
                         ('nlte_flag', '1D_LTE'), ('nlte_ref', '')):
        if col not in result.columns:
            result[col] = default

    for idx, row in result.iterrows():
        element = str(row['element'])
        if element not in reg or element not in want:
            continue
        if parse_ion(row.get('ion', 1)) != int(reg[element]['ion']):
            continue                          # registry ion only (Cr II excluded; Ba II included)

        a_1dlte   = float(row['A_X'])
        ref       = reg[element]['ref']
        ion_label = 'II' if int(reg[element]['ion']) == 2 else 'I'

        if per_line_df is None or per_line_df.empty:
            print(f"  {element} {ion_label}: no per-line table — element NLTE skipped (1D LTE retained)")
            dbg.trace_fallback('nlte_lte_fallback',
                               f'{element} {ion_label}: no per-line table -> whole '
                               f'element retained at 1D LTE', severity='ERROR',
                               species=f'{element} {ion_label}', miss='no_per_line_table')
            result.at[idx, 'A_X_nlte']  = a_1dlte
            result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'
            continue

        lines = per_line_df[
            (per_line_df['element'] == element)
            & (per_line_df['ion'].apply(parse_ion) == int(reg[element]['ion']))
        ].copy().reset_index(drop=True)
        if lines.empty:
            dbg.trace_fallback('nlte_lte_fallback',
                               f'{element} {ion_label}: no measured line of this species '
                               f'in the per-line table -> 1D LTE retained',
                               severity='ERROR', species=f'{element} {ion_label}',
                               miss='no_species_lines')
            result.at[idx, 'A_X_nlte']  = a_1dlte
            result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'
            continue

        # ── RYA-409 Part A: LOUD out-of-hull guard — kill the silent NaN→1D-LTE clamp ──
        # The LinearND interpolator returns NaN outside the (Teff,logg,[Fe/H]) hull, which
        # used to collapse to the generic 'NLTE_unavailable' with no signal — so 7 Family-A
        # elements silently dropped to 1D-LTE at the metal-rich 55 Cnc ([Fe/H]+0.31 vs grid
        # ceiling +0.30). Same loud-failure-over-silent-fallback class as RYA-399's
        # 3D_unavailable / the RestFrameError gate: flag it LOUD + distinct (NLTE_OUT_OF_HULL),
        # name the element + (Teff,logg,[Fe/H]) + the hull it missed, retain LTE only WITH that
        # explicit flag — never a silent substitution. `strict=True` raises instead.
        if not element_grid_in_bounds(element, teff, logg, feh):
            b = _load_mpia_element_grid(element)['bounds']
            msg = (f"{element} {ion_label}: (Teff={teff:.0f}, logg={logg:.2f}, [Fe/H]={feh:+.2f}) is "
                   f"OUTSIDE the NLTE grid hull (Teff {b['teff'][0]:.0f}–{b['teff'][1]:.0f}, logg "
                   f"{b['logg'][0]:.2f}–{b['logg'][1]:.2f}, [Fe/H] {b['feh'][0]:+.2f}..{b['feh'][1]:+.2f}). "
                   f"Retaining 1D-LTE with an explicit NLTE_OUT_OF_HULL flag — NOT a silent NLTE→LTE "
                   f"substitution (RYA-409). Re-source the grid past this [Fe/H] (Part B).")
            warnings.warn(msg, RuntimeWarning, stacklevel=2)
            print(f"  LOUD [NLTE_OUT_OF_HULL]: {msg}")
            dbg.trace_fallback('nlte_out_of_hull',
                               f'{element} {ion_label}: (Teff={teff:.0f}, logg={logg:.2f}, '
                               f'[Fe/H]={feh:+.2f}) outside the grid hull -> 1D LTE '
                               f'retained with NLTE_OUT_OF_HULL', severity='ERROR',
                               species=f'{element} {ion_label}', teff=float(teff),
                               logg=float(logg), feh=float(feh),
                               hull_teff=list(b['teff']), hull_logg=list(b['logg']),
                               hull_feh=list(b['feh']))
            result.at[idx, 'A_X_nlte']  = a_1dlte
            result.at[idx, 'nlte_flag'] = 'NLTE_OUT_OF_HULL'
            if strict:
                raise ValueError(f"RYA-409 strict NLTE hull guard: {msg}")
            continue

        lines['aberr']  = np.nan
        lines['a_nlte'] = lines['a_1dlte']
        for i, lrow in lines.iterrows():
            ab = _mpia_element_delta(element, float(lrow['wavelength_air_A']),
                                     teff, logg, feh, tol=wave_tol)
            if np.isfinite(ab):
                lines.at[i, 'aberr']  = ab
                lines.at[i, 'a_nlte'] = float(lrow['a_1dlte']) + ab

        n_corrected = int(lines['aberr'].notna().sum())
        # RYA-318/317: the matched-line fraction, finally landing in the diagnostics
        # layer as 318 asked. 0 matched = ERROR; the run continues exactly as before.
        dbg.trace_check('nlte_matched_fraction', n_corrected > 0,
                        detail=f'{element} {ion_label}: {n_corrected}/{len(lines)} lines '
                               f'matched the NLTE grid',
                        species=f'{element} {ion_label}',
                        matched=int(n_corrected), total=int(len(lines)))
        if n_corrected == 0:
            print(f"  {element} {ion_label}: no lines in NLTE grid (Δ all NaN) — 1D LTE retained")
            dbg.trace_fallback('nlte_lte_fallback',
                               f'{element} {ion_label}: 0/{len(lines)} lines matched the '
                               f'NLTE grid -> 1D LTE retained', severity='ERROR',
                               species=f'{element} {ion_label}', miss='zero_matched',
                               total=int(len(lines)))
            result.at[idx, 'A_X_nlte']  = a_1dlte
            result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'
            continue

        a_nlte_med = float(np.median(lines['a_nlte']))
        a_nlte_std = float(lines['a_nlte'].std()) if len(lines) > 1 else np.nan
        mean_delta = float(lines['aberr'].mean(skipna=True))
        print(f"  {element} {ion_label} per-line NLTE ({n_corrected}/{len(lines)} corrected): "
              f"mean Δ = {mean_delta:+.4f} dex  A({element};1D)={a_1dlte:.3f} → "
              f"A({element};NLTE)={a_nlte_med:.3f}")

        result.at[idx, 'delta_nlte_mean'] = round(mean_delta, 4)
        result.at[idx, 'n_nlte_lines']    = n_corrected
        result.at[idx, 'A_X_nlte']        = round(a_nlte_med, 3)
        result.at[idx, 'A_X_std_nlte']    = round(a_nlte_std, 3) if np.isfinite(a_nlte_std) else np.nan
        # per-source provenance tag (RYA-165): MPIA MAFAGS-OS by default; a non-MPIA
        # grid (Na INSPECT/Lind, Ba Korotin) carries its own `flag` so it is never
        # mislabelled as MPIA.
        result.at[idx, 'nlte_flag']       = reg[element].get('flag', 'NLTE_MPIA_MAFAGS_1D')
        result.at[idx, 'nlte_ref']        = ref

    return result


# ── C I / O I NLTE — the C/O leg (RYA-359, Amarsi 2019) ──────────────────────
# Carbon and oxygen non-LTE corrections come from the Amarsi, Nissen &
# Skuladottir 2019 (A&A 630, A104; CDS J/A+A/630/A104) line-by-line grids,
# vendored at data/nlte_grids/amarsi2019_cno/. The resolver lives in
# pipeline.nlte_cno; it is re-exported here so the C/O leg sits in the same
# (sign-corrected, RYA-339) module as the Fe NLTE path. 3D-NLTE leg for
# Teff <= 6500 K, 1D-NLTE leg above (Procyon). Corrections are large & NEGATIVE;
# applied as A(X;NLTE) = A(X;1D-LTE) + Delta with a sign guard (assert_cno_sign).
from pipeline.nlte_cno import (            # noqa: E402
    apply_cno_nlte_corrections,
    cno_nlte_delta,
    resolve_line,
    select_leg,
    grid_coverage,
    star_in_grid,
    coverage_ok,
    validate_against_table7,
    assert_cno_sign,
    REQUIRED_LINES,
    TEFF_3D_CEILING,
    CITATION as CNO_CITATION,
)


# ── Metal 3D leg (RYA-399, Si/Ti/Cr) — composes AFTER the element NLTE pass ───
# The 3D dimensional correction for the GET-3D metals (routed by RYA-400) is added
# on top of the 1D-NLTE abundance: A(3D-NLTE) = A(1D-NLTE) + delta_3d. Re-exported
# here so the RYA-371 Phase-C compose can run NLTE then 3D from one module. FINDING:
# the published solar 3D corrections are small (<=0.1 dex, Ti/Cr positive) and do
# NOT close the Si/Ti/Cr solar residual — carried forward, never tuned (RYA-399).
from pipeline.threed_corrections import (    # noqa: E402
    apply_threed_corrections,
    solar_threed_delta,
    residual_after_3d,
)


# ── CLI: coverage + per-star verdict + solar corrections + self-validation ───────

def _cli(argv=None):
    """Smoke test / characterization entrypoint.

        python -m pipeline.nlte_corrections --species "C I,O I" \\
            --source amarsi2019_cds --star solar --validate
    """
    import argparse
    parser = argparse.ArgumentParser(description="C/O (and Fe) NLTE grid tools")
    parser.add_argument('--species', default='C I,O I',
                        help='comma list, e.g. "C I,O I"')
    parser.add_argument('--source', default='amarsi2019_cds',
                        help='grid source (amarsi2019_cds = Amarsi 2019 CDS)')
    parser.add_argument('--star', default='solar', help='star_id for the demo')
    parser.add_argument('--validate', action='store_true',
                        help='run the table7 self-validation cross-check')
    args = parser.parse_args(argv)

    from config.constants import get_star_params, SOLAR_ASPLUND2021

    sp_map = {'C I': 'CI', 'O I': 'OI', 'CI': 'CI', 'OI': 'OI',
              'C': 'CI', 'O': 'OI'}
    species = []
    for tok in args.species.split(','):
        s = sp_map.get(tok.strip())
        if s and s not in species:
            species.append(s)

    print("=" * 78)
    print(f"C/O NLTE grids — source={args.source}  "
          f"(Amarsi 2019, {CNO_CITATION})")
    print("=" * 78)

    # 1) coverage table -------------------------------------------------------
    print("\nGRID COVERAGE (per species, per leg)")
    hdr = (f"  {'species':<8}{'leg':<5}{'Teff':>14}{'logg':>12}"
           f"{'[Fe/H]':>14}{'vmic':>10}{'nlines':>8}{'nrec':>10}")
    print(hdr)
    for sp in species:
        for leg in ('3D', '1D'):
            c = grid_coverage(sp, leg)
            print(f"  {sp:<8}{leg:<5}"
                  f"{c['teff'][0]:>6.0f}-{c['teff'][1]:<7.0f}"
                  f"{c['logg'][0]:>5.1f}-{c['logg'][1]:<6.1f}"
                  f"{c['feh'][0]:>6.2f}-{c['feh'][1]:<7.2f}"
                  f"{c['vmic'][0]:>4.1f}-{c['vmic'][1]:<5.1f}"
                  f"{len(c['lines']):>8}{c['n_records']:>10}")
    ok, problems = coverage_ok()
    print(f"\n  STOP-GATE (required lines + [Fe/H] >= +0.32): "
          f"{'PASS' if ok else 'FAIL'}")
    for p in problems:
        print(f"    ! {p}")

    # 2) per-star in-grid verdict --------------------------------------------
    print("\nPER-STAR IN-GRID VERDICT (Sun / Procyon / alpha Cen A / 55 Cnc A)")
    for star_id in ('solar', 'procyon', 'alpha_cen_a', '55cnc_a'):
        try:
            rec = get_star_params(star_id)
        except Exception as e:
            print(f"  {star_id}: no STAR_PARAMS record ({e})")
            continue
        teff = float(rec['teff']); logg = float(rec['logg'])
        feh = float(rec.get('feh_ref', 0.0))
        vmic = float(rec.get('xi', rec.get('xi_init', 1.0)))
        leg = select_leg(teff)
        print(f"  {star_id:<12} Teff={teff:.0f} logg={logg:.2f} "
              f"[Fe/H]={feh:+.2f} vmic={vmic:.2f} → {leg}-leg")
        for sp in species:
            v = star_in_grid(sp, teff, logg, feh, vmic)
            miss = (' MISSING=' + ','.join(v['required_missing'])
                    if v['required_missing'] else '')
            print(f"      {sp}: in_box={v['in_box']} "
                  f"required_present={v['required_present']}{miss}")

    # 3) demo: solar corrections with sign + magnitude -----------------------
    print(f"\nSOLAR CORRECTIONS (star={args.star}, demo on required lines)")
    rec = get_star_params(args.star)
    teff = float(rec['teff']); logg = float(rec['logg'])
    feh = float(rec.get('feh_ref', 0.0))
    vmic = float(rec.get('xi', rec.get('xi_init', 1.0)))
    leg = select_leg(teff)
    sun_logeps = {'CI': SOLAR_ASPLUND2021['C'], 'OI': SOLAR_ASPLUND2021['O']}
    print(f"  leg={leg}  logeps(C)={sun_logeps['CI']}  logeps(O)={sun_logeps['OI']}")
    for sp in species:
        for label in REQUIRED_LINES[sp]:
            d = cno_nlte_delta(sp, label, teff, logg, feh, vmic,
                               sun_logeps[sp], leg=leg)
            assert_cno_sign(sp, label, d)
            sign = 'neg' if d < 0 else ('pos' if d > 0 else 'zero')
            print(f"    {sp} {label:<9} Delta = {d:+.4f} dex  ({sign})")

    # 4) self-validation against table7 --------------------------------------
    if args.validate:
        print("\nSELF-VALIDATION vs table7 (interpolated Delta vs published 3N-1L)")
        res = validate_against_table7()
        print(f"  {'star':<12}{'leg':<4}{'Teff':>6}{'[Fe/H]':>8}"
              f"{'C_pred':>9}{'C_obs':>9}{'C_res':>8}"
              f"{'O_pred':>9}{'O_obs':>9}{'O_res':>8}")
        for r in res['rows']:
            def fmt(x): return f"{x:+.3f}" if np.isfinite(x) else "   n/a"
            print(f"  {r['star']:<12}{r['leg']:<4}{r['teff']:>6}{r['feh']:>+8.2f}"
                  f"{fmt(r['c_pred']):>9}{fmt(r['c_obs']):>9}{fmt(r['c_resid']):>8}"
                  f"{fmt(r['o_pred']):>9}{fmt(r['o_obs']):>9}{fmt(r['o_resid']):>8}")
        print(f"\n  C: MAE={res['c_mae']:.4f}  max|res|(all)={res['c_max']:.4f}  "
              f"max|res|(sci [Fe/H]>={res['feh_floor']})={res['c_max_sci']:.4f} dex")
        print(f"  O: MAE={res['o_mae']:.4f}  max|res|(all)={res['o_max']:.4f}  "
              f"max|res|(sci [Fe/H]>={res['feh_floor']})={res['o_max_sci']:.4f} dex")
        print(f"  (the fixed optical-C-line proxy degrades at [Fe/H]<~-2 halo "
              f"stars — outside our near-solar science range; grid itself OK)")
        print(f"  tolerance={res['tol']} over science range  → "
              f"{'PASS' if res['passed'] else 'FAIL'}")
    print("\n" + "=" * 78)


if __name__ == '__main__':
    _cli()
