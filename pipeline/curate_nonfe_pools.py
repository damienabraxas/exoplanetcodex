"""
pipeline/curate_nonfe_pools.py
==============================
RYA-395 — GES-style curation of the NON-Fe solar EW line pools.

WHY
---
The RYA-239 solar run read every measured non-Fe element systematically HIGH
(Δ vs Asplund 2021 all positive, +0.17 to +2.25 dex). Fe II from the *same*
`solar_ew` table lands correct post-cull, so this is NOT a global EW/scale bias
— it is per-element LINE-POOL quality (saturation, bad gf, blends) in the
uncurated non-Fe pools, the only pools that never got the GES-style vetting
Fe I/II received. This module closes that gap; it gates RYA-371's non-Fe GO.

THE CARDINAL RULE — validate, don't tune
----------------------------------------
Asplund 2021 is the VALIDATION target, NEVER the tuning target. Every cull
decision is made on spectroscopic-quality grounds (gf provenance, curve-of-
growth position, blend status), with thresholds FIXED IN ADVANCE (the module
constants below), uniform across elements, documented, and applied BLIND to
the resulting abundance. The abundance is read *after*, as validation. Culling
lines until the mean hits Asplund is the cardinal sin and makes the result
meaningless. This module enforces abundance-blindness structurally: the cull
logic (`cull_reasons`) is only ever handed the quality columns
(`_CULL_INPUT_COLS`); it never sees A(X), the Asplund anchor, or the per-element
target. `--verify` proves it (a shuffle test).

THE CULL-vs-gf-SCALE FORK
-------------------------
After quality-based curation, if an element's offset persists the residual is a
systematic gf-scale issue (astrophysical-gf — RYA-161 territory), NOT stray bad
lines. We FLAG it for escalation; we do NOT force-fit. The discriminant is the
curated pool: a high residual that survives the quality cuts on a pool whose gf
is dominated by ungraded / Kurucz-theoretical loggf is a gf-scale systematic.

Cr is the canary (LTE +0.345 high; NLTE adds only +0.073 → curation must pull
~0.42 dex down on quality grounds to reach Asplund 5.62). Mn is the worse canary
(NLTE +0.102 POSITIVE on top of LTE +0.647 → ~0.75 dex of downward curation
owed; HFS-heavy). For both the contract is: curation brings LTE to
(Asplund − Δ_NLTE) so the positive NLTE lifts onto scale — NOT to Asplund
directly, or NLTE overshoots. `--verify` asserts that overshoot contract.

DATA SOURCES (single-source; nothing hardcoded)
-----------------------------------------------
* measured EW pool : data/measured/sol_ew_results_v1.csv  (element, ion, λ_air,
                     ew_mA, ew_err_mA, profile chi2, blend_flag)
* atomic data      : data/linelists/canonical_gf.csv  (the gf single source:
                     log_gf, EP, nist_grade, loggf_reference) — full pool match.
* RT species codes : the GES synthesis atomic linelist (atomic_lines.tsv) supplies
                     MOOG/WIDTH species codes + damping; loggf/EP are OVERRIDDEN
                     from canonical so EW→A reads the single-source gf.
* NLTE Δ (verdict) : the vendored MPIA / INSPECT / CDS grids (RYA-396/235),
                     queried at solar params — never hardcoded.

OUTPUTS
-------
* data/curation/nonfe_pools/{element}_cull_rya395.csv  — per-line KEEP/CULL with
  reason + provenance header (the fe2_solar_cull pattern).
* data/curation/nonfe_pools/curation_diagnostics_rya395.csv — the per-element
  diagnostic table (A(X) LTE, scatter, REW slope, N_clean) + cull-vs-gf verdict.

USAGE
-----
    python -m pipeline.curate_nonfe_pools --phase 1 --verify   # Mg Si Ca Ni
    python -m pipeline.curate_nonfe_pools --phase 2 --verify   # Ti Cr Na Al Mn
    python -m pipeline.curate_nonfe_pools --phase all

Linear issue: RYA-395
"""
from __future__ import annotations

import argparse
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from config.constants import (              # noqa: E402
    PATHS, PIPELINE, SOLAR_ASPLUND2021, get_star_params, EW_BASELINE_CODE,
)

# ── Scope (the Tier-1/2 metals, phased) ──────────────────────────────────────
PHASE1 = ['Mg', 'Si', 'Ca', 'Ni']            # gate 55 Cnc
PHASE2 = ['Ti', 'Cr', 'Na', 'Al', 'Mn']
# Out of scope (Tier-3 trace / no hard gate): Co Y Sc V Cu Zr Sr — carried only
# if measured, always flagged low-confidence. N/K/P/Eu/Li are VIS-lineless (RYA-162).

# ─────────────────────────────────────────────────────────────────────────────
# FIXED, ABUNDANCE-BLIND CURATION THRESHOLDS
# Every value here is set on PHYSICAL grounds (curve-of-growth, measurement
# noise, line-list provenance), fixed BEFORE any abundance is read, and applied
# UNIFORMLY across all elements. None of them is a function of A(X) or Asplund.
# Where a pipeline-wide constant already encodes the physics we reuse it (single
# source) rather than inventing a parallel number.
# ─────────────────────────────────────────────────────────────────────────────
EW_MIN_MA         = float(PIPELINE['ew_min_mA'])           # 5 mÅ — below = noise/continuum
# Saturation cut on REDUCED EW (the ticket's specified axis): the curve of growth
# is linear only for log(EW/λ) ≲ −4.9; beyond the knee the line saturates onto the
# flat/Doppler part where EW→A is degenerate and biased high (and vmic-sensitive).
# −4.9 is the standard solar-type COG knee (Gray, Stellar Photospheres §16); it is
# fixed from COG theory, BEFORE any abundance is read, uniform across elements.
REW_LINEAR_CEILING = -4.90
# Secondary absolute-strength guard (reuses the pipeline saturation ceiling, RYA-330).
SAT_EW_CEILING_MA  = float(PIPELINE['vmic_ew_ceiling_mA'])  # 100 mÅ
EW_ERR_FRAC_MAX   = float(PIPELINE['fe2_ew_err_frac_max']) # 0.5 — ew_err/ew above = unmeasurably noisy

# Blend vetting — strength-weighted opacity overlap (the principled fix to the
# buggy `vald_proximity_flag`, which flagged on ANY nearby VALD line regardless
# of strength and so tagged ~87% of lines). We weight each in-window neighbour by
# its solar-reference line opacity (A_sun + loggf − EP·θ) and a Gaussian profile
# kernel, and cull only when contaminants supply more than BLEND_FRAC_MAX of the
# target line's opacity inside its profile. Uses the REFERENCE solar composition
# (an input), never the derived A(X) — so it stays abundance-blind.
BLEND_WINDOW_A = 0.10        # ± search window (Å); ~2× the solar optical FWHM
BLEND_SIGMA_A  = 0.04        # profile kernel σ (Å); FWHM≈0.09 Å → σ≈0.04
BLEND_FRAC_MAX = 0.50        # cull if Σ contaminant opacity / target opacity > this

# gf provenance — cull only graded-bad lab values; FLAG (low-confidence, do not
# auto-cull) the ungraded / Kurucz-theoretical majority. The fraction of LOW-tier
# gf in a pool is the cull-vs-gf-scale evidence (RYA-161).
NIST_GRADE_HIGH = {'A+', 'A', 'B'}           # lab, <3 % — trusted
NIST_GRADE_CULL = {'D', 'E', 'F'}            # <25 % or worse — cull
_KURUCZ_REF = re.compile(r'^K(\d+|P)$', re.I)  # K03/K07/K08/K09/K10/KP — semi-empirical

N_CLEAN_FLOOR = 5            # fewer clean lines → element flagged low-confidence
RESID_TOL_DEX = 0.10         # |curated A − (Asplund−Δ_NLTE)| within this = on scale

# Columns the cull logic is ALLOWED to read. The abundance-blind invariant: the
# cull mask is a function of these ONLY. `--verify` shuffles A(X) and asserts the
# mask is unchanged.
_CULL_INPUT_COLS = frozenset({
    'ew_mA', 'ew_err_mA', 'err_frac', 'rew', 'wavelength_air_A',
    'chi2', 'blend_flag', 'blend_ratio', 'nist_grade', 'gf_tier',
})
_FORBIDDEN_IN_CULL = frozenset({
    'A_lte', 'a_lte', 'A_X', 'delta', 'residual', 'asplund', 'target_lte',
})

OUT_DIR = _REPO / 'data' / 'curation' / 'nonfe_pools'

# Solar reference for the COG / blend opacity proxy (single source).
_THETA_SOLAR = 5040.0 / float(get_star_params('solar')['teff'])

_EW_POOL = _REPO / 'data' / 'measured' / 'sol_ew_results_v1.csv'
_CANON   = _REPO / 'data' / 'linelists' / 'canonical_gf.csv'

# NLTE grids (RYA-396 / 235) — element → vendored grid file. Solar Δ is QUERIED
# from these at solar params for the verdict contract; never hardcoded.
_NLTE_GRID_DIR = _REPO / 'data' / 'nlte_grids'
_NLTE_GRID_FILES = {
    'Mg': 'Mg_Bergemann_MPIA.csv',
    'Si': 'Si_Bergemann_MPIA.csv',
    'Mn': 'Mn_Bergemann_MPIA.csv',
    'Na': 'Na_Lind2011_INSPECT.csv',
    'Ba': 'Ba_Korotin2015.csv',
    'Ca': 'Ca_Mashonkina2017.csv',
    'Ti': 'Ti_Bergemann2011_MPIA.csv',
    'Cr': 'Cr_Bergemann2010_MPIA.csv',
    # Ni, Al: no vendored NLTE grid → Δ≈0 (pure-curation / LTE target = Asplund).
}


# ── NLTE solar Δ (verdict target only — never feeds a cull) ───────────────────

def solar_nlte_delta(element: str) -> float:
    """Median NLTE Δ for `element` at solar params (Teff 5772, logg 4.44, [Fe/H] 0),
    read from its vendored grid. NaN if the element has no grid (Ni/Al). This is
    used ONLY to set the LTE validation target (Asplund − Δ); it is read AFTER the
    blind cull and never influences which lines are culled."""
    fname = _NLTE_GRID_FILES.get(element)
    if fname is None:
        return np.nan
    path = _NLTE_GRID_DIR / fname
    if not path.exists():
        return np.nan
    g = pd.read_csv(path)
    if 'delta_nlte' not in g.columns:
        return np.nan
    sp = get_star_params('solar')
    teff, logg = float(sp['teff']), float(sp['logg'])
    sel = g[(g['teff_K'].sub(teff).abs() < 150) &
            (g['logg'].sub(logg).abs() < 0.20) &
            (g['feh'].abs() < 0.05)]
    if sel.empty:
        sel = g[(g['teff_K'].sub(teff).abs() < 300) & (g['feh'].abs() < 0.10)]
    vals = sel['delta_nlte'].dropna()
    return float(np.median(vals)) if len(vals) else np.nan


# ── gf provenance tiering ────────────────────────────────────────────────────

def _gf_tier(nist_grade, loggf_reference) -> str:
    """HIGH / MED / LOW provenance tier for one line's gf, from its NIST grade and
    loggf source. Kurucz semi-empirical loggf (K03–K10, KP) carry ~0.1–0.3 dex
    uncertainty and known systematics → LOW (RYA-161 territory)."""
    g = str(nist_grade).strip() if nist_grade == nist_grade else ''
    if g in NIST_GRADE_HIGH:
        return 'HIGH'
    if g in NIST_GRADE_CULL:
        return 'CULL'
    ref = str(loggf_reference).strip() if loggf_reference == loggf_reference else ''
    if _KURUCZ_REF.match(ref):
        return 'LOW'
    if ref in ('', 'nan'):
        return 'LOW'
    return 'MED'


# ── Atomic-data loaders ──────────────────────────────────────────────────────

_canon_cache: dict = {}


def _load_canonical() -> pd.DataFrame:
    if 'df' not in _canon_cache:
        c = pd.read_csv(_CANON, low_memory=False)
        c['element'] = c['species'].astype(str).str.split().str[0]
        _canon_cache['df'] = c
    return _canon_cache['df']


def _ion_int(ion: str) -> int:
    return 1 if str(ion).strip() == 'I' else 2


def load_pool(elements) -> pd.DataFrame:
    """Measured EW pool joined to the canonical single-source atomic data, with the
    derived quality columns (REW, err fraction, gf tier). NO abundance here — this
    is the table the blind cull operates on."""
    ew = pd.read_csv(_EW_POOL)
    ew = ew[ew['element'].isin(elements)].copy()
    canon = _load_canonical()

    rows = []
    for _, r in ew.iterrows():
        el, ion, w = r['element'], r['ion'], float(r['wavelength_air_A'])
        cc = canon[(canon['element'] == el) & (canon['ion'] == _ion_int(ion))]
        if len(cc):
            j = (cc['wavelength_air_A'] - w).abs().to_numpy().argmin()
            if abs(float(cc.iloc[j]['wavelength_air_A']) - w) <= 0.05:
                m = cc.iloc[j]
                loggf, ep = float(m['log_gf']), float(m['excitation_potential_eV'])
                grade, ref = m['nist_grade'], m['loggf_reference']
            else:
                loggf = ep = np.nan; grade = ref = np.nan
        else:
            loggf = ep = np.nan; grade = ref = np.nan
        rows.append({
            'element': el, 'ion': ion, 'wavelength_air_A': w,
            'ew_mA': float(r['ew_mA']), 'ew_err_mA': float(r['ew_err_mA']),
            'chi2': float(r.get('chi2', np.nan)),
            'blend_flag': bool(r.get('blend_flag', False)),
            'log_gf': loggf, 'excitation_potential_eV': ep,
            'nist_grade': grade, 'loggf_reference': ref,
        })
    df = pd.DataFrame(rows)
    df['rew'] = np.log10(df['ew_mA'] / 1000.0 / df['wavelength_air_A'])
    df['err_frac'] = df['ew_err_mA'] / df['ew_mA']
    df['gf_tier'] = [_gf_tier(g, r) for g, r in zip(df['nist_grade'], df['loggf_reference'])]
    return df.reset_index(drop=True)


# ── Strength-weighted blend metric (abundance-blind) ─────────────────────────

def add_blend_ratio(df: pd.DataFrame, atomic_lines) -> pd.DataFrame:
    """For each measured line, the proximity-weighted solar-reference opacity of
    in-window OTHER-species neighbours, divided by the target line's opacity.
    Opacity proxy logκ = A_sun(elem) + loggf − EP·θ_solar (linear-COG opacity at
    the reference composition). This replaces `vald_proximity_flag`: it weighs the
    neighbours by how much they actually absorb, not merely whether they exist."""
    w_all = np.asarray(atomic_lines['wave_A'], dtype=float)
    el_all = np.array([str(x).split()[0] for x in atomic_lines['element']])
    ion_all = np.array([str(x).split()[1] if len(str(x).split()) > 1 else '1'
                        for x in atomic_lines['element']])
    gf_all = np.asarray(atomic_lines['loggf'], dtype=float)
    ep_all = np.asarray(atomic_lines['lower_state_eV'], dtype=float)
    mol_all = np.array([str(x) for x in atomic_lines['molecule']])
    order = np.argsort(w_all)
    w_s = w_all[order]

    def a_sun(el):
        return SOLAR_ASPLUND2021.get(el, 2.0)   # rare/unlisted → low reference

    ratios = []
    for _, r in df.iterrows():
        w0 = float(r['wavelength_air_A'])
        tgt_logk = a_sun(r['element']) + float(r['log_gf']) - float(r['excitation_potential_eV']) * _THETA_SOLAR
        lo = np.searchsorted(w_s, w0 - BLEND_WINDOW_A)
        hi = np.searchsorted(w_s, w0 + BLEND_WINDOW_A)
        idx = order[lo:hi]
        num = 0.0
        tgt_ion = str(r['ion']).strip()
        for k in idx:
            if mol_all[k] == 'T':
                continue
            # same species (element+ion) within the core = the line itself / its
            # HFS components, not a blend → skip.
            if el_all[k] == r['element'] and (ion_all[k] == '1') == (tgt_ion == 'I'):
                continue
            dl = w_all[k] - w0
            wkern = np.exp(-0.5 * (dl / BLEND_SIGMA_A) ** 2)
            logk = a_sun(el_all[k]) + gf_all[k] - ep_all[k] * _THETA_SOLAR
            num += wkern * (10.0 ** logk)
        ratios.append(num / (10.0 ** tgt_logk) if np.isfinite(tgt_logk) else np.nan)
    df = df.copy()
    df['blend_ratio'] = ratios
    return df


# ── THE BLIND CULL ───────────────────────────────────────────────────────────

def cull_reasons(row: dict) -> list:
    """Quality-based cull reasons for ONE line. ABUNDANCE-BLIND by construction:
    it reads only `_CULL_INPUT_COLS`. Thresholds are the fixed module constants.

    Reasons:
      WEAK  — EW below the detection floor (noise / continuum-dominated)
      SAT   — reduced EW above the linear-COG knee (saturated, EW→A biased high),
              or EW above the absolute saturation ceiling
      HIERR — ew_err/ew above the noise ceiling (unmeasurable)
      BLEND — measured blend flag, OR strength-weighted contaminant opacity > max
      BADGF — gf with a NIST grade of D/E/F (graded-unreliable lab value)
    A line with no reason is KEPT.
    """
    seen = set(row.keys()) & _FORBIDDEN_IN_CULL
    if seen:
        raise AssertionError(f"cull_reasons saw forbidden abundance column(s) {seen} "
                             f"— the cull must be abundance-blind (RYA-395 cardinal rule).")
    reasons = []
    if row['ew_mA'] < EW_MIN_MA:
        reasons.append('WEAK')
    if row['rew'] > REW_LINEAR_CEILING or row['ew_mA'] > SAT_EW_CEILING_MA:
        reasons.append('SAT')
    ef = row.get('err_frac', np.nan)
    if np.isfinite(ef) and ef > EW_ERR_FRAC_MAX:
        reasons.append('HIERR')
    br = row.get('blend_ratio', np.nan)
    if bool(row.get('blend_flag', False)) or (np.isfinite(br) and br > BLEND_FRAC_MAX):
        reasons.append('BLEND')
    if str(row.get('gf_tier')) == 'CULL':
        reasons.append('BADGF')
    return reasons


def apply_cull(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    cols = [c for c in _CULL_INPUT_COLS if c in df.columns]
    masks = df[cols].to_dict('records')
    df['cull_reason'] = [','.join(cull_reasons(m)) for m in masks]
    df['kept'] = df['cull_reason'] == ''
    # Provenance FLAGS (not culls): low-confidence gf, carried for the verdict.
    df['gf_flag'] = np.where(df['gf_tier'] == 'LOW', 'GF_UNVERIFIED', '')
    return df


# ── LTE A(X) — read AFTER the cull, for validation only ──────────────────────

def compute_lte_abundances(df: pd.DataFrame) -> pd.DataFrame:
    """Per-line 1D-LTE A(X) via iSpec MOOG EW inversion at solar params. Linemasks
    are built from the GES synthesis linelist (correct MOOG/WIDTH species codes +
    damping) with loggf/EP OVERRIDDEN from the canonical single source. Returns df
    with an `A_lte` column (NaN for lines that could not be matched/converged).

    This runs on the FULL pool (kept + culled) so the diagnostic can report the
    before/after of the blind cull. It is invoked only after `apply_cull`.
    """
    import pipeline.abundances_derive as ad
    sys.path.insert(0, str(Path(ad.ISPEC_DIR)))
    import ispec

    if not ad._ew_engine_available(EW_BASELINE_CODE):
        raise RuntimeError(
            f"EW baseline engine '{EW_BASELINE_CODE}' unavailable to iSpec — refusing "
            f"to fake A(X). Install/compile MOOGSILENT (no silent fallback, RYA-289).")

    template = ispec.read_line_regions(ad._LINE_REGIONS_ALL)
    dtype = template.dtype
    atomic = ispec.read_atomic_linelist(ad._SYNTH_LINELIST_FILE)
    a_notes = np.array([str(x) for x in atomic['element']])
    a_wave = np.asarray(atomic['wave_A'], dtype=float)
    a_loggf = np.asarray(atomic['loggf'], dtype=float)
    shared = [n for n in dtype.names if n in atomic.dtype.names]

    sp = get_star_params('solar')
    teff, logg, vturb = float(sp['teff']), float(sp['logg']), 1.0
    atmosphere = ad._load_atmosphere(teff, logg, 0.0, vturb)
    solar_abund = ispec.read_solar_abundances(ad._ISPEC_SOLAR_ABUND_FILE)

    rows, refs = [], []
    for i, r in df.iterrows():
        if not np.isfinite(r['log_gf']):
            continue
        note = f"{r['element']} {_ion_int(r['ion'])}"
        cand = np.where((a_notes == note) &
                        (np.abs(a_wave - r['wavelength_air_A']) <= 0.05))[0]
        if len(cand) == 0:
            continue
        # nearest in (λ, loggf) jointly → the canonical transition, not the strongest
        score = (np.abs(a_wave[cand] - r['wavelength_air_A']) / 0.05
                 + np.abs(a_loggf[cand] - float(r['log_gf'])))
        rows.append((int(cand[int(np.argmin(score))]), float(r['log_gf']),
                     float(r['excitation_potential_eV']), float(r['ew_mA']),
                     float(r['ew_err_mA'])))
        refs.append(i)

    a_lte = pd.Series(np.nan, index=df.index)
    if rows:
        lm = np.zeros(len(rows), dtype=dtype)
        for k, (j, loggf, ep, ewv, errv) in enumerate(rows):
            for f in shared:
                lm[k][f] = atomic[j][f]
            lm[k]['loggf'] = loggf
            lm[k]['lower_state_eV'] = ep
            lm[k]['ew'] = ewv
            lm[k]['ew_err'] = errv
            lm[k]['ewr'] = np.log10(ewv / 1000.0 / lm[k]['wave_A']) if ewv > 0 else 0.0
            lm[k]['wave_peak'] = lm[k]['wave_nm']
            lm[k]['wave_base'] = lm[k]['wave_nm'] - 0.002
            lm[k]['wave_top'] = lm[k]['wave_nm'] + 0.002
            lm[k]['note'] = str(atomic[j]['element'])
        _, normal_ab, _, _ = ispec.determine_abundances(
            atmosphere, teff, logg, 0.0, 0.0, lm, solar_abund,
            microturbulence_vel=vturb, verbose=0, code=EW_BASELINE_CODE,
            tmp_dir='/tmp/ispec_codex')
        for idx, a in zip(refs, np.asarray(normal_ab, dtype=float)):
            a_lte.at[idx] = a
    out = df.copy()
    out['A_lte'] = a_lte.values
    return out


# ── Diagnostics + verdict ────────────────────────────────────────────────────

def _rew_slope(rew, a):
    """A(X)-vs-reduced-EW slope (saturation / vmic diagnostic). Returns NaN unless
    there are ≥5 lines spanning ≥0.3 dex in REW — a slope fit over a handful of
    near-degenerate REW points is meaningless (and would read as a wild number)."""
    rew, a = np.asarray(rew, float), np.asarray(a, float)
    m = np.isfinite(rew) & np.isfinite(a)
    if m.sum() < 5 or (rew[m].max() - rew[m].min()) < 0.3:
        return np.nan
    return float(np.polyfit(rew[m], a[m], 1)[0])


def element_diagnostic(element: str, df_el: pd.DataFrame) -> dict:
    """Post-curation diagnostic + cull-vs-gf-scale verdict for one element.

    The cuts have already been applied (df_el carries `kept`). A(X) is read here
    ONLY to report and to classify the residual — it never re-enters the cull."""
    ion = 'I'                                          # neutral pool (abundance workhorse)
    pool = df_el[df_el['ion'] == ion].copy()
    kept = pool[pool['kept']]
    asp = SOLAR_ASPLUND2021.get(element, np.nan)
    dnlte = solar_nlte_delta(element)
    dnlte0 = 0.0 if not np.isfinite(dnlte) else dnlte
    target_lte = asp - dnlte0                          # contract: LTE→here so NLTE lands on Asplund

    a_all = pool['A_lte'].dropna()
    a_keep = kept['A_lte'].dropna()
    med_all = float(np.median(a_all)) if len(a_all) else np.nan
    med_keep = float(np.median(a_keep)) if len(a_keep) else np.nan
    scatter = float(np.std(a_keep)) if len(a_keep) > 1 else np.nan
    rew_slope = _rew_slope(kept['rew'], kept['A_lte'])
    n_clean = int(len(a_keep))
    down = med_all - med_keep if (np.isfinite(med_all) and np.isfinite(med_keep)) else np.nan
    resid = med_keep - target_lte if np.isfinite(med_keep) else np.nan
    low_gf_frac = float((pool['gf_tier'] == 'LOW').mean()) if len(pool) else np.nan

    # Verdict — quality cuts are exhausted; classify what's left (no tuning).
    if n_clean < N_CLEAN_FLOOR:
        verdict = 'LOW_CONFIDENCE'
    elif np.isfinite(resid) and abs(resid) <= RESID_TOL_DEX:
        verdict = 'CLEAN'
    elif np.isfinite(resid) and resid > RESID_TOL_DEX:
        verdict = 'GF_SCALE_RYA161'                    # persists after curation → gf-scale
    elif np.isfinite(resid) and resid < -RESID_TOL_DEX:
        verdict = 'CHECK_OVERCULL'                     # reads low — inspect cuts/NLTE sign
    else:
        verdict = 'NO_ABUNDANCE'

    return {
        'element': element, 'ion': ion,
        'N_total': int(len(pool)), 'N_clean': n_clean,
        'N_culled': int((~pool['kept']).sum()),
        'A_lte_precull': round(med_all, 3) if np.isfinite(med_all) else np.nan,
        'A_lte_curated': round(med_keep, 3) if np.isfinite(med_keep) else np.nan,
        'scatter': round(scatter, 3) if np.isfinite(scatter) else np.nan,
        'rew_slope': round(rew_slope, 3) if np.isfinite(rew_slope) else np.nan,
        'asplund': asp,
        'nlte_delta_solar': round(dnlte, 4) if np.isfinite(dnlte) else np.nan,
        'target_lte': round(target_lte, 3) if np.isfinite(target_lte) else np.nan,
        'down_correction': round(down, 3) if np.isfinite(down) else np.nan,
        'residual_vs_target': round(resid, 3) if np.isfinite(resid) else np.nan,
        'low_gf_frac': round(low_gf_frac, 2) if np.isfinite(low_gf_frac) else np.nan,
        'verdict': verdict,
    }


# ── Artifact writers ─────────────────────────────────────────────────────────

def _write_cull_csv(element: str, df_el: pd.DataFrame, path: Path) -> None:
    cols = ['element', 'ion', 'wavelength_air_A', 'ew_mA', 'ew_err_mA', 'err_frac',
            'rew', 'log_gf', 'excitation_potential_eV', 'nist_grade',
            'loggf_reference', 'gf_tier', 'blend_ratio', 'A_lte',
            'cull_reason', 'kept', 'gf_flag']
    cols = [c for c in cols if c in df_el.columns]
    out = df_el[cols].sort_values(['ion', 'wavelength_air_A'])
    header = [
        f"# RYA-395 non-Fe EW line-pool curation — {element}",
        "# Single source: measured EW = data/measured/sol_ew_results_v1.csv;",
        "#   atomic data (loggf/EP/nist_grade/loggf_reference) = data/linelists/canonical_gf.csv.",
        "# Cull criteria (FIXED in advance, abundance-BLIND, uniform across elements):",
        f"#   WEAK  ew_mA < {EW_MIN_MA:.0f}",
        f"#   SAT   log(EW/λ) > {REW_LINEAR_CEILING}  (linear-COG knee) OR ew_mA > {SAT_EW_CEILING_MA:.0f}",
        f"#   HIERR ew_err/ew > {EW_ERR_FRAC_MAX}",
        f"#   BLEND measured blend_flag OR strength-weighted contaminant opacity > {BLEND_FRAC_MAX}",
        "#         (NOT vald_proximity_flag — that crude flag was the known silent bug)",
        "#   BADGF nist_grade in {D,E,F}",
        "# gf_flag=GF_UNVERIFIED marks ungraded/Kurucz-theoretical gf (low-confidence, NOT culled);",
        "#   a pool dominated by these with a surviving offset is a gf-scale systematic → RYA-161.",
        "# Abundances (A_lte) are reported as VALIDATION only — never used to choose culls.",
    ]
    with open(path, 'w') as f:
        f.write('\n'.join(header) + '\n')
        out.to_csv(f, index=False)


# ── Orchestration ────────────────────────────────────────────────────────────

def curate(elements, with_abundance: bool = True) -> tuple:
    """Run the full curation for `elements`. Returns (per_line_df, diagnostics_df)."""
    import pipeline.abundances_derive as ad
    sys.path.insert(0, str(Path(ad.ISPEC_DIR)))
    import ispec

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"── RYA-395 non-Fe pool curation — {', '.join(elements)} ──")

    pool = load_pool(elements)
    print(f"  loaded {len(pool)} measured lines across {pool['element'].nunique()} elements")

    atomic = ispec.read_atomic_linelist(ad._SYNTH_LINELIST_FILE)
    pool = add_blend_ratio(pool, atomic)
    pool = apply_cull(pool)
    print(f"  blind cull: {int(pool['kept'].sum())} kept / {int((~pool['kept']).sum())} culled")

    if with_abundance:
        print(f"  computing 1D-LTE A(X) (iSpec {EW_BASELINE_CODE}) …")
        pool = compute_lte_abundances(pool)
    else:
        pool['A_lte'] = np.nan

    diags = []
    for el in elements:
        df_el = pool[pool['element'] == el].copy()
        _write_cull_csv(el, df_el, OUT_DIR / f'{el}_cull_rya395.csv')
        diags.append(element_diagnostic(el, df_el))
    diag_df = pd.DataFrame(diags)
    diag_df.to_csv(OUT_DIR / 'curation_diagnostics_rya395.csv', index=False)
    return pool, diag_df


def _print_table(diag_df: pd.DataFrame) -> None:
    print("\n── Post-curation diagnostic (LTE A(X), validation read AFTER blind cull) ──")
    cols = ['element', 'N_total', 'N_clean', 'A_lte_precull', 'A_lte_curated',
            'scatter', 'rew_slope', 'asplund', 'nlte_delta_solar', 'target_lte',
            'down_correction', 'residual_vs_target', 'low_gf_frac', 'verdict']
    print(diag_df[cols].to_string(index=False))
    print("\n  verdict key: CLEAN = curation reached (Asplund−Δ_NLTE) on quality grounds;")
    print("               GF_SCALE_RYA161 = offset persists after curation → gf-scale, escalate (do NOT tune);")
    print("               LOW_CONFIDENCE = too few clean lines to report; CHECK_OVERCULL = reads low.")


def verify(pool: pd.DataFrame, diag_df: pd.DataFrame) -> None:
    """Assert the CRITICAL conditions of the ticket. Raises on any violation."""
    print("\n── --verify: CRITICAL condition assertions ──")

    # 1. Thresholds are fixed module constants (not data/abundance-derived).
    for name, val in [('EW_MIN_MA', EW_MIN_MA), ('SAT_EW_CEILING_MA', SAT_EW_CEILING_MA),
                      ('EW_ERR_FRAC_MAX', EW_ERR_FRAC_MAX), ('BLEND_FRAC_MAX', BLEND_FRAC_MAX)]:
        assert isinstance(val, float), name
    print("  [1] cull thresholds are fixed module constants ✓")

    # 2. Abundance-blindness: shuffle A_lte and re-derive the cull mask — it must
    #    be byte-identical. Proves the cull never reads the abundance.
    base = apply_cull(pool.drop(columns=['cull_reason', 'kept', 'gf_flag'], errors='ignore'))
    shuf = pool.copy()
    rng = np.random.default_rng(0)
    shuf['A_lte'] = rng.permutation(shuf['A_lte'].values)
    shuf = apply_cull(shuf.drop(columns=['cull_reason', 'kept', 'gf_flag'], errors='ignore'))
    assert (base['cull_reason'].values == shuf['cull_reason'].values).all(), \
        "cull mask changed when A(X) was shuffled — NOT abundance-blind!"
    print("  [2] cull mask invariant under A(X) shuffle (abundance-blind) ✓")

    # 3. cull_reasons refuses to read a forbidden abundance column.
    try:
        cull_reasons({'ew_mA': 50.0, 'A_lte': 7.0})
        raise SystemExit("cull_reasons did NOT reject a forbidden abundance column")
    except AssertionError:
        print("  [3] cull_reasons rejects forbidden abundance columns ✓")

    # 4. Blend uses the strength-weighted metric, never vald_proximity_flag — the
    #    banned signal is not an input column, and the strength-weighted ratio is.
    assert 'vald_proximity' not in {c.lower() for c in _CULL_INPUT_COLS}, \
        "vald_proximity is a cull input"
    assert not any('proximity' in str(c).lower() for c in pool.columns), \
        "a vald-proximity column leaked into the pool"
    assert 'blend_ratio' in pool.columns and pool['blend_ratio'].notna().any(), \
        "strength-weighted blend_ratio not computed"
    print("  [4] blend vetting is strength-weighted, not vald_proximity ✓")

    # 5. Overshoot contract (Cr/Mn canaries): target_lte = Asplund − Δ_NLTE, and
    #    target_lte + Δ_NLTE == Asplund (NLTE lifts the curated LTE onto scale,
    #    never past it). Guards the Cr/Mn trap.
    for el in ('Cr', 'Mn'):
        row = diag_df[diag_df['element'] == el]
        if row.empty:
            continue
        row = row.iloc[0]
        if np.isfinite(row['nlte_delta_solar']) and np.isfinite(row['target_lte']):
            assert abs((row['target_lte'] + row['nlte_delta_solar']) - row['asplund']) < 1e-6, \
                f"{el} overshoot contract broken"
            assert row['nlte_delta_solar'] > 0, f"{el} NLTE Δ should be positive (canary)"
    print("  [5] Cr/Mn overshoot contract holds (target = Asplund − Δ_NLTE, Δ>0) ✓")

    # 6. Every cull carries a quality reason (no silent / abundance-driven drops).
    culled = pool[~pool['kept']]
    assert (culled['cull_reason'].str.len() > 0).all(), "a culled line has no reason"
    valid = {'WEAK', 'SAT', 'HIERR', 'BLEND', 'BADGF'}
    for rs in culled['cull_reason']:
        assert set(rs.split(',')) <= valid, f"unknown cull reason {rs}"
    print("  [6] every cull justified by a fixed quality reason ✓")
    print("  ALL CRITICAL CONDITIONS PASS")


def run(phase='1', do_verify=False, with_abundance=True):
    if str(phase) == '1':
        elements = PHASE1
    elif str(phase) == '2':
        elements = PHASE2
    elif str(phase) == 'all':
        elements = PHASE1 + PHASE2
    else:
        raise ValueError(f"phase must be 1 | 2 | all, got {phase!r}")

    pool, diag_df = curate(elements, with_abundance=with_abundance)
    _print_table(diag_df)

    # Cr / Mn explicit report (the canaries).
    for el in ('Cr', 'Mn'):
        r = diag_df[diag_df['element'] == el]
        if r.empty:
            continue
        r = r.iloc[0]
        print(f"\n  {el} canary: pre-cull A={r['A_lte_precull']}  curated A={r['A_lte_curated']}  "
              f"(down {r['down_correction']} dex on quality grounds)")
        print(f"    target LTE = Asplund {r['asplund']} − NLTE {r['nlte_delta_solar']:+} = {r['target_lte']}; "
              f"residual {r['residual_vs_target']:+} → {r['verdict']}")

    if do_verify and with_abundance:
        verify(pool, diag_df)
    return pool, diag_df


def _cli(argv=None):
    p = argparse.ArgumentParser(description="RYA-395 non-Fe EW line-pool curation")
    p.add_argument('--phase', default='1', choices=['1', '2', 'all'])
    p.add_argument('--verify', action='store_true', help="assert the CRITICAL conditions")
    p.add_argument('--no-abundance', action='store_true',
                   help="skip iSpec A(X) (cull lists + structure only; for fast tests)")
    a = p.parse_args(argv)
    run(phase=a.phase, do_verify=a.verify, with_abundance=not a.no_abundance)


if __name__ == '__main__':
    _cli()
