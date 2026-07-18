"""
pipeline/engine_selection.py — the two-engine floor selector (RYA-525)
======================================================================
The ratified law (RYA-525 description §2–3 + the three ratification comments +
the 2026-07-10 appendix): every element is measured on BOTH engines —
**Engine-A** (1D-NLTE = EW + grid delta) and **Engine-B** (synthesis) — and the
reported value is the quality-selected best, chosen PER LINE on line physics ALONE
and aggregated to the element. This module is the pure, reference-blind selector.

The pre-declared criterion (encode-don't-tune; NEVER keys on proximity to a
reference value — tuning firewall RYA-161). All thresholds come from
`config.constants.TWO_ENGINE`; none are inline here.

Per line:
  clause 1  validity gates — an engine is ELIGIBLE only if it produced a value and
            passed its own gate (Engine-A: grid in-hull; Engine-B: med_red_chi2 ≤
            synth_chi2_gate, an eligibility floor, NOT a quality selector).
  clause 2  exactly one eligible → report it.
  clause 3  both eligible → LINE REGIME decides:
              CLEAN-WEAK (unsaturated AND unblended AND not a problem-child) → Engine-A
              HARD       (blended OR saturated OR problem-child/HFS)          → Engine-B
  clause 4  INDETERMINATE regime only (neither clearly clean-weak nor clearly hard)
            → lower line-scatter σ; exact tie → 1D-NLTE (the anchor / differential
            zero-point, so ties stay on one scale). Clause 4 governs ONLY the border;
            clause 3 governs every clear regime — they do not overlap.
  clause 5  neither eligible → no value (caller records a disposition; the loud-fail
            guard forbids a silent single-engine PASS).
  clause 6  the rejected engine per line is recorded + shown but EXCLUDED from the
            value AND the uncertainty budget.

Per element:
  aggregation — the reported value = inverse-variance combine of the per-LINE WINNERS
            (each line contributes only its winning engine's value + error). The
            budget carries winning-line uncertainties only (clause 6, per line).
  cross-engine-mix guard (the Ti lesson) — if the winners span BOTH engines AND the
            mean cross-engine Δ exceeds `cross_engine_mix_gate`, FLAG + adjudicate;
            never silently average two disagreeing scales.
  cross-engine Δ — always recorded as a SEPARATE diagnostic, never folded into σ.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from config.constants import TWO_ENGINE, NLTE_CORRECTION_ELEMENTS

# ── engine + regime labels ────────────────────────────────────────────────────
ENGINE_A = 'engineA_1dnlte'   # 1D-NLTE = EW + grid delta
ENGINE_B = 'engineB_synth'    # synthesis (Turbospectrum LTE + TS-native NLTE)

CLEAN_WEAK = 'clean-weak'
HARD = 'hard'
INDETERMINATE = 'indeterminate'


class TwoEngineError(RuntimeError):
    """A two-engine-floor law was violated (missing grid, silent single engine, …)."""


class ReferenceProximityError(TwoEngineError):
    """A selection input referenced proximity to a reference value — the tuning firewall."""


# ── inputs / outputs ──────────────────────────────────────────────────────────
@dataclass(frozen=True)
class LineEngines:
    """Per-line inputs from BOTH engines + reference-blind line physics.

    Reference-blind BY CONSTRUCTION: there is no reference-value field here, so the
    selector cannot key on |A − reference|. `a_value`/`b_value` is None when that
    engine produced nothing for this line."""
    wavelength: float
    species: str
    # Engine-A (1D-NLTE)
    a_value: Optional[float]
    a_err: Optional[float]
    a_in_hull: bool = True         # the NLTE grid covers the star's params
    # Engine-B (synthesis)
    b_value: Optional[float] = None
    b_err: Optional[float] = None
    b_chi2: Optional[float] = None  # med_red_chi2 of the synth fit
    # reference-blind line physics
    ew_mA: Optional[float] = None
    blend_flag: bool = False
    is_problem_child: bool = False  # on the RYA-463 registry (blend/HFS/special)


@dataclass(frozen=True)
class LineWinner:
    wavelength: float
    species: str
    engine: str
    value: float
    err: Optional[float]
    regime: str
    reason: str
    rejected_engine: Optional[str]
    rejected_value: Optional[float]
    cross_engine_delta: Optional[float]   # b − a where both present (diagnostic only)


@dataclass(frozen=True)
class ElementRecord:
    species: str
    value: float                 # the reported (recommended) element value
    err: float
    n_lines: int
    selected_engines: tuple      # engines that won ≥1 line
    cross_engine_mix: bool       # winners span both engines
    mix_flagged: bool            # mix AND |mean Δ| > gate → adjudicate, do not trust the mean
    mean_cross_engine_delta: Optional[float]
    engineA_value: Optional[float]   # element-level A-only aggregate (DIAGNOSTIC, not reported)
    engineB_value: Optional[float]   # element-level B-only aggregate (DIAGNOSTIC, not reported)
    per_line: tuple = field(default_factory=tuple)


# ── clause 1: validity gates ──────────────────────────────────────────────────
def _eligible_A(l: LineEngines) -> bool:
    return l.a_value is not None and l.a_in_hull


def _eligible_B(l: LineEngines, cfg=TWO_ENGINE) -> bool:
    if l.b_value is None:
        return False
    # med_red_chi2 gate is an ELIGIBILITY floor (synth didn't catastrophically fail),
    # NOT a quality selector — a missing chi2 is treated as "ran, no failure signal".
    return l.b_chi2 is None or l.b_chi2 <= cfg['synth_chi2_gate']


# ── clause 3/4: reference-blind line-regime classifier ────────────────────────
def classify_regime(l: LineEngines, cfg=TWO_ENGINE) -> str:
    """CLEAN-WEAK / HARD / INDETERMINATE from line physics ONLY (never a reference)."""
    elem = l.species.split()[0] if l.species else ''
    saturated = l.ew_mA is not None and l.ew_mA > cfg['saturation_knee_mA']
    hfs = elem in cfg['hfs_elements']
    if l.blend_flag or saturated or l.is_problem_child or hfs:
        return HARD
    # clean-weak requires a POSITIVE clean signal (a known-unsaturated EW), not merely
    # the absence of hard flags — an unknown EW is INDETERMINATE, not clean.
    if l.ew_mA is not None and l.ew_mA <= cfg['saturation_knee_mA'] and not l.blend_flag \
            and not l.is_problem_child and not hfs:
        return CLEAN_WEAK
    return INDETERMINATE


def _mk(l, engine, regime, reason, delta) -> LineWinner:
    if engine == ENGINE_A:
        val, err, rej, rejval = l.a_value, l.a_err, ENGINE_B, l.b_value
    else:
        val, err, rej, rejval = l.b_value, l.b_err, ENGINE_A, l.a_value
    return LineWinner(l.wavelength, l.species, engine, float(val),
                      (float(err) if err is not None else None), regime, reason,
                      (rej if rejval is not None else None),
                      (float(rejval) if rejval is not None else None), delta)


def select_line(l: LineEngines, cfg=TWO_ENGINE) -> Optional[LineWinner]:
    """The per-line winner (or None if neither engine is eligible → clause 5)."""
    eA, eB = _eligible_A(l), _eligible_B(l, cfg)
    delta = (float(l.b_value) - float(l.a_value)) if (l.a_value is not None and l.b_value is not None) else None
    if not eA and not eB:
        return None
    regime = classify_regime(l, cfg)
    if eA and not eB:                                   # clause 2
        return _mk(l, ENGINE_A, regime, 'only Engine-A eligible', delta)
    if eB and not eA:                                   # clause 2
        return _mk(l, ENGINE_B, regime, 'only Engine-B eligible', delta)
    # both eligible — clause 3 (clear regimes)
    if regime == CLEAN_WEAK:
        return _mk(l, ENGINE_A, regime, 'clean-weak line → 1D-NLTE (cleanest for weak lines)', delta)
    if regime == HARD:
        return _mk(l, ENGINE_B, regime, 'hard line (blend/saturation/HFS) → synthesis', delta)
    # clause 4 — INDETERMINATE regime only: lower line-scatter σ; exact tie → anchor default
    a_s = l.a_err if (l.a_err is not None) else np.inf
    b_s = l.b_err if (l.b_err is not None) else np.inf
    if b_s < a_s:
        return _mk(l, ENGINE_B, regime, 'indeterminate regime → lower line-scatter σ (Engine-B)', delta)
    if a_s < b_s:
        return _mk(l, ENGINE_A, regime, 'indeterminate regime → lower line-scatter σ (Engine-A)', delta)
    tie = ENGINE_A if cfg['tie_default_engine'] == 'engineA_1dnlte' else ENGINE_B
    return _mk(l, tie, regime, 'indeterminate tie → 1D-NLTE anchor-scale default', delta)


# ── element aggregation ───────────────────────────────────────────────────────
def _inverse_variance(values, errs):
    values = np.asarray(values, float)
    errs = np.asarray([e if (e is not None and np.isfinite(e) and e > 0) else np.nan for e in errs], float)
    if np.all(np.isfinite(errs)):
        w = 1.0 / errs ** 2
        val = float(np.sum(w * values) / np.sum(w))
        err = float(np.sqrt(1.0 / np.sum(w)))
        return val, err
    # fall back to a scatter-weighted (median + spread) combine when any err is missing
    val = float(np.median(values))
    err = float(np.std(values, ddof=1)) if len(values) > 1 else float('nan')
    return val, err


def aggregate_element(species: str, winners, cfg=TWO_ENGINE) -> ElementRecord:
    """Combine the per-LINE winners into the reported element value (clause: aggregation)."""
    winners = [w for w in winners if w is not None]
    if not winners:
        raise TwoEngineError(
            f"{species}: no line eligible on either engine — record a disposition, "
            "never a silent value (RYA-525 §2.5)")
    value, err = _inverse_variance([w.value for w in winners], [w.err for w in winners])
    engines = tuple(sorted({w.engine for w in winners}))
    mix = len(engines) > 1
    deltas = [w.cross_engine_delta for w in winners if w.cross_engine_delta is not None]
    mean_delta = float(np.mean(deltas)) if deltas else None
    mix_flagged = bool(mix and mean_delta is not None and abs(mean_delta) > cfg['cross_engine_mix_gate'])
    # element-level per-engine diagnostic aggregates (BOTH values per line; not reported)
    a_vals = [(w.value if w.engine == ENGINE_A else w.rejected_value) for w in winners]
    b_vals = [(w.value if w.engine == ENGINE_B else w.rejected_value) for w in winners]
    a_vals = [v for v in a_vals if v is not None]
    b_vals = [v for v in b_vals if v is not None]
    engA = float(np.median(a_vals)) if a_vals else None
    engB = float(np.median(b_vals)) if b_vals else None
    return ElementRecord(
        species=species, value=value, err=err, n_lines=len(winners),
        selected_engines=engines, cross_engine_mix=mix, mix_flagged=mix_flagged,
        mean_cross_engine_delta=mean_delta, engineA_value=engA, engineB_value=engB,
        per_line=tuple(winners))


def select_element(species: str, lines, cfg=TWO_ENGINE) -> ElementRecord:
    """End-to-end: classify+select every line, then aggregate to the element value."""
    return aggregate_element(species, [select_line(l, cfg) for l in lines], cfg)


# ── loud-fail guards (RYA-525 §3; siblings of RYA-409/518) ────────────────────
# The RYA-526 two-engine coverage ledger is the pre-declared EXCEPTION LIST. A
# disposition of 'acquire-task'/'build-task' means a genuinely-absent grid; 'wired-both'
# means both engines MUST run; 'wired-one'/'LTE-only-by-design' are documented, owned
# single-engine states (recorded, never a silent PASS).
_LEDGER_TASK = {'acquire-task', 'build-task'}
_LEDGER_ONE_OK = {'wired-one', 'LTE-only-by-design'}


def require_grid_or_raise(species: str, ledger_disposition: str, produced_engines: set) -> None:
    """Guard (a): a missing synthesis grid RAISES — never a silent EW-1D PASS."""
    if ledger_disposition in _LEDGER_TASK:
        raise TwoEngineError(
            f"{species}: synthesis grid genuinely absent (ledger disposition="
            f"'{ledger_disposition}') — acquire/build it (RYA-526/540), do NOT report a "
            "silent single-engine value (RYA-525 §3)")
    if ledger_disposition == 'wired-both' and produced_engines != {ENGINE_A, ENGINE_B}:
        raise TwoEngineError(
            f"{species}: ledger says wired-both but the run produced only {produced_engines or '∅'} "
            "— a wired grid failed to run; RAISE, do not downgrade (RYA-525 §3)")


def assert_cross_engine_recorded(species: str, record: ElementRecord,
                                 ledger_disposition: str) -> None:
    """Guard (b): a single-engine reported value with NO cross-engine record RAISES,
    unless the ledger documents the element as a one-engine / LTE-by-design state."""
    if len(record.selected_engines) == 1 and record.mean_cross_engine_delta is None \
            and ledger_disposition not in _LEDGER_ONE_OK:
        raise TwoEngineError(
            f"{species}: reported on one engine with no cross-engine record and no "
            f"documented disposition (ledger='{ledger_disposition}') — RAISE (RYA-525 §3)")


def assert_reference_blind(selector_inputs) -> None:
    """Guard (c) / tuning firewall (RYA-161): the selector must never see a reference
    value. `LineEngines` has no reference field by construction; this asserts a caller
    did not smuggle one in (e.g. an `a_ref`/`reference`/`asplund`/`delta_vs_*` key)."""
    banned = ('ref', 'reference', 'asplund', 'delta_vs', 'proximity', 'closer', 'target_value')
    for name in selector_inputs:
        low = str(name).lower()
        if any(b in low for b in banned):
            raise ReferenceProximityError(
                f"selection input '{name}' references a reference value — selection must key "
                "on line physics only (RYA-525 §2.2, tuning firewall RYA-161)")


# ── Guard (d): ratified-excluded species (RYA-558) ───────────────────────────
# The reference-blind floor selects on line physics alone — but it must not be able to
# report a species the science has RATIFIED as excluded. The archetype is Cr II: the
# RYA-232 −0.777 dex was two saturated lines in the COG damping wing (a line-matching
# artifact, not a real systematic), so Cr II was deliberately excluded (RYA-240) and Cr
# is reported as Cr I. Left unguarded, the floor picked Cr II 5.676 in the RYA-527 re-emit
# only because it looked better than a worse raw-EW artifact (8.354) — the exact failure
# this guard prevents. A ratified-excluded species stays a cross-engine DIAGNOSTIC; it is
# NEVER the reported value. Promotion (e.g. Cr II) is a separate decision gated on clean
# unsaturated weak lines — not something the blind floor may do implicitly.
#
# Single source of truth = the NLTE registry: a registered element's ratified ionisation
# stage is fixed (NLTE_CORRECTION_ELEMENTS[el]['ion']; Cr=1, Sr/Ba=2), so a species on a
# different ion is excluded-from-value. Explicit, cited exclusions are listed below.
_ION_NUM = {'I': 1, 'II': 2, 'III': 3}
RATIFIED_EXCLUDED_SPECIES = {
    'Cr II': 'RYA-240 — COG/damping-wing saturation artifact (the RYA-232 −0.777 dex was '
             '2 saturated lines); Cr is reported as Cr I. Bergemann & Cescutti 2010.',
}


def ratified_reported_ion(element: str):
    """The ratified reported ionisation stage for an element, from the NLTE registry
    (Cr=1 Cr I, Sr/Ba=2). None if the element is not registry-ion-locked."""
    e = NLTE_CORRECTION_ELEMENTS.get(element)
    return e.get('ion') if e else None


def is_ratified_excluded_species(species: str) -> bool:
    """True if `species` (e.g. 'Cr II') is a ratified EXCLUSION — a cross-engine
    diagnostic only, NEVER the reported value of the reference-blind floor. Sourced from
    the registry (a registered element's ratified ion is fixed; a different ion is
    excluded) plus the explicit RATIFIED_EXCLUDED_SPECIES list."""
    parts = str(species).split()
    if len(parts) < 2:
        return False
    if species in RATIFIED_EXCLUDED_SPECIES:
        return True
    ratified = ratified_reported_ion(parts[0])
    ion_num = _ION_NUM.get(parts[1])
    return ratified is not None and ion_num is not None and ion_num != ratified


def exclusion_reason(species: str) -> str:
    """Cited reason a species is ratified-excluded (for the diagnostic label / log)."""
    if species in RATIFIED_EXCLUDED_SPECIES:
        return RATIFIED_EXCLUDED_SPECIES[species]
    return (f"not the registry-ratified ion for {str(species).split()[0]} "
            f"(NLTE_CORRECTION_ELEMENTS ion={ratified_reported_ion(str(species).split()[0])})")


def assert_not_excluded_value(species: str) -> None:
    """Guard (d): a ratified-excluded species must NEVER be the reported floor value.
    Loud-fail (record it as a diagnostic instead) — the blind floor may not implicitly
    promote a ratified exclusion (RYA-558/240)."""
    if is_ratified_excluded_species(species):
        raise TwoEngineError(
            f"{species}: ratified-excluded species cannot be the reported value of the "
            f"reference-blind floor ({exclusion_reason(species)}) — keep it as a "
            f"cross-engine DIAGNOSTIC, report the ratified ion (RYA-558/240).")
