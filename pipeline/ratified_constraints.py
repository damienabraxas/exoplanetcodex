#!/usr/bin/env python3
"""
pipeline/ratified_constraints.py — RYA-674 §2C
==============================================
A RATIFIED CONSTRAINT THAT A DOWNSTREAM MODULE CAN SILENTLY VIOLATE IS NOT RATIFIED —
IT IS A SUGGESTION.

The defect class
----------------
RYA-669's Phase 2 execution surfaced three violations of ratified decisions, all in one
re-emit, all the same shape: *a constraint Ryan ratified, checked by the wrong signal
or by no signal at all, in a module that never knew the constraint existed.*

  * **Fe method_scale desync (RYA-553/665/681).** The 1D→3D correction's idempotency was
    keyed on a prose label a *different* process wrote. The label could not follow the
    value, so the correction re-armed and A(Fe I) landed at 7.416 — inside `FE_GATE` by
    0.006 dex, because the correction magnitude equals the gate half-width.
  * **Li 1.409 veto leak (RYA-563/103/458).** Li carries the registry `upper_limit`
    disposition: the reference-blind two-engine floor may never emit a synthesis point
    value for it. The July re-emit adopted 1.409 anyway.
  * **Cr II 5.676 exclusion leak (RYA-240/558).** Cr II is a ratified exclusion — a
    cross-engine diagnostic, never the reported value. The floor picked it because it
    looked better than a worse artifact, and the verdict proposed it.

Each was fixed once, in the module where it bit. None of the fixes stopped the *next*
module from making the same mistake, because each fix lived where the bug was rather
than where the emission is.

What this module does instead
-----------------------------
It makes the check a property of the ACT OF EMITTING, not of any one emitter. Every
path that writes a per-element result — the phase_c verdict generator, the gold
reference builder, the disposition report, the two-engine record emitter — calls
`assert_ratified_constraints_satisfied()` on its candidate rows before writing.
Violation raises, naming the module, the constraint, the ratifying ticket and the
offending row. `pipeline/provenance_honesty.py` (RYA-596/653) is the template: make the
contradiction unrepresentable, not merely detectable.

Two row kinds, and why the distinction is load-bearing
------------------------------------------------------
`RowKind.ELEMENT_VALUE` is an element-level assertion — "element X's reported / frozen /
proposed value is V". `RowKind.SPECIES_RECORD` is a per-species measurement record in a
diagnostic table, where by design a table may legitimately carry species we would never
report (the two-engine artifact records Fe II, Ti II, Si II alongside the reported ions).

Collapsing the two would be wrong in both directions: applied to species records, the
Cr II exclusion would forbid the very diagnostic RYA-558 ratified keeping; applied only
to element rows, the floor could keep emitting Li 1.409 as a species record for
something downstream to adopt. So a species record may carry an excluded/vetoed species
only if it is marked `diagnostic_only` — the demotion RYA-558 and RYA-563 both specify —
and an element-level row may never source its value from one.

Single-sourcing
---------------
No constraint here re-implements a ratified rule. `is_upper_limit_disposition` and
`RATIFIED_EXCLUDED_SPECIES` come from `pipeline.engine_selection` (registry-sourced);
the scale classification comes from `pipeline.solar_scale_provenance` (RYA-681); the
correction identifiers come from `config/corrections_registry.yaml`. This module is
where they are *invoked at emission time*, not a second copy of any of them.

Adding a constraint
-------------------
See `docs/SCIENCE_STANDARDS.md` § "Ratified Constraints". It needs a Ryan decision, an
entry here citing the ratifying RYA-#, and a test in
`tests/test_ratified_constraints.py`. The registry is append-only + revoked (with a
revocation ticket), never silent-delete.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Mapping, Optional

__all__ = [
    'RatifiedConstraintViolation',
    'RowKind',
    'ConstraintType',
    'RatifiedConstraint',
    'EmissionRow',
    'RATIFIED_CONSTRAINTS',
    'constraint_by_id',
    'check_ratified_constraints',
    'assert_ratified_constraints_satisfied',
]


class RatifiedConstraintViolation(RuntimeError):
    """An emission path was about to write a row that violates a ratified constraint.

    Raised, never warned (RYA-518 discipline). The message names the emitting module,
    the constraint, the ratifying ticket and the row, because "something violated a
    constraint somewhere" is not actionable at 2am.
    """


class RowKind:
    #: An element-level assertion: "element X's reported/frozen/proposed value is V".
    ELEMENT_VALUE = 'ELEMENT_VALUE'
    #: A per-species record in a diagnostic table (the two-engine floor's own output).
    SPECIES_RECORD = 'SPECIES_RECORD'


class ConstraintType:
    FORBIDDEN_VALUE = 'FORBIDDEN_VALUE'
    EXCLUDED_SPECIES = 'EXCLUDED_SPECIES'
    REQUIRED_CORRECTION = 'REQUIRED_CORRECTION'
    EXCLUDED_LINE = 'EXCLUDED_LINE'


# ── the normalised view of one candidate row ─────────────────────────────────
# Emission paths have four different row schemas (verdict rows, gold CSV rows,
# disposition dataclasses, two-engine records) and unifying them is not this ticket's
# job. So the checks read a normalised VIEW, and each field names every schema's
# spelling of it in one place — rather than every check learning all four.

_VALUE_KEYS = ('A_measured', 'A_X_nlte', 'A_X', 'reported', 'value', 'v3')
_SOURCE_KEYS = ('source', 'value_source', 'channel', 'provenance', 'selected_engines',
                'method_scale', 'note', 'owed', 'promotion_reason')
#: Substrings that mark a value as having come from the two-engine synthesis floor.
#: These are the strings the emitters actually write (`selected_engines` entries are
#: `engineA_1dnlte` / `engineB_synth`; the RYA-527 ladder writes "two-engine ... floor").
_FLOOR_SOURCE_MARKERS = ('two-engine', 'two_engine', 'engineb_synth', 'enginea_',
                         'synthesis floor', 'synth floor', 'b_synth')


@dataclass(frozen=True)
class EmissionRow:
    """One candidate output row, normalised. `raw` is always the original mapping."""
    element: str
    kind: str
    value: Optional[float]
    ion: Optional[str]
    species: Optional[str]
    diagnostic_only: bool
    source_text: str
    verdict: str
    raw: Mapping[str, Any] = field(repr=False, default_factory=dict)

    def summary(self) -> str:
        v = 'None' if self.value is None else f'{self.value:.3f}'
        return (f"{self.species or self.element} value={v} kind={self.kind}"
                + (f" [diagnostic_only]" if self.diagnostic_only else ""))


def _as_float(x) -> Optional[float]:
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def _normalise(row: Any, kind: str) -> EmissionRow:
    if hasattr(row, 'as_dict'):                       # ElementDisposition & friends
        row = row.as_dict()
    if not isinstance(row, Mapping):
        raise RatifiedConstraintViolation(
            f"cannot check a candidate row of type {type(row).__name__} — emission paths "
            f"must hand this module mappings (or objects with .as_dict()), so the check "
            f"can never silently skip a row it did not understand (RYA-674)")
    element = str(row.get('element') or '').strip()
    ion = row.get('ion')
    ion = None if ion is None or str(ion).strip().lower() in ('', 'nan', 'none') else str(ion).strip()
    species = row.get('species')
    species = str(species).strip() if species else (f"{element} {ion}" if element and ion else None)
    value = None
    for k in _VALUE_KEYS:
        if k in row:
            value = _as_float(row.get(k))
            if value is not None:
                break
    bits = []
    for k in _SOURCE_KEYS:
        v = row.get(k)
        if isinstance(v, (list, tuple)):
            bits.extend(str(x) for x in v)
        elif v is not None:
            bits.append(str(v))
    return EmissionRow(
        element=element, kind=kind, value=value, ion=ion, species=species,
        diagnostic_only=bool(row.get('diagnostic_only')),
        source_text=' | '.join(bits).lower(),
        verdict=str(row.get('verdict') or '').strip().upper(),
        raw=row)


def _looks_floor_sourced(row: EmissionRow) -> bool:
    return any(m in row.source_text for m in _FLOOR_SOURCE_MARKERS)


# ── the three ratified constraints ───────────────────────────────────────────
# Each check_fn returns a reason string (VIOLATION) or None (PASS). None of them
# re-derives a ratified rule; each invokes the module that already owns it.

def _check_upper_limit_veto(row: EmissionRow) -> Optional[str]:
    """`Li_6707_veto_1_409` — RYA-563 (with RYA-103/458).

    Ratified: Li I 6707.84 carries the registry disposition
    `required_treatment=upper_limit` — "CN-blended; a clean low value is a RED FLAG.
    Carried as UPPER_LIMIT, never a point value." The reference-blind two-engine floor
    may NOT emit a synthesis point value for such an element; the floor's 1.409 is a
    DIAGNOSTIC, and the phase_c upper limit governs.

    Membership is read from the registry via
    `engine_selection.is_upper_limit_disposition` — never a hardcoded `{'Li'}`, so an
    element acquiring or losing the disposition needs no edit here (RYA-563's own rule).
    """
    from pipeline.engine_selection import is_upper_limit_disposition
    if not row.element or not is_upper_limit_disposition(row.element):
        return None
    if row.diagnostic_only:
        return None                       # the ratified demotion — a diagnostic is fine
    if row.kind == RowKind.SPECIES_RECORD:
        if row.value is not None:
            return (f"the two-engine floor reports a point value A({row.element}) = "
                    f"{row.value:.3f} for an element whose registry disposition is "
                    f"upper_limit. RYA-563: record it under `diagnostic_only`, never as "
                    f"the reported value.")
        return None
    # ELEMENT_VALUE: the upper limit itself is legitimate (phase_c reports 0.727 and the
    # gold freezes it at tier `upper_limit`). What is forbidden is ADOPTING the floor's
    # synthesis point value in its place — the concrete RYA-669 leak.
    if row.value is not None and _looks_floor_sourced(row):
        return (f"A({row.element}) = {row.value:.3f} is sourced from the two-engine "
                f"synthesis floor ({row.source_text[:120]}), but {row.element} carries the "
                f"registry upper_limit disposition. RYA-563: the floor value is a "
                f"DIAGNOSTIC; the phase_c upper limit governs the reported value.")
    if row.verdict == 'PASS':
        return (f"{row.element} is emitted with verdict PASS, but its registry disposition "
                f"is upper_limit — an upper limit is structurally never a PASS point value "
                f"(RYA-563/103/458).")
    return None


def _check_excluded_species(row: EmissionRow) -> Optional[str]:
    """`Cr_II_species_exclusion` — RYA-240 / RYA-558.

    Ratified: Cr II is excluded from the reported value. The RYA-232 −0.777 dex was two
    saturated lines in the COG damping wing (a line-matching artifact, not a real
    systematic), so Cr is reported as Cr I and Cr II stays a cross-engine DIAGNOSTIC.

    Scope is the EXPLICIT, cited exclusion list
    `engine_selection.RATIFIED_EXCLUDED_SPECIES` — single-sourced, and deliberately NOT
    `is_ratified_excluded_species()`, whose registry-ion rule additionally excludes
    Ti II / Si II. Those are real exclusions in the floor's own selector, but Ryan has
    not ratified them as emission-time constraints and RYA-674 adds the three known
    ones with no interpretation. A future ratification widens this scope by editing the
    single source, not this function.
    """
    from pipeline.engine_selection import RATIFIED_EXCLUDED_SPECIES, exclusion_reason
    species = row.species
    if not species or species not in RATIFIED_EXCLUDED_SPECIES:
        return None
    if row.diagnostic_only:
        return None                       # the ratified demotion — a diagnostic is fine
    if row.value is None:
        return None                       # no value asserted, nothing to exclude
    where = ('the two-engine floor reports' if row.kind == RowKind.SPECIES_RECORD
             else 'this emission reports')
    return (f"{where} a value A({species}) = {row.value:.3f} for a RATIFIED-EXCLUDED "
            f"species. {exclusion_reason(species)} A ratified-excluded species is a "
            f"cross-engine DIAGNOSTIC (mark the row `diagnostic_only`), never the "
            f"reported value (RYA-240/558).")


def _check_fe_reported_scale(row: EmissionRow) -> Optional[str]:
    """`Fe_1D_3D_correction_required_on_solar_report` — RYA-553 (hardened RYA-681/674).

    Ratified: the reported solar A(Fe I) is quoted on the true 3D scale — the tabulated
    Magic-2013 1D→3D correction is applied exactly once, at the reported-value layer,
    after NLTE. Three distinct ways an emission can break that, all checked here:

      1. the value is on NEITHER scale — the doubled-correction signature (7.416);
      2. the row's declaration and its value disagree — gold v3's shape, the state that
         re-arms the correction;
      3. the value is self-consistently on the 1D-NLTE scale — the correction is simply
         missing from a reported solar Fe number.

    Applies to ELEMENT_VALUE rows only. A per-species floor record is a pre-report
    engine measurement (Fe I 7.58 on the raw EW leg) and is not a claim about the
    reported anchor.
    """
    from pipeline.solar_scale_provenance import (
        REPORTED_SCALE_CORRECTED_ELEMENTS, SCALE_3D_NLTE, ScaleProvenanceError,
        declared_scale, scale_centres, scale_discrimination_halfwidth, scale_from_value)
    if row.kind != RowKind.ELEMENT_VALUE:
        return None
    if row.element not in REPORTED_SCALE_CORRECTED_ELEMENTS or row.value is None:
        return None
    el, a_x = row.element, row.value
    from_value = scale_from_value(el, a_x)
    centres = scale_centres(el)
    half = scale_discrimination_halfwidth(el)
    if from_value is None:
        return (f"A({el}) = {a_x:.3f} is on NEITHER recognised abundance scale "
                f"(1D-NLTE centre {centres['1D-NLTE']:.3f}, 3D-NLTE centre "
                f"{centres['3D-NLTE']:.3f}, discrimination half-width ±{half:.3f}). The "
                f"signature of the RYA-553 correction applied more than once — RYA-669 "
                f"measured exactly this at 7.416, and FE_GATE cannot catch it because the "
                f"correction magnitude equals the gate half-width.")
    try:
        declared, source = declared_scale(row.raw)
    except ScaleProvenanceError as exc:
        return f"the {el} row's scale declaration is unreadable: {exc}"
    if declared is not None and declared != from_value:
        return (f"the {el} row CONTRADICTS ITSELF: declared '{declared}' (from {source}) "
                f"but A({el}) = {a_x:.3f} sits on the '{from_value}' scale. This is the "
                f"RYA-669 defect — a post-correction value under a pre-correction "
                f"declaration re-arms the RYA-553 correction on the next trip round the "
                f"gold → verdict → candidate → freeze loop.")
    if from_value != SCALE_3D_NLTE:
        return (f"the reported A({el}) = {a_x:.3f} is on the '{from_value}' scale, but a "
                f"reported solar {el} anchor must carry the tabulated Magic-2013 1D→3D "
                f"correction (RYA-553) and sit on '{SCALE_3D_NLTE}' (centre "
                f"{centres[SCALE_3D_NLTE]:.3f} ±{half:.3f}).")
    return None


@dataclass(frozen=True)
class RatifiedConstraint:
    constraint_id: str
    scope: str
    constraint_type: str
    check_fn: Callable[[EmissionRow], Optional[str]]
    provenance_ticket: str
    #: Set ONLY by an explicit Ryan-ratified reversal, naming the revoking ticket. The
    #: registry is append-only + revoked, never silent-delete (RYA-674 §2D).
    revoked_by: Optional[str] = None
    applies_to: tuple[str, ...] = (RowKind.ELEMENT_VALUE, RowKind.SPECIES_RECORD)

    @property
    def active(self) -> bool:
        return self.revoked_by is None


RATIFIED_CONSTRAINTS: tuple[RatifiedConstraint, ...] = (
    RatifiedConstraint(
        constraint_id='Li_6707_veto_1_409',
        scope="any element carrying the registry `required_treatment=upper_limit` "
              "disposition (today: Li, from Li I 6707.84) — its reported value may not "
              "be the two-engine synthesis point value (the leak was A(Li) 1.409)",
        constraint_type=ConstraintType.FORBIDDEN_VALUE,
        check_fn=_check_upper_limit_veto,
        provenance_ticket='RYA-563 (with RYA-103 / RYA-458)',
    ),
    RatifiedConstraint(
        constraint_id='Cr_II_species_exclusion',
        scope="the explicitly-ratified excluded species in "
              "`engine_selection.RATIFIED_EXCLUDED_SPECIES` (today: Cr II) — never the "
              "reported value of any emission; cross-engine diagnostic only",
        constraint_type=ConstraintType.EXCLUDED_SPECIES,
        check_fn=_check_excluded_species,
        provenance_ticket='RYA-240 / RYA-558',
    ),
    RatifiedConstraint(
        constraint_id='Fe_1D_3D_correction_required_on_solar_report',
        scope="every element with a reported-layer 1D→3D correction registered in "
              "`config/corrections_registry.yaml` (today: solar Fe) — the reported "
              "anchor must be on the post-correction scale, declared, and applied once",
        constraint_type=ConstraintType.REQUIRED_CORRECTION,
        check_fn=_check_fe_reported_scale,
        provenance_ticket='RYA-553 (hardened RYA-681 / RYA-674)',
        applies_to=(RowKind.ELEMENT_VALUE,),
    ),
)


def constraint_by_id(constraint_id: str) -> RatifiedConstraint:
    for c in RATIFIED_CONSTRAINTS:
        if c.constraint_id == constraint_id:
            return c
    raise KeyError(f"no ratified constraint {constraint_id!r}; registered: "
                   f"{[c.constraint_id for c in RATIFIED_CONSTRAINTS]}")


def check_ratified_constraints(candidate_output: Iterable[Any], module_name: str,
                               row_kind: str = RowKind.ELEMENT_VALUE) -> list[str]:
    """Every violation in `candidate_output`, as human-readable strings. Never raises
    on a violation — that is `assert_ratified_constraints_satisfied`'s job. Exposed so
    a report can ENUMERATE violations (e.g. auditing a historical artifact) without
    halting on the first."""
    if row_kind not in (RowKind.ELEMENT_VALUE, RowKind.SPECIES_RECORD):
        raise ValueError(f"unknown row_kind {row_kind!r}")
    violations: list[str] = []
    for raw in candidate_output:
        row = _normalise(raw, row_kind)
        if not row.element:
            continue                       # a summary/header row carries no element
        for c in RATIFIED_CONSTRAINTS:
            if not c.active or row_kind not in c.applies_to:
                continue
            reason = c.check_fn(row)
            if reason:
                violations.append(
                    f"[{module_name}] {c.constraint_id} ({c.constraint_type}, ratified "
                    f"{c.provenance_ticket}) VIOLATED by row {row.summary()}: {reason}")
    return violations


def assert_ratified_constraints_satisfied(candidate_output: Iterable[Any],
                                          module_name: str,
                                          row_kind: str = RowKind.ELEMENT_VALUE) -> None:
    """THE emission-time gate. Call before writing any per-element result.

    `candidate_output` is the rows about to be written; `module_name` is the emitting
    stage, so the failure says who was about to write it. Loud-fail on any violation —
    never a warning, never a filtered row (RYA-518). Fixing a violation means fixing the
    upstream decision or obtaining a ratified reversal, not exempting the caller.
    """
    violations = check_ratified_constraints(candidate_output, module_name, row_kind)
    if not violations:
        return
    raise RatifiedConstraintViolation(
        f"RYA-674: {len(violations)} ratified constraint violation(s) in the candidate "
        f"output of {module_name} — REFUSING TO WRITE.\n  - "
        + "\n  - ".join(violations)
        + "\n\nA ratified constraint is not a suggestion: it was ratified by the ticket "
          "named on each line. Do not exempt this call site — fix the value, mark the "
          "row `diagnostic_only` where the ratifying ticket specifies that demotion, or "
          "obtain an explicit Ryan-ratified reversal in a new ticket "
          "(docs/SCIENCE_STANDARDS.md § Ratified Constraints).")
