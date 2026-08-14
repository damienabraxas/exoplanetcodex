#!/usr/bin/env python3
"""
pipeline/solar_scale_provenance.py — RYA-681
============================================
THE ABUNDANCE SCALE A FROZEN VALUE IS ON IS *DATA*, NOT PROSE.

Why this module exists
----------------------
RYA-553 applies the tabulated Magic-2013 1D→3D solar Fe correction (−0.05 dex) at
the REPORTED-value layer, inside `scripts/phase_c_verdict_rya371.py`. That script
reads its Fe anchor from the FROZEN gold reference (RYA-469 `read_solar_reference`)
and the gold reference is itself REBUILT FROM the verdict channel
(`scripts/build_solar_reference_v2_rya522.py`). So the correction sits inside a
closed loop:

    gold vN  ──►  phase_c verdict  ──►  gold candidate  ──►  freeze  ──►  gold vN+1

Every trip round that loop must apply the correction EXACTLY ONCE, at the trip where
the anchor is still on the 1D-NLTE scale, and NEVER again. RYA-553's guard asserted
that idempotency against a PROSE LABEL — `method_scale` — which a *different* process
writes:

    if '3D' not in scale.upper():            # scripts/phase_c_verdict_rya371.py
        a_meas = round(a_meas + dex, 3)

and the builder that writes that label hardcoded it:

    if el == "Fe":                           # scripts/build_solar_reference_v2_rya522.py
        return "1D-NLTE (Fe I)", "..."

so the label could never follow the value. RYA-665 froze gold v3 with the POST-
correction number `A_X 7.466` under the PRE-correction label `1D-NLTE (Fe I)`, the
guard read the label, saw no `3D`, and re-armed. RYA-669 measured the result: a fresh
phase_c run lands A(Fe I) at **7.416**, and every gate stays green
(`test_solar_calibration_gate.py` 9 passed, RYA-632 ledger guard 0 undocumented).
That is not a one-off — it is a RATCHET: left alone, each freeze cycle subtracts
another 0.05 dex.

The fix (three parts, all here)
-------------------------------
1. **An explicit, machine-readable declaration.** `SCALE_STATE_COLUMN` is a real
   column on the gold row carrying a canonical token (`1D-NLTE` / `3D-NLTE`). It is
   written by the gold builder FROM the verdict's own `fe_1d3d_correction` record —
   data copied from data, never a hand-maintained string. It is the authority.

2. **A value-side cross-check that needs no column at all.** The two scale states are
   separated by a KNOWN, tabulated distance (`FE_1D3D_SOLAR_OFFSET`), and both scale
   centres are independently published quantities (Asplund 2021 3D-true 7.46, and
   7.46 + 0.05 for 1D-NLTE — the very relationship gold's own Fe note asserts). So
   "which scale is this number actually on" is DECIDABLE from the number: the
   hypotheses are non-overlapping, and the decision boundary sits at their midpoint.
   No new free parameter — the discrimination half-width is |correction| / 2,
   single-sourced from the correction itself.

3. **Refuse to proceed on contradiction.** When declaration and value disagree — gold
   v3's exact shape — this module RAISES. A silent re-arm becomes a loud load failure.
   It generalises: any element that ever acquires a reported-layer scale correction
   gets the same treatment for free.

The correction is now idempotent with respect to THE VALUE ITSELF, which no separate
process can desynchronise.

What this deliberately does NOT do
----------------------------------
It does not "fix" gold v3. Gold v3 is write-once and frozen (RYA-469; sha256
47ad869e…83421) and the wrong cell is its `method_scale` label, not its value —
7.466 IS the correct 3D-scale number. Re-freezing a v4 is a ratification call
(RYA-669). Until that lands, `resolve_gold_scale('Fe', <gold v3 row>)` raises, and
phase_c cannot regenerate. That is the intended, honest outcome: the channel is
blocked LOUDLY instead of quietly emitting 7.416.

RYA-674 — the declaration generalises from one string to a LIST
--------------------------------------------------------------
`scale_state` answers exactly one question ("which of two Fe scales is this on?"). The
question the project actually has is "which tabulated corrections does this number
already carry?", and an element may one day carry more than one. So the authoritative
declaration is now `corrections_applied` — a JSON list of identifiers from
`config/corrections_registry.yaml` — and `scale_state` is DERIVED from it rather than
stored independently. One stored fact, two views: two stored facts is how RYA-669
happened. The reader below prefers, in order:

    corrections_applied  (RYA-674, authoritative)
  > scale_state          (RYA-681, still written; corroborating)
  > method_scale prose   (legacy frozen rows only)

and RAISES if any two of them disagree.
"""
from __future__ import annotations

import json
from typing import Any, Mapping, NamedTuple

import numpy as np

from config.constants import CORRECTIONS_3D, SOLAR_ASPLUND2021

__all__ = [
    'ScaleProvenanceError',
    'SCALE_STATE_COLUMN',
    'CORRECTIONS_APPLIED_COLUMN',
    'FE_1D3D_CORRECTION_ID',
    'classify_gold_scale',
    'try_apply_reported_scale_correction',
    'ScaleContradiction',
    'SCALE_SELF_CONTRADICTION',
    'SCALE_OFF_SCALE',
    'SCALE_UNDECLARED',
    'SCALE_1D_NLTE',
    'SCALE_3D_NLTE',
    'REPORTED_SCALE_CORRECTED_ELEMENTS',
    'correction_dex',
    'scale_centres',
    'scale_discrimination_halfwidth',
    'scale_from_value',
    'read_corrections_applied',
    'encode_corrections_applied',
    'corrections_applied_for_state',
    'scale_state_from_corrections',
    'declared_scale',
    'resolve_gold_scale',
    'apply_reported_scale_correction',
    'method_scale_label',
]


class ScaleProvenanceError(RuntimeError):
    """A frozen row's abundance-scale provenance is missing, ambiguous or
    self-contradictory. Raised INSTEAD of guessing — guessing is what produced
    the RYA-669 double correction."""


# The explicit machine-readable declaration. Written by the gold builder from the
# verdict's own correction record; read here in preference to any prose label.
SCALE_STATE_COLUMN = 'scale_state'

# RYA-674 — THE authoritative declaration: a JSON list of correction identifiers from
# config/corrections_registry.yaml naming every tabulated correction this number
# already carries. `[]` means "none applied" and is a POSITIVE statement; a missing
# column means "this row predates the schema" and is not the same thing (the old ''
# fallback silently meant "apply", which is how RYA-669 happened).
CORRECTIONS_APPLIED_COLUMN = 'corrections_applied'

# The one correction that currently exists at the reported layer. Named here so the
# gold builder and the guards agree on the token without re-typing it; the token's
# MEANING lives in config/corrections_registry.yaml, not here.
FE_1D3D_CORRECTION_ID = '1D_3D_solar_Fe_Magic2013'

SCALE_1D_NLTE = '1D-NLTE'
SCALE_3D_NLTE = '3D-NLTE'
_CANONICAL_STATES = (SCALE_1D_NLTE, SCALE_3D_NLTE)

# Elements carrying a tabulated 1D→3D correction at the REPORTED-value layer.
# Keyed off CORRECTIONS_3D so registering another element needs no edit here.
_REPORTED_KEY_SUFFIX = '_1D3D_solar_dex'
REPORTED_SCALE_CORRECTED_ELEMENTS = tuple(sorted(
    k[: -len(_REPORTED_KEY_SUFFIX)] for k in CORRECTIONS_3D
    if k.endswith(_REPORTED_KEY_SUFFIX)))


def correction_dex(element: str) -> float | None:
    """The tabulated reported-layer 1D→3D correction for `element`, or None."""
    return (float(CORRECTIONS_3D[f'{element}{_REPORTED_KEY_SUFFIX}'])
            if f'{element}{_REPORTED_KEY_SUFFIX}' in CORRECTIONS_3D else None)


def scale_centres(element: str) -> dict[str, float]:
    """Where each scale's anchor SHOULD sit, from published quantities only.

    3D-NLTE centre = Asplund 2021 3D-true A(X).
    1D-NLTE centre = that, MINUS the (negative) tabulated correction, i.e. the
    3D-true value plus the granulation term. Both terms are external literature —
    never our own output — so the centres are not circular.
    """
    dex = correction_dex(element)
    if dex is None:
        raise ScaleProvenanceError(
            f"{element} carries no reported-layer 1D→3D correction in CORRECTIONS_3D")
    three_d = float(SOLAR_ASPLUND2021[element])
    return {SCALE_3D_NLTE: three_d, SCALE_1D_NLTE: round(three_d - dex, 6)}


def scale_discrimination_halfwidth(element: str) -> float:
    """Half the correction magnitude — the midpoint between the two scale centres.

    NOT a tuned tolerance and NOT a gate: it is the geometric decision boundary of a
    two-hypothesis classification. Bands of this half-width around the two centres are
    exactly adjacent and never overlap, so a value can never be assigned to both, and
    the boundary moves automatically if the tabulated correction is ever revised.
    """
    dex = correction_dex(element)
    if dex is None:
        raise ScaleProvenanceError(
            f"{element} carries no reported-layer 1D→3D correction in CORRECTIONS_3D")
    return abs(dex) / 2.0


def scale_from_value(element: str, a_x: float) -> str | None:
    """Which scale the NUMBER itself says it is on, or None if it is on neither.

    Returns SCALE_1D_NLTE / SCALE_3D_NLTE when `a_x` lies within the discrimination
    half-width of that scale's centre. Returns None when the value sits more than
    that from BOTH centres — an anchor that has drifted off every recognised scale,
    which callers must treat as loud-fail territory rather than assume a state.
    """
    if a_x is None or not np.isfinite(float(a_x)):
        return None
    half = scale_discrimination_halfwidth(element)
    best, best_d = None, float('inf')
    for state, centre in scale_centres(element).items():
        d = abs(float(a_x) - centre)
        if d <= half and d < best_d:
            best, best_d = state, d
    return best


def _canonicalise(token: Any) -> str | None:
    """Map a declaration/label to a canonical state token, or None if it says nothing.

    Accepts the explicit column's canonical tokens and, for LEGACY frozen rows that
    predate the column, the `method_scale` prose (e.g. '1D-NLTE (Fe I)'). Prose is a
    FALLBACK read, never the authority — it is exactly what desynchronised.
    """
    if token is None:
        return None
    s = str(token).strip()
    if not s or s.lower() in ('nan', 'none'):
        return None
    u = s.upper()
    has_3d, has_1d = '3D' in u, '1D' in u
    if has_3d and not has_1d:
        return SCALE_3D_NLTE
    if has_1d and not has_3d:
        return SCALE_1D_NLTE
    return None


def read_corrections_applied(row: Mapping[str, Any]) -> list[str] | None:
    """The `corrections_applied` declaration on one row, or None if it has none.

    RYA-674. Accepts the JSON-encoded string a CSV round-trips and the real list a
    freshly-built candidate row carries. Distinguishes THREE states, deliberately:

      * ``None``  — the row carries no `corrections_applied` column at all (a legacy
        frozen row). "Undeclared" is not "none applied".
      * ``[]``    — the row positively declares that no correction is baked in.
      * ``[...]`` — the identifiers whose corrections this number already carries.

    A present-but-unparseable cell RAISES: a declaration we cannot read is worse than
    no declaration, because it looks like one.
    """
    if row is None or not hasattr(row, 'get'):
        return None
    raw = row.get(CORRECTIONS_APPLIED_COLUMN)
    if raw is None:
        return None
    if isinstance(raw, (list, tuple)):
        return [str(x) for x in raw]
    if isinstance(raw, float) and np.isnan(raw):
        return None
    s = str(raw).strip()
    if not s or s.lower() in ('nan', 'none'):
        return None
    try:
        parsed = json.loads(s)
    except ValueError as exc:
        raise ScaleProvenanceError(
            f"{CORRECTIONS_APPLIED_COLUMN}={raw!r} is not JSON ({exc}). A correction "
            f"declaration that cannot be read is not a declaration (RYA-674).") from exc
    if not isinstance(parsed, list):
        raise ScaleProvenanceError(
            f"{CORRECTIONS_APPLIED_COLUMN}={raw!r} decoded to {type(parsed).__name__}, "
            f"expected a JSON list of correction identifiers (RYA-674)")
    return [str(x) for x in parsed]


def encode_corrections_applied(identifiers) -> str:
    """The canonical CSV cell for a `corrections_applied` list — sorted JSON, so two
    rows carrying the same corrections are byte-identical."""
    return json.dumps(sorted({str(x) for x in (identifiers or [])}))


def corrections_applied_for_state(element: str, state: str) -> list[str]:
    """The correction identifiers a row on `state` must declare, for `element`.

    Derived from `config/corrections_registry.yaml` (scope + post_scale), never a
    literal list here: registering a second reported-layer correction for an element
    needs no edit to this module.
    """
    if state not in _CANONICAL_STATES:
        raise ScaleProvenanceError(f"{state!r} is not a canonical scale token")
    from pipeline.corrections_registry import correction, corrections_for_element
    return sorted(cid for cid in corrections_for_element(element)
                  if str(correction(cid).get('post_scale')) == state)


def scale_state_from_corrections(element: str, applied) -> str | None:
    """The canonical scale state implied by a `corrections_applied` list, or None.

    The list is the stored fact; the state is a VIEW of it. `[]` for a
    reported-layer-corrected element positively means the pre-correction scale.
    """
    if applied is None:
        return None
    have = set(applied)
    post = set(corrections_applied_for_state(element, SCALE_3D_NLTE))
    if not post:
        return None                     # element carries no reported-layer correction
    if post & have:
        return SCALE_3D_NLTE
    return SCALE_1D_NLTE


def declared_scale(row: Mapping[str, Any]) -> tuple[str | None, str]:
    """(canonical state, where it came from) for one gold row.

    Reads, in order of authority: the RYA-674 `corrections_applied` list, the RYA-681
    `scale_state` column, and — only for frozen rows that predate both — the
    `method_scale` prose. Where a row carries more than one of these they must AGREE:
    a disagreement is the RYA-669 shape one layer in, and RAISES rather than resolving
    to whichever the reader happens to prefer.
    """
    if row is None:
        return None, 'absent'
    element = str(row.get('element') or '') if hasattr(row, 'get') else ''
    readings: list[tuple[str, str]] = []          # (state, source)

    applied = read_corrections_applied(row)
    if applied is not None and element:
        state = scale_state_from_corrections(element, applied)
        if state is not None:
            readings.append((state, CORRECTIONS_APPLIED_COLUMN))

    explicit = row.get(SCALE_STATE_COLUMN) if hasattr(row, 'get') else None
    state = _canonicalise(explicit)
    if state is not None:
        if state not in _CANONICAL_STATES:                      # pragma: no cover
            raise ScaleProvenanceError(
                f"{SCALE_STATE_COLUMN}={explicit!r} is not a canonical scale token")
        readings.append((state, SCALE_STATE_COLUMN))

    label = row.get('method_scale') if hasattr(row, 'get') else None
    state = _canonicalise(label)
    if state is not None:
        readings.append((state, 'method_scale (legacy prose fallback)'))

    if not readings:
        return None, 'absent'
    distinct = {s for s, _ in readings}
    if len(distinct) > 1:
        detail = '; '.join(f"{src} says {st}" for st, src in readings)
        raise ScaleProvenanceError(
            f"the {element or '?'} row's own scale declarations DISAGREE — {detail}. "
            f"`{CORRECTIONS_APPLIED_COLUMN}` is authoritative and the others are views of "
            f"it, so they cannot legitimately differ; a row that states two scales states "
            f"none (RYA-674).")
    return readings[0][0], readings[0][1]


#: RYA-815 — the contradiction kinds, named so a caller can localise on them.
SCALE_OFF_SCALE = "off-scale"
SCALE_UNDECLARED = "undeclared"
SCALE_SELF_CONTRADICTION = "self-contradiction"


class ScaleContradiction(NamedTuple):
    """A reference self-contradiction, described rather than raised.

    RYA-681 made these LOUD, which was right — a silent 0.05 dex/cycle ratchet is
    far worse. RYA-815 makes the loudness LOCAL: the caller can turn ONE bad row
    into ONE INDETERMINATE element instead of vetoing all 28. Nothing is softened;
    `resolve_gold_scale` still raises, and this carries the same message verbatim.
    """
    kind: str
    element: str
    message: str
    cell: str            # the offending cell, NAMED (e.g. "Fe.method_scale")
    declared: str = ""
    from_value: str = ""


def classify_gold_scale(element: str, row: Mapping[str, Any], a_x: float):
    """(state, None) when the row is self-consistent, else (None, ScaleContradiction).

    The non-raising twin of `resolve_gold_scale`. Both share this one body so the
    two can never drift apart — a second copy of the contradiction logic is exactly
    how a guard comes to disagree with itself.
    """
    from_value = scale_from_value(element, a_x)
    centres = scale_centres(element)
    half = scale_discrimination_halfwidth(element)
    if from_value is None:
        return None, ScaleContradiction(
            kind=SCALE_OFF_SCALE, element=element, cell=f"{element}.A_X",
            message=(
                f"A({element}) = {a_x} is on NEITHER recognised abundance scale "
                f"(1D-NLTE centre {centres[SCALE_1D_NLTE]:.3f}, 3D-NLTE centre "
                f"{centres[SCALE_3D_NLTE]:.3f}, discrimination half-width ±{half:.3f}). "
                f"Refusing to guess whether the RYA-553 1D→3D correction has been applied. "
                f"Either the anchor has moved materially or the correction has been applied "
                f"more than once (RYA-681/669)."))
    declared, source = declared_scale(row)
    if declared is None:
        return None, ScaleContradiction(
            kind=SCALE_UNDECLARED, element=element,
            cell=f"{element}.{SCALE_STATE_COLUMN}|method_scale",
            from_value=from_value,
            message=(
                f"the frozen {element} row declares NO abundance scale: it carries neither a "
                f"'{SCALE_STATE_COLUMN}' column nor a scale-bearing 'method_scale' label. "
                f"A frozen anchor must state the scale it is on (RYA-681)."))
    if declared != from_value:
        return None, ScaleContradiction(
            kind=SCALE_SELF_CONTRADICTION, element=element,
            cell=f"{element}.{source}", declared=declared, from_value=from_value,
            message=(
                f"the frozen {element} row CONTRADICTS ITSELF: it is labelled '{declared}' "
                f"(from {source}) but its value A({element}) = {a_x} sits on the "
                f"'{from_value}' scale (centre {centres[from_value]:.3f}, ±{half:.3f}; the "
                f"'{declared}' centre is {centres[declared]:.3f}). This is the RYA-669 "
                f"defect: a post-correction value frozen under a pre-correction label, which "
                f"re-arms the RYA-553 1D→3D correction and double-applies it. Refusing to "
                f"load. Fixing the label requires a re-freeze (RYA-669 / RYA-527), NOT a "
                f"code exemption — gold is write-once (RYA-469)."))
    return declared, None


def resolve_gold_scale(element: str, row: Mapping[str, Any], a_x: float) -> str:
    """THE authoritative scale state of a frozen gold row. RAISES on contradiction.

    Unchanged behaviour (RYA-681): three ways this fails LOUDLY rather than guessing —
    the value on neither scale, nothing declaring a scale, or declaration and value
    disagreeing (gold v3's exact shape, the state that used to silently re-arm the
    RYA-553 correction and produce 7.416).

    RYA-815 note: this now DELEGATES to `classify_gold_scale` and raises what it
    reports. The refusal is identical; a caller that wants to localise the failure to
    one element calls the classifier instead. Sharing one body is deliberate — two
    copies of this logic would eventually disagree, and a guard that disagrees with
    itself is worse than no guard.
    """
    state, contradiction = classify_gold_scale(element, row, a_x)
    if contradiction is not None:
        raise ScaleProvenanceError(contradiction.message)
    return state


def try_apply_reported_scale_correction(element: str, a_x: float,
                                        row: Mapping[str, Any]):
    """Non-raising twin. Returns `(a_out, record, contradiction_or_None)`.

    RYA-815: lets a caller keep going for the other 27 elements when ONE reference
    row is self-contradictory. On contradiction the value is returned UNCHANGED and
    UNCORRECTED — deliberately: the whole point of RYA-681 is that we must not emit
    a number we cannot vouch for, so the caller is expected to withhold the verdict
    (INDETERMINATE), not to publish `a_x` as though it were fine.
    """
    state, contradiction = classify_gold_scale(element, a_x=a_x, row=row)
    if contradiction is not None:
        return round(float(a_x), 3), {
            'applied': False,
            'reason': f'reference self-contradiction ({contradiction.kind})',
            'contradiction_kind': contradiction.kind,
            'contradiction_cell': contradiction.cell,
            'contradiction_message': contradiction.message,
            'correction_dex': 0.0,
            'scale': None,
            'gold_scale_state': None,
            'gold_scale_source': None,
        }, contradiction
    a_out, record = apply_reported_scale_correction(element, a_x, row)
    return a_out, record, None


def apply_reported_scale_correction(element: str, a_x: float,
                                    row: Mapping[str, Any]) -> tuple[float, dict]:
    """Apply the reported-layer 1D→3D correction to a gold anchor, exactly once.

    Returns `(a_out, record)`. Idempotency is keyed on the VALUE (cross-checked
    against the row's explicit declaration), so no separate label-maintaining
    process can desynchronise it. Raises ScaleProvenanceError rather than emitting
    a number it cannot vouch for.
    """
    dex = correction_dex(element)
    if dex is None:
        raise ScaleProvenanceError(
            f"{element} carries no reported-layer 1D→3D correction in CORRECTIONS_3D")
    state = resolve_gold_scale(element, row, a_x)
    a_pre = round(float(a_x), 3)
    _, source = declared_scale(row)
    if state == SCALE_1D_NLTE:
        a_out = round(float(a_x) + dex, 3)
        record = {
            'applied': True,
            'source': 'Magic et al. 2013 (Stagger 3D, A&A 557 A26)',
            'correction_dex': dex,
            'a_1dnlte_pre': a_pre,
            'a_3dnlte_post': a_out,
            'scale': SCALE_3D_NLTE,
            'gold_scale_state': state,
            'gold_scale_source': source,
            # RYA-674 §2B outcome 3: the value was on the pre-correction scale, so the
            # correction is applied AND its identifier is added to the declaration the
            # freeze will carry. The identifier list is derived from the registry.
            'corrections_applied': corrections_applied_for_state(element, SCALE_3D_NLTE),
        }
        return a_out, record
    record = {
        'applied': False,
        'reason': f'gold anchor already on the {state} scale '
                  f'(declared via {source}, corroborated by the value)',
        'correction_dex': 0.0,
        'a_1dnlte_pre': a_pre,
        'a_3dnlte_post': a_pre,
        'scale': state,
        'gold_scale_state': state,
        'gold_scale_source': source,
        # RYA-674 §2B outcome 1 (idempotent skip) — the declaration the value already
        # satisfies, carried forward unchanged so the next freeze restates it.
        'corrections_applied': corrections_applied_for_state(element, state),
    }
    return a_pre, record


def method_scale_label(element: str, state: str) -> str:
    """The `method_scale` prose for a gold row, DERIVED from the canonical state.

    The gold builder must call this instead of hardcoding a string, so the prose can
    never disagree with the value it labels (RYA-681). Prose stays human-readable;
    the machine reads `scale_state`.
    """
    if state not in _CANONICAL_STATES:
        raise ScaleProvenanceError(f"{state!r} is not a canonical scale token")
    ion = 'I'
    if state == SCALE_3D_NLTE:
        return f'3D-NLTE ({element} {ion}, Magic 2013)'
    return f'1D-NLTE ({element} {ion})'
