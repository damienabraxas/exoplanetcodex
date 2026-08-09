#!/usr/bin/env python3
"""
pipeline/reliability_contract.py — the `reliable` flag's CONSUMER-side vocabulary
================================================================================
RYA-691 (rule ratified RYA-679; emission-time enforcement ratified RYA-699).

The rule itself lives in `scripts/solar_profile_fit.assess_reliability`:

    reliable = (not railed) AND dEW_dA >= RELIABLE_DEWDA          # RYA-679

`red_chi2` is reported and review-flagged but does NOT gate — RYA-679 measured that the
full-window statistic tracks how crowded the window is, not how good the element's
measurement is. **Nothing in this module re-derives that rule.** This is the consumer
side: what a downstream reader must do with a flag somebody else already computed.

Why the vocabulary is a module and not three f-strings
-----------------------------------------------------
RYA-691 made `scripts/rya527_two_engine_run.py` honour the flag at every read, raising
where it is present-and-false. That fixed the one consumer that had eight ungated reads.
It could not fix the *next* consumer, because the contract lived inside a script's
private helper — a module that never imported it could not violate it loudly, only
quietly.

Ryan ratified this as a RYA-674 constraint on 2026-08-08, which requires the check to
run at EMISSION time in any module. So the two basis strings the producer writes and the
classifier the emission gate reads have to share one definition. That is all this file
is: `reliability_basis()` builds the string, `classify_reliability_basis()` reads it, and
they cannot drift because they are eight lines apart.

Three states, and why "no flag" is not a failure
------------------------------------------------
`GATED`   — the artifact carried `reliable=True`. The measurement cleared RYA-679.
`UNGATED` — the artifact carries no reliability flag **and the record says why**. This is
            a real and legitimate state: the RYA-491/237 CNO cross-arm artifact is a
            multi-indicator reconciliation, not a profile fit, so it has no `dEW_dA` and
            no `railed` to test. Forcing a uniform check over artifacts with genuinely
            different semantics would fabricate agreement rather than find it (RYA-691
            §3A said so explicitly).
`DEMOTED` — the artifact carried `reliable=False`. Never emitted, ever.

The `absent_reason` is mandatory, not optional. An UNGATED basis with no reason is itself
a violation: "no flag" is legitimate only when the record states in words which artifact
lacked the flag and why that is correct for its shape. Without the reason, `UNGATED` is
indistinguishable from an emitter that simply never looked — which is the RYA-691 defect
restated.
"""
from __future__ import annotations

__all__ = [
    'ReliabilityState', 'RELIABILITY_RULE_TICKET', 'UNGATED', 'GATED_PREFIX',
    'reliability_basis', 'classify_reliability_basis', 'RELIABILITY_BASIS_KEYS',
]

#: The ticket that ratified the rule itself. Named in every message so a reader lands on
#: the measured argument (RYA-679 §"Why no red_chi2 ceiling") rather than on this file.
RELIABILITY_RULE_TICKET = 'RYA-679'

UNGATED = 'UNGATED'
GATED_PREFIX = f'{RELIABILITY_RULE_TICKET} reliability-gated'

#: The keys an emitted row may use to carry its basis. Listed here rather than in the
#: check, for the same reason `_VALUE_KEYS` is: emission paths have different schemas and
#: unifying them is not this ticket's job, so every spelling is named in ONE place.
RELIABILITY_BASIS_KEYS = ('engineB_reliability', 'reliability_basis', 'reliability')


class ReliabilityState:
    GATED = 'GATED'
    UNGATED = 'UNGATED'
    DEMOTED = 'DEMOTED'
    #: An UNGATED basis carrying no reason — see the module docstring.
    UNREASONED = 'UNREASONED'


def reliability_basis(gated: bool | None, *, key: str = 'reliable',
                      absent_reason: str | None = None,
                      detail: str | None = None) -> str:
    """Build the basis string recorded alongside a value.

    `gated=None` means the artifact carried no flag; `absent_reason` is then REQUIRED.
    Producers should not hand-write these strings — round-tripping through
    `classify_reliability_basis` is what makes the emission gate meaningful.

    `detail` replaces the default `key=True` tail for a gated basis whose evidence is
    not one boolean. Co I and Zr II are gated on a COUNT ("3 of 5 fitted lines cleared")
    and Mg I on a conjunction ("target_reliable=True, lines_concordant=True") — all
    genuinely RYA-679-gated, none expressible as a single flag. The distinction the
    classifier cares about is gated-vs-not, so the tail is free text; what it must not
    be is a string outside the vocabulary entirely.
    """
    if gated is None:
        if not (absent_reason or '').strip():
            raise ValueError(
                "an UNGATED reliability basis requires an absent_reason: 'no flag' is a "
                "legitimate state only when the record says which artifact lacked the "
                "flag and why that is correct for its shape (RYA-691)")
        return f"{UNGATED} — {absent_reason}"
    if not gated:
        raise ValueError(
            "a demoted measurement has no basis string — it must not be emitted at all "
            f"({RELIABILITY_RULE_TICKET}/RYA-691)")
    return f"{GATED_PREFIX}: {detail or f'{key}=True'}"


def classify_reliability_basis(basis) -> str:
    """Read a recorded basis back to a `ReliabilityState`.

    Deliberately NOT a substring sniff over prose — RYA-695 measured what that costs
    (a row whose text *explained* that it declined the floor value was failed for
    containing the phrase 'two-engine'). Classification keys on the leading token the
    builder writes, so a basis that merely *mentions* a demotion in its reason text is
    read as what it is.
    """
    text = str(basis or '').strip()
    if not text:
        return ReliabilityState.UNREASONED
    if text.startswith(UNGATED):
        tail = text[len(UNGATED):].lstrip(' —-:').strip()
        return ReliabilityState.UNGATED if tail else ReliabilityState.UNREASONED
    if text.startswith(GATED_PREFIX):
        return ReliabilityState.GATED
    # Anything else is an unrecognised vocabulary: a hand-written string, or a producer
    # that predates this module. Treated as UNREASONED — refused rather than assumed
    # good, because assuming good is the whole defect.
    return ReliabilityState.UNREASONED
