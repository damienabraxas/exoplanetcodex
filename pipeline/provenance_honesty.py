"""
pipeline/provenance_honesty.py  (RYA-653 — generalizes the RYA-596 tripwire)
===========================================================================
THE INVARIANT: **no stage emits or inherits a causal/provenance claim it did not
itself establish.**

RYA-596 found `scripts/phase_c_verdict_rya371.py` asserting ONE cause for every
blank row — "no independent-gf line survives the graded cull" — from a
fallthrough that never checked it. The claim was false for Ca/Ti/Ni/Na/Al/Sr,
which all carry graded survivors and DID produce a curated A(X); their blank is
the ratified RYA-522 `owed`-tier HOLD round-tripping through the frozen gold
(RYA-469), not a cull.

The RYA-632 audit then found the IDENTICAL unchecked fallthrough one file over,
in `scripts/build_solar_reference_v2_rya522.py::_scale_and_note` — where it is
actively lying in gold v2's **Ba** row against a Ba that RYA-559 has since
measured (A(Ba) 2.410).

This module is the SINGLE SOURCE OF TRUTH for that tripwire. Both emitters import
it; neither owns a copy. (A forked second copy would be the very duplication the
invariant forbids — the same defect class as the four-way rest-frame/gsig fork
RYA-643 had to single-source.)

The tell is always on the row itself: a zero-survivor claim sitting next to
`n_lines > 0`, or next to a measured value. Both are contradictions on their
face, and both are made unrepresentable here.
"""
from __future__ import annotations

import math

# The claim itself — the thing that may only be said about a genuinely empty row.
ZERO_SURVIVOR_CLAIM = 'no independent-gf line survives the graded cull'
# The phase_c verdict channel that carries it (RYA-596 spelling, preserved verbatim:
# it is written into every generation of the verdict artifact).
ZERO_SURVIVOR_CHANNEL = f'EW present; {ZERO_SURVIVOR_CLAIM}'


def claims_zero_survivors(cause) -> bool:
    """True if `cause` asserts that the graded cull left nothing.

    Substring, not equality: the gold builder writes the bare claim while the
    verdict writes it behind an "EW present; " prefix. One claim, two spellings —
    the check must not be fooled by either.
    """
    return isinstance(cause, str) and ZERO_SURVIVOR_CLAIM in cause


def _survivor_count(n_lines) -> int:
    """n_lines as an int. Missing/NaN (pandas, or a verdict row without the key)
    counts as 0 — an absent count is not evidence of survivors, and this check
    must never manufacture a red from a blank cell."""
    if n_lines is None:
        return 0
    try:
        n = float(n_lines)
    except (TypeError, ValueError):
        return 0
    return 0 if math.isnan(n) else int(n)


def _is_measured(a_measured) -> bool:
    if a_measured is None:
        return False
    try:
        return not math.isnan(float(a_measured))
    except (TypeError, ValueError):
        return False


def assert_blank_cause_is_honest(element, cause, n_lines, a_measured=None,
                                 site='phase_c verdict'):
    """RYA-596/653 tripwire: refuse to emit a zero-survivor cause on a row that
    contradicts it.

    A row may claim the graded cull left nothing ONLY if it actually carries
    nothing: no graded survivors AND no measured value. Anything else is a cause
    the emitting stage did not establish.

    `a_measured` is optional so the phase_c call site keeps its original,
    n_lines-only contract (RYA-596) unchanged; the gold builder passes it,
    because the Ba defect is exactly the measured-value shape (n_lines=1 AND
    A(Ba)=2.410 sitting next to "no line survives").
    """
    if not claims_zero_survivors(cause):
        return
    n = _survivor_count(n_lines)
    measured = _is_measured(a_measured)
    if not (n > 0 or measured):
        return
    contradiction = []
    if n > 0:
        contradiction.append(f"n_lines={n} graded survivor(s)")
    if measured:
        contradiction.append(f"a measured A({element})={float(a_measured):.3f}")
    raise AssertionError(
        f"RYA-596/653 [{site}]: {element} row claims '{ZERO_SURVIVOR_CLAIM}' but "
        f"carries {' and '.join(contradiction)}. A blank on such a row is a "
        f"gold-tier HOLD (RYA-522) or an unadopted measurement — never a cull. "
        f"The emitting stage must state the cause it actually established.")
