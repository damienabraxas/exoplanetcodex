"""
tests/gold_scale_blocker.py — RYA-681 / RYA-669
===============================================
THE ONE PLACE that states why the phase_c verdict cannot currently be regenerated.

Gold v3 is self-contradictory: its Fe row holds the POST-correction value 7.466
under the PRE-correction label `1D-NLTE (Fe I)` (RYA-665 moved `A_X` without moving
`method_scale`). Before RYA-681 the RYA-553 idempotency guard read that label, saw
no `3D`, and re-applied the Magic-2013 −0.05 dex correction — regenerating the
verdict at A(Fe I) **7.416**, with every gate green (RYA-669).

RYA-681 re-keys that guard on the VALUE, so the contradiction is now a loud refusal
to load (`ScaleProvenanceError`) instead of a silent second subtraction. The correct
consequence is that **anything which re-runs `build_verdicts` against the live frozen
gold cannot run at all** until the label is fixed — and fixing it means re-freezing a
gold v4, which is a ratification call (RYA-669 / RYA-527), not something code may do:
gold is write-once (RYA-469).

So the tests below are blocked BY THE DATA, not broken by the change. They are
marked xfail with this reason WHILE the contradiction stands, and the marking is
CONDITIONAL on the live gold row — the moment a self-consistent v4 is promoted, the
condition evaporates and every one of them runs for real again. Nothing is
suppressed unconditionally, and no test is allowed to pass on a doubled anchor: the
verdict-artifact assertions in `tests/test_fe_1d3d_idempotency_rya681.py` and the
RYA-166 gate both run unblocked and both fail on 7.416.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import data_namespace as ns                      # noqa: E402
from pipeline.solar_scale_provenance import (                  # noqa: E402
    REPORTED_SCALE_CORRECTED_ELEMENTS, declared_scale, scale_from_value)


def regeneration_blocked() -> str | None:
    """Why `build_verdicts` cannot run against the live frozen gold, or None.

    Returns a human-readable reason iff some reported-layer scale-corrected row in
    the CURRENT gold reference contradicts itself (label vs value). Any failure to
    READ the gold returns None — this helper exists to explain a known blocker, not
    to swallow unrelated breakage.
    """
    try:
        ab, version = ns.read_solar_reference('CURRENT')
    except Exception:
        return None
    for el in REPORTED_SCALE_CORRECTED_ELEMENTS:
        sub = ab[ab['element'] == el]
        if sub.empty:
            continue
        row = sub.iloc[0]
        try:
            a_x = float(row['A_X'])
        except (TypeError, ValueError):
            continue
        declared, source = declared_scale(row)
        from_value = scale_from_value(el, a_x)
        if declared is not None and from_value is not None and declared != from_value:
            return (
                f"BLOCKED by the frozen gold reference, not by this test: gold {version}'s "
                f"{el} row is self-contradictory — A({el}) = {a_x} is on the {from_value} "
                f"scale but the row is labelled {declared} (via {source}). phase_c now "
                f"REFUSES to regenerate rather than re-apply the RYA-553 1D→3D correction "
                f"to an already-corrected anchor (RYA-669 measured that at 7.416). "
                f"Clearing this needs a ratified gold v4 whose label matches its value "
                f"(RYA-669/RYA-527) — gold is write-once (RYA-469), so no code change can "
                f"do it. This mark is CONDITIONAL on the live row and retires itself the "
                f"moment a self-consistent gold is promoted.")
    return None


def xfail_if_regeneration_blocked() -> None:
    """Call at the top of any test that re-runs the verdict off the live gold."""
    reason = regeneration_blocked()
    if reason:
        pytest.xfail(reason)
