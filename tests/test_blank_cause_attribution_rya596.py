"""
tests/test_blank_cause_attribution_rya596.py
============================================
RYA-596 — a blank A(X) in the solar verdict must name its ACTUAL cause.

THE BUG THIS LOCKS OUT
----------------------
`_classify`'s in-pool fallthrough used to assert ONE cause for every blank row:
"EW present; no independent-gf line survives the graded cull" — blaming the
RYA-398/456 graded-gf firewall. That claim was FALSE for Ca/Ti/Ni/Na/Al/Sr, all
of which DO have graded survivors (2/10/2/2/1/1 lines) and DID produce a curated
A(X) — visible in gold v1 (Ca 6.324, Ti 5.471, Ni 6.946, Na 6.264, Al 7.406).

Their real cause is the RYA-522 gold tiering: an `owed`-tier row freezes NO value
BY RATIFIED DESIGN (Ryan 2026-07-05 — "suspect → held, not immortalised"), and
this verdict reads the frozen gold back in (RYA-469). The blank is a deliberate
HOLD round-tripping through the gold, not a cull.

The tell was on the row itself the whole time: `n_lines > 0` next to "no line
survives". These tests make that contradiction unrepresentable.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as P  # noqa: E402

# The five the ticket names, plus Sr — same mechanism, same phantom cause.
HELD_WITH_SURVIVORS = ('Ca', 'Ti', 'Ni', 'Na', 'Al', 'Sr')


def _rows():
    ab, ew, phase_a = P._load()
    rows = P.build_verdicts(ab, ew, phase_a)
    P._apply_kittpeak(rows, P._load_kittpeak())
    return {r['element']: r for r in rows}


def test_no_row_claims_zero_survivors_while_carrying_survivors():
    """THE INVARIANT. A row may claim the graded cull left nothing only if it
    actually carries nothing. build_verdicts also raises on this, so a regression
    fails loudly at generation time rather than shipping a false verdict."""
    for el, r in _rows().items():
        if r['channel'] == P.ZERO_SURVIVOR_CHANNEL:
            assert r['n_lines'] == 0, (
                f"{el} claims '{P.ZERO_SURVIVOR_CHANNEL}' but carries "
                f"n_lines={r['n_lines']} graded survivor(s) — that blank is a "
                f"gold-tier HOLD (RYA-522), not a cull.")


def test_the_five_are_held_not_culled():
    """Ca/Ti/Ni/Na/Al (+Sr) are blank because the ratified `owed` tier withholds
    their value — NOT because the graded-gf firewall emptied their pools."""
    rows = _rows()
    for el in HELD_WITH_SURVIVORS:
        r = rows[el]
        assert r['A_measured'] is None, f"{el} unexpectedly carries a frozen value"
        assert r['n_lines'] > 0, f"{el} should have graded survivors"
        assert P.ZERO_SURVIVOR_CHANNEL not in r['channel'], (
            f"{el} still blames the graded cull for a tier hold")
        assert 'HELD' in r['channel'] and 'RYA-522' in r['channel'], (
            f"{el} channel must name the gold-tier hold: {r['channel']!r}")


def test_genuine_zero_survivor_rows_keep_the_cull_explanation():
    """The counter-case must NOT be over-corrected: Mg/Y/Zr/Eu really do end the
    graded cull with zero surviving lines, and must still say so."""
    rows = _rows()
    for el in ('Mg', 'Y', 'Zr', 'Eu'):
        r = rows[el]
        assert r['n_lines'] == 0, f"{el} is expected to have zero graded survivors"
        assert r['channel'] == P.ZERO_SURVIVOR_CHANNEL, (
            f"{el} lost its honest zero-survivor explanation: {r['channel']!r}")


def test_tripwire_raises_on_a_contradictory_row():
    """The generation-time guard is real, not decorative."""
    with pytest.raises(AssertionError, match='RYA-596'):
        P._assert_blank_cause_is_honest('Ca', P.ZERO_SURVIVOR_CHANNEL, n_lines=2)
    # and stays quiet on the two legitimate shapes
    P._assert_blank_cause_is_honest('Mg', P.ZERO_SURVIVOR_CHANNEL, n_lines=0)
    P._assert_blank_cause_is_honest('Ca', "EW: 2 curated line(s); value HELD", n_lines=2)


def test_blank_is_never_silent():
    """RYA-596 smoke test: every blank row states an explicit, non-empty cause —
    no silent blanks in either direction."""
    for el, r in _rows().items():
        if r['A_measured'] is None:
            assert r['channel'] and r['owed'], f"{el} is blank with no stated cause"
            assert r['verdict'] != 'PASS', f"{el} cannot PASS with no value"
