"""
tests/test_fe_1d3d_idempotency_rya681.py
========================================
THE HOLE RYA-681 CLOSES: nothing in this repo distinguished A(Fe I) 7.466 from
7.416.

RYA-669 substituted a double-corrected verdict for the committed one and ran the
gates. `tests/test_solar_calibration_gate.py` returned **9 passed**; the RYA-632
ledger consistency guard reported **0 undocumented failures**. Both green on a
verdict whose Fe anchor is 0.05 dex wrong, because:

  * `FE_GATE = [7.410, 7.510]` is a symmetric ±0.05 window on the 3D-true reference
    and the 1D→3D correction is also 0.05 — so a doubled correction lands at
    `A_correct − 0.05`, inside the window for any correct anchor ≥ 7.46;
  * RYA-632 compares verdicts, tiers and counts, never values (by design);
  * the RYA-553 idempotency guard keyed on the prose `method_scale` label, which the
    gold BUILDER writes, so a freeze that moved `A_X` without the label re-armed it.

Every test below is a direction of that hole, and the file is written so the
**failing** direction is exercised explicitly rather than argued for in a comment.
"""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from config.constants import CORRECTIONS_3D, SOLAR_ASPLUND2021          # noqa: E402
from pipeline.solar_scale_provenance import (                            # noqa: E402
    SCALE_1D_NLTE, SCALE_3D_NLTE, SCALE_STATE_COLUMN, ScaleProvenanceError,
    apply_reported_scale_correction, correction_dex, declared_scale,
    method_scale_label, resolve_gold_scale, scale_centres,
    scale_discrimination_halfwidth, scale_from_value)

DEX = CORRECTIONS_3D['Fe_1D3D_solar_dex']          # −0.05
A_CORRECT = 7.466        # the reported 3D-scale anchor (gold v3's VALUE — correct)
A_DOUBLED = 7.416        # RYA-669's double-corrected anchor
A_UNCORRECTED = 7.516    # the 1D-NLTE measurement (gold v2's value)

VERDICT = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'


# ── the gold-row shapes this defect is made of ───────────────────────────────
# v2: 1D value under a 1D label       → self-consistent, correction DUE
GOLD_V2_FE = {'element': 'Fe', 'A_X': A_UNCORRECTED, 'method_scale': '1D-NLTE (Fe I)'}
# v3: 3D value under a 1D label       → SELF-CONTRADICTORY (RYA-665 froze this)
GOLD_V3_FE = {'element': 'Fe', 'A_X': A_CORRECT, 'method_scale': '1D-NLTE (Fe I)'}
# v4-shape: 3D value, 3D label, and the explicit machine-readable declaration
GOLD_V4_FE = {'element': 'Fe', 'A_X': A_CORRECT,
              'method_scale': method_scale_label('Fe', SCALE_3D_NLTE),
              SCALE_STATE_COLUMN: SCALE_3D_NLTE}


def _fe_verdict_row(a_measured, *, applied=True, pre=None):
    """A verdict-shaped Fe row at a chosen anchor, as phase_c would emit it."""
    pre = pre if pre is not None else round(a_measured - DEX, 3)
    return {'element': 'Fe', 'asplund2021': SOLAR_ASPLUND2021['Fe'],
            'A_measured': a_measured,
            'fe_1d3d_correction': {
                'applied': applied, 'correction_dex': DEX if applied else 0.0,
                'a_1dnlte_pre': pre, 'a_3dnlte_post': a_measured,
                'scale': SCALE_3D_NLTE}}


# ── §2C — THE REGRESSION TEST, BOTH DIRECTIONS ───────────────────────────────

def _gate_a1_on(fe_row):
    """Run the RYA-166 A1 gate logic against an arbitrary Fe verdict row.

    Imports the real gate module and calls its assertion with `_fe_verdict`
    monkeypatched, so this exercises the SHIPPED gate, not a copy of it.
    """
    import tests.test_solar_calibration_gate as gate
    orig = gate._fe_verdict
    gate._fe_verdict = lambda: fe_row
    try:
        gate.test_A1_a_fe_in_fe_gate()
    finally:
        gate._fe_verdict = orig


def test_gate_PASSES_on_the_correct_anchor_7466():
    """Direction 1 — the correct anchor still clears, comfortably."""
    _gate_a1_on(_fe_verdict_row(A_CORRECT))


def test_gate_FAILS_on_the_double_corrected_anchor_7416():
    """Direction 2 — THE HOLE. This exact substitution returned 9 passed before
    RYA-681. It must now fail, and fail on the SCALE, not on the range."""
    with pytest.raises(AssertionError, match=r'not on the 3D-NLTE scale'):
        _gate_a1_on(_fe_verdict_row(A_DOUBLED))


def test_gate_still_FAILS_on_the_uncorrected_anchor_7516():
    """The dropped-correction teeth FE_GATE already had are not lost."""
    with pytest.raises(AssertionError):
        _gate_a1_on(_fe_verdict_row(A_UNCORRECTED, applied=False, pre=A_UNCORRECTED))


def test_the_doubled_anchor_is_INSIDE_fe_gate_which_is_why_range_cannot_catch_it():
    """The premise of the fix, asserted rather than assumed: FE_GATE genuinely
    cannot see 7.416, so the added check had to be orthogonal to it."""
    from config.constants import FE_GATE_LOWER, FE_GATE_UPPER
    lo = SOLAR_ASPLUND2021['Fe'] + FE_GATE_LOWER
    hi = SOLAR_ASPLUND2021['Fe'] + FE_GATE_UPPER
    assert lo <= A_DOUBLED <= hi, "premise broken: 7.416 is meant to be inside FE_GATE"
    assert abs(DEX) == pytest.approx(hi - SOLAR_ASPLUND2021['Fe']), (
        "the correction size equals the gate half-width — that identity is the whole "
        "reason a symmetric window cannot distinguish correct from double-corrected")


def test_the_committed_verdict_is_the_correct_anchor():
    """The shipped artifact must be the 7.466 one — a 7.416 verdict must never land."""
    if not VERDICT.exists():
        pytest.skip('solar verdict artifact absent')
    fe = next(v for v in json.loads(VERDICT.read_text())['verdicts']
              if v['element'] == 'Fe')
    assert float(fe['A_measured']) == pytest.approx(A_CORRECT, abs=1e-9)
    assert scale_from_value('Fe', float(fe['A_measured'])) == SCALE_3D_NLTE


# ── §2A — the re-keyed idempotency guard ─────────────────────────────────────

def test_correction_applies_exactly_once_on_a_self_consistent_1d_row():
    a, rec = apply_reported_scale_correction('Fe', A_UNCORRECTED, GOLD_V2_FE)
    assert a == pytest.approx(A_CORRECT)
    assert rec['applied'] is True and rec['scale'] == SCALE_3D_NLTE


def test_correction_is_SKIPPED_on_a_self_consistent_3d_row():
    """The v4 shape: value and declaration both say 3D → no second subtraction."""
    a, rec = apply_reported_scale_correction('Fe', A_CORRECT, GOLD_V4_FE)
    assert a == pytest.approx(A_CORRECT)
    assert rec['applied'] is False
    assert rec['gold_scale_source'] == SCALE_STATE_COLUMN   # data, not prose


def test_guard_is_idempotent_under_repeated_application():
    """Feed the output back in, as the gold↔verdict loop does. It must converge,
    not ratchet. Under the OLD label-keyed guard this walked 7.516 → 7.466 → 7.416
    → 7.366 …"""
    row = dict(GOLD_V2_FE)
    a = A_UNCORRECTED
    a, rec = apply_reported_scale_correction('Fe', a, row)
    # the builder now writes the state back as DATA, from the correction record
    row = {'element': 'Fe', 'A_X': a,
           'method_scale': method_scale_label('Fe', rec['scale']),
           SCALE_STATE_COLUMN: rec['scale']}
    for _ in range(5):
        a, rec = apply_reported_scale_correction('Fe', a, row)
        assert a == pytest.approx(A_CORRECT), "the correction ratcheted"
        assert rec['applied'] is False


def test_gold_v3_shape_RAISES_instead_of_double_correcting():
    """THE DEFECT. A 3D value under a 1D label used to silently re-arm the
    correction and emit 7.416; it must now be a loud refusal to load."""
    with pytest.raises(ScaleProvenanceError, match=r'CONTRADICTS ITSELF'):
        apply_reported_scale_correction('Fe', A_CORRECT, GOLD_V3_FE)


def test_an_already_doubled_anchor_RAISES_as_on_neither_scale():
    with pytest.raises(ScaleProvenanceError, match=r'NEITHER recognised abundance scale'):
        apply_reported_scale_correction(
            'Fe', A_DOUBLED, {'element': 'Fe', 'A_X': A_DOUBLED,
                              'method_scale': '3D-NLTE (Fe I, Magic 2013)'})


def test_an_undeclared_scale_RAISES_rather_than_defaulting_to_1d():
    """The old guard's fallback was `''`, and `'3D' not in ''` is True — an absent
    label silently meant "apply the correction". Now it means "stop"."""
    with pytest.raises(ScaleProvenanceError, match=r'declares NO abundance scale'):
        apply_reported_scale_correction('Fe', A_CORRECT,
                                        {'element': 'Fe', 'A_X': A_CORRECT})


def test_explicit_scale_state_beats_the_prose_label():
    state, source = declared_scale(GOLD_V4_FE)
    assert (state, source) == (SCALE_3D_NLTE, SCALE_STATE_COLUMN)
    state, source = declared_scale(GOLD_V2_FE)
    assert state == SCALE_1D_NLTE and 'method_scale' in source


def test_explicit_column_contradicting_the_value_also_RAISES():
    """The explicit column is authoritative but NOT unchecked — a mis-stamped
    scale_state is caught by the same value cross-check."""
    bad = {'element': 'Fe', 'A_X': A_CORRECT, 'method_scale': '3D-NLTE',
           SCALE_STATE_COLUMN: SCALE_1D_NLTE}
    with pytest.raises(ScaleProvenanceError, match=r'CONTRADICTS ITSELF'):
        resolve_gold_scale('Fe', bad, A_CORRECT)


# ── the classifier's own properties (no tuned tolerance) ─────────────────────

def test_discrimination_boundary_is_derived_from_the_correction_not_tuned():
    half = scale_discrimination_halfwidth('Fe')
    assert half == pytest.approx(abs(DEX) / 2.0)
    c = scale_centres('Fe')
    assert c[SCALE_3D_NLTE] == pytest.approx(SOLAR_ASPLUND2021['Fe'])
    assert c[SCALE_1D_NLTE] == pytest.approx(SOLAR_ASPLUND2021['Fe'] + abs(DEX))
    # the two bands are exactly adjacent: no value is ever on both scales
    assert c[SCALE_1D_NLTE] - c[SCALE_3D_NLTE] == pytest.approx(2 * half)


@pytest.mark.parametrize('value,expected', [
    (A_UNCORRECTED, SCALE_1D_NLTE),
    (A_CORRECT, SCALE_3D_NLTE),
    (A_DOUBLED, None),
    (7.510, SCALE_1D_NLTE),
    (7.460, SCALE_3D_NLTE),
    (7.760, None),            # gross loggf slip
])
def test_scale_from_value_classification(value, expected):
    assert scale_from_value('Fe', value) == expected


def test_correct_anchor_clears_the_scale_check_comfortably():
    """Confirm the fix does not leave 7.466 on a knife edge: it sits 0.006 dex from
    the 3D centre with a ±0.025 discrimination half-width — ~4× margin — while the
    doubled anchor misses by 0.019."""
    half = scale_discrimination_halfwidth('Fe')
    margin = half - abs(A_CORRECT - SOLAR_ASPLUND2021['Fe'])
    assert margin == pytest.approx(0.019, abs=1e-9) and margin > 0
    assert abs(A_DOUBLED - SOLAR_ASPLUND2021['Fe']) - half == pytest.approx(0.019, abs=1e-9)


# ── the builder end of the loop (the desync SOURCE) ──────────────────────────

def test_builder_derives_the_fe_label_from_the_verdict_not_a_literal():
    """`build_solar_reference_v2_rya522._scale_and_note` hardcoded '1D-NLTE (Fe I)'
    for Fe regardless of the verdict — that literal is why the label could never
    follow the value. It must now track the verdict's own correction record."""
    import build_solar_reference_v2_rya522 as B
    v3d = _fe_verdict_row(A_CORRECT)
    scale, _note = B._scale_and_note('Fe', 'EW', 'PASS', A_CORRECT,
                                     SOLAR_ASPLUND2021['Fe'], 'gold', vrow=v3d)
    assert '3D' in scale.upper() and '1D' not in scale.upper(), (
        f"builder still stamps a 1D label on a 3D anchor: {scale!r}")
    assert B._scale_state_from_verdict('Fe', v3d) == SCALE_3D_NLTE


def test_builder_REFUSES_a_verdict_whose_record_contradicts_its_value():
    import build_solar_reference_v2_rya522 as B
    lying = {'element': 'Fe', 'A_measured': A_DOUBLED,
             'fe_1d3d_correction': {'applied': True, 'correction_dex': DEX,
                                    'a_1dnlte_pre': A_CORRECT,
                                    'a_3dnlte_post': A_DOUBLED, 'scale': SCALE_3D_NLTE}}
    with pytest.raises(ScaleProvenanceError):
        B._scale_state_from_verdict('Fe', lying)


def test_builder_REFUSES_a_verdict_with_no_scale_record():
    import build_solar_reference_v2_rya522 as B
    with pytest.raises(ScaleProvenanceError, match=r'UNDECLARED'):
        B._scale_state_from_verdict('Fe', {'element': 'Fe', 'A_measured': A_CORRECT})


# ── phase_c against the ACTUAL frozen gold (the standing blocker) ────────────

def test_phase_c_loud_fails_on_the_live_frozen_gold_row():
    """Gold v3 as frozen is self-contradictory (3D value, 1D label). phase_c must
    now REFUSE to regenerate rather than silently emit 7.416. This test is expected
    to keep passing until a ratified v4 re-freeze fixes the label (RYA-669) — at
    which point it flips to the skip branch below, which is the intended resolution.
    """
    from pipeline import data_namespace as ns
    try:
        ab, version = ns.read_solar_reference('CURRENT')
    except Exception as e:                                    # pragma: no cover
        pytest.skip(f'gold reference unreadable: {e}')
    fe = ab[ab['element'] == 'Fe'].iloc[0]
    a_x = float(fe['A_X'])
    declared, _src = declared_scale(fe)
    from_value = scale_from_value('Fe', a_x)
    if declared == from_value:
        pytest.skip(f'gold {version} Fe row is self-consistent ({declared}) — '
                    f'the RYA-669 re-freeze has landed')
    with pytest.raises(ScaleProvenanceError, match=r'CONTRADICTS ITSELF'):
        resolve_gold_scale('Fe', fe, a_x)


def test_correction_magnitude_stays_single_sourced():
    from config.constants import FE_1D3D_SOLAR_OFFSET
    assert correction_dex('Fe') == pytest.approx(-FE_1D3D_SOLAR_OFFSET)
