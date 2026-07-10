"""
tests/test_two_engine_floor_rya525.py
=====================================
The two-engine floor selector (pipeline/engine_selection.py) — the ratified,
reference-BLIND, per-line → per-element criterion (RYA-525 §2–3 + the three
ratification comments + the 2026-07-10 appendix). Encode-don't-tune: these pin the
LAW on constructed lines, never tuned on real solar data.
"""
import numpy as np
import pytest

from config.constants import TWO_ENGINE
from pipeline import engine_selection as es
from pipeline.engine_selection import (
    LineEngines, ENGINE_A, ENGINE_B, CLEAN_WEAK, HARD, INDETERMINATE,
    classify_regime, select_line, aggregate_element, select_element,
    TwoEngineError, ReferenceProximityError,
)

KNEE = TWO_ENGINE['saturation_knee_mA']


def _line(**kw):
    base = dict(wavelength=5000.0, species='Si I', a_value=7.50, a_err=0.05,
                b_value=7.55, b_err=0.05, b_chi2=1.0, ew_mA=30.0,
                blend_flag=False, is_problem_child=False, a_in_hull=True)
    base.update(kw)
    return LineEngines(**base)


# ── clause 3/4 regime classification (reference-blind) ────────────────────────
class TestRegime:
    def test_clean_weak(self):
        assert classify_regime(_line(ew_mA=KNEE - 10)) == CLEAN_WEAK

    def test_saturated_is_hard(self):
        assert classify_regime(_line(ew_mA=KNEE + 10)) == HARD

    def test_blended_is_hard(self):
        assert classify_regime(_line(ew_mA=20.0, blend_flag=True)) == HARD

    def test_problem_child_is_hard(self):
        assert classify_regime(_line(ew_mA=20.0, is_problem_child=True)) == HARD

    def test_hfs_element_is_hard(self):
        assert classify_regime(_line(species='Mn I', ew_mA=20.0)) == HARD  # Mn ∈ hfs_elements

    def test_unknown_ew_is_indeterminate_not_clean(self):
        # absence of hard flags is NOT clean — a clean call needs a positive unsaturated EW
        assert classify_regime(_line(ew_mA=None)) == INDETERMINATE


# ── clause 1 validity gates ───────────────────────────────────────────────────
class TestEligibility:
    def test_synth_chi2_gate_disqualifies_engineB(self):
        # B ineligible (bad fit) → A wins even on a HARD line
        w = select_line(_line(ew_mA=KNEE + 50, b_chi2=TWO_ENGINE['synth_chi2_gate'] + 5))
        assert w.engine == ENGINE_A and 'only Engine-A eligible' in w.reason

    def test_out_of_hull_disqualifies_engineA(self):
        w = select_line(_line(ew_mA=10.0, a_in_hull=False))   # clean-weak but A out of hull
        assert w.engine == ENGINE_B and 'only Engine-B eligible' in w.reason

    def test_neither_eligible_returns_none(self):
        assert select_line(_line(a_value=None, b_value=None)) is None


# ── clause 2/3/4 selection ────────────────────────────────────────────────────
class TestSelection:
    def test_clause3_clean_weak_picks_engineA(self):
        w = select_line(_line(ew_mA=20.0))
        assert w.engine == ENGINE_A and w.regime == CLEAN_WEAK

    def test_clause3_hard_picks_engineB(self):
        w = select_line(_line(ew_mA=KNEE + 40))
        assert w.engine == ENGINE_B and w.regime == HARD

    def test_clause4_indeterminate_picks_lower_sigma(self):
        w = select_line(_line(ew_mA=None, a_err=0.09, b_err=0.03))
        assert w.regime == INDETERMINATE and w.engine == ENGINE_B and 'lower line-scatter' in w.reason

    def test_clause4_exact_tie_defaults_to_anchor_engineA(self):
        w = select_line(_line(ew_mA=None, a_err=0.05, b_err=0.05))
        assert w.engine == ENGINE_A and 'anchor-scale default' in w.reason

    def test_rejected_engine_and_delta_recorded(self):
        w = select_line(_line(ew_mA=20.0, a_value=7.50, b_value=7.55))
        assert w.engine == ENGINE_A and w.rejected_engine == ENGINE_B
        assert abs(w.rejected_value - 7.55) < 1e-9
        assert abs(w.cross_engine_delta - (7.55 - 7.50)) < 1e-9   # b − a


# ── the tuning firewall (adversarial: reference-closer ≠ quality) ─────────────
class TestTuningFirewall:
    def test_reference_closer_engine_does_not_win_when_it_is_lower_quality(self):
        # Construct: a CLEAN-WEAK line. Engine-B happens to sit closer to a hypothetical
        # reference (say 7.60), but the line is clean-weak so the PHYSICS winner is
        # Engine-A. The selector must pick Engine-A — it never sees the reference.
        w = select_line(_line(ew_mA=20.0, a_value=7.48, b_value=7.60))
        assert w.engine == ENGINE_A   # physics (clean-weak), not proximity to 7.60

    def test_hard_line_picks_synthesis_even_if_EW_is_reference_closer(self):
        # A blended/saturated line where the EW (Engine-A) value is closer to a reference:
        # physics says HARD → synthesis; proximity is irrelevant.
        w = select_line(_line(ew_mA=KNEE + 30, blend_flag=True, a_value=7.50, b_value=7.30))
        assert w.engine == ENGINE_B

    def test_assert_reference_blind_raises_on_smuggled_reference_key(self):
        with pytest.raises(ReferenceProximityError):
            es.assert_reference_blind(['ew_mA', 'blend_flag', 'delta_vs_asplund'])
        es.assert_reference_blind(['ew_mA', 'blend_flag', 'b_chi2'])   # clean set: no raise


# ── element aggregation + cross-engine-mix guard ──────────────────────────────
class TestAggregation:
    def test_inverse_variance_combine_of_winners(self):
        # two clean-weak lines, both won by Engine-A; equal errs → simple mean
        lines = [_line(wavelength=5000, ew_mA=20, a_value=7.40, a_err=0.10),
                 _line(wavelength=5001, ew_mA=20, a_value=7.60, a_err=0.10)]
        rec = select_element('Si I', lines)
        assert rec.selected_engines == (ENGINE_A,)
        assert abs(rec.value - 7.50) < 1e-9
        assert not rec.cross_engine_mix and not rec.mix_flagged

    def test_budget_uses_only_winning_line_errors(self):
        # inverse-variance err of two 0.10 measurements = 0.10/sqrt(2)
        lines = [_line(wavelength=5000, ew_mA=20, a_value=7.5, a_err=0.10),
                 _line(wavelength=5001, ew_mA=20, a_value=7.5, a_err=0.10)]
        rec = select_element('Si I', lines)
        assert abs(rec.err - 0.10 / np.sqrt(2)) < 1e-6

    def test_mix_agreeing_engines_not_flagged(self):
        # one clean-weak (A wins) + one hard (B wins); engines AGREE within the gate → no flag
        lines = [_line(wavelength=5000, ew_mA=20, a_value=7.50, b_value=7.52),
                 _line(wavelength=5001, ew_mA=KNEE + 20, blend_flag=True, a_value=7.51, b_value=7.53)]
        rec = select_element('Si I', lines)
        assert rec.cross_engine_mix and not rec.mix_flagged

    def test_mix_disagreeing_engines_flagged_the_ti_lesson(self):
        # winners span both engines AND |mean Δ| > gate (≈0.11, the Ti split) → adjudicate
        d = TWO_ENGINE['cross_engine_mix_gate'] + 0.05
        lines = [_line(wavelength=5000, ew_mA=20, a_value=7.50, b_value=7.50 + d),
                 _line(wavelength=5001, ew_mA=KNEE + 20, blend_flag=True, a_value=7.50, b_value=7.50 + d)]
        rec = select_element('Ti I', lines)
        assert rec.cross_engine_mix and rec.mix_flagged

    def test_empty_raises_never_silent(self):
        with pytest.raises(TwoEngineError):
            aggregate_element('Si I', [None, None])


# ── loud-fail guards (RYA-525 §3) ─────────────────────────────────────────────
class TestGuards:
    def test_missing_grid_disposition_raises(self):
        for disp in ('acquire-task', 'build-task'):
            with pytest.raises(TwoEngineError):
                es.require_grid_or_raise('V', disp, {ENGINE_A})

    def test_wired_both_missing_an_engine_raises(self):
        with pytest.raises(TwoEngineError):
            es.require_grid_or_raise('Si', 'wired-both', {ENGINE_A})   # B failed to run

    def test_wired_both_with_both_engines_ok(self):
        es.require_grid_or_raise('Si', 'wired-both', {ENGINE_A, ENGINE_B})   # no raise

    def test_wired_one_single_engine_ok(self):
        es.require_grid_or_raise('Fe', 'wired-one', {ENGINE_A})   # documented one-engine → no raise

    def test_single_engine_no_record_undocumented_raises(self):
        rec = select_element('Zz I', [_line(ew_mA=20, a_value=7.5, b_value=None)])
        assert rec.selected_engines == (ENGINE_A,) and rec.mean_cross_engine_delta is None
        with pytest.raises(TwoEngineError):
            es.assert_cross_engine_recorded('Zz', rec, 'wired-both')   # undocumented single → raise

    def test_single_engine_documented_disposition_ok(self):
        rec = select_element('Zr I', [_line(ew_mA=20, a_value=7.5, b_value=None)])
        es.assert_cross_engine_recorded('Zr', rec, 'LTE-only-by-design')   # documented → no raise


# ── all thresholds come from config ───────────────────────────────────────────
def test_thresholds_sourced_from_config_not_inline():
    # moving the saturation knee reclassifies a borderline line — proves config-driven
    cfg = dict(TWO_ENGINE)
    cfg['saturation_knee_mA'] = 25.0
    assert classify_regime(_line(ew_mA=30.0), cfg) == HARD          # 30 > 25 → saturated
    cfg['saturation_knee_mA'] = 50.0
    assert classify_regime(_line(ew_mA=30.0), cfg) == CLEAN_WEAK    # 30 ≤ 50 → clean
