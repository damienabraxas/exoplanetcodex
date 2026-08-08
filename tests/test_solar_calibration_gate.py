"""
tests/test_solar_calibration_gate.py
====================================
THE Beta-blocking solar calibration gate (RYA-166 Layer 3; gates RYA-527).

ONE command answers "is Solar Beta-ready":
    pytest tests/test_solar_calibration_gate.py -v

Each assertion below is a single pass/fail check, sourced from COMMITTED artifacts
+ ratified, cited thresholds (never a threshold loosened to fit a run — validate-
don't-tune). It aggregates checks that otherwise live in separate scripts:

  A1  3D-corrected A(Fe I) in the real FE_GATE [7.41, 7.51] (RYA-553: Magic-2013
      1D→3D solar Fe correction applied to the reported anchor; the old scale-aware
      diagnostic band [7.44,7.58] is retired here)
  A2  Fe I-vs-Fe II ionization balance <= FE_IONISATION_GATE (RYA-406 synth arbiter)
  A3  Fe I reduced-EW (vmic) slope |slope| < FE_REW_SLOPE_GATE (RYA-330)
  A4  Fe I RAW line-to-line scatter <= 0.1398 honest floor (RYA-407, wired RYA-446)
  A5  Fe I sigma_REPORTED <= 0.05 dex — the SE+systematics budget (RYA-158/166 Step
      2.5), a DIFFERENT, smaller quantity than A4's raw scatter (guarded below)
  A6  RYA-400 physics-regime audit run() == 0 (zero failures), wired IN-process
  A7  RYA-429 no-silent-drop invariant holds for the full solar linelist
  A8  Every TARGET element carries an honest verdict class (no silent absence)

The two 'sigma' checks (A4 raw scatter vs A5 reported) are DELIBERATELY separate and
must not be conflated — that conflation is the exact bug the RYA-166 rewrite fixes.
"""
import contextlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from config.constants import (ACCEPTANCE_PROFILES, FE_IONISATION_GATE,          # noqa: E402
                              FE_IONIZATION_SYNTH_ARBITER, FE_REW_SLOPE_GATE,
                              FE_GATE_LOWER, FE_GATE_UPPER, CORRECTIONS_3D,
                              TARGET_ELEMENTS)
from pipeline.solar_scale_provenance import (  # noqa: E402  RYA-681 scale identity
    SCALE_1D_NLTE, SCALE_3D_NLTE, scale_discrimination_halfwidth, scale_from_value)

# ── committed artifacts (guard these, never re-run the pipeline in the gate) ──
VERDICT = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'
FE1_RESID = ROOT / 'data' / 'audit' / 'fe1_scatter' / 'fe1_per_line_residuals_rya407.csv'
UNCERT = ROOT / 'data' / 'audit' / 'uncertainty' / 'solar_uncertainty_rya158.json'
LEDGER = ROOT / 'data' / 'results' / 'rejection_ledger_solar_rya429.json'

ASPLUND_FE = 7.46
# RYA-553: the tabulated Magic-2013 1D→3D solar Fe correction is now APPLIED to the
# reported Fe I anchor (in the verdict), so the anchor sits on the true 3D scale and
# FE_GATE [7.41,7.51] is again the real solar reported-value gate — the scale-aware
# diagnostic band [7.44,7.58] is retired here (it still describes the un-corrected
# 1D-NLTE measurement layer / off-solar stars, RYA-550).
FE_GATE_LO = ASPLUND_FE + FE_GATE_LOWER    # 7.41
FE_GATE_HI = ASPLUND_FE + FE_GATE_UPPER    # 7.51
SIGMA_REPORTED_GATE = 0.05     # Asplund 2021 reported-uncertainty target (RYA-158/166)
VALID_VERDICTS = {'PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP'}


def _verdict():
    if not VERDICT.exists():
        pytest.skip(f"solar verdict artifact absent: {VERDICT.name}")
    return json.loads(VERDICT.read_text())


def _fe_verdict():
    return next(v for v in _verdict()['verdicts'] if v['element'] == 'Fe')


def _fe1_pool():
    if not FE1_RESID.exists():
        pytest.skip(f"Fe I residual pool absent: {FE1_RESID.name}")
    return pd.read_csv(FE1_RESID)


def _uncertainty_fe():
    if not UNCERT.exists():
        pytest.skip(f"uncertainty budget absent ({UNCERT.name}); run "
                    "`python -m pipeline.uncertainty_stack --star solar` first")
    d = json.loads(UNCERT.read_text())
    return next(r for r in d['per_element'] if r['element'] == 'Fe')


# ── A1 — 3D-corrected A(Fe I) on the real FE_GATE (RYA-553) ──────────────────

def test_A1_a_fe_in_fe_gate():
    """The reported Fe I anchor is the 3D-corrected value (Magic-2013 1D→3D applied,
    RYA-553); it must sit inside the real FE_GATE [7.41,7.51]. A gross zero-point error
    (e.g. a loggf slip ~+0.3 → 7.76) still falls outside, so the gate keeps its teeth.

    RYA-681 — FE_GATE ALONE IS NOT ENOUGH, and narrowing it is not the answer.
    FE_GATE is a symmetric ±0.05 window on the 3D-true reference: a *physical*
    acceptance tolerance. The 1D→3D correction is also 0.05, i.e. exactly the gate's
    half-width, so a DOUBLED correction lands at `A_correct − 0.05`, which is inside
    the window for any correct anchor at or above 7.46. That is structural, not a
    coincidence of today's numbers — and it is why RYA-669's 7.416 passed all nine of
    these assertions. Narrowing the window to h < 0.044 would catch today's instance
    while leaving the correct anchor only 0.038 dex of headroom, and would still miss
    a smaller doubled correction: it fits the coincidence, not the defect.

    The defect is a PROVENANCE error, so the added check is a provenance check. The
    two scale states are separated by a known, tabulated distance, so "which scale is
    this number on" is decidable without any tuned tolerance (see
    pipeline.solar_scale_provenance). The reported anchor must be identifiably on the
    3D scale: 7.466 → 3D ✓, 7.416 → NEITHER ✗, 7.516 → 1D ✗. FE_GATE is unchanged.
    """
    fe = _fe_verdict()
    a_fe = float(fe['A_measured'])
    corr = fe.get('fe_1d3d_correction')
    assert corr is not None, (
        "Fe verdict carries no fe_1d3d_correction record — the reported anchor's scale "
        "is undeclared (RYA-553/681)")
    assert corr.get('correction_dex') in (CORRECTIONS_3D['Fe_1D3D_solar_dex'], 0.0), (
        "Fe 1D→3D correction magnitude drifted from CORRECTIONS_3D (single-source broken)")
    assert corr.get('scale') == SCALE_3D_NLTE, (
        f"Fe verdict declares scale {corr.get('scale')!r}, not {SCALE_3D_NLTE} — the "
        f"reported anchor is not on the 3D scale (RYA-553)")

    # RYA-681 scale-identity: the NUMBER itself must land on the 3D scale.
    on_scale = scale_from_value('Fe', a_fe)
    assert on_scale == SCALE_3D_NLTE, (
        f"reported A(Fe I)={a_fe} is not on the 3D-NLTE scale (classified {on_scale!r}; "
        f"3D centre {ASPLUND_FE}, 1D centre {ASPLUND_FE + abs(CORRECTIONS_3D['Fe_1D3D_solar_dex']):.3f}, "
        f"discrimination ±{scale_discrimination_halfwidth('Fe'):.3f}). A doubled RYA-553 "
        f"correction lands here (7.416, RYA-669) and FE_GATE cannot see it.")

    # …and the correction record must be arithmetically self-consistent: the value it
    # was applied TO must itself be on the 1D scale. With a doubled correction the
    # recorded pre-value is 7.466, which is a 3D-scale number — caught here too.
    if corr.get('applied') is True:
        pre, post = float(corr['a_1dnlte_pre']), float(corr['a_3dnlte_post'])
        assert round(pre + CORRECTIONS_3D['Fe_1D3D_solar_dex'], 3) == post, (
            f"fe_1d3d_correction is not self-consistent: {pre} {CORRECTIONS_3D['Fe_1D3D_solar_dex']:+} "
            f"!= {post}")
        assert scale_from_value('Fe', pre) == SCALE_1D_NLTE, (
            f"the RYA-553 correction was applied to {pre}, which is NOT on the 1D-NLTE "
            f"scale (classified {scale_from_value('Fe', pre)!r}) — the correction has been "
            f"applied to an already-corrected anchor (RYA-669/681)")
        assert post == pytest.approx(a_fe, abs=1e-9), (
            f"the verdict's reported A(Fe I)={a_fe} is not the correction's own post-value {post}")

    assert FE_GATE_LO <= a_fe <= FE_GATE_HI, (
        f"3D-corrected A(Fe I)={a_fe} outside FE_GATE [{FE_GATE_LO:.2f},{FE_GATE_HI:.2f}]")


# ── A2 — Fe I/II ionization balance on the synth arbiter (RYA-406) ───────────

def test_A2_ionization_balance_within_gate():
    dfe = FE_IONIZATION_SYNTH_ARBITER['solar']['dFe']       # ratified, RYA-341/406
    assert abs(dfe) <= FE_IONISATION_GATE, (
        f"|dFe(I-II)|={abs(dfe)} > FE_IONISATION_GATE {FE_IONISATION_GATE}")


# ── A3 — Fe I reduced-EW (vmic) slope (RYA-330) ──────────────────────────────

def test_A3_rew_slope_flat_within_tolerance():
    d = _fe1_pool()
    slope = float(np.polyfit(d['REW'], d['a_1dlte'], 1)[0])
    assert abs(slope) < FE_REW_SLOPE_GATE, (
        f"|REW slope|={abs(slope):.4f} >= FE_REW_SLOPE_GATE {FE_REW_SLOPE_GATE}")


# ── A4 — Fe I RAW scatter vs the RYA-407 honest floor (wired RYA-446) ─────────

def test_A4_raw_scatter_within_honest_floor():
    floor = float(ACCEPTANCE_PROFILES['G']['fe1_scatter_max'])
    assert floor == 0.1398, "RYA-446 wiring drifted: G fe1_scatter_max != 0.1398 (RYA-407)"
    scatter = float(_fe_verdict()['sigma'])                 # committed live pool scatter
    assert scatter <= floor, f"Fe I raw scatter {scatter} > honest floor {floor}"


# ── A5 — Fe I sigma_REPORTED (SE + Type B) <= 0.05, distinct from raw scatter ─

def test_A5_sigma_reported_within_gold_gate():
    fe = _uncertainty_fe()
    sr = float(fe['sigma_reported'])
    # guard the conflation bug (RYA-166 rewrite): reported != raw scatter, and is
    # the SE+systematics combination, not a stub/placeholder.
    assert sr < float(fe['raw_sigma']), (
        f"sigma_reported {sr} not smaller than raw scatter {fe['raw_sigma']} "
        "— reported uncertainty must be the SE+systematics quantity, not raw scatter")
    assert np.isfinite(fe['sigma_SE']) and fe['sigma_SE'] > 0, "sigma_SE missing/stub"
    assert abs(sr - float(fe['raw_sigma'])) > 1e-6, "sigma_reported == raw scatter (conflation)"
    assert sr <= SIGMA_REPORTED_GATE, (
        f"Fe I sigma_reported {sr} > gold gate {SIGMA_REPORTED_GATE} dex")


# ── A6 — RYA-400 physics-regime audit, wired in-process ──────────────────────

def test_A6_physics_regime_audit_zero_failures():
    import audit_physics_regime_rya400 as A
    with contextlib.redirect_stdout(io.StringIO()):
        rc = A.run()
    assert rc == 0, "RYA-400 physics-regime audit reports failures (run() != 0)"


# ── A7 — RYA-429 no-silent-drop invariant on the full solar linelist ─────────

def test_A7_no_silent_drop_invariant_holds():
    if not LEDGER.exists():
        pytest.skip(f"rejection ledger absent: {LEDGER.name}")
    a = json.loads(LEDGER.read_text())
    assert a['invariant_holds'] is True
    assert a['n_unaccounted'] == 0 and a['n_unclassified'] == 0
    assert a['n_measured'] + a['n_rejected'] == a['n_in_coverage']


# ── A8 — every target element honestly classified, none silently absent ──────

def test_A8_every_target_element_has_honest_verdict():
    by_el = {v['element']: v for v in _verdict()['verdicts']}
    missing = [el for el in TARGET_ELEMENTS if el not in by_el]
    assert not missing, f"target elements silently absent from the verdict: {missing}"
    bad = {el: by_el[el].get('verdict') for el in TARGET_ELEMENTS
           if by_el[el].get('verdict') not in VALID_VERDICTS}
    assert not bad, f"elements with no honest verdict class: {bad}"


# ── the single aggregate: ONE assertion = Solar Beta-ready ───────────────────

def test_solar_is_beta_ready():
    """Aggregate PASS/FAIL over every gate assertion — the RYA-527 pass condition."""
    checks = []

    def _run(name, fn):
        try:
            fn()
            checks.append((name, True, ''))
        except pytest.skip.Exception as e:      # missing committed artifact
            checks.append((name, None, f'SKIP: {e}'))
        except AssertionError as e:
            checks.append((name, False, str(e).splitlines()[0]))

    _run('A1 3D-corrected A(Fe) in FE_GATE', test_A1_a_fe_in_fe_gate)
    _run('A2 ionization balance', test_A2_ionization_balance_within_gate)
    _run('A3 REW slope', test_A3_rew_slope_flat_within_tolerance)
    _run('A4 raw scatter floor', test_A4_raw_scatter_within_honest_floor)
    _run('A5 sigma_reported <= 0.05', test_A5_sigma_reported_within_gold_gate)
    _run('A6 RYA-400 audit', test_A6_physics_regime_audit_zero_failures)
    _run('A7 RYA-429 invariant', test_A7_no_silent_drop_invariant_holds)
    _run('A8 verdict coverage', test_A8_every_target_element_has_honest_verdict)

    print("\n  SOLAR CALIBRATION GATE (RYA-166 Layer 3)")
    for name, ok, note in checks:
        tag = 'PASS' if ok else ('SKIP' if ok is None else 'FAIL')
        print(f"    [{tag}] {name}{('  -> ' + note) if note else ''}")
    failed = [n for n, ok, _ in checks if ok is False]
    assert not failed, f"Solar NOT Beta-ready — failing gate assertions: {failed}"


if __name__ == '__main__':
    test_solar_is_beta_ready()
    print("\n  Solar Beta-ready: ALL gate assertions PASS")
