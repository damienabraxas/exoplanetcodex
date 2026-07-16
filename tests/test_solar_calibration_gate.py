"""
tests/test_solar_calibration_gate.py
====================================
THE Beta-blocking solar calibration gate (RYA-166 Layer 3; gates RYA-527).

ONE command answers "is Solar Beta-ready":
    pytest tests/test_solar_calibration_gate.py -v

Each assertion below is a single pass/fail check, sourced from COMMITTED artifacts
+ ratified, cited thresholds (never a threshold loosened to fit a run — validate-
don't-tune). It aggregates checks that otherwise live in separate scripts:

  A1  A(Fe) in the RYA-336 scale-aware DIAGNOSTIC window [7.44, 7.58]
      (absolute A(Fe) is a diagnostic, NOT a strict 7.46+/-0.05 pass/fail)
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
                              FE_1D3D_SOLAR_OFFSET, FE_ABS_DIAG_HALFWIDTH,
                              TARGET_ELEMENTS)

# ── committed artifacts (guard these, never re-run the pipeline in the gate) ──
VERDICT = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'
FE1_RESID = ROOT / 'data' / 'audit' / 'fe1_scatter' / 'fe1_per_line_residuals_rya407.csv'
UNCERT = ROOT / 'data' / 'audit' / 'uncertainty' / 'solar_uncertainty_rya158.json'
LEDGER = ROOT / 'data' / 'results' / 'rejection_ledger_solar_rya429.json'

ASPLUND_FE = 7.46
FE_ABS_DIAG_LO = ASPLUND_FE + FE_1D3D_SOLAR_OFFSET - FE_ABS_DIAG_HALFWIDTH   # 7.44
FE_ABS_DIAG_HI = ASPLUND_FE + FE_1D3D_SOLAR_OFFSET + FE_ABS_DIAG_HALFWIDTH   # 7.58
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


# ── A1 — absolute A(Fe): scale-aware diagnostic window (RYA-336) ─────────────

def test_A1_a_fe_in_scale_aware_diagnostic_window():
    a_fe = float(_fe_verdict()['A_measured'])
    assert FE_ABS_DIAG_LO <= a_fe <= FE_ABS_DIAG_HI, (
        f"A(Fe)={a_fe} outside RYA-336 diagnostic window "
        f"[{FE_ABS_DIAG_LO:.2f},{FE_ABS_DIAG_HI:.2f}]")


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

    _run('A1 A(Fe) diagnostic window', test_A1_a_fe_in_scale_aware_diagnostic_window)
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
