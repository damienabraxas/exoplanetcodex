"""
tests/test_uncertainty_stack_rya158.py
======================================
RYA-158 (folded into RYA-166 Step 2.5) — the reported-uncertainty budget.

Guards the DISCIPLINE + the committed budget, NOT a re-run (Type B re-derives
abundances at perturbed params, which needs the EW->abundance pipeline; the
committed audit JSON is the stable artifact — house style).

The cardinal invariant: sigma_reported is the standard-error-of-the-mean combined
with stellar-parameter systematics, a DIFFERENT and SMALLER number than the raw
line-to-line scatter (RYA-407 floor). Conflating them was the month-long bug the
RYA-166 rewrite exists to prevent.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import ACCEPTANCE_PROFILES        # noqa: E402
import pipeline.uncertainty_stack as U                  # noqa: E402

BUDGET = ROOT / 'data' / 'audit' / 'uncertainty' / 'solar_uncertainty_rya158.json'


def _budget():
    if not BUDGET.exists():
        pytest.skip(f"budget absent ({BUDGET.name}); run "
                    "`python -m pipeline.uncertainty_stack --star solar`")
    return json.loads(BUDGET.read_text())


def _fe():
    return next(r for r in _budget()['per_element'] if r['element'] == 'Fe')


# ── the module is implemented (no longer a stub) ─────────────────────────────

def test_module_is_implemented_not_a_stub():
    assert hasattr(U, 'run') and callable(U.run)
    # the stub raised NotImplementedError unconditionally; the real run() has the
    # Type B derivative-step + solar delta_p machinery.
    assert isinstance(U.TYPE_B_STEPS, dict) and 'teff_K' in U.TYPE_B_STEPS


# ── raw_sigma is CITED from the single source, not recomputed ─────────────────

def test_fe_raw_sigma_is_the_ratified_floor():
    raw = U._raw_sigma_and_N('solar')
    fe_raw_sigma, n, _ = raw['Fe']
    assert fe_raw_sigma == float(ACCEPTANCE_PROFILES['G']['fe1_scatter_max'])  # 0.1398, RYA-407
    assert n == 62


# ── solar delta_p: the Sun is the zero-point (logg exact, [Fe/H]=0) ──────────

def test_solar_delta_p_are_the_suns_own_uncertainties():
    _, deltas = U._solar_params_and_deltas()
    assert deltas['teff_K'] == 1.0        # solar Teff ~1 K
    assert deltas['logg'] == 0.0          # solar logg effectively exact -> no logg term
    assert deltas['feh'] == 0.0           # Sun DEFINES [Fe/H]=0 -> no FeH term
    assert deltas['vturb_kms'] > 0        # vmic carries the only non-trivial Type B term


# ── the committed budget: arithmetic + the raw-vs-reported distinction ───────

def test_fe_sigma_SE_is_scatter_over_sqrt_N():
    fe = _fe()
    expect = float(fe['raw_sigma']) / np.sqrt(int(fe['n_lines']))
    assert abs(float(fe['sigma_SE']) - expect) < 1e-3


def test_fe_sigma_reported_is_quadrature_sum():
    fe = _fe()
    terms = [fe['sigma_SE'], fe['sigma_B_Teff'], fe['sigma_B_logg'],
             fe['sigma_B_vmic'], fe['sigma_B_FeH']]
    expect = float(np.sqrt(sum(float(t) ** 2 for t in terms)))
    assert abs(float(fe['sigma_reported']) - expect) < 1e-3


def test_fe_reported_is_NOT_raw_scatter_the_conflation_guard():
    fe = _fe()
    assert float(fe['sigma_reported']) < float(fe['raw_sigma'])
    assert abs(float(fe['sigma_reported']) - float(fe['raw_sigma'])) > 0.05


def test_fe_solar_terms_with_zero_delta_are_zero():
    fe = _fe()
    assert float(fe['sigma_B_logg']) == 0.0     # delta_logg = 0
    assert float(fe['sigma_B_FeH']) == 0.0      # delta_FeH = 0


def test_fe_sigma_reported_in_expected_ballpark():
    # RYA-166 rewrite predicted ~0.02-0.03 once genuinely computed; landing far off
    # is the STOP-and-report condition (not a reason to move the 0.05 target).
    sr = float(_fe()['sigma_reported'])
    assert 0.01 <= sr <= 0.05, f"Fe I sigma_reported {sr} outside the predicted ballpark"
