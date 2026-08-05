"""
tests/test_co_validation_rya390.py
==================================
RYA-390 Part B — guard the three-way CO validation engine: velocity recovery on a
synthetic shifted spectrum, the CRIRES smoothing, and the real provisional-product
run (telluric-dominated headline + telluric-reference cross-check).
"""
import sys
from pathlib import Path

import numpy as np
from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz (RYA-313)
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline import co_validation_rya390 as v  # noqa: E402

COND = ROOT / 'data' / 'audit' / 'crires_co_conditioned'


def _synthetic(wave, centers, depth=0.5, sigma=0.1):
    f = np.ones_like(wave)
    for c in centers:
        f -= depth * np.exp(-0.5 * ((wave - c) / sigma) ** 2)
    return f


def test_measure_velocity_recovers_known_shift():
    # reference at rest; "observed" = reference shifted by a known velocity
    w = np.linspace(22900.0, 23000.0, 4000)
    centers = np.array([22930.0, 22955.0, 22978.0, 22995.0])
    v_true = 30.0
    f_ref = _synthetic(w, centers)                                   # rest
    f_obs = _synthetic(w, centers * (1.0 + v_true / v.C_KMS))        # Doppler-shifted lines
    out = v.measure_velocity(w, f_obs, w, f_ref, smooth=False)
    assert abs(out['v_kms'] - v_true) <= 1.0
    assert out['peak_xcorr'] > 0.9


def test_smooth_to_crires_broadens():
    w = np.linspace(22900.0, 23000.0, 8000)
    f = _synthetic(w, [22950.0], depth=0.8, sigma=0.02)   # very narrow line
    fs = v._smooth_to_crires(w, f)
    assert fs.min() > f.min()           # smoothing makes the core shallower
    assert _trapezoid(1 - fs, w) == pytest.approx(_trapezoid(1 - f, w), rel=0.05)  # ~flux conserved


@pytest.mark.skipif(not (COND / v.CONDITIONED['K2192']).exists(),
                    reason="conditioned CO product not present")
def test_run_telluric_dominated_and_refs_agree():
    rep = v.run()
    # the two INDEPENDENT telluric atlases must agree (validates the reference side)
    tc = rep['telluric_reference_crosscheck']
    assert tc['status'] == 'RAN' and tc['telluric_depth_corr'] > 0.8
    for s, r in rep['settings'].items():
        # provisional product is telluric-dominated: high xcorr with telluric at v≈0,
        # low solar recovery → the decisive Part B finding (regression guard)
        rt = r['diagnostic_residual_telluric_vs_photatl_atm']
        assert abs(rt['v_kms']) <= v.V_TELLURIC_KMS
        assert rt['peak_xcorr'] > 0.7
        assert 'TELLURIC-DOMINATED' in r['verdict']
        # telluric-model check 2 now RUNS — RYA-380 persists molecfit mtrans (was
        # BLOCKED). It compares mtrans vs the independent telluric atlases; on the
        # on-chip bandhead order the overlap is photatl-atmospheric (Wallace segment
        # misses it). The corrected product being telluric-dominated corroborates a
        # weak telluric MODEL: depth_corr is low, not a PASS.
        c2 = r['check2_telluric_model_vs_Wallace']
        assert c2['status'] == 'RAN' and c2['verdict'] != 'BLOCKED'
        assert 'depth_corr' in c2 and 'per_reference' in c2
