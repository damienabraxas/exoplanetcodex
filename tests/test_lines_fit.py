# tests/test_lines_fit.py
import pytest
import numpy as np
from pipeline.lines_fit import (
    _voigt_abs, _gauss_abs, _integrate_profile,
    _fit_profile, _local_renorm, _ew_error, _ew_scores
)


# ── _ew_scores unit tests (RYA-166: the numeric line subscores had only bounds
#    checks via the loader; here we pin their behaviour directly) ───────────────

def test_ew_scores_bounded_unit_interval():
    for ew, err, chi2 in [(50.0, 2.0, 0.5), (5.0, 5.0, 50.0), (250.0, 1.0, 0.01)]:
        snr, chi2s, sat = _ew_scores(ew, err, chi2)
        for s in (snr, chi2s, sat):
            assert 0.0 <= s <= 1.0

def test_ew_scores_snr_monotonic_in_precision():
    """Higher EW/err (better SNR) -> higher snr subscore, saturating at 1."""
    lo, _, _ = _ew_scores(30.0, 6.0, 0.5)     # SNR 5
    hi, _, _ = _ew_scores(120.0, 2.0, 0.5)    # SNR 60 -> saturates
    assert hi >= lo and hi == pytest.approx(1.0, abs=1e-6)

def test_ew_scores_chi2_penalises_bad_fits():
    """chi2 subscore = 1 for a clean fit, decreases as chi2 rises, floors at 0."""
    clean, mid, bad = (_ew_scores(50.0, 2.0, c)[1] for c in (0.5, 5.0, 100.0))
    assert clean == pytest.approx(1.0, abs=1e-6)
    assert clean > mid > 0.0
    assert bad == pytest.approx(0.0, abs=1e-6)

def test_ew_scores_saturation_rises_with_ew():
    """saturation subscore ramps from 0 (weak) to 1 (strong) across the band."""
    weak, mid, strong = (_ew_scores(e, 2.0, 0.5)[2] for e in (80.0, 150.0, 250.0))
    assert weak == pytest.approx(0.0, abs=1e-6)
    assert 0.0 < mid < 1.0
    assert strong == pytest.approx(1.0, abs=1e-6)

def test_voigt_abs_peak():
    """Voigt profile should reach minimum at line centre."""
    x = np.linspace(6299, 6301, 1000)
    y = _voigt_abs(x, 6300.0, 0.5, 0.025, 0.010)
    assert y.min() == pytest.approx(0.5, abs=0.01)
    assert y[500] < y[0]  # core deeper than wing

def test_gauss_abs_peak():
    x = np.linspace(6299, 6301, 1000)
    y = _gauss_abs(x, 6300.0, 0.5, 0.025)
    assert y.min() == pytest.approx(0.5, abs=0.01)

def test_integrate_profile_synthetic():
    """EW of a known synthetic Gaussian should be recoverable."""
    x = np.linspace(5575, 5577, 500)
    # Gaussian: depth=0.3, sigma=0.05 → EW ≈ depth * sigma * sqrt(2π) * 1000 mÅ
    popt = [5576.0, 0.3, 0.05]
    ew = _integrate_profile(x, popt, 'gaussian')
    expected = 0.3 * 0.05 * np.sqrt(2 * np.pi) * 1000  # ~75 mÅ
    assert abs(ew - expected) < 5.0, f"EW={ew:.1f}, expected~{expected:.1f}"

def test_fit_profile_recovers_synthetic_gaussian():
    """Fitter should recover known parameters on clean synthetic data."""
    x = np.linspace(5575, 5577, 200)
    true_popt = [5576.0, 0.4, 0.030]
    y = _gauss_abs(x, *true_popt) + np.random.normal(0, 0.001, len(x))
    popt, pcov, profile_type, chi2 = _fit_profile(x, y, 5576.0)
    assert popt is not None, "Fit should succeed on clean synthetic data"
    assert abs(popt[0] - 5576.0) < 0.02, "Line centre should be recovered"
    assert abs(popt[1] - 0.4) < 0.05, "Depth should be recovered"

def test_local_renorm_flat_continuum():
    """Renorm of flat continuum should return unity."""
    x = np.linspace(5574, 5578, 400)
    y = np.ones_like(x)
    y_norm, cont = _local_renorm(x, y, 5576.0)
    assert np.allclose(y_norm, 1.0, atol=0.005)
