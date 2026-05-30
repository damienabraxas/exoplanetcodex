"""
tests/test_acquire.py
=====================
Tests for pipeline/spectra_acquire.py.
"""

import importlib.util
import sys
from pathlib import Path

import pytest
import numpy as np

# pipeline/spectra_acquire.py is a valid dotted import now but we keep
# the explicit load for consistency with how the test suite was set up.
_spec = importlib.util.spec_from_file_location(
    "spectra_acquire",
    Path(__file__).parent.parent / "pipeline" / "spectra_acquire.py",
)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["spectra_acquire"] = _mod
_spec.loader.exec_module(_mod)

correct_radial_velocity      = _mod.correct_radial_velocity
estimate_snr                 = _mod.estimate_snr
make_synthetic_demo_spectrum = _mod.make_synthetic_demo_spectrum
spectrum_summary             = _mod.spectrum_summary


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def synthetic_spectrum():
    wav  = np.linspace(3780, 6910, 5000)
    flux = np.ones_like(wav) * 10_000 + np.random.default_rng(0).normal(0, 100, 5000)
    return wav, np.clip(flux, 0, None)


# ── correct_radial_velocity ───────────────────────────────────────────────────

class TestCorrectRadialVelocity:
    def test_returns_array_same_shape(self, synthetic_spectrum):
        wav, _ = synthetic_spectrum
        assert correct_radial_velocity(wav, rv_kms=27.58).shape == wav.shape

    def test_positive_rv_shifts_blueward(self, synthetic_spectrum):
        wav, _ = synthetic_spectrum
        assert np.all(correct_radial_velocity(wav, rv_kms=27.58) < wav)

    def test_zero_rv_no_change(self, synthetic_spectrum):
        wav, _ = synthetic_spectrum
        np.testing.assert_allclose(correct_radial_velocity(wav, rv_kms=0.0), wav, rtol=1e-10)

    def test_doppler_formula(self):
        from config.constants import PHYSICS
        wav = np.array([5000.0])
        rv  = 29979.2458
        expected = 5000.0 / (1.0 + rv / PHYSICS['c_kms'])
        np.testing.assert_allclose(correct_radial_velocity(wav, rv), [expected], rtol=1e-9)


# ── estimate_snr ──────────────────────────────────────────────────────────────

class TestEstimateSnr:
    def test_returns_float(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        assert isinstance(estimate_snr(wav, flux), float)

    def test_snr_positive(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        assert estimate_snr(wav, flux) > 0

    def test_higher_snr_for_lower_noise(self):
        wav = np.linspace(5900, 6100, 1000)
        snr_clean = estimate_snr(wav, np.ones(1000)*10_000 + np.random.default_rng(2).normal(0, 10, 1000))
        snr_noisy = estimate_snr(wav, np.ones(1000)*10_000 + np.random.default_rng(1).normal(0, 500, 1000))
        assert snr_clean > snr_noisy

    def test_zero_noise_returns_zero(self):
        wav = np.linspace(5900, 6100, 200)
        assert estimate_snr(wav, np.ones(200) * 5000.0) == 0.0


# ── make_synthetic_demo_spectrum ──────────────────────────────────────────────

class TestMakeSyntheticDemoSpectrum:
    def test_returns_two_arrays(self):
        wav, flux = make_synthetic_demo_spectrum()
        assert isinstance(wav, np.ndarray) and isinstance(flux, np.ndarray)

    def test_same_length(self):
        wav, flux = make_synthetic_demo_spectrum()
        assert len(wav) == len(flux)

    def test_wavelength_range(self):
        from config.constants import PIPELINE
        wav, _ = make_synthetic_demo_spectrum()
        assert wav.min() >= PIPELINE['wav_min_A'] - 1
        assert wav.max() <= PIPELINE['wav_max_A'] + 1

    def test_flux_non_negative(self):
        _, flux = make_synthetic_demo_spectrum()
        assert np.all(flux >= 0)

    def test_reproducible(self):
        _, flux1 = make_synthetic_demo_spectrum()
        _, flux2 = make_synthetic_demo_spectrum()
        np.testing.assert_array_equal(flux1, flux2)


# ── spectrum_summary ──────────────────────────────────────────────────────────

class TestSpectrumSummary:
    def test_returns_dict(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        assert isinstance(spectrum_summary(wav, flux, {}), dict)

    def test_required_keys(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        result = spectrum_summary(wav, flux, {})
        for key in ('n_pixels', 'wav_min_A', 'wav_max_A', 'snr_estimate'):
            assert key in result

    def test_n_pixels(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        assert spectrum_summary(wav, flux, {})['n_pixels'] == len(wav)

    def test_wavelength_bounds(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        result = spectrum_summary(wav, flux, {})
        assert result['wav_min_A'] == pytest.approx(wav.min())
        assert result['wav_max_A'] == pytest.approx(wav.max())
