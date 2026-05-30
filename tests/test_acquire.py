"""
tests/test_acquire.py
=====================
Tests for pipeline/01_acquire.py.
"""

import importlib.util
import sys
from pathlib import Path

import pytest
import numpy as np

# pipeline/01_acquire.py can't be imported via normal dotted path because
# Python identifiers can't start with a digit. Load it explicitly.
_spec = importlib.util.spec_from_file_location(
    "acquire_01",
    Path(__file__).parent.parent / "pipeline" / "01_acquire.py",
)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["acquire_01"] = _mod
_spec.loader.exec_module(_mod)

correct_radial_velocity   = _mod.correct_radial_velocity
estimate_snr              = _mod.estimate_snr
make_synthetic_demo_spectrum = _mod.make_synthetic_demo_spectrum
spectrum_summary          = _mod.spectrum_summary


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def synthetic_spectrum():
    """A small synthetic spectrum for fast testing."""
    wav = np.linspace(3780, 6910, 5000)
    flux = np.ones_like(wav) * 10_000 + np.random.default_rng(0).normal(0, 100, 5000)
    return wav, np.clip(flux, 0, None)


# ── correct_radial_velocity ───────────────────────────────────────────────────

class TestCorrectRadialVelocity:
    def test_returns_array_same_shape(self, synthetic_spectrum):
        wav, _ = synthetic_spectrum
        result = correct_radial_velocity(wav, rv_kms=27.58)
        assert result.shape == wav.shape

    def test_positive_rv_shifts_wavelengths_blueward(self, synthetic_spectrum):
        # Star moving away → observed λ longer → rest λ shorter after correction
        wav, _ = synthetic_spectrum
        corrected = correct_radial_velocity(wav, rv_kms=27.58)
        assert np.all(corrected < wav)

    def test_zero_rv_no_change(self, synthetic_spectrum):
        wav, _ = synthetic_spectrum
        corrected = correct_radial_velocity(wav, rv_kms=0.0)
        np.testing.assert_allclose(corrected, wav, rtol=1e-10)

    def test_doppler_formula_correctness(self):
        from config.constants import PHYSICS
        wav = np.array([5000.0])
        rv = 29979.2458  # 10% of c
        corrected = correct_radial_velocity(wav, rv_kms=rv)
        expected = 5000.0 / (1.0 + rv / PHYSICS['c_kms'])
        np.testing.assert_allclose(corrected, [expected], rtol=1e-9)


# ── estimate_snr ──────────────────────────────────────────────────────────────

class TestEstimateSnr:
    def test_returns_float(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        snr = estimate_snr(wav, flux)
        assert isinstance(snr, float)

    def test_snr_positive(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        snr = estimate_snr(wav, flux)
        assert snr > 0

    def test_higher_snr_for_lower_noise(self):
        wav = np.linspace(5900, 6100, 1000)
        flux_clean = np.ones(1000) * 10_000
        flux_noisy = flux_clean + np.random.default_rng(1).normal(0, 500, 1000)
        snr_clean = estimate_snr(wav, flux_clean + np.random.default_rng(2).normal(0, 10, 1000))
        snr_noisy = estimate_snr(wav, flux_noisy)
        assert snr_clean > snr_noisy

    def test_zero_noise_returns_zero(self):
        wav = np.linspace(5900, 6100, 200)
        flux = np.ones(200) * 5000.0  # std = 0 → SNR = 0
        snr = estimate_snr(wav, flux)
        assert snr == 0.0


# ── make_synthetic_demo_spectrum ──────────────────────────────────────────────

class TestMakeSyntheticDemoSpectrum:
    def test_returns_two_arrays(self):
        wav, flux = make_synthetic_demo_spectrum()
        assert isinstance(wav, np.ndarray)
        assert isinstance(flux, np.ndarray)

    def test_arrays_same_length(self):
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
        wav1, flux1 = make_synthetic_demo_spectrum()
        wav2, flux2 = make_synthetic_demo_spectrum()
        np.testing.assert_array_equal(flux1, flux2)


# ── spectrum_summary ──────────────────────────────────────────────────────────

class TestSpectrumSummary:
    def test_returns_dict(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        result = spectrum_summary(wav, flux, {})
        assert isinstance(result, dict)

    def test_required_keys(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        result = spectrum_summary(wav, flux, {})
        for key in ('n_pixels', 'wav_min_A', 'wav_max_A', 'snr_estimate'):
            assert key in result

    def test_n_pixels_matches_input(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        result = spectrum_summary(wav, flux, {})
        assert result['n_pixels'] == len(wav)

    def test_wavelength_bounds(self, synthetic_spectrum):
        wav, flux = synthetic_spectrum
        result = spectrum_summary(wav, flux, {})
        assert result['wav_min_A'] == pytest.approx(wav.min())
        assert result['wav_max_A'] == pytest.approx(wav.max())
