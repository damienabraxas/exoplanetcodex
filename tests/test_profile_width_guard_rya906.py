"""RYA-906/911 — the width guard tests the TOTAL profile width, not one component of a
degenerate pair.

The guard this replaces rejected a line when its fitted Gaussian sigma sat on the
instrumental bound, reasoning that "no observed feature can be narrower than the
instrument profile". Measured on the committed HARPS Fe II pool, it fired on 25 lines and
EVERY ONE was a Voigt fit; it never fired on a Gaussian. In a Voigt profile sigma and
gamma are degenerate, so a railed sigma records where the optimiser resolved that
degeneracy, not how wide the line is.
"""
from __future__ import annotations

import math

import pytest

from scripts.measure_band_profilefit import (
    _irreducible_sigma_kms, _physical_floor_fwhm, _total_width_below_physical_floor,
    voigt_fwhm)
from config.constants import STAR_PARAMS

C_KMS = 299792.458


def test_the_irreducible_width_is_derived_from_the_stars_own_parameters():
    """Not a literal. Thermal + microturbulent, from the ratified solar row."""
    sun = STAR_PARAMS["solar"]
    k_b, m_fe = 1.380649e-23, 55.845 * 1.66054e-27
    expect = math.sqrt(k_b * float(sun["teff"]) / m_fe
                       + (float(sun["xi"]) * 1000.0) ** 2 / 2.0) / 1000.0
    assert _irreducible_sigma_kms() == pytest.approx(expect, rel=1e-12)


def test_the_floor_EXCLUDES_macroturbulence_so_it_stays_a_lower_bound():
    """vmac is 3.8 km/s — several times the thermal width. Including it would turn a
    floor into an expectation and reject perfectly good lines."""
    assert _irreducible_sigma_kms() < float(STAR_PARAMS["solar"]["vmac"])
    # and it must sit below the fit's seed (instrumental (+) 1.7 km/s stellar), which
    # DOES include macroturbulence
    assert _irreducible_sigma_kms() < 1.7


def test_voigt_fwhm_reduces_to_the_gaussian_case():
    sigma = 0.03
    assert voigt_fwhm(sigma, None) == pytest.approx(2 * sigma * math.sqrt(2 * math.log(2)))
    assert voigt_fwhm(sigma, 0.0) == pytest.approx(2 * sigma * math.sqrt(2 * math.log(2)))


def test_voigt_fwhm_grows_with_the_lorentzian_component():
    """The point of the whole fix: gamma carries width that sigma does not see."""
    g0 = voigt_fwhm(0.03, None)
    assert voigt_fwhm(0.03, 0.02) > g0
    assert voigt_fwhm(0.03, 0.06) > voigt_fwhm(0.03, 0.02)


# ── the correction itself, as a positive control ─────────────────────────────

def test_a_railed_sigma_with_a_broad_lorentzian_is_NOT_rejected():
    """🔴 THE CORRECTION. Fe II 6084.102 on HARPS: sigma 0.0225 sat exactly on the
    instrumental bound and the OLD guard rejected it — while gamma was 0.0591 and the
    total FWHM 0.1395 A, squarely among the lines that were kept. The width was in the
    Lorentzian. It must survive."""
    lam, sigma_inst = 6084.102, 0.0225
    popt = [lam, 0.5, 0.0225, 0.0591]
    assert not _total_width_below_physical_floor(popt, "voigt", lam, sigma_inst)
    # the old test would have rejected it: sigma IS on the bound
    assert abs(popt[2] - sigma_inst) < 1e-4


def test_a_genuinely_narrow_line_IS_still_rejected():
    """The guard must not become a no-op. A profile with no Lorentzian and a sigma at the
    instrumental bound really is narrower than thermal+microturbulent broadening allows."""
    lam, sigma_inst = 6084.102, 0.0225
    assert _total_width_below_physical_floor([lam, 0.5, sigma_inst, 0.0], "gaussian",
                                             lam, sigma_inst)


def test_the_floor_scales_with_wavelength():
    """Both the instrumental and the Doppler terms are proportional to lambda, so a fixed
    Angstrom threshold would be wrong at one end of any band."""
    lo = _physical_floor_fwhm(4000.0, 4000.0 / 115000 / 2.35482)
    hi = _physical_floor_fwhm(6900.0, 6900.0 / 115000 / 2.35482)
    assert hi > lo
    assert hi / lo == pytest.approx(6900.0 / 4000.0, rel=1e-6)
