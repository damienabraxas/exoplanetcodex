"""The width guard at the SECOND profile-fit site — RYA-911 follow-on.

RYA-906/911 (PR #315) replaced the sigma-on-the-floor test with a total-Voigt-FWHM test
in `scripts/measure_band_profilefit.py`. It did NOT reach `pipeline/measure/profile_fit.py`
— the registered `ProfileFitHandler` that `resolve_handler()` returns for VIS and
red-optical — which kept the refuted test.

⚠️ THESE TESTS DELIBERATELY IMPORT NOTHING FROM `scripts/measure_band_profilefit.py`.
That module's import chain loads the Kitt Peak atlas, so it is unimportable off Sirius —
which is *why* the corrected physics was never shared, and why a test that could only run
in one place would be reproducing the original mistake.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from pipeline.band_policy import resolve
from pipeline.band_products import REW_SATURATION_CEILING
from pipeline.line_width import (irreducible_sigma_kms, physical_floor_fwhm, voigt_fwhm,
                                 total_width_below_physical_floor)
from pipeline.lines_fit import _gauss_abs, _voigt_abs
from pipeline.measure.profile_fit import ProfileFitHandler
from config.constants import STAR_PARAMS

C = 5500.0
R = 115000.0            # HARPS


def _measure(flux):
    h = ProfileFitHandler()
    wav = np.arange(C - 2.0, C + 2.0, 0.01)
    return h.measure_line(wav, flux, element="Fe", ion="II", wavelength_A=C,
                          instrument="harps", policy=resolve(C), pre_normalised=True,
                          context={"resolving_power": R})


def _floor():
    return ProfileFitHandler().widths(R, C)[1]


def test_the_handler_is_what_vis_actually_resolves_to():
    """If this stops being true the guard below is guarding nothing."""
    from pipeline.measure import resolve_handler
    assert isinstance(resolve_handler(C), ProfileFitHandler)
    assert "profile-fit" in resolve(C).permitted_methods


def test_a_railed_sigma_with_a_healthy_gamma_is_KEPT():
    """🔴 THE REGRESSION. This line is 2.9x WIDER than the instrument profile and the old
    guard threw it away for being too narrow.

    The width is in the Lorentzian: sigma sits exactly on its bound while gamma carries
    the line. Judging the fit by sigma alone measures where the degeneracy landed.
    """
    s_min = _floor()
    gamma = 0.060
    lm = _measure(_voigt_abs(np.arange(C - 2.0, C + 2.0, 0.01), C, 0.10, s_min, gamma))

    # The premise: this really is a comfortably wide line.
    total = voigt_fwhm(s_min, gamma)
    assert total > 2.5 * (C / R)
    assert total > physical_floor_fwhm(C, s_min)
    # ...and it must not be excluded for some OTHER reason, or the test proves nothing.
    assert lm.rew < REW_SATURATION_CEILING, "test line drifted into saturation"

    assert lm.in_aggregate, f"kept-line regression: {lm.excluded_reason}"
    assert not lm.excluded_reason


def test_a_genuinely_sub_physical_line_is_still_QUARANTINED():
    """The guard must not have become a no-op — a green guard proves nothing until it can
    fail. A pure Gaussian narrower than the instrument has no Lorentzian to hide in."""
    s_min = _floor()
    lm = _measure(_gauss_abs(np.arange(C - 2.0, C + 2.0, 0.01), C, 0.10, s_min * 0.6))
    assert not lm.in_aggregate
    assert lm.excluded_reason.startswith("UNDER-PHYSICAL-WIDTH")


def test_the_retired_verdict_is_not_emitted_here_any_more():
    """FIT-PINNED was refuted, not renamed. It must not reappear as a reason string."""
    import inspect

    import pipeline.measure.profile_fit as pf
    src = inspect.getsource(pf)
    emitted = [ln for ln in src.splitlines()
               if "FIT-PINNED" in ln and "#" not in ln.split("FIT-PINNED")[0]]
    assert not emitted, f"FIT-PINNED still emitted as a verdict: {emitted}"


def test_the_fit_numbers_are_recorded_not_just_used():
    """RYA-911: a quantity that decides a verdict and is then dropped leaves a sentence a
    human can read and nothing a program can. This site dropped all four."""
    s_min = _floor()
    lm = _measure(_voigt_abs(np.arange(C - 2.0, C + 2.0, 0.01), C, 0.10, s_min, 0.060))
    assert lm.profile_sigma_A == pytest.approx(s_min, abs=1e-4)
    assert lm.profile_sigma_floor_A == pytest.approx(s_min)
    assert lm.profile_gamma_A is not None and lm.profile_gamma_A > 0
    assert lm.red_chi2 is not None


def test_profile_sigma_is_an_angstrom_and_is_not_sigma_A():
    """`sigma_A` is one sigma on A in DEX (RYA-847). `profile_sigma_A` is a width in
    ANGSTROM. One character apart, no shared units — RYA-911 nearly conflated them."""
    s_min = _floor()
    lm = _measure(_voigt_abs(np.arange(C - 2.0, C + 2.0, 0.01), C, 0.10, s_min, 0.060))
    assert lm.profile_sigma_A < 1.0          # an Angstrom width at 5500 A
    assert lm.sigma_A is None                # the profile fit constrains no chi2 curvature


def test_the_floor_derives_from_config_and_excludes_macroturbulence():
    """A true lower bound, not an expectation — and never a literal in the module."""
    sun = STAR_PARAMS["solar"]
    k_b, m_fe = 1.380649e-23, 55.845 * 1.66054e-27
    expect = math.sqrt(k_b * float(sun["teff"]) / m_fe
                       + (float(sun["xi"]) * 1000.0) ** 2 / 2.0) / 1000.0
    assert irreducible_sigma_kms() == pytest.approx(expect, rel=1e-12)
    assert irreducible_sigma_kms() < float(sun["vmac"])


def test_the_floor_scales_with_wavelength():
    """sigma_inst is lambda/R, so the floor is a per-line quantity, never a constant."""
    assert physical_floor_fwhm(8000.0, 8000.0 / R * (1 / 2.35482)) > \
           physical_floor_fwhm(4000.0, 4000.0 / R * (1 / 2.35482))


def test_a_gaussian_fit_is_judged_on_its_gaussian_width_alone():
    """`gamma_of` returns None for a Gaussian — not 0.0, which would read as measured."""
    s_min = _floor()
    popt_g = (C, 0.1, s_min * 0.5)
    assert total_width_below_physical_floor(popt_g, "gaussian", C, s_min)
    popt_wide = (C, 0.1, s_min * 4.0)
    assert not total_width_below_physical_floor(popt_wide, "gaussian", C, s_min)


def test_both_fitters_share_one_implementation():
    """The whole point. Two routes may not carry two copies of one physical verdict."""
    import pipeline.line_width as lw
    import inspect
    src = inspect.getsource(inspect.getmodule(ProfileFitHandler))
    assert "from pipeline.line_width import" in src
    # the physics itself appears once, here
    assert "0.5346" in inspect.getsource(lw.voigt_fwhm)
