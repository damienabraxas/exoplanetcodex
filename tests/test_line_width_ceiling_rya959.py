"""The width CEILING — RYA-959.

RYA-906/911 gave `pipeline.line_width` a floor and no mirror. A fit that failed for any
reason therefore escaped UPWARD: `sigma_max` is 0.40 A, ~5x the widest Gaussian solar
Doppler physics permits, and the optimiser used all of it. On the first fresh RYA-959
Fe I VIS run 960 of 2342 HARPS fits sat pinned at exactly 0.400 with `red_chi2` ~ 0.02 —
they looked converged, and they integrated to 400-1000 mA "Fe I lines". That is how
RYA-958's contaminated pool was built.

These tests pin the RELATIONSHIPS, not the numbers (RYA-870): every assertion below stays
true if `STAR_PARAMS['solar']` is re-ratified, because each one is derived from config the
same way the code is. A test that hardcoded 0.0826 A would have to be edited the next time
vmac moves, and would then be pinning the edit rather than the physics.
"""
from __future__ import annotations

import math

import pytest

from config.constants import STAR_PARAMS
from pipeline import line_width as LW


HARPS_SIGMA_INST_5500 = 5500.0 / 115000.0 / 2.35482


def test_ceiling_is_above_the_floor_everywhere():
    """The two bounds must not cross — a crossing would quarantine every line at once."""
    for lam in (3800.0, 5000.0, 5500.0, 6910.0, 9200.0, 13000.0):
        for sig_inst in (0.005, 0.0203, 0.05):
            floor = LW.physical_floor_fwhm(lam, sig_inst)
            ceil = LW.voigt_fwhm(LW.physical_ceiling_sigma_A(lam, sig_inst), None)
            assert ceil > floor, f"ceiling {ceil} <= floor {floor} at {lam} A"


def test_ceiling_contains_every_doppler_term_the_floor_omits():
    """The floor is a floor BECAUSE it drops vmac and vsini; the ceiling must add them
    back. Asserted as a relationship so re-ratifying either parameter cannot silently
    turn the ceiling back into the floor."""
    assert LW.max_stellar_sigma_kms() > LW.irreducible_sigma_kms()
    sun = STAR_PARAMS["solar"]
    expected = math.hypot(LW.irreducible_sigma_kms(),
                          math.hypot(float(sun["vmac"]), float(sun["vsini"])))
    assert LW.max_stellar_sigma_kms() == pytest.approx(expected)


def test_ceiling_scales_with_wavelength_and_never_below_the_instrument():
    """Doppler broadening is a VELOCITY, so the width in Angstrom must scale with lambda —
    the RYA-837 mistake of holding an Angstrom value fixed across a wavelength lever."""
    a = LW.physical_ceiling_sigma_A(4000.0, 0.0)
    b = LW.physical_ceiling_sigma_A(8000.0, 0.0)
    assert b == pytest.approx(2.0 * a, rel=1e-9)
    # An instrument wider than every Doppler term still sets the ceiling.
    assert LW.physical_ceiling_sigma_A(5500.0, 5.0) > 5.0


def test_the_040_optimiser_bound_is_far_outside_physics():
    """The defect, stated as a test. If someone tightens `sigma_max` toward physics this
    fails and they must come here and say so — which is the point."""
    ceil = LW.physical_ceiling_sigma_A(5500.0, HARPS_SIGMA_INST_5500)
    assert ceil < 0.40 / 3.0, (
        "sigma_max=0.40 A is supposed to be far outside the physical ceiling; if it is "
        "no longer, the fitters' bound changed and this guard's premise changed with it")


def test_pinned_fit_is_convicted_and_the_reason_says_pinned():
    popt = [5500.0, 0.5, 0.40, 0.40]
    assert LW.gaussian_sigma_above_physical_ceiling(popt, 5500.0, HARPS_SIGMA_INST_5500)
    why = LW.over_physical_width_reason(popt, "voigt", 5500.0, HARPS_SIGMA_INST_5500, 0.40)
    assert why.startswith("OVER-PHYSICAL-WIDTH:")
    assert "PINNED" in why, "a fit ON the bound must be diagnosed as pinned, not merely wide"


def test_a_physical_solar_line_is_not_convicted_by_either_bound():
    """sigma 0.031 A at 5500 A is the observed solar Fe I width quoted by both fitters."""
    popt = [5500.0, 0.4, 0.031, 0.02]
    assert not LW.gaussian_sigma_above_physical_ceiling(
        popt, 5500.0, HARPS_SIGMA_INST_5500)
    assert not LW.total_width_below_physical_floor(
        popt, "voigt", 5500.0, HARPS_SIGMA_INST_5500)


def test_gamma_is_deliberately_unbounded():
    """Pressure broadening is LORENTZIAN and real. Bounding gamma here would quarantine
    the strongest lines in the band — which is precisely where RYA-945's laboratory gf
    backbone lives — so the sigma ceiling must ignore gamma entirely."""
    narrow = [5500.0, 0.4, 0.031, 0.001]
    damped = [5500.0, 0.4, 0.031, 0.900]
    assert not LW.gaussian_sigma_above_physical_ceiling(
        narrow, 5500.0, HARPS_SIGMA_INST_5500)
    assert not LW.gaussian_sigma_above_physical_ceiling(
        damped, 5500.0, HARPS_SIGMA_INST_5500)


# ── the implied-width backstop ───────────────────────────────────────────────

def test_equiv_rect_factor_is_derived_not_typed():
    """RYA-958's hand diagnosis wrote this as the bare literal 1.06."""
    assert LW.EQUIV_RECT_FACTOR == pytest.approx(1.06447, abs=1e-5)


def test_implied_width_reproduces_the_rya958_worst_case():
    """Fe I 5602.541 in the stale pool: 269 mA integrated onto a depth-0.008 core."""
    assert LW.implied_width_A(269.0, 0.008) == pytest.approx(31.6, abs=0.2)
    assert LW.implied_width_exceeds_ceiling(269.0, 0.008, 1.2)


def test_implied_ceiling_is_the_integration_window_not_a_multiplier():
    """🔴 THE FIRST VERSION OF THIS WAS TUNED AND IS REFUTED.

    It was `3.0 x the Doppler FWHM`. Measured on the RYA-959 run, the implied/Doppler
    ratio of sigma-clean fits is a smooth continuum (50/90/95/99th percentiles
    1.22/2.90/3.34/4.04) with no valley, so any multiple cuts a populated distribution.
    The bound is now a fact about the harness: `_integrate_profile` integrates over the
    fit window only, so an equivalent rectangle wider than that window claims more
    absorption than the interval can hold at the observed depth.
    """
    assert LW.implied_width_ceiling_A(1.2) == 1.2
    assert LW.implied_width_ceiling_A(2.4) == 2.4
    assert not hasattr(LW, "IMPLIED_WIDTH_ALLOWANCE"), (
        "the tuned multiplier was removed; re-introducing one needs a measured split "
        "to sit inside, and this population does not have one")


def test_implied_width_is_nan_safe_and_never_convicts_on_a_missing_depth():
    for bad in (0.0, -0.1, float("nan")):
        assert math.isnan(LW.implied_width_A(50.0, bad))
        assert not LW.implied_width_exceeds_ceiling(50.0, bad, 1.2)


def test_both_fitters_import_the_same_ceiling():
    """RYA-906/911's floor landed on ONE of the two profile fitters and the other kept a
    refuted test for a whole release. The module exists so that cannot recur; this
    asserts both sites still resolve to these functions and not to local copies."""
    import pipeline.measure.profile_fit as pf
    assert pf.gaussian_sigma_above_physical_ceiling is LW.gaussian_sigma_above_physical_ceiling
    assert pf.implied_width_exceeds_ceiling is LW.implied_width_exceeds_ceiling
