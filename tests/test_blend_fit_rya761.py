"""RYA-761 — the N-component blend fitter.

The first test is the one that matters for the production path: at n = 1 the general
fitter must reproduce the single-component fitter BIT-FOR-BIT. If it does not, every
banked EW that was ever produced by `_fit_profile` is at risk, and the refactor that
removes the duplicate optimiser call is not safe to make.
"""
from __future__ import annotations

import numpy as np
import pytest

from pipeline.lines_fit import _fit_blend, _fit_profile, _gauss_abs, _voigt_abs

RNG = np.random.default_rng(20260811)


def _spectrum(centres, depths, sigma=0.030, gamma=0.008, noise=0.0,
              half_width=0.60, n_pix=300, seed=None):
    """A synthetic window: sum of Voigt absorbers on a unit continuum."""
    wav = np.linspace(centres[0] - half_width, centres[0] + half_width, n_pix)
    flux = np.ones_like(wav)
    for c, d in zip(centres, depths):
        flux = flux - (1.0 - _voigt_abs(wav, c, d, sigma, gamma))
    if noise:
        rng = np.random.default_rng(seed) if seed is not None else RNG
        flux = flux + rng.normal(0.0, noise, size=wav.size)
    return wav, flux


# ── the reduction guarantee ──────────────────────────────────────────────────
#
# `_fit_profile` now DELEGATES to `_fit_blend`, so comparing the two against each
# other would be tautological. These are the values the single-component fitter
# produced BEFORE the delegation, captured from the pre-refactor code and pinned here.
# They are the real regression barrier: every banked EW in the project came out of
# that code path, so if the delegation moved a digit, these fail.
REFERENCE = {
    'd030_s030': ('voigt', 5.795028591905097e-18,
                  [6800.0, 0.30000000286659784, 0.02999999553651971,
                   0.008000007615795654]),
    'd005_s020': ('voigt', 7.871412585223308e-24,
                  [6800.0, 0.04999999999659146, 0.020000000003935613,
                   0.008000000009708824]),
    'd070_s045': ('voigt', 9.598784960327637e-22,
                  [6800.0, 0.700000000035413, 0.04499999996947458,
                   0.008000000051867389]),
    'kp_s008':   ('voigt', 1.1597040740970366e-17,
                  [8000.000000000003, 0.20000001223266795, 0.007999992061359925,
                   0.008000009036982298]),
}


@pytest.mark.parametrize("key,depth,sigma", [
    ('d030_s030', 0.30, 0.030),
    ('d005_s020', 0.05, 0.020),
    ('d070_s045', 0.70, 0.045),
])
def test_single_component_matches_pre_refactor_values_exactly(key, depth, sigma):
    """The production single-component path must be unchanged to the last bit."""
    wav, flux = _spectrum([6800.0], [depth], sigma=sigma)
    kind_ref, chi_ref, p_ref = REFERENCE[key]

    p, _, kind, chi = _fit_profile(wav, flux, 6800.0)

    assert kind == kind_ref
    assert chi == chi_ref, "chi2 moved — the delegation changed the optimiser call"
    np.testing.assert_array_equal(p, np.array(p_ref))


def test_kitt_peak_resolution_path_matches_pre_refactor_values():
    """The RYA-713 instrument-resolution arguments must survive the delegation."""
    wav, flux = _spectrum([8000.0], [0.20], sigma=0.008)
    kind_ref, chi_ref, p_ref = REFERENCE['kp_s008']

    p, _, kind, chi = _fit_profile(
        wav, flux, 8000.0, sigma_init=0.008, sigma_min=0.005, sigma_max=0.40)

    assert kind == kind_ref and chi == chi_ref
    np.testing.assert_array_equal(p, np.array(p_ref))


def test_fit_profile_and_fit_blend_agree_at_n1():
    """Belt and braces: the delegation is real, so these are the same call."""
    wav, flux = _spectrum([7200.0], [0.25], noise=0.002, seed=7)
    p_a, _, k_a, chi_a = _fit_profile(wav, flux, 7200.0)
    p_b, _, k_b, chi_b = _fit_blend(wav, flux, [7200.0])
    assert k_a == k_b and chi_a == chi_b
    np.testing.assert_array_equal(p_a, p_b)


# ── what the second component actually buys ──────────────────────────────────

def test_two_component_recovers_a_blend_the_single_fit_cannot():
    """The RYA-761 case: a real neighbour 0.25 A away.

    The single-component fit must absorb the neighbour into an inflated depth/width;
    the two-component fit must recover the target's true depth.
    """
    true_d, sep = 0.30, 0.25
    wav, flux = _spectrum([7500.0, 7500.0 + sep], [true_d, 0.22], sigma=0.030)

    p1, _, _, _ = _fit_profile(wav, flux, 7500.0)
    p2, _, _, _ = _fit_blend(wav, flux, [7500.0, 7500.0 + sep])

    err1 = abs(p1[1] - true_d)
    err2 = abs(p2[1] - true_d)
    assert err2 < err1, f"two-component no better: {err2:.4f} vs {err1:.4f}"
    assert err2 < 0.02, f"target depth not recovered: {p2[1]:.4f} vs {true_d}"


def test_interloper_centre_is_not_free():
    """RYA-761 forbids freeing the neighbour's wavelength.

    Fitting with a deliberately WRONG interloper centre must not silently slide that
    component onto the real line: the parameter simply does not exist in the vector.
    Shared width => 4 target params + (n-1) depths.
    """
    wav, flux = _spectrum([7500.0, 7500.30], [0.30, 0.20])
    popt, _, kind, _ = _fit_blend(wav, flux, [7500.0, 7500.30])
    expected = (4 if kind == 'voigt' else 3) + 1
    assert len(popt) == expected, (
        f"parameter vector length {len(popt)} != {expected}; a free interloper "
        f"centre would add a parameter")


def test_three_components_fit():
    """The generalisation is not special-cased to two."""
    wav, flux = _spectrum([6900.0, 6900.28, 6899.70], [0.28, 0.18, 0.15])
    popt, _, kind, chi2 = _fit_blend(wav, flux, [6900.0, 6900.28, 6899.70])
    assert popt is not None and np.isfinite(chi2)
    assert len(popt) == (4 if kind == 'voigt' else 3) + 2
    assert abs(popt[1] - 0.28) < 0.03


def test_blend_fit_does_not_invent_a_neighbour_that_is_absent():
    """A clean isolated line fitted as a 2-component blend must give the extra
    component ~zero depth. If it does not, the model is manufacturing structure and
    every 'recovery' it reports would be suspect."""
    wav, flux = _spectrum([7100.0], [0.30])
    popt, _, kind, _ = _fit_blend(wav, flux, [7100.0, 7100.30])
    interloper_depth = popt[-1]
    assert interloper_depth < 0.02, (
        f"invented a neighbour of depth {interloper_depth:.4f} where none exists")
