"""
RYA-311 — the both-axis (FITEXY) straight-line fit, and why an OLS slope will not do.

Mucciarelli 2011 (A&A 528 A44) requires the A(Fe I)-vs-reduced-EW relation to be fitted
with uncertainties in BOTH coordinates, because sigma_EW enters the abscissa directly and
the ordinate through the curve of growth. The first test below MEASURES that claim on
synthetic data rather than restating it: with abscissa error present, OLS lands well below
the true slope and `fitexy` does not. A test that merely re-fitted clean data would pass
under either method and demonstrate nothing.
"""
from __future__ import annotations

import numpy as np
import pytest

from pipeline.fitexy import chi2_at_slope, fitexy


def _dataset(seed: int, n: int = 120, slope: float = 0.60, sx: float = 0.25,
             sy: float = 0.05):
    """A line observed with error on BOTH axes. `sx` is deliberately large enough that
    regression dilution is visible above the sampling noise."""
    rng = np.random.default_rng(seed)
    x_true = rng.uniform(-1.0, 1.0, n)
    y_true = 0.25 + slope * x_true
    x = x_true + rng.normal(0.0, sx, n)
    y = y_true + rng.normal(0.0, sy, n)
    return x, y, np.full(n, sx), np.full(n, sy), slope


def test_ols_is_biased_towards_zero_and_fitexy_is_not():
    """The Mucciarelli point, measured. Not a restatement of it."""
    ols_err, exy_err = [], []
    for seed in range(25):
        x, y, sx, sy, truth = _dataset(seed)
        ols = float(np.polyfit(x, y, 1)[0])
        exy = fitexy(x, y, sx, sy).slope
        ols_err.append(ols - truth)
        exy_err.append(exy - truth)
    ols_bias, exy_bias = float(np.mean(ols_err)), float(np.mean(exy_err))
    # The bound is not invented: classical regression dilution attenuates an OLS slope by
    # var(x_true) / (var(x_true) + sigma_x^2). x_true ~ U(-1,1) so var = 1/3.
    sx = 0.25
    predicted = 0.60 * ((1 / 3) / (1 / 3 + sx ** 2)) - 0.60
    assert ols_bias == pytest.approx(predicted, abs=0.02), (
        f"OLS bias {ols_bias:+.4f} vs the dilution prediction {predicted:+.4f}")
    # FITEXY removes essentially all of it. The bound asserts "unbiased", not a particular
    # residual from 25 seeds.
    assert abs(exy_bias) < 0.01, f"fitexy bias {exy_bias:+.4f}"
    assert abs(exy_bias) < abs(ols_bias) / 5.0


def test_recovers_the_true_slope_within_its_own_error():
    x, y, sx, sy, truth = _dataset(seed=101)
    r = fitexy(x, y, sx, sy)
    assert abs(r.slope - truth) < 3.0 * r.sigma_slope


def test_sigma_slope_is_the_delta_chi2_equals_one_halfwidth():
    """The quoted error must BE the stated confidence construction, not a proxy for it.

    The two sides are checked separately because they are not equal: chi2(b) stops being
    a parabola once b^2 sigma_x^2 sits in the denominator. `sigma_slope` is their mean, so
    asserting Delta-chi2 = 1 at slope +/- sigma_slope would be asserting a symmetry the
    fit does not have.
    """
    x, y, sx, sy, _ = _dataset(seed=7)
    r = fitexy(x, y, sx, sy)
    c0 = chi2_at_slope(r.slope, x, y, sx, sy)
    assert chi2_at_slope(r.slope + r.sigma_slope_hi, x, y, sx, sy) - c0 \
        == pytest.approx(1.0, abs=1e-6)
    assert chi2_at_slope(r.slope - r.sigma_slope_lo, x, y, sx, sy) - c0 \
        == pytest.approx(1.0, abs=1e-6)
    assert r.sigma_slope == pytest.approx(0.5 * (r.sigma_slope_lo + r.sigma_slope_hi))


def test_scaled_error_is_the_formal_one_times_sqrt_chi2_red():
    """Both numbers are reported so neither can be quietly preferred (RYA-161). They must
    also be the two things they claim to be."""
    rng = np.random.default_rng(3)
    n = 80
    x = rng.uniform(-1, 1, n)
    # Understate the errors so chi2_red is comfortably above 1 and the two differ.
    y = 0.3 * x + rng.normal(0, 0.20, n)
    r = fitexy(x, y, np.full(n, 0.01), np.full(n, 0.05))
    assert r.chi2_red > 3.0
    assert r.sigma_slope_scaled == pytest.approx(r.sigma_slope * np.sqrt(r.chi2_red),
                                                 rel=1e-6)
    assert r.sigma_slope_scaled > r.sigma_slope


def test_a_zero_error_is_refused_rather_than_floored():
    """A zero error is infinite weight. Silently flooring it would let one point dictate
    the fit while the caller still believed every point was weighted by its own error."""
    x, y, sx, sy, _ = _dataset(seed=11, n=20)
    sy = sy.copy()
    sy[3] = 0.0
    with pytest.raises(ValueError, match="strictly positive"):
        fitexy(x, y, sx, sy)


def test_too_few_points_is_refused():
    with pytest.raises(ValueError, match="at least 3"):
        fitexy([0.0, 1.0], [0.0, 1.0], [0.1, 0.1], [0.1, 0.1])


# ── RYA-311: one EW inversion, one call site ─────────────────────────────────────
# `invert_linemasks` was split out of `_ew_to_abundance` so the xi sweep could hold the
# line pool fixed across trial microturbulences. The whole value of that split is that it
# did NOT create a second definition of "our EW inversion" -- writing a fresh
# `determine_abundances` call in the sweep script is the RYA-845 double-definition defect,
# where two copies of one quantity are free to drift apart. grep would miss a call written
# as `getattr(ispec, "determine_abundances")`; the AST does not.

def test_exactly_one_determine_abundances_call_in_the_ew_route():
    import ast as _ast
    from pathlib import Path as _Path

    root = _Path(__file__).resolve().parents[1]
    src = (root / "pipeline" / "abundances_derive.py").read_text()
    tree = _ast.parse(src)
    calls = [n for n in _ast.walk(tree)
             if isinstance(n, _ast.Call)
             and isinstance(n.func, _ast.Attribute)
             and n.func.attr == "determine_abundances"]
    assert len(calls) == 1, (
        f"{len(calls)} ispec.determine_abundances call sites in abundances_derive "
        f"(lines {[c.lineno for c in calls]}). There must be exactly one: every EW "
        f"inversion in this project goes through invert_linemasks (RYA-311).")

    sweep = (root / "scripts" / "rya311_solar_xi_fitexy.py").read_text()
    assert "determine_abundances" not in sweep, (
        "the RYA-311 sweep calls determine_abundances directly -- it must route through "
        "abundances_derive.invert_linemasks so the sweep and production cannot disagree "
        "about what the EW inversion is.")
