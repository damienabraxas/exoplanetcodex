"""
pipeline/fitexy.py — RYA-311
============================
Straight-line fit with uncertainties in BOTH coordinates (Press et al. 1992,
*Numerical Recipes* §15.3, `fitexy`), and the formal error on its slope.

🔴 WHY THIS IS NOT `np.polyfit`.

The microturbulence solve zeroes the slope of A(Fe I) against reduced equivalent width
log(EW/λ). Both axes carry the SAME measurement error: σ_EW propagates into the abscissa
directly (σ_REW = σ_EW / (EW ln10)) and into the ordinate through the curve of growth
(σ_A from re-inverting at EW ± 1σ_EW). An ordinary least-squares fit assumes the abscissa
is exact. It is not, the two errors are perfectly correlated in origin, and the resulting
slope is biased towards zero (regression dilution) while its quoted error is too small.

Mucciarelli 2011 (A&A 528 A44) makes precisely this point for the vmic solve and
prescribes the Press both-axis fit. That is the method this module implements, so
`delta_xi` is a propagated formal error rather than an OLS number that understates it.

THE FIT
-------
Minimise, over intercept `a` and slope `b`,

    chi2(a, b) = SUM_i  (y_i - a - b x_i)^2 / (sigma_y_i^2 + b^2 sigma_x_i^2)

The `b^2 sigma_x^2` in the denominator is the whole difference from OLS: an abscissa
error acts on the residual scaled by the slope itself, which is what makes the problem
non-linear in `b`. For any fixed `b` the intercept is closed-form,

    a(b) = SUM_i w_i (y_i - b x_i) / SUM_i w_i ,    w_i = 1 / (sigma_y_i^2 + b^2 sigma_x_i^2)

so the minimisation is one-dimensional in `b`.

THE SLOPE ERROR
---------------
`sigma_b` is the Δchi2 = 1 half-width about the minimum, with `a` re-minimised at every
trial `b` (the standard one-interesting-parameter confidence interval, and what NR's
`fitexy` returns). It is a FORMAL error: it says how well the slope is determined GIVEN
the stated point errors, and nothing else.

⚠️ IT IS ONLY THE WHOLE ANSWER WHEN THE POINT ERRORS EXPLAIN THE SCATTER. Fe I line
abundances scatter by ~0.14 dex in the Sun (RYA-407) for reasons that are not EW
measurement error at all -- gf errors, blends, continuum placement. When that dominates,
chi2_red >> 1 and the Δchi2 = 1 error is an underestimate of how well the slope is really
known. `FitExyResult` therefore carries BOTH numbers and neither is hidden:

  * `sigma_slope`        -- the formal Δchi2 = 1 error. The method's own answer.
  * `sigma_slope_scaled` -- the same error times sqrt(chi2_red), the usual inflation for
                            a fit whose residuals exceed its stated errors.

Report the formal one as the result and the scaled one beside it; do not silently swap
whichever is more convenient (RYA-161).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FitExyResult:
    """One both-axis straight-line fit."""
    slope: float
    intercept: float
    sigma_slope: float          #: formal Delta-chi2 = 1 error on the slope
    sigma_slope_scaled: float   #: sigma_slope * sqrt(chi2_red) -- the inflated alternative
    chi2: float
    dof: int

    @property
    def chi2_red(self) -> float:
        return self.chi2 / self.dof if self.dof > 0 else float("nan")


def _weights(b: float, sx: np.ndarray, sy: np.ndarray) -> np.ndarray:
    return 1.0 / (sy ** 2 + (b * sx) ** 2)


def _intercept(b: float, x: np.ndarray, y: np.ndarray,
               sx: np.ndarray, sy: np.ndarray) -> float:
    w = _weights(b, sx, sy)
    return float(np.sum(w * (y - b * x)) / np.sum(w))


def chi2_at_slope(b: float, x: np.ndarray, y: np.ndarray,
                  sx: np.ndarray, sy: np.ndarray) -> float:
    """chi2 minimised over the intercept at a FIXED slope `b`."""
    a = _intercept(b, x, y, sx, sy)
    w = _weights(b, sx, sy)
    return float(np.sum(w * (y - a - b * x) ** 2))


def fitexy(x, y, sigma_x, sigma_y, *,
           n_scan: int = 4001, span: float = 8.0) -> FitExyResult:
    """Fit y = a + b x with errors on both axes; return the slope and its formal error.

    `sigma_x` and `sigma_y` must be strictly positive -- a zero error would give a point
    infinite weight, which is not an approximation this fit makes. Callers with a genuinely
    error-free axis should say so and use OLS rather than pass zeros here.

    The 1-D minimisation is a dense scan followed by a parabolic refinement rather than a
    library optimiser: chi2(b) is smooth but not convex for badly-conditioned data, and a
    local optimiser started at the OLS slope can settle in the wrong basin. The scan span
    is `span` times the OLS slope error about the OLS slope, which brackets the minimum for
    any data where the both-axis correction is a correction rather than a different answer.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sx = np.asarray(sigma_x, dtype=float)
    sy = np.asarray(sigma_y, dtype=float)
    if not (x.shape == y.shape == sx.shape == sy.shape):
        raise ValueError("fitexy: x, y, sigma_x, sigma_y must have the same shape")
    n = x.size
    if n < 3:
        raise ValueError(f"fitexy: need at least 3 points to fit 2 parameters and quote "
                         f"an error; got {n}")
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)) \
            or np.any(~np.isfinite(sx)) or np.any(~np.isfinite(sy)):
        raise ValueError("fitexy: non-finite input")
    if np.any(sx <= 0) or np.any(sy <= 0):
        raise ValueError("fitexy: sigma_x and sigma_y must be strictly positive")

    # Bracket from the OLS solution (y-errors only), which is the b -> small limit.
    w0 = 1.0 / sy ** 2
    sw = np.sum(w0)
    sxw, syw = np.sum(w0 * x), np.sum(w0 * y)
    sxx, sxy = np.sum(w0 * x * x), np.sum(w0 * x * y)
    denom = sw * sxx - sxw ** 2
    if denom <= 0:
        raise ValueError("fitexy: degenerate abscissa (all x identical?)")
    b_ols = (sw * sxy - sxw * syw) / denom
    sb_ols = np.sqrt(sw / denom)
    # A pure-OLS error can be tiny; keep the scan wide enough to see the both-axis minimum,
    # which moves by roughly the regression-dilution factor.
    half = max(span * sb_ols, 4.0 * abs(b_ols), 1e-6)

    grid = np.linspace(b_ols - half, b_ols + half, int(n_scan))
    c2 = np.array([chi2_at_slope(b, x, y, sx, sy) for b in grid])
    k = int(np.argmin(c2))
    if k in (0, len(grid) - 1):
        raise ValueError("fitexy: chi2 minimum sits on the scan boundary -- widen `span`")
    # Parabolic refinement through the three points about the discrete minimum.
    b0, b1, b2 = grid[k - 1], grid[k], grid[k + 1]
    f0, f1, f2 = c2[k - 1], c2[k], c2[k + 1]
    d = (f0 - 2 * f1 + f2)
    b_hat = b1 + 0.5 * (b0 - b2) * (f0 - f2) / (4 * d) if d > 0 else b1
    c2_min = chi2_at_slope(b_hat, x, y, sx, sy)

    # Delta chi2 = 1 half-width, found by bisection on each side of the minimum.
    def _edge(direction: int) -> float:
        lo, hi = b_hat, b_hat + direction * max(half, 1e-9)
        # Grow until chi2 exceeds the +1 contour, then bisect.
        for _ in range(60):
            if chi2_at_slope(hi, x, y, sx, sy) - c2_min >= 1.0:
                break
            hi = b_hat + (hi - b_hat) * 2.0
        else:
            return float("nan")
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if chi2_at_slope(mid, x, y, sx, sy) - c2_min >= 1.0:
                hi = mid
            else:
                lo = mid
        return abs(0.5 * (lo + hi) - b_hat)

    up, dn = _edge(+1), _edge(-1)
    sigma_b = float(np.nanmean([up, dn]))
    dof = n - 2
    chi2_red = c2_min / dof if dof > 0 else float("nan")
    scaled = sigma_b * np.sqrt(chi2_red) if np.isfinite(chi2_red) and chi2_red > 0 else sigma_b
    return FitExyResult(slope=float(b_hat),
                        intercept=_intercept(b_hat, x, y, sx, sy),
                        sigma_slope=sigma_b,
                        sigma_slope_scaled=float(scaled),
                        chi2=float(c2_min),
                        dof=int(dof))
