"""RYA-848 — the curvature σ that IS the published C/N/O σ_stat.

Every test here asserts a RELATIONSHIP or an analytically known answer, never a magic
constant. RYA-845's lesson: `assert X == 0.100` stayed green for the entire life of a
double-count, because a constant agrees with itself no matter what the code does.

`test_the_old_arithmetic_is_red_*` are the ones that matter. A regression test that has
never been shown to FAIL on the defect it describes is decoration — so the pre-RYA-848
expression is reproduced verbatim below and each fix is asserted against it.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from pipeline.cno_synthesis import CURVATURE_PROBE_STEP_DEX as STEP
from pipeline.cno_synthesis import curvature_sigma


def _parabola(a0: float, curv: float):
    """χ²(a) = curv · (a − a0)²  — a minimum whose true σ is known in closed form."""
    return lambda a: curv * (a - a0) ** 2


def _old_sigma(chi2_fn, a_best, chi2_min, a_hi, step=0.05):
    """The pre-RYA-848 expression, VERBATIM, so the fixes can be shown to change it."""
    a_hi_probe = min(a_best + step, a_hi)
    dchi = max(chi2_fn(a_hi_probe) - chi2_min, 1e-6)
    return float(np.clip(abs(a_hi_probe - a_best) / np.sqrt(dchi), 0.0, 1.0))


# ── (a) the rescale ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("red_chi2", [1.0, 4.0, 16.0, 33.358, 66.553])
def test_sigma_scales_as_sqrt_red_chi2(red_chi2):
    """THE RELATIONSHIP: σ ∝ √red_chi2, exactly.

    Not a value check. The same curvature evaluated against a χ² that is N times too
    large must give a σ that is √N times larger, because rescaling to red_chi2 = 1 is
    the whole correction. This is what makes the solar O σ move 0.049 → 0.400.
    """
    f = _parabola(8.0, curv=1000.0)
    kw = dict(a_best=8.0, chi2_min=0.0, a_lo=5.0, a_hi=11.0, edge=False)
    unit = curvature_sigma(f, red_chi2=1.0, **kw)
    scaled = curvature_sigma(f, red_chi2=red_chi2, **kw)
    assert scaled == pytest.approx(unit * math.sqrt(red_chi2), rel=1e-12)


def test_sigma_matches_the_closed_form_for_a_known_parabola():
    """For χ² = c·(a−a₀)² rescaled by R, Δχ²_resc(step) = c·step²/R, so
    σ = step / √(c·step²/R) = √(R/c). Checked against the algebra, not a fixture."""
    for curv, red in ((1000.0, 1.0), (250.0, 9.0), (40.0, 66.553)):
        got = curvature_sigma(_parabola(8.0, curv), a_best=8.0, chi2_min=0.0,
                              red_chi2=red, a_lo=5.0, a_hi=11.0, edge=False)
        assert got == pytest.approx(math.sqrt(red / curv), rel=1e-12)


def test_a_constant_baseline_offset_cancels():
    """Δχ² is a DIFFERENCE, so an additive model-inadequacy floor must not move σ.

    This is the property that makes the rescaled curvature comparable across bands
    where an ABSOLUTE red_chi2 is not (RYA-843: the good NIR line 10145.561 sits at
    red_chi2 = 72 while a genuinely good one sits at 2.6).
    """
    curv, a0 = 500.0, 8.0
    bare = curvature_sigma(_parabola(a0, curv), a_best=a0, chi2_min=0.0, red_chi2=4.0,
                           a_lo=5.0, a_hi=11.0, edge=False)
    OFFSET = 1.0e6
    shifted = curvature_sigma(lambda a: _parabola(a0, curv)(a) + OFFSET, a_best=a0,
                              chi2_min=OFFSET, red_chi2=4.0, a_lo=5.0, a_hi=11.0,
                              edge=False)
    assert shifted == pytest.approx(bare, rel=1e-12)


# ── (b) the railed probe ──────────────────────────────────────────────────────

def test_a_railed_fit_reports_not_measurable_never_a_tight_bar():
    """A fit on the bound is unbounded on the rail side: σ is NaN, not a number."""
    a_hi = 6.61
    got = curvature_sigma(lambda a: 1000.0 - 100.0 * a, a_best=a_hi, chi2_min=1000.0 - 661.0,
                          red_chi2=50.0, a_lo=3.0, a_hi=a_hi, edge=True)
    assert math.isnan(got)


def test_the_old_arithmetic_is_red_on_the_railed_case():
    """🔴 THE DEFECT, DEMONSTRATED. The old one-sided probe clamps at a_hi, so a fit
    sitting ON a_hi gets a zero-length step and reports the TIGHTEST POSSIBLE error bar
    for the LEAST constrained possible fit. If this ever stops being 0.0, the test is
    no longer describing the bug it was written for."""
    a_hi = 6.61
    mono = lambda a: 1000.0 - 100.0 * a           # noqa: E731 — falls onto the bound
    assert _old_sigma(mono, a_best=a_hi, chi2_min=mono(a_hi), a_hi=a_hi) == 0.0
    new = curvature_sigma(mono, a_best=a_hi, chi2_min=mono(a_hi), red_chi2=50.0,
                          a_lo=3.0, a_hi=a_hi, edge=True)
    assert math.isnan(new), "the fix must not report a number where the old code reported 0.000"


def test_the_old_arithmetic_is_red_on_the_rescale():
    """🔴 The old expression has no red_chi2 in it at all, so it cannot respond to an
    uncalibrated χ². Same curvature, red_chi2 = 66.553: the old value is unchanged and
    the new one is √66.553 ≈ 8.16x larger."""
    f = _parabola(8.0, curv=1000.0)
    old = _old_sigma(f, a_best=8.0, chi2_min=0.0, a_hi=11.0)
    new = curvature_sigma(f, a_best=8.0, chi2_min=0.0, red_chi2=66.553,
                          a_lo=5.0, a_hi=11.0, edge=False)
    assert new == pytest.approx(old * math.sqrt(66.553), rel=1e-9)
    assert new > old


def test_the_probe_is_two_sided_and_takes_the_worse_side():
    """Asymmetric χ²: shallow below the minimum, steep above. The one-sided probe saw
    only the steep side and reported the tighter σ; the honest answer is the shallow
    side's, which is larger."""
    a0 = 8.0

    def asym(a):
        return 4000.0 * (a - a0) ** 2 if a >= a0 else 100.0 * (a - a0) ** 2

    got = curvature_sigma(asym, a_best=a0, chi2_min=0.0, red_chi2=1.0,
                          a_lo=5.0, a_hi=11.0, edge=False)
    shallow = math.sqrt(1.0 / 100.0)
    steep = math.sqrt(1.0 / 4000.0)
    assert got == pytest.approx(shallow, rel=1e-12)
    assert got > steep


def test_a_probe_outside_the_bracket_is_skipped_not_clamped():
    """Clamping is what produced the zero-length step. Near (but not on) the upper
    bound the upper probe is unavailable, so the answer must come from the LOWER side
    alone — and must equal the lower side computed on its own."""
    a_hi, a0 = 8.02, 8.0                       # a0 + 0.05 lies outside
    f = _parabola(a0, curv=900.0)
    got = curvature_sigma(f, a_best=a0, chi2_min=0.0, red_chi2=1.0,
                          a_lo=5.0, a_hi=a_hi, edge=False)
    assert got == pytest.approx(math.sqrt(1.0 / 900.0), rel=1e-12)
    assert STEP == pytest.approx(0.05)


def test_no_measurable_curvature_is_nan_not_a_number():
    """A perfectly flat objective has no curvature to invert. NaN says so; the old code
    said 1.000, which reads like a measured 1 dex. (Removing the clip that produces
    that 1.000 in the live path is RYA-847's, not this ticket's.)"""
    got = curvature_sigma(lambda a: 1234.5, a_best=8.0, chi2_min=1234.5, red_chi2=10.0,
                          a_lo=5.0, a_hi=11.0, edge=False)
    assert math.isnan(got)


# ── the published quantity ────────────────────────────────────────────────────

def test_the_consumer_reads_this_value_so_it_cannot_be_silently_renamed():
    """`_uncertainty_budget` publishes `sigma_fit` verbatim as stat['N'] and stat['O'].

    Pinned because the blast radius is the point: if the key stops being read, or the
    budget stops sourcing σ_stat from the fit, this test should be revisited
    deliberately rather than the coupling being rediscovered later.
    """
    import inspect

    from pipeline import cno_synthesis
    src = inspect.getsource(cno_synthesis._uncertainty_budget)
    assert "sigma_fit" in src
    for el in ("'N'", "'O'"):
        assert f"stat[{el}]" in src
