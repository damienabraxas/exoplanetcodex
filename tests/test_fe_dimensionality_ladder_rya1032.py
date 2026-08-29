"""RYA-1032 — the Fe 1D→⟨3D⟩ climb is DIMENSIONALITY, and the axis contract proves it.

🔴 WHAT THIS PROTECTS. The headline "the atmosphere step is ~2.6x the model-family
spread" is only meaningful if each difference holds the axes it claims to hold. An
early version of the deriver asserted "scale fixed" across `1D-LTE` vs `<3D>-LTE` —
which holds nothing, because the registry's `scale` is two axes (dimensionality x
LTE/NLTE) wedged into one string. The guard caught it. These tests keep it caught.

The other thing protected here is a REFUSAL: the ⟨3D⟩ NLTE effect must never be
computed by differencing the two published cells. On this holding both are medians
that collided on exactly 7.552, so that subtraction returns 0.000 for a real +0.032
per-line shift (RYA-1099 Finding 3; the correct statistic is RYA-1083's paired
differential). A test that only checked numbers would not notice the deriver quietly
starting to publish that zero as physics.
"""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))

from rya1032_fe_dimensionality_ladder import (  # noqa: E402
    LadderError, analyse, load_registry, load_rungs, step)


@pytest.fixture(scope="module")
def a():
    return analyse()


def test_every_rung_is_present_and_live(a):
    """A partial ladder must raise, not silently report fewer rungs."""
    assert len(a["ladder"]) == 7
    assert all(r["A"] is not None for r in a["ladder"])
    with pytest.raises(LadderError, match="missing rungs"):
        analyse("a_holding_that_does_not_exist")


def test_dimensionality_step_holds_family_and_lte_ness(a):
    s = a["dimensionality_step"]
    assert s["axes_hi"]["family"] == s["axes_lo"]["family"] == "gerber"
    assert s["axes_hi"]["nlte"] == s["axes_lo"]["nlte"] == "LTE"
    assert s["axes_hi"]["dim"] != s["axes_lo"]["dim"]
    assert s["delta"] == pytest.approx(0.101, abs=5e-4)


def test_family_spread_is_deconfounded_by_a_measured_nuisance(a):
    """The family spread crosses atlas9 -> marcs-ges, so that must be measured.

    Without this control the ratio is not defensible: the spread could have been
    mostly atmosphere. It is not — the LTE nuisance is +0.004.
    """
    assert a["atmosphere_nuisance"]["delta"] == pytest.approx(0.004, abs=5e-4)
    assert a["model_family_spread"]["delta"] == pytest.approx(0.043, abs=5e-4)
    assert a["model_family_spread_deconfounded"] == pytest.approx(0.039, abs=5e-4)
    # the nuisance must be small compared to the spread it is correcting,
    # otherwise "model family" is the wrong name for that difference
    assert abs(a["atmosphere_nuisance"]["delta"]) < 0.2 * abs(
        a["model_family_spread"]["delta"])


def test_the_verdict_dimensionality_beats_model_family(a):
    assert a["ratio_dimensionality_over_family"] >= 2.0
    assert a["ratio_dimensionality_over_family_deconfounded"] >= 2.0
    assert a["dimensionality_step"]["delta"] > 2 * abs(
        a["model_family_spread"]["delta"])


def test_the_ladder_is_non_monotonic_mean3d_overshoots_full_3d(a):
    """⟨3D⟩ lands ABOVE Amarsi full-3D — the mean-3D approximation over-corrects."""
    assert a["overshoot_vs_amarsi_full3d"] == pytest.approx(0.040, abs=1e-3)
    assert a["overshoot_vs_amarsi_full3d"] > 0
    assert a["ladder_monotonic_toward_amarsi"] is False


def test_pool_size_differences_are_declared_not_hidden(a):
    """Bergemann is n=61 against everything else's n=67 — say so."""
    assert a["model_family_spread"]["same_pool"] is False
    assert a["dimensionality_step"]["same_pool"] is True


def test_it_refuses_to_compute_the_nlte_effect_from_two_medians(a):
    """The published pair collided; differencing it is the RYA-1099 defect."""
    e = a["mean3d_nlte_effect"]
    assert e["computed_here"] is False
    assert e["median_collision"] is True
    assert e["published_difference_of_medians"] == 0.0
    assert "paired_differential" in e["why"]


def test_the_axis_guard_actually_fires():
    """NEGATIVE CONTROL — a guard that never fails proves nothing.

    Assert a contract that is false for this pair and require the raise. Without
    this, every assertion above could pass against a `step()` whose checks had
    been silently removed.
    """
    reg, (_, rungs) = load_registry(), load_rungs("solar_kpno_molecfit_corrected")

    # dimensionality DOES change between these two, so "fixed" must raise
    with pytest.raises(LadderError, match="must be FIXED"):
        step(reg, rungs, "synth-mean3D-LTE-gerber-stagger", "synth-1D-LTE-gerber",
             fixed=("dim",), label="bogus")

    # family does NOT change between these two, so "varies" must raise
    with pytest.raises(LadderError, match="must VARY"):
        step(reg, rungs, "synth-mean3D-LTE-gerber-stagger", "synth-1D-LTE-gerber",
             varies="family", label="bogus")
