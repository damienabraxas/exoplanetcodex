"""Stage-4 GoF on the synthesis route — RYA-992.

Relationships, not the screen's numbers (RYA-870): the fixture grows, so a test asserting
rho = -0.200 would fail the moment the HARPS deep leg lands. The one number pinned here is
`SYNTH_GOF_REF_CUT`, which is a statement about the objective (chi2 doubles) rather than a
measurement — and the per-arm scales are asserted as RELATIONSHIPS, never as 1.357.
"""
from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from pipeline.synth_gof import (ARM_SCALE, SYNTH_GOF_REF_ARM, SYNTH_GOF_REF_CUT,
                                measure_arm_scale, synth_gof_cut, synth_gof_ok,
                                synth_gof_reason)

KP = SYNTH_GOF_REF_ARM


def test_the_reference_threshold_is_the_point_where_chi2_doubles():
    """🔴 DERIVED, NOT A QUOTA. `frac_rise_weaker` is a FRACTIONAL rise, so 1.0 is the
    objective merely doubling. A percentile cut would move with the sample; this does not."""
    assert SYNTH_GOF_REF_CUT == 1.0
    assert synth_gof_cut(SYNTH_GOF_REF_ARM) == SYNTH_GOF_REF_CUT


def test_a_fit_that_barely_resists_a_change_in_A_is_refused():
    assert not synth_gof_ok(0.0, KP)
    assert not synth_gof_ok(0.99, KP)
    assert synth_gof_ok(1.0, KP)
    assert synth_gof_ok(50.0, KP)


def test_a_missing_statistic_is_NOT_a_pass():
    """An unscored fit is not a passing one. Admitting `None` would let exactly the
    population stage 4 exists to score walk through it (RYA-833)."""
    assert not synth_gof_ok(None, KP)
    assert not synth_gof_ok(float("nan"), KP)
    assert "UNSCORED" in synth_gof_reason(None, KP)


def test_every_refusal_carries_a_coded_reason():
    """Quarantine, never cull (RYA-711): a rejected line says why, in a parseable code."""
    assert synth_gof_reason(0.4, KP).startswith("SYNTH-GOF:")
    assert synth_gof_reason(None, KP).startswith("SYNTH-GOF-UNSCORED:")
    assert "0.40" in synth_gof_reason(0.4, KP)


def test_the_refusal_names_the_arm_and_its_own_cut():
    """A line rejected on HARPS must not be told it failed the Kitt Peak bar."""
    msg = synth_gof_reason(1.1, "harps")
    assert "harps" in msg
    assert f"{synth_gof_cut('harps'):.3f}" in msg


def test_the_criterion_is_offset_invariant():
    """F2. `frac_rise_weaker` describes the SHAPE of the chi2 surface under a change in A,
    so shifting every abundance by a constant cannot change it. `instrument` is a label,
    not a measurement, so it cannot carry an offset in either."""
    params = set(inspect.signature(synth_gof_ok).parameters)
    assert params == {"frac_rise_weaker", "instrument"}, (
        "the criterion must depend on the chi2-shape statistic and the ARM ALONE — an "
        "abundance or a reference reaching this signature would make stage 4 "
        "offset-dependent")


def test_no_reference_abundance_is_reachable():
    """🔒 RYA-161, enforced by signature — including on the arm-scale derivation."""
    for fn in (synth_gof_ok, synth_gof_reason, synth_gof_cut, measure_arm_scale):
        for banned in ("reference", "gold", "expected", "target", "truth", "abundance"):
            assert not any(banned in p for p in inspect.signature(fn).parameters), fn


def test_it_rejects_a_minority_not_a_quota():
    """Sanity on the shape of the cut: a stage-4 gate that rejected most of a GRADED anchor
    would be measuring something other than fit quality."""
    graded_like = [0.4, 0.9, 1.2, 2.5, 4.8, 9.0, 20.0, 33.0]
    kept = sum(synth_gof_ok(v, KP) for v in graded_like)
    assert kept >= len(graded_like) * 0.6


# --------------------------------------------------------------------------- per-arm cut

def test_the_cut_is_not_one_number_across_arms():
    """🔴 THE WHOLE POINT. A flat 1.0 for every instrument is the defect this replaces."""
    cuts = {a: synth_gof_cut(a) for a in ARM_SCALE}
    assert len(set(cuts.values())) > 1, (
        f"every arm derived the same cut {cuts} — a fixed constant wearing a registry")
    assert cuts["harps"] != pytest.approx(SYNTH_GOF_REF_CUT, rel=0.02)


def test_harps_bar_went_UP_not_down():
    """RYA-995 expected the HARPS cut to FALL because its frac_rise distribution sits
    lower. The RYA-986 two-arm pool says the opposite: at a given frac_rise a HARPS line
    is never better than a Kitt Peak line, so its bar rises. Pinning the DIRECTION guards
    the finding; the magnitude is free to be re-measured."""
    assert synth_gof_cut("harps") > synth_gof_cut(SYNTH_GOF_REF_ARM)


def test_an_unmeasured_arm_raises_instead_of_defaulting():
    """Falling back to the reference cut would be the transferred constant coming back,
    silently, for every instrument nobody remembered to measure (RYA-833)."""
    with pytest.raises(KeyError) as e:
        synth_gof_cut("espresso")
    assert "espresso" in str(e.value)
    with pytest.raises(KeyError):
        synth_gof_ok(99.0, "espresso")       # even a value that would obviously pass


def test_the_reference_arm_is_unity_by_definition():
    assert ARM_SCALE[SYNTH_GOF_REF_ARM].scale == 1.0


def test_every_registry_row_says_what_it_was_measured_on():
    """A scale multiplied into every graded line owes its provenance (RYA-873)."""
    for name, row in ARM_SCALE.items():
        assert row.instrument == name
        assert row.scale > 0
        assert row.measured_on.strip() and row.note.strip()


# ------------------------------------------------ the derivation, on a synthetic fixture

def _synthetic_arm(n_pix: int, *, n=400, seed=0):
    """One arm's (frac_rise, sigma_A) columns, from the algebra the module documents.

    frac_rise = B**2 / (sigma_A**2 * dof). Only `n_pix` — the pixels in the fit window —
    differs between the two arms here: the abundance constraint `sigma_A` is drawn from
    the SAME distribution, so anything the derivation reports is the geometry alone.
    """
    rng = np.random.default_rng(seed)
    sigma = 10 ** rng.normal(-1.1, 0.25, n)          # ~0.08 dex, lognormal
    B = rng.normal(2.95, 0.15, n)                    # the bracket, near-constant
    frac = B ** 2 / (sigma ** 2 * n_pix)
    return frac, sigma


def test_two_resolutions_derive_two_different_cuts():
    """🔴 THE TEST THE BRIEF ASKS FOR. Same fits, same sigma_A distribution, only the
    window sampling differs — and a fixed 1.0 cannot tell them apart."""
    coarse = _synthetic_arm(1200, seed=1)            # a lower-resolution arm
    fine = _synthetic_arm(2400, seed=1)              # twice the pixels in the window
    scale = measure_arm_scale(*coarse, *fine)

    # dof halved => k = B**2/dof doubled => the coarse arm's cut is ~2x the fine arm's.
    assert scale == pytest.approx(2.0, rel=0.05), scale
    assert scale != pytest.approx(1.0, rel=0.10), (
        "a fixed cut assumes scale == 1; this fixture must refuse that")


def test_the_derivation_ignores_a_pure_quality_difference():
    """A scale is window GEOMETRY, not fit quality. Two arms with identical sampling but
    genuinely different constraint must derive the SAME cut — the worse arm then fails it
    more often, which is the gate working, not the gate needing loosening."""
    good_f, good_s = _synthetic_arm(1800, seed=2)
    bad_f, bad_s = good_f / 4.0, good_s * 2.0        # 4x flatter chi2, 2x worse sigma_A
    assert measure_arm_scale(bad_f, bad_s, good_f, good_s) == pytest.approx(1.0, rel=0.02)


def test_the_derivation_refuses_a_pool_too_small_to_measure_on():
    f, s = _synthetic_arm(1500, n=12, seed=3)
    ref = _synthetic_arm(1500, seed=3)
    with pytest.raises(ValueError, match="noise"):
        measure_arm_scale(f, s, *ref)


def test_the_derivation_refuses_mismatched_columns():
    f, s = _synthetic_arm(1500, seed=4)
    with pytest.raises(ValueError, match="differ in length"):
        measure_arm_scale(f, s[:-1], f, s)


def test_the_derivation_is_a_ratio_not_a_percentile():
    """A quota would put a fixed FRACTION of each arm below the cut. This does not: feed it
    an arm whose lines are uniformly better-constrained and the cut stays put, so the
    better arm simply keeps more."""
    f, s = _synthetic_arm(1600, seed=5)
    better_f, better_s = f * 4.0, s / 2.0            # same geometry, 4x steeper chi2
    assert measure_arm_scale(better_f, better_s, f, s) == pytest.approx(1.0, rel=0.02)
    cut = SYNTH_GOF_REF_CUT
    assert (better_f >= cut).mean() > (f >= cut).mean()
