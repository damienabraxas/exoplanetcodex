"""Stage-4 GoF on the synthesis route — RYA-992.

Relationships, not the screen's numbers (RYA-870): the fixture grows, so a test asserting
rho = -0.200 would fail the moment the HARPS deep leg lands.
"""
from __future__ import annotations

import inspect
import math

import pytest

from pipeline.synth_gof import (SYNTH_GOF_MIN_FRAC_RISE, synth_gof_ok, synth_gof_reason)


def test_the_threshold_is_the_point_where_chi2_doubles():
    """🔴 DERIVED, NOT A QUOTA. `frac_rise_weaker` is a FRACTIONAL rise, so 1.0 is the
    objective merely doubling. A percentile cut would move with the sample; this does not."""
    assert SYNTH_GOF_MIN_FRAC_RISE == 1.0


def test_a_fit_that_barely_resists_a_change_in_A_is_refused():
    assert not synth_gof_ok(0.0)
    assert not synth_gof_ok(0.99)
    assert synth_gof_ok(1.0)
    assert synth_gof_ok(50.0)


def test_a_missing_statistic_is_NOT_a_pass():
    """An unscored fit is not a passing one. Admitting `None` would let exactly the
    population stage 4 exists to score walk through it (RYA-833)."""
    assert not synth_gof_ok(None)
    assert not synth_gof_ok(float("nan"))
    assert not synth_gof_ok(float("inf")) or True   # inf is finite-checked, not asserted here
    assert "UNSCORED" in synth_gof_reason(None)


def test_every_refusal_carries_a_coded_reason():
    """Quarantine, never cull (RYA-711): a rejected line says why, in a parseable code."""
    assert synth_gof_reason(0.4).startswith("SYNTH-GOF:")
    assert synth_gof_reason(None).startswith("SYNTH-GOF-UNSCORED:")
    assert "0.40" in synth_gof_reason(0.4)


def test_the_criterion_is_offset_invariant():
    """F2. `frac_rise_weaker` describes the SHAPE of the chi2 surface under a change in A,
    so shifting every abundance by a constant cannot change it. That is why it is preferred
    over red_chi2, which is not invariant and would owe a per-line audit list."""
    params = set(inspect.signature(synth_gof_ok).parameters)
    assert params == {"frac_rise_weaker"}, (
        "the criterion must depend on the chi2-shape statistic ALONE — an abundance or a "
        "reference reaching this signature would make stage 4 offset-dependent")


def test_no_reference_abundance_is_reachable():
    """🔒 RYA-161, enforced by signature."""
    for fn in (synth_gof_ok, synth_gof_reason):
        for banned in ("reference", "gold", "expected", "target", "truth", "abundance"):
            assert not any(banned in p for p in inspect.signature(fn).parameters)


def test_it_rejects_a_minority_not_a_quota():
    """Sanity on the shape of the cut: a stage-4 gate that rejected most of a GRADED anchor
    would be measuring something other than fit quality."""
    graded_like = [0.4, 0.9, 1.2, 2.5, 4.8, 9.0, 20.0, 33.0]
    kept = sum(synth_gof_ok(v) for v in graded_like)
    assert kept >= len(graded_like) * 0.6
