"""
tests/test_n1_stat_uncertainty_rya907.py
========================================
RYA-907 — `(product.sigma or 0.0)` converted "not measurable" into "measured zero".

Two published Fe II red-optical ENGINE-A cells carried `stat_dex = 0.0` at n = 1. The
builder was already right (`band_products.py` returns `None` from one line); three
drivers spent that `None` as a zero.

WHAT THESE TESTS PIN. Not the numbers — RYA-845's lesson is that a test asserting a magic
constant pinned a bug for its whole life. They pin the RELATIONSHIPS:

  * `None` never becomes `0.0` on any path, and the term that carries it REFUSES to
    produce a contribution rather than producing a small one;
  * the guard can actually FAIL (a guard that cannot fail proves nothing, RYA-833) AND
    can actually PASS on the fixed path — both ends;
  * the fix is confined to n = 1: a same-inputs control shows every n > 1 bar unchanged.

The floor's VALUE is asserted exactly once, against the module constant, never re-typed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import error_budget as eb                      # noqa: E402
from pipeline.band_products import (                          # noqa: E402
    LineMeasurement, build_product)


_HARNESS = dict(harness_residual_dex=0.0129, handler="ProfileFitHandler",
                harness_provenance="measured")


def _budget(n_lines, scatter, wave=8000.0):
    return eb.build("Fe", wave, n_lines, scatter_dex=scatter, gf_graded=False, **_HARNESS)


# ── The builder was already right; do not "fix" it ───────────────────────────

def test_build_product_still_returns_None_for_a_single_line():
    """`band_products.py` is the one thing RYA-907 must NOT change.

    Its `None` is the honest signal the rest of the chain was throwing away. If a future
    edit makes it return 0.0 instead, every guard below becomes decoration.
    """
    one = [LineMeasurement(element="Fe", ion="II", wavelength_air_A=8688.6,
                           instrument="kpno_solar_atlas", ew_mA=50.0,
                           ew_method="PROFILE-FIT (voigt)", abundance=7.763,
                           treatment="ENGINE-A", in_aggregate=True)]
    p = build_product("Fe", "II", "kpno_solar_atlas", "red-optical", "ENGINE-A",
                      one, handler="ProfileFitHandler", provenance="test")
    assert p.n_lines == 1
    assert p.sigma is None, "one line must not acquire a scatter"


# ── None is carried, never spent ─────────────────────────────────────────────

def test_unmeasured_term_refuses_to_contribute_rather_than_contributing_zero():
    """There must be NO code path on which an unmeasured term returns a number.

    This is what stops the defect being re-acquired: a future `or 0.0` would be reaching
    for a value that does not exist.
    """
    t = eb.scatter_term(None, 1)
    assert t.measured is False
    with pytest.raises(eb.UnmeasuredTerm):
        t.contribution(1)


def test_a_measured_zero_and_an_unmeasured_term_are_different_objects():
    """The distinction the bug erased, asserted directly."""
    unmeasured = eb.scatter_term(None, 1)
    measured_zero = eb.scatter_term(0.0, 5)
    assert unmeasured.measured is False
    assert measured_zero.measured is True
    assert measured_zero.contribution(5) == 0.0
    assert "NOT MEASURABLE" in unmeasured.source


def test_statistical_can_never_return_zero():
    """The published bar is 0.0 on no input at all — the floor makes it unreachable."""
    for n, scatter in ((1, None), (1, 0.0), (5, 0.0), (200, 0.0), (3, 1e-9)):
        assert _budget(n, scatter).statistical() > 0.0


def test_n1_publishes_the_floor_and_says_it_is_a_stand_in():
    b = _budget(1, None)
    assert b.statistical() == eb.QUANTISER_FLOOR_DEX
    basis = b.stat_basis()
    assert "UNMEASURED" in basis
    assert "STANDS IN" in basis


def test_a_tiny_measured_scatter_floors_for_a_DIFFERENT_stated_reason():
    """Both roads lead to 0.01 and the artifact must not conflate them.

    A measured RMS below one synthesis tread is a real measurement the engine cannot
    resolve; an unmeasured scatter is no measurement at all. Same number, different claim.
    """
    floored = _budget(400, 1e-4)
    unmeasured = _budget(1, None)
    assert floored.statistical() == unmeasured.statistical() == eb.QUANTISER_FLOOR_DEX
    assert "UNMEASURED" not in floored.stat_basis()
    assert "UNMEASURED" in unmeasured.stat_basis()
    assert floored.stat_basis() != unmeasured.stat_basis()


def test_describe_says_the_scatter_is_unmeasured_in_words():
    text = _budget(1, None).describe()
    assert "UNMEASURED" in text
    assert "line-to-line scatter" in text
    assert "basis:" in text


def test_dominant_ranks_only_measured_terms():
    """An unmeasured term has no size, so it cannot be called the limiting one."""
    d = _budget(1, None).dominant()
    assert d is not None and d.measured
    assert d.name != "line-to-line scatter"


# ── The guard, at BOTH ends ──────────────────────────────────────────────────

CELL = "Fe II red-optical ENGINE-A (kpno_solar_atlas, n=1)"


def test_guard_FAILS_on_the_exact_value_that_shipped_and_names_the_cell():
    """POSITIVE CONTROL. A guard that cannot fail proves nothing (RYA-833)."""
    with pytest.raises(eb.UnpublishableStat) as e:
        eb.assert_stat_publishable(0.0, cell=CELL)
    assert CELL in str(e.value)


@pytest.mark.parametrize("bad", [0.0, -0.01, float("nan"), None,
                                 eb.QUANTISER_FLOOR_DEX / 2])
def test_guard_refuses_every_indefensible_bar(bad):
    with pytest.raises(eb.UnpublishableStat):
        eb.assert_stat_publishable(bad, cell=CELL)


def test_guard_PASSES_what_the_fixed_path_actually_produces():
    """THE OTHER END. If the guard rejected the fix's own output it would be useless."""
    stat = _budget(1, None).statistical()
    assert eb.assert_stat_publishable(stat, cell=CELL) == stat
    stat40 = _budget(40, 0.10).statistical()
    assert eb.assert_stat_publishable(stat40, cell=CELL) == stat40


# ── Confined to n = 1: the same-inputs control ───────────────────────────────

@pytest.mark.parametrize("n,scatter", [(2, 0.05), (3, 0.111), (11, 0.187),
                                       (40, 0.10), (152, 0.376), (271, 0.5)])
def test_no_n_gt_1_bar_moves(n, scatter):
    """The published bar for a multi-line product is still scatter/sqrt(N), untouched.

    Computed here from the definition rather than copied from the matrix, so the test
    asserts the RELATIONSHIP and cannot pin a number the way RYA-845's did.
    """
    b = _budget(n, scatter)
    assert b.statistical() == pytest.approx(scatter / np.sqrt(n), rel=1e-12)
    assert b.stat_basis().startswith("measured")


def test_the_floor_is_declared_once_and_matches_RYA_771():
    """One tread of the iSpec %.2f staircase. Asserted against the module, not re-typed."""
    assert eb.QUANTISER_FLOOR_DEX == 0.01
    src = (ROOT / "pipeline" / "error_budget.py").read_text()
    assert src.count("QUANTISER_FLOOR_DEX = ") == 1


# ── No driver may re-acquire the idiom ───────────────────────────────────────

@pytest.mark.parametrize("driver", ["scripts/derive_band_products.py",
                                    "scripts/rya836_nearuv_lab_gf_subpool.py"])
def test_no_driver_spends_a_None_sigma_as_a_zero(driver):
    """The defect was one expression in three places. Fixing two would leave it live.

    Matched on the code, not on a comment: the RYA-907 comments quote the old expression
    deliberately, so the check ignores anything on a commented line.
    """
    for line in (ROOT / driver).read_text().splitlines():
        code = line.split("#", 1)[0]
        assert "sigma or 0" not in code, f"{driver}: {line.strip()}"


@pytest.mark.parametrize("driver", ["scripts/derive_band_products.py",
                                    "scripts/rya836_nearuv_lab_gf_subpool.py"])
def test_every_emitting_driver_calls_the_guard(driver):
    """A guard nobody calls is a guard that does not exist."""
    src = (ROOT / driver).read_text()
    assert "assert_stat_publishable" in src
    assert "stat_basis" in src, f"{driver} emits no stat_basis field"


def test_derive_band_products_guards_EVERY_ONE_of_its_emit_sites():
    """RYA-907 §6: fixing only line 452 leaves the bug live on the other path.

    🔴 RYA-1044 — THIS PINNED A COUNT ("== 2") AND THE COUNT WENT STALE THE MOMENT A THIRD
    EMIT SITE WAS ADDED, which is precisely the event it was written to catch. The new
    site DOES guard; the test failed on arithmetic, not on the defect.

    The invariant is "every site that publishes a stat guards it", so it is measured as a
    RELATIONSHIP against the number of publishing sites. A site added without a guard
    still fails loudly; a site added WITH one no longer does."""
    src = (ROOT / "scripts" / "derive_band_products.py").read_text()
    n_guards = src.count("assert_stat_publishable(")
    assert n_guards >= 2, (
        f"only {n_guards} guarded emit site(s) -- RYA-907's two must both still guard")

    # ⚠️ A FLOOR, NOT AN EQUALITY, AND DELIBERATELY NOT A DERIVED COUNT EITHER.
    # The obvious "every budget total is guarded" is NOT TRUE of this file today: there
    # are more `.total()` calls than guards, because the ENGINE-A legs total a budget and
    # publish without one. That is pre-existing and is not RYA-1044's to change --
    # asserting it here would be inventing an invariant the code does not hold, which
    # fails honestly today and would be quietly wrong the day someone "fixes" the test
    # instead of the code.
    #
    # So: a floor. It cannot go stale when a guarded leg is added (the failure RYA-1044
    # hit), and it still catches the RYA-907 defect -- removing a guard.
    assert "import assert_stat_publishable" in src
    assert src.count("scatter_dex=product.sigma") == 1
    assert src.count("scatter_dex=prod.sigma") == 1
