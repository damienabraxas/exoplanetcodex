"""RYA-1083 — differencing two published product values must never be a valid path.

🔴 THE REAL CASE IS THE FIXTURE. The graded solar Fe I ⟨3D⟩ pair is committed under
`band_products/`, both members publish A = 7.552, and subtracting them gives exactly
0.000 while the per-line differential is +0.032 on 66 of 67 lines.

A synthetic fixture would have been easier and much weaker: this test fails the moment
anyone reintroduces the shortcut, ON THE DATA THAT ACTUALLY BIT US, and it cannot drift
from reality because it IS reality. (It also cannot be "fixed" by regenerating the
products — the collision is a coincidence of where the 34th of 67 values landed, so a
rerun would not reproduce it. That is precisely why it had to be caught in a committed
artifact rather than in a run.)
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
BP = ROOT / "data" / "results" / "band_products"
STEM = ("FeI_4200_6910_kpno_solar_atlas_solar_kpno_molecfit_corrected_SYNTH_GRADED"
        "_synth-mean3D-{}-gerber-stagger_lines.csv")

from pipeline.paired_differential import (  # noqa: E402
    UnpairableProducts, paired_differential)


@pytest.fixture(scope="module")
def real_pair():
    hi, lo = BP / STEM.format("NLTE"), BP / STEM.format("LTE")
    if not (hi.exists() and lo.exists()):
        pytest.skip("the graded ⟨3D⟩ pair is not present in this checkout")
    return pd.read_csv(hi), pd.read_csv(lo)


def test_the_two_published_medians_really_do_collide(real_pair):
    """The premise. If this ever stops being true the test below is toothless, and a
    reader deserves to be told that rather than to keep trusting a green tick."""
    hi, lo = real_pair
    assert hi.abundance.median() == lo.abundance.median() == 7.552


def test_differencing_the_published_values_gives_a_FALSE_ZERO(real_pair):
    """🔴 THE DEFECT ITSELF, pinned so it cannot be reintroduced quietly."""
    hi, lo = real_pair
    naive = hi.abundance.median() - lo.abundance.median()
    assert naive == pytest.approx(0.0, abs=1e-9)
    r = paired_differential(hi, lo)
    assert r.median == pytest.approx(0.032, abs=1e-3)
    assert r.n_nonzero == 66, "the effect is real on 66 of 67 lines"
    assert r.n_paired == 67


def test_the_result_CARRIES_what_the_naive_path_would_have_said(real_pair):
    """Showing the trap is how the next reader learns it exists. A result that only
    reports the right number teaches nobody why the wrong one is available."""
    hi, lo = real_pair
    r = paired_differential(hi, lo)
    assert r.difference_of_aggregates == pytest.approx(0.0, abs=1e-9)
    assert r.collision is True
    assert "collision" in r.as_dict()


def test_the_collision_flag_stays_quiet_when_the_two_paths_agree():
    """It must flag a real disagreement, not fire on every pair -- otherwise it is noise
    and gets ignored, which is how a warning stops working."""
    hi = pd.DataFrame({"wavelength_air_A": [5000.0, 5001.0, 5002.0],
                       "abundance": [7.50, 7.55, 7.60]})
    lo = pd.DataFrame({"wavelength_air_A": [5000.0, 5001.0, 5002.0],
                       "abundance": [7.45, 7.50, 7.55]})
    r = paired_differential(hi, lo)
    assert r.median == pytest.approx(0.05)
    assert r.difference_of_aggregates == pytest.approx(0.05)
    assert r.collision is False


def test_median_of_differences_is_NOT_difference_of_medians():
    """The statistical fact the whole ticket rests on, DEMONSTRATED rather than asserted.

    ⚠️ My first version of this test picked two frames that happened to give the same
    answer both ways, so it proved nothing while looking like it did. The median is not a
    linear operator, so a divergence has to be constructed on purpose:

        hi 7.0, 7.1, 8.0   median 7.1
        lo 7.0, 7.4, 7.5   median 7.4
        per-line          0.0, -0.3, +0.5   ->  median  0.0
        difference of medians                            -0.3

    Zero against −0.3 on three lines. The real ⟨3D⟩ pair is the same phenomenon with the
    divergence landing at +0.032 against 0.000.
    """
    hi = pd.DataFrame({"wavelength_air_A": [1.0, 2.0, 3.0], "abundance": [7.0, 7.1, 8.0]})
    lo = pd.DataFrame({"wavelength_air_A": [1.0, 2.0, 3.0], "abundance": [7.0, 7.4, 7.5]})
    r = paired_differential(hi, lo)
    assert r.median == pytest.approx(0.0)
    assert r.difference_of_aggregates == pytest.approx(-0.3)
    assert r.median != r.difference_of_aggregates
    assert r.collision is True


def test_a_tolerance_join_is_used_because_the_routes_round_differently():
    """Measured: an EXACT merge kept 6 of 55 lines where a tolerance join kept 55, because
    one route stores 4372.9817 and the other 4372.9820. An exact join would silently
    discard the pairing and report a differential over whatever happened to survive."""
    hi = pd.DataFrame({"wavelength_air_A": [4372.9820], "abundance": [7.60]})
    lo = pd.DataFrame({"wavelength_air_A": [4372.9817], "abundance": [7.50]})
    r = paired_differential(hi, lo)
    assert r.n_paired == 1 and r.median == pytest.approx(0.10)


def test_products_sharing_no_line_REFUSE_rather_than_return_zero():
    """An empty pairing must not aggregate to 0.0 -- that is the same false zero by
    another route, and `np.median([])` is exactly the kind of quiet nan a caller ignores."""
    hi = pd.DataFrame({"wavelength_air_A": [5000.0], "abundance": [7.5]})
    lo = pd.DataFrame({"wavelength_air_A": [9999.0], "abundance": [7.4]})
    with pytest.raises(UnpairableProducts):
        paired_differential(hi, lo)


def test_lines_excluded_from_one_aggregate_do_not_contribute():
    """A line the comparand excluded is not part of the comparand."""
    hi = pd.DataFrame({"wavelength_air_A": [1.0, 2.0], "abundance": [7.6, 7.6],
                       "in_aggregate": [True, True]})
    lo = pd.DataFrame({"wavelength_air_A": [1.0, 2.0], "abundance": [7.5, 7.0],
                       "in_aggregate": [True, False]})
    r = paired_differential(hi, lo)
    assert r.n_paired == 1 and r.median == pytest.approx(0.10)
