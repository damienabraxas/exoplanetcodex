"""Ryan 2026-11-19: "We do all uncertainty the Asplund way." Asplund+2021 Sect. 2.1 prices one
or two lines from the goodness of the profile fit, three or more from the standard error of the
mean -- a two-line product is never priced from the scatter of two numbers."""
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]

import rya1230_cno_budget as b  # noqa: E402


def test_one_line_is_its_own_fit_sigma():
    assert b.fit_stat([0.19]) == pytest.approx(0.19)


def test_two_lines_are_the_sigma_of_their_mean():
    assert b.fit_stat([0.19, 0.76]) == pytest.approx(math.hypot(0.19, 0.76) / 2)


def test_two_identical_lines_do_not_get_a_zero_bar():
    # The scatter of two equal numbers is 0.000 -- the bug this rule closes.
    assert b.fit_stat([0.1, 0.1]) > 0.07


def test_three_lines_are_refused_they_are_priced_by_scatter():
    assert b.MAX_FIT_PRICED_FEATURES == 2
    with pytest.raises(ValueError):
        b.fit_stat([0.1, 0.1, 0.1])


def test_the_assembler_switches_to_scatter_exactly_above_two():
    src = (ROOT / "scripts/rya1230_cno_budget.py").read_text()
    assert "len(_feats) <= MAX_FIT_PRICED_FEATURES" in src
    assert "nf > MAX_FIT_PRICED_FEATURES" in src
