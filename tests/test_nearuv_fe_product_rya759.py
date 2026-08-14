"""
tests/test_nearuv_fe_product_rya759.py
RYA-759 — the near-UV Fe product's line selection.

These cover the two things that actually went wrong while building it, both of which
produced a WRONG-BUT-PLAUSIBLE product rather than an error until the code was made to
fail loudly:

  1. iSpec spells neutral iron 'Fe 1' (space, arabic), NOT 'Fe I'. Matching the roman
     form selects ZERO rows — and an empty selection aggregates to an empty product, not
     to a crash, unless something refuses.
  2. `theoretical_ew` is ALL ZERO in this list (our VALD converter never populates it).
     Ranking or windowing on it silently ranks everything equal. `theoretical_depth` is
     the populated column, and selection must use it.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from rya759_nearuv_fe_product import (       # noqa: E402
    DEPTH_CEIL, DEPTH_FLOOR, FE_I, select_lines)

_FIELDS = [("element", "U8"), ("wave_A", "f8"), ("loggf", "f8"),
           ("lower_state_eV", "f8"), ("theoretical_depth", "f8"),
           ("theoretical_ew", "f8")]


def _list(rows):
    a = np.zeros(len(rows), dtype=_FIELDS)
    for i, (el, wv, gf, ep, depth) in enumerate(rows):
        a[i] = (el, wv, gf, ep, depth, 0.0)      # theoretical_ew is 0, as in the real list
    return a


def test_selects_only_neutral_iron_as_ispec_spells_it():
    ll = _list([
        (FE_I,   3100.0, -1.0, 2.0, 0.5),
        ("Fe 2", 3200.0, -1.0, 2.0, 0.5),        # ionised — not this product
        ("Co 1", 3300.0, -1.0, 2.0, 0.5),        # the most common species in the band
    ])
    got = select_lines(ll, lo_A=3000, hi_A=3780, n=10, teff=5772, min_sep_A=1.0)
    assert list(got["wave_A"]) == [3100.0]


def test_roman_numeral_spelling_is_not_silently_empty():
    """'Fe I' selects nothing; that must RAISE and name the real values, not return []."""
    ll = _list([("Fe I", 3100.0, -1.0, 2.0, 0.5)])     # roman — the wrong spelling
    with pytest.raises(SystemExit) as e:
        select_lines(ll, lo_A=3000, hi_A=3780, n=10, teff=5772, min_sep_A=1.0)
    assert "Fe 1" in str(e.value)                      # tells you what it DID find


def test_depth_window_drops_invisible_and_black_cores():
    ll = _list([
        (FE_I, 3100.0, -1.0, 2.0, DEPTH_FLOOR - 0.01),   # too shallow to see
        (FE_I, 3200.0, -1.0, 2.0, 0.50),                 # keep
        (FE_I, 3300.0, -1.0, 2.0, DEPTH_CEIL + 0.01),    # chi2 ~flat in A(X)
    ])
    got = select_lines(ll, lo_A=3000, hi_A=3780, n=10, teff=5772, min_sep_A=1.0)
    assert list(got["wave_A"]) == [3200.0]


def test_ranks_by_depth_not_by_the_all_zero_theoretical_ew():
    """The deeper line must win. If ranking used theoretical_ew (all 0.0) it could not."""
    ll = _list([(FE_I, 3100.0, -1.0, 2.0, 0.20),
                (FE_I, 3500.0, -1.0, 2.0, 0.80)])
    got = select_lines(ll, lo_A=3000, hi_A=3780, n=1, teff=5772, min_sep_A=1.0)
    assert list(got["wave_A"]) == [3500.0]


def test_minimum_separation_spreads_the_set_across_the_band():
    """Two lines 0.5 A apart are the same photons twice; only the deeper survives."""
    ll = _list([(FE_I, 3200.0, -1.0, 2.0, 0.80),
                (FE_I, 3200.5, -1.0, 2.0, 0.70),
                (FE_I, 3600.0, -1.0, 2.0, 0.60)])
    got = select_lines(ll, lo_A=3000, hi_A=3780, n=10, teff=5772, min_sep_A=4.0)
    assert list(got["wave_A"]) == [3200.0, 3600.0]


def test_band_edges_are_respected():
    ll = _list([(FE_I, 2999.0, -1.0, 2.0, 0.5),
                (FE_I, 3400.0, -1.0, 2.0, 0.5),
                (FE_I, 3781.0, -1.0, 2.0, 0.5)])
    got = select_lines(ll, lo_A=3000, hi_A=3780, n=10, teff=5772, min_sep_A=1.0)
    assert list(got["wave_A"]) == [3400.0]


def test_all_zero_depth_raises_rather_than_returning_nothing():
    ll = _list([(FE_I, 3100.0, -1.0, 2.0, 0.0),
                (FE_I, 3400.0, -1.0, 2.0, 0.0)])
    with pytest.raises(SystemExit) as e:
        select_lines(ll, lo_A=3000, hi_A=3780, n=10, teff=5772, min_sep_A=1.0)
    assert "theoretical depth" in str(e.value)
