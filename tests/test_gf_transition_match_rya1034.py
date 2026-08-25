"""A tier is earned by a verified transition, never by wavelength proximity — RYA-1034.

The gf grader matched a laboratory entry to a canonical line on wavelength AND lower-state
energy, then took `m.iloc[argmin]` -- the NEAREST of however many candidates fell inside
the tolerance. Both the laboratory tables and the canonical list carry `lower_level` /
`upper_level`, and the promotion path wrote them into its audit row without ever comparing
them. So the one piece of evidence that identifies a transition was collected and
discarded.

Two properties are pinned here:

  1. An AMBIGUOUS match raises instead of resolving by proximity. Two candidates inside
     the tolerance means the tolerance cannot separate them; breaking that tie by rounding
     is how a blend acquires a laboratory sigma it never earned (RYA-161).
  2. `transition_matches` returns None -- NOT True -- when either side is silent about its
     levels. "Could not ask the question" must never read as "confirmed".

Property 2 is the one worth guarding hardest. A boolean helper that returns True on
missing data turns an absent check into a passing one, which is exactly the RYA-833 shape.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pipeline.gf_grades import AmbiguousLineMatch, WAVE_TOL_A, _nearest, transition_matches


def test_silence_is_not_confirmation():
    """🔴 The rule that keeps an unasked question from reading as a passed one."""
    row = {"lower_level": "4s ^2S_1/2", "upper_level": "5p ^2P_3/2"}
    assert transition_matches(row, None, "5p ^2P_3/2") is None
    assert transition_matches(row, "4s ^2S_1/2", None) is None
    assert transition_matches({}, "4s ^2S_1/2", "5p ^2P_3/2") is None
    assert transition_matches({"lower_level": "nan", "upper_level": "nan"},
                              "4s ^2S_1/2", "5p ^2P_3/2") is None


def test_same_transition_confirms_across_notation_dialects():
    """Burheim ships LaTeX; NIST and VALD do not. The same transition must confirm."""
    tex = {"lower_level": "4s $^2$S$_{1/2}$", "upper_level": "5p $^2$P$_{3/2}$"}
    assert transition_matches(tex, "4s ^2S_1/2", "5p ^2P_3/2") is True


def test_a_different_upper_level_does_not_confirm():
    """The Al 11253/11254 pair differs ONLY in the J of both levels -- 3d 2D3/2 - 4f 2F5/2
    against 3d 2D5/2 - 4f 2F7/2, 1.7 A apart with lower-state energies 0.0002 eV apart.
    Wavelength and ep cannot separate them; the term symbols can."""
    row = {"lower_level": "3d $^2$D$_{3/2}$", "upper_level": "4f $^2$F$_{5/2}$"}
    assert transition_matches(row, "3d ^2D_5/2", "4f ^2F_7/2") is False
    assert transition_matches(row, "3d ^2D_3/2", "4f ^2F_5/2") is True


def test_ambiguous_match_raises_instead_of_picking_the_nearest():
    """🔴 The regression: `m.iloc[argmin]` silently chose among candidates."""
    df = pd.DataFrame({"wavelength_air_A": [5000.000, 5000.010],
                       "elo_eV": [3.0, 3.0],
                       "log_gf": [-1.0, -2.0]})
    with pytest.raises(AmbiguousLineMatch) as e:
        _nearest(df, 5000.004, 3.0, "wavelength_air_A", "elo_eV")
    msg = str(e.value)
    assert "5000.0000" in msg and "5000.0100" in msg, (
        f"the refusal must NAME the candidates it will not choose between: {msg}")


def test_an_unambiguous_match_still_resolves():
    """Guard the guard in the other direction: the check must not refuse good matches.
    A grade check tested only on the case that motivated it is how I broke the
    red-optical cell with the first version of the tier-vs-gf guard."""
    df = pd.DataFrame({"wavelength_air_A": [5000.000, 5000.500],
                       "elo_eV": [3.0, 3.0],
                       "log_gf": [-1.0, -2.0]})
    got = _nearest(df, 5000.001, 3.0, "wavelength_air_A", "elo_eV")
    assert got is not None and np.isclose(got["log_gf"], -1.0)


def test_no_match_is_still_none():
    df = pd.DataFrame({"wavelength_air_A": [6000.0], "elo_eV": [3.0], "log_gf": [-1.0]})
    assert _nearest(df, 5000.0, 3.0, "wavelength_air_A", "elo_eV") is None
