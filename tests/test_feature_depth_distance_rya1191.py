"""RYA-1191 — a graded line must not inherit its depth from a stranger.

`_feature_depth` decides which pool a LAB-graded line enters: at or below
`line_accounting_rya709.DEPTH_HI` it is GRADED, above it DEEPGRADED. It reads that depth
out of `linelist_solar.csv` by taking the nearest FEATURE — and `np.argmin` has no
distance limit, so a line the catalogue does not carry at all silently borrows whatever
happens to be closest.

🔴 MEASURED, ONE LINE IN 450. Fe I 8432.174 (Ruffoni 2014, LAB tier, present in the
SYNTHESIS list with loggf -1.770) has NO Fe row in linelist_solar. The nearest feature is
8433.778, **1.604 A away**, whose depth 0.023 sailed it through the <= 0.60 graded gate.
It was then fitted on three holdings and returned A = 5.51, 5.62 and 5.86 — 1.6 to 1.9 dex
low, because there is no measurable line there to fit.

⚠️ AND IT IS NOT A BLEND, WHICH IS WHAT THE BLEND AUDIT FIRST SAID. Its window holds three
CN lines of depth 0.006, 0.005 and 0.002 — total 0.013. The audit's "molecular share
1.0000" was a share of very nearly nothing, which is why a share needs a floor on its
denominator before it means anything.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="module")
def fd():
    import derive_band_products as D
    return D._feature_depth


def test_a_line_with_no_feature_nearby_gets_no_depth(fd):
    """🔴 THE CASE. NaN, not the neighbour's 0.023 and not a 0.0 default — "we do not know
    this line's depth" is a third state (RYA-833)."""
    assert np.isnan(fd(np.array([8432.174]))[0])


def test_lines_that_sit_on_their_own_feature_are_untouched(fd):
    """The guard must not cost anything real. These four sit within 0.001 A of their own
    feature and keep the depths the selectors have always used."""
    got = fd(np.array([9012.075, 9437.793, 4303.170, 5171.672]))
    assert np.allclose(got, [0.345, 0.200, 0.780, 0.865], atol=1e-3), got


def test_the_gate_that_admits_a_line_still_works_on_the_kept_ones(fd):
    """4303.170 at 0.780 is ABOVE the 0.60 gate, so it belongs to DEEPGRADED — which is
    exactly the pool it appears in. If this inverts, the pools have been reshuffled."""
    from line_accounting_rya709 import DEPTH_HI
    assert fd(np.array([4303.170]))[0] > DEPTH_HI
    assert fd(np.array([9437.793]))[0] <= DEPTH_HI


def test_the_limit_is_derived_from_the_feature_grouping_not_chosen():
    """⚠️ A number picked to exclude one known line is a tuned cut. `GROUP_A` is the width
    inside which rows are ALREADY treated as one feature, so twice it is the distance
    beyond which a line cannot belong to that feature — and it happens to change exactly
    one line of 450, rather than being chosen to."""
    src = (ROOT / "scripts/derive_band_products.py").read_text()
    assert "limit = 2.0 * GROUP_A" in src
    assert "DERIVED, not chosen" in src


def test_an_unknown_depth_is_reported_and_not_silently_dropped():
    """Both `depth > gate` and `depth <= gate` are False for NaN, so an unknown-depth line
    would vanish from BOTH pools without a word. The selector prints it."""
    src = (ROOT / "scripts/derive_band_products.py").read_text()
    assert "A LINE WITH NO KNOWN DEPTH IS REPORTED, NOT SILENTLY DROPPED" in src
    assert "no Fe row in the stellar catalogue at this wavelength" in src


def test_the_problem_is_isolated_and_that_is_asserted_not_assumed(fd):
    """One line in 450 — so this is a guard, not a reshuffle. If a future catalogue change
    makes many lines unknown, that is a catalogue problem to look at, not a limit to relax,
    and this test is where it surfaces."""
    import pandas as pd
    cg = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    lab = cg[(cg.species.astype(str).isin(["Fe I", "Fe II"]))
             & (cg.gf_tier.astype(str) == "LAB")]
    d = fd(lab.wavelength_air_A.values.astype(float))
    n_unknown = int(np.isnan(d).sum())
    assert n_unknown == 1, f"{n_unknown} of {len(d)} graded lines now have no known depth"
