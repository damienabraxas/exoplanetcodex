"""The zero-point cap, and the property that makes it worth splitting — RYA-987.

`zero_point_term` charged sigma/sqrt(n) over the anchor. While the anchor was 7 lines that
was also the whole story. RYA-984 grew it to 163, the statistical term fell ~4x, and a
second term that had been hiding under it became the larger one.

These tests pin the RELATIONSHIPS, not the run's numbers (RYA-870). The measured cap is a
snapshot of a growing anchor — the HARPS deep leg is owed, RYA-977 adds ~65 lines — so a
test asserting 0.0343 would fail the moment the project made progress, and would be
pinning the snapshot rather than the physics.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from pipeline.error_budget import zero_point_cap


def _anchor(n_per_source: dict[str, int], sigma: dict[str, float], scatter=0.20, seed=0):
    rng = np.random.default_rng(seed)
    a, s, g = [], [], []
    for src, n in n_per_source.items():
        a.extend(rng.normal(7.45, scatter, n))
        s.extend([src] * n)
        g.extend([sigma[src]] * n)
    return a, s, g


def test_statistical_shrinks_as_one_over_root_n():
    """The half that MORE LINES buys."""
    small = zero_point_cap(*_anchor({"L": 25}, {"L": 0.04}))
    large = zero_point_cap(*_anchor({"L": 400}, {"L": 0.04}))
    ratio = small.statistical_dex / large.statistical_dex
    # 16x the lines -> ~4x the precision, allowing for the sampled scatter
    assert 3.0 < ratio < 5.5, ratio


def test_systematic_does_NOT_shrink_when_lines_come_from_the_same_source():
    """🔴 THE WHOLE POINT. Adding lines from a laboratory cannot beat that laboratory's own
    normalisation — every line it reports carries the same shift. Treating the cited sigmas
    as independent would let this average down as 1/sqrt(n), which is the manufactured
    precision RYA-968 §3.2 exists to prevent."""
    few = zero_point_cap(*_anchor({"L": 20}, {"L": 0.04}))
    many = zero_point_cap(*_anchor({"L": 2000}, {"L": 0.04}))
    assert few.systematic_dex == pytest.approx(many.systematic_dex, abs=1e-9)
    assert many.systematic_dex == pytest.approx(0.04, abs=1e-9)


def test_spreading_the_anchor_across_INDEPENDENT_sources_does_reduce_it():
    """The converse, and the reason the term is source-keyed rather than a flat constant:
    independent laboratories DO average against each other."""
    one = zero_point_cap(*_anchor({"A": 300}, {"A": 0.06}))
    three = zero_point_cap(*_anchor({"A": 100, "B": 100, "C": 100},
                                    {"A": 0.06, "B": 0.06, "C": 0.06}))
    assert three.systematic_dex < one.systematic_dex
    # equal thirds of equal sigma -> sqrt(3) * (sigma/3)
    assert three.systematic_dex == pytest.approx(0.06 / math.sqrt(3), abs=1e-9)


def test_the_limited_by_verdict_and_its_consequence():
    stat_limited = zero_point_cap(*_anchor({"L": 9}, {"L": 0.01}, scatter=0.30))
    assert stat_limited.limited_by == "statistical"
    assert stat_limited.lines_to_halve() is not None, (
        "when the statistical term dominates, more lines DO help and the function must "
        "say how many")

    sys_limited = zero_point_cap(*_anchor({"L": 400}, {"L": 0.06}, scatter=0.10))
    assert sys_limited.limited_by == "systematic"
    assert sys_limited.lines_to_halve() is None, (
        "once the systematic exceeds half the cap, NO finite anchor halves it — the "
        "honest way to say 'more lines will not help'")


def test_combined_is_quadrature_and_never_below_either_term():
    c = zero_point_cap(*_anchor({"A": 50, "B": 50}, {"A": 0.03, "B": 0.05}))
    assert c.combined_dex == pytest.approx(math.hypot(c.statistical_dex, c.systematic_dex))
    assert c.combined_dex >= c.statistical_dex
    assert c.combined_dex >= c.systematic_dex


def test_no_reference_abundance_can_enter_the_computation():
    """🔒 RYA-161. The firewall is enforced by the SIGNATURE, not by a convention: there is
    no parameter to pass a published solar abundance into. A cap that could be computed
    against gold could be asked 'how close are we', and gold carries its own hidden zero
    point — answering that question would launder one unknown into another."""
    import inspect
    params = set(inspect.signature(zero_point_cap).parameters)
    assert params == {"abundances", "sources", "cited_sigmas"}
    for banned in ("reference", "gold", "expected", "target", "truth", "literature"):
        assert not any(banned in p for p in params)


def test_a_missing_cited_sigma_is_refused_not_treated_as_zero():
    """A blank sigma silently read as 0.0 would SHRINK the systematic — the RYA-873 defect
    shape, where a zero nobody measured was printed beside the word MEASURED."""
    a, s, g = _anchor({"L": 10}, {"L": 0.04})
    g[3] = float("nan")
    with pytest.raises(ValueError, match="cited laboratory sigma"):
        zero_point_cap(a, s, g)


def test_an_anchor_of_one_line_is_refused():
    with pytest.raises(ValueError, match="at least two"):
        zero_point_cap([7.45], ["L"], [0.04])


def test_the_registered_anchor_matches_its_parts():
    """The registry names products, and `load` must refuse a partial anchor rather than
    quietly reporting a floor better than the one that was measured."""
    from pipeline.anchor_pools import ANCHORS
    pool = ANCHORS["rya984_graded_163"]
    assert len(pool.parts) == 2, "163 = RYA-967's 55 shallow + RYA-984's 108 deep"
    assert any("DEEPGRADED" in p for p in pool.parts), (
        "the deep set must be named by its own selector tag — RYA-984 found the two "
        "products collided on one filename before the selector entered the stem")
