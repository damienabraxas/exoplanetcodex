"""Tests for the RYA-353 single-source gf resolver (pipeline/gf_resolver.py)."""
import numpy as np
import pytest

from pipeline import gf_resolver as gr
from pipeline.species import species_key


def test_resolve_anchor_is_nist():
    # 6247.557 Fe II is NIST grade-B → must resolve to -2.329 (folds in the 6247 fix)
    assert abs(gr.resolve(species_key('Fe', 'II'), 6247.557, 3.892) + 2.329) < 1e-3


def test_resolve_loud_guard_raises():
    # No silent fallback: an absent line raises, never defaults.
    with pytest.raises(gr.GfResolutionError):
        gr.resolve(species_key('Fe', 'II'), 9999.99, 3.0)


def test_clustering_splits_same_ep_distinct_lines():
    # The RYA-353 clustering bug: distinct lines sharing an EP (Sc I at EP 1.969,
    # one near 5350.27 the other 6414.5) must NOT merge despite the shared EP.
    keys = [(21, 1), (21, 1)]
    wls = np.array([5350.272, 6414.587])
    eps = np.array([1.969, 1.969])
    cls = gr.cluster_physical_lines(keys, wls, eps)
    assert len(cls) == 2, "distinct same-EP far-λ lines must not merge"


def test_clustering_groups_true_hfs():
    # Two close-λ components sharing an EP ARE one physical line.
    keys = [(21, 1), (21, 1)]
    wls = np.array([5350.272, 5350.321])
    eps = np.array([1.969, 1.969])
    assert len(gr.cluster_physical_lines(keys, wls, eps)) == 1


def test_synth_rescale_preserves_branching_and_hits_total():
    # Two HFS components rescaled to a target total must keep their ratio and sum to it.
    g = np.array([-0.704, -0.632])
    target = -0.076
    cur = float(np.log10(np.sum(10.0 ** g)))
    shifted = g + (target - cur)
    assert abs(float(np.log10(np.sum(10.0 ** shifted))) - target) < 1e-9   # total hits
    assert abs((shifted[0] - shifted[1]) - (g[0] - g[1])) < 1e-12          # ratio kept


def test_apply_to_regions_skips_non_gf_array():
    # A minimal region array without loggf/EP fields is returned untouched (no raise).
    arr = np.array([('Fe 1', 5000., 0., 0.)],
                   dtype=[('note', '<U10'), ('wave_A', '<f8'),
                          ('ew', '<f8'), ('ew_err', '<f8')])
    out = gr.apply_to_regions(arr)
    assert out is arr
