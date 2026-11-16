import numpy as np
import pytest
from pipeline.radiation_probe_grid import gaussian_probe_grid


def test_probe_preserves_every_node_and_ranks_only_positive_indicators():
    wave=np.array([1.,1.5,2.,2.5,3.]);scores=np.array([0.,9.,2.,1.])
    expanded,points,weights,panels=gaussian_probe_grid(wave,scores,8)
    np.testing.assert_array_equal(panels,[1,2]);np.testing.assert_array_equal(expanded[np.searchsorted(expanded,wave)],wave)
    assert len(expanded)==len(wave)+8 and np.all(weights>0)
    assert np.all((points[:4]>1.5)&(points[:4]<2.))and np.all((points[4:]>2.)&(points[4:]<2.5))
    np.testing.assert_allclose(weights.sum(),1.,rtol=2e-16)
    # Shared score ties resolve by wavelength order, without a reference target.
    np.testing.assert_array_equal(gaussian_probe_grid(wave,np.ones(4),4)[3],[0])


def test_probe_refuses_cap_partial_groups_missing_indicators_and_bad_grid():
    with pytest.raises(ValueError,match='20000'):gaussian_probe_grid(np.arange(1.,19971.),np.ones(19969),32)
    for scores,budget in [([1.,1.],3),([0.,0.],4),([1.,np.nan],4),([-1.,1.],4)]:
        with pytest.raises(ValueError):gaussian_probe_grid([1.,2.,3.],scores,budget)
    with pytest.raises(ValueError):gaussian_probe_grid([1.,1.,2.],[1.,1.],4)
