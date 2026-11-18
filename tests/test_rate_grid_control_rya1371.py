import numpy as np
import pytest
from pipeline.rate_grid_control import midpoint_grid


def test_midpoint_refinement_preserves_endpoints_and_original_nodes():
    coarse=np.array([90.,130.2,844.8568,6876.57]);fine=midpoint_grid(coarse,2)
    np.testing.assert_array_equal(fine[::4],coarse)
    assert len(fine)==13 and fine[0]==coarse[0]and fine[-1]==coarse[-1]


@pytest.mark.parametrize('wave,passes',[([2.,1.],1),([1.,1.],1),([1.,float('nan')],1),([1.,2.],True),([1.,2.],3)])
def test_invalid_or_unbounded_controls_refuse(wave,passes):
    with pytest.raises(ValueError):midpoint_grid(wave,passes)


def test_local_control_refines_only_gaps_inside_declared_interval():
    coarse=np.array([1.,2.,3.,4.]);fine=midpoint_grid(coarse,2,interval_vac_nm=[2.,3.])
    np.testing.assert_array_equal(fine,[1.,2.,2.25,2.5,2.75,3.,4.])
    with pytest.raises(ValueError,match='20000'):midpoint_grid(np.linspace(1.,2.,6795),2)
