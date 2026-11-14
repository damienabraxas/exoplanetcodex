import numpy as np
import pytest
from pipeline.adaptive_continuum_grid import adaptive_grid,panel_indicators
from pipeline.adaptive_continuum_grid import nested_parabolic_measure


def test_adaptive_indicator_allocates_curved_region_and_preserves_all_nodes():
    base=np.array([1.,2.,3.]);fine=np.array([1.,1.5,2.,2.5,3.]);F=np.array([[1.],[1.],[1.],[4.],[1.]])
    score=panel_indicators(base,fine,F);assert score[0]==0 and score[1]>0
    out,selected=adaptive_grid(base,fine,score,2)
    np.testing.assert_array_equal(selected,[1]);np.testing.assert_array_equal(out,[1.,1.5,2.,2.25,2.5,2.75,3.])
    # Linear positive integrands have no trapezoid error, at any depth scale.
    np.testing.assert_array_equal(panel_indicators(base,fine,fine[:,None]*[1.,17.]),0)


def test_adaptive_refuses_missing_nodes_invalid_fields_and_size_bypass():
    with pytest.raises(ValueError):panel_indicators([1.,2.],[1.,1.6,2.1],np.ones((3,1)))
    with pytest.raises(ValueError):panel_indicators([1.,2.],[1.,1.5,2.],-np.ones((3,1)))
    with pytest.raises(ValueError,match='20000'):adaptive_grid([1.,2.],[1.,1.5,2.],[1.],20000)


def test_nested_parabolic_measure_integrates_cubic_without_negative_weights():
    base=np.array([1.,2.,4.]);fine=np.array([1.,1.25,1.5,1.75,2.,3.,4.])
    w=nested_parabolic_measure(base,fine)
    assert np.all(w>0)
    for power in range(4):np.testing.assert_allclose(w@fine**power,(4**(power+1)-1)/(power+1),rtol=1e-15)
