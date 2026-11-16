import numpy as np
import pytest
from pipeline.thermal_profile_contract import background_balance,collision_summary


def test_scattering_cannot_mask_absorption_defect():
    result=background_balance([[1e6]],[[2.]],[[1e6-1]],[[1.]])
    assert result['thermal_absorption']['max_absolute_signed_defect']==pytest.approx(1/3)
    assert result['total_at_J_equals_B']['max_absolute_signed_defect']<1e-6


def test_scattering_counted_once_and_absent_absorption_not_pass():
    result=background_balance([[2.,3.]],[[0.,2.]],[[2.,2.]],[[1.,2.]])
    assert result['thermal_absorption']['present_cells']==1
    assert result['total_at_J_equals_B']['max_absolute_signed_defect']==0
    assert result['pure_scattering_cells']==1
    empty=background_balance([[0.]],[[0.]],[[0.]],[[1.]])
    assert empty['thermal_absorption']['max_absolute_signed_defect'] is None


@pytest.mark.parametrize('chi,eta,sca,B',[([[1.]],[[1.]],[[2.]],[[1.]]),([[1.]],[[float('nan')]],[[0.]],[[1.]]),([[1.]],[[1.]],[[0.]],[[0.]])])
def test_invalid_background_refuses(chi,eta,sca,B):
    with pytest.raises(ValueError):background_balance(chi,eta,sca,B)


def test_collision_one_way_and_disconnected_are_explicit():
    C=np.zeros((3,3,1));C[1,0,0]=1
    result=collision_summary(C,np.ones((3,1)))
    assert result==dict(pair_depths=1,one_way_pair_depths=1,max_absolute_signed_defect=1.,levels_without_collision_edges_1based=[3])
