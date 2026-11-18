import numpy as np
import pytest
from pipeline.joint_oxygen_lte_control import source_direction_inverse_control,validate_joint_plan
from pipeline.native_oxygen_collisions import collision_matrix
from pipeline.thermal_profile_contract import collision_summary


def test_mixed_channels_preserve_original_source_direction_not_all_upward():
    records=[[str(i),kind,'0','1']for i,kind in enumerate(['CE','CI','CH','CH0'])]
    q=np.full((4,2),1e-9);T=np.array([5000.,6000.]);ne=np.array([1e10,2e10]);H=np.array([1e14,2e14]);HI=H*2
    star=np.ones((41,2));candidate=star.copy();candidate[0]*=1.2;candidate[1]*=.8
    E=np.arange(41,dtype=float);g=np.ones(41);ek=11604.
    native=collision_matrix(records,q,T,ne,H,HI,star,E,g,ek)
    result,audit=source_direction_inverse_control(records,q,T,ne,H,HI,star,candidate,E,g,ek)
    assert all(x['source_forward_exact']for x in audit.values())
    assert audit['CE']['preserved_direction']==audit['CH']['preserved_direction']=='downward'
    assert audit['CI']['preserved_direction']==audit['CH0']['preserved_direction']=='upward'
    assert np.any(result[1,0]!=native[1,0])and np.any(result[0,1]!=native[0,1])
    assert collision_summary(result,candidate)['max_absolute_signed_defect']<3e-16
    assert collision_summary(native,candidate)['max_absolute_signed_defect']>.1
    with pytest.raises(ValueError,match='total'):source_direction_inverse_control(records,q,T,ne,H,HI,star,candidate*2,E,g,ek)


def test_joint_plan_refuses_partial_thermal_or_unmatched_controls():
    plan=dict(joint_lte_control=dict(method='held-LW-finite-level-source-forward'),emission_profile_control=dict(profile='held-LW'),boundary_control='native-zero-ghost-linear',native_physical_control={},continuum_interpolation_control={})
    validate_joint_plan(plan)
    for key in list(plan):
        wrong=dict(plan);wrong.pop(key)
        with pytest.raises(ValueError):validate_joint_plan(wrong)
    with pytest.raises(ValueError):validate_joint_plan(dict(plan,fixed_population_thermal_control={}))
    with pytest.raises(ValueError):validate_joint_plan(dict(plan,emission_profile_control=dict(profile='source-native')))
