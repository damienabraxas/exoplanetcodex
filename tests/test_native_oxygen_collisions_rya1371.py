"""Discriminate source directions, reverse channel and explicit collider density roles."""
import numpy as np
import pytest
from pipeline.native_oxygen_collisions import collision_matrix

def fixture():
    return dict(records=[['0','CH','0','1'],['1','CH0','2','39']],coefficients_cm3=np.array([[2.,3.],[4.,5.]]),temperature_K=np.array([5000.,6000.]),ne_cm3=np.array([7.,8.]),H_ground_raw_cm3=np.array([11.,12.]),H_I_total_cm3=np.array([13.,14.]),nstar_cm3=np.broadcast_to(np.arange(1.,42.)[:,None],(41,2)).copy(),energy_eV=np.arange(41.),g=np.ones(41),e_k=11604.)

def test_native_directions_lte_inverse_and_total_H_control_are_distinct():
    a=fixture();C=collision_matrix(**a)
    np.testing.assert_array_equal(C[0,1],np.array([22.,36.]))
    np.testing.assert_array_equal(C[1,0],np.array([44.,72.]))
    np.testing.assert_array_equal(C[39,2],np.array([44.,60.]))
    np.testing.assert_allclose(C[2,39],np.array([44.,60.])*3/40)
    D=collision_matrix(**a,control='directed_CH0_only');np.testing.assert_array_equal(D[2,39],0);np.testing.assert_array_equal(D[39,2],C[39,2])
    H=collision_matrix(**a,control='total_H_substitution');np.testing.assert_allclose(H[0,1]/C[0,1],a['H_I_total_cm3']/a['H_ground_raw_cm3'])
    R=collision_matrix(**a,control='LW_CH_direction');np.testing.assert_array_equal(R[1,0],C[0,1]);assert np.any(R!=C)
    np.testing.assert_array_equal(C[np.arange(41),np.arange(41)],0)

def test_ne_is_electron_only_and_bad_inputs_refuse():
    a=fixture();C=collision_matrix(**a);b=dict(a,ne_cm3=a['ne_cm3']*100);np.testing.assert_array_equal(collision_matrix(**b),C)
    for bad in [dict(a,H_ground_raw_cm3=np.array([-1.,12.])),dict(a,coefficients_cm3=np.ones((1,2))),dict(a,records=[['0','CP','0','1'],['1','CH0','2','39']])]:
        with pytest.raises(ValueError):collision_matrix(**bad)
    with pytest.raises(ValueError):collision_matrix(**a,control='implicit-LTE-fallback')


def test_electron_ce_ci_rates_and_lte_ratios_use_no_hydrogen_density():
    a=fixture();a['records']=[['0','CE','0','1'],['1','CI','2','39']]
    C=collision_matrix(**a);T=a['temperature_K'];ne=a['ne_cm3'];q=a['coefficients_cm3']
    np.testing.assert_allclose(C[0,1],q[0]*ne*np.sqrt(T))
    np.testing.assert_allclose(C[1,0],C[0,1]*2)
    np.testing.assert_allclose(C[39,2],q[1]*ne*np.sqrt(T)*np.exp(-37*a['e_k']/T))
    np.testing.assert_allclose(C[2,39],C[39,2]*3/40)
    b=dict(a,H_ground_raw_cm3=a['H_ground_raw_cm3']*100,H_I_total_cm3=a['H_I_total_cm3']*100)
    np.testing.assert_array_equal(collision_matrix(**b),C)
    b=dict(a,ne_cm3=ne*2);np.testing.assert_array_equal(collision_matrix(**b),2*C)
