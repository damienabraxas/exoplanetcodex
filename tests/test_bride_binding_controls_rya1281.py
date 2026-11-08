import numpy as np
import pytest
from pipeline.bride_binding_controls import control_state


def fixture():
    records=np.array([(1,5000.,1.,1e5,1e20,1.01e20),(2,6000.,2.,1e5,2e20,1.98e20)],dtype=[(x,float)for x in ['depth_index','temperature_K','electron_pressure_dyn_cm2','xi_cm_s','MARCS_LW_nH_m3','nH_network_m3']])
    return dict(ndep=2,height_m=np.array([100.,0.])),records


def test_only_named_H_vector_changes_and_no_source_aliasing():
    s,r=fixture();a=control_state(s,r,'model_mass_H',kBoltzmann_SI=1.38e-23);b=control_state(s,r,'native_network_H',kBoltzmann_SI=1.38e-23)
    assert a.ne_m3==b.ne_m3 and a.temperature_K==b.temperature_K and a.height_m==b.height_m
    np.testing.assert_array_equal(a.nH_m3,r['MARCS_LW_nH_m3']);np.testing.assert_array_equal(b.nH_m3,r['nH_network_m3'])
    assert a.electron_pressure_Pa==(0.1,0.2)and a.vturb_m_s==(1000.,1000.)
    r['temperature_K'][0]=9000;s['height_m'][0]=200
    assert a.temperature_K[0]==5000 and a.height_m[0]==100


def test_invalid_choice_depth_and_nonfinite_state_refused():
    s,r=fixture()
    with pytest.raises(ValueError):control_state(s,r,'preferred',kBoltzmann_SI=1.38e-23)
    r['depth_index'][0]=2
    with pytest.raises(ValueError):control_state(s,r,'model_mass_H',kBoltzmann_SI=1.38e-23)
    s,r=fixture();r['temperature_K'][0]=np.nan
    with pytest.raises(ValueError):control_state(s,r,'native_network_H',kBoltzmann_SI=1.38e-23)
