import numpy as np
import pytest
from pipeline.bride_cn_diagnostic import frequency_offsets,external_state_extinction
from pipeline.wavelength_util import air_to_vac


def coords(wave,mode,medium):
    return frequency_offsets(wave,np.array([7771.957]),np.array([2e9]),np.array([2.9979e10]),np.array([7771.957e-8]),input_medium=medium,mode=mode)


def test_explicit_frames_and_zero_centres():
    air=np.array([7771.,7771.957,7773.]);vac=air_to_vac(air)
    for mode,grid,medium in [('native_linear_air',air,'air'),('exact_air_frequency',air,'air'),('exact_vacuum_frequency',vac,'vacuum')]:
        v=coords(grid,mode,medium);assert v.shape==(3,1);assert v[1,0]==0.;assert v[0,0]<0<v[2,0]
    for mode,medium in [('native_linear_air','vacuum'),('exact_vacuum_frequency','air'),('unknown','air')]:
        with pytest.raises(ValueError):coords(air,mode,medium)
    with pytest.raises(ValueError):coords([1900.],'native_linear_air','air')


def test_linear_limit_and_two_dimensional_kernel_contract():
    wave=np.array([7771.956999,7771.957001]);linear=coords(wave,'native_linear_air','air');exact=coords(wave,'exact_air_frequency','air')
    np.testing.assert_allclose(linear,exact,rtol=2e-10)
    v=np.array([-1.,0.,1.]);fields=([0.,0.],[1.,1.],[2.,2.],[1.,1.],[1.,1.],[3.,3.])
    np.testing.assert_array_equal(external_state_extinction(v,*fields),external_state_extinction(np.column_stack([v,v]),*fields))
    with pytest.raises(ValueError):external_state_extinction(np.ones((2,3)),*fields)
