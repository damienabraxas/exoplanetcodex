import numpy as np
import pytest
from pipeline.bride_cn_diagnostic import external_state_extinction


def evaluate(**override):
    args=dict(offset_doppler=np.array([-1.,0.,1.]),damping_a=np.array([0.]),doppler_width_Hz=np.array([2.]),lower_population_m3=np.array([3.]),gf=np.array([1.]),stimulated_emission=np.array([1.]),integrated_strength_m2_Hz=np.array([5.]))
    args.update(override);return external_state_extinction(**args)


def test_gaussian_limit_symmetry_and_zero_population():
    chi=evaluate();assert chi.shape==(3,1);np.testing.assert_array_equal(chi,chi[::-1])
    assert chi[0,0]/chi[1,0]==pytest.approx(np.exp(-1.))
    np.testing.assert_array_equal(evaluate(lower_population_m3=np.array([0.])),np.zeros((3,1)))


def test_frequency_integral_and_population_scaling():
    v=np.linspace(-12,12,20001);chi=evaluate(offset_doppler=v)[:,0]
    assert np.trapezoid(chi,v)*2==pytest.approx(15.,rel=1e-12)
    np.testing.assert_allclose(evaluate(lower_population_m3=np.array([6.])),2*evaluate())


@pytest.mark.parametrize('override',[dict(damping_a=[-1.]),dict(doppler_width_Hz=[0.]),dict(stimulated_emission=[1.1]),dict(lower_population_m3=[float('nan')]),dict(gf=[1.,2.]),dict(offset_doppler=[[0.,1.]])])
def test_invalid_external_state_refused(override):
    with pytest.raises(ValueError):evaluate(**override)
