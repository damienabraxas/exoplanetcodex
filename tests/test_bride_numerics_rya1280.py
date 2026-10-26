"""Analytic slab controls with source-derived coefficient truncation bounds.

Initial 4096-eps-only experiment failed four cases; retain its receipt in the
validation report. It omitted documented Taylor truncation, so the subsequent
bound includes independently derived series remainders and conditioning, rather
than choosing a constant from observed errors. No abundance target is involved.
"""
import numpy as np
import pytest
pytest.importorskip('lightweaver')
from pipeline.bride_numerical_validation import constant_slab, homogeneous_error_bound

TOL = 4096*np.finfo(np.float64).eps

@pytest.mark.parametrize('solver',['piecewise_linear_1d','piecewise_besser_1d','piecewise_bezier3_1d'])
@pytest.mark.parametrize('tau',[1e-10, .01, 1., 20.])
@pytest.mark.parametrize('mu',[.2,1.])
@pytest.mark.parametrize('depths',[9,65])
def test_formal_solution_matches_independent_slab_and_thin_limit(solver,tau,mu,depths):
    actual,expected=constant_slab(tau_vertical=tau,source=4.,lower_intensity=2.,
        upper_intensity=3.,mu=mu,depths=depths,solver=solver)
    bound=homogeneous_error_bound(tau_vertical=tau,mu=mu,depths=depths,
        source=4.,lower_intensity=2.,solver=solver)
    np.testing.assert_allclose(actual,expected,rtol=0,atol=bound)
    if tau==1e-10:
        np.testing.assert_allclose(actual,2.+2.*tau/mu,rtol=TOL,atol=TOL)


def test_opposite_boundary_cannot_contaminate_absorption_only_outgoing_ray():
    kw=dict(tau_vertical=.5,source=4.,lower_intensity=2.,mu=.7,depths=17,
            solver='piecewise_besser_1d')
    a,e=constant_slab(**kw,upper_intensity=3.)
    b,_=constant_slab(**kw,upper_intensity=300.)
    np.testing.assert_array_equal(a,b)
    wrong,_=constant_slab(**dict(kw,lower_intensity=3.),upper_intensity=2.)
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(wrong,e,rtol=TOL,atol=TOL)


def test_exact_vacuum_is_not_silently_relabelled_a_successful_backend_control():
    with pytest.raises(ValueError,match='positive'):
        constant_slab(tau_vertical=0,source=4.,lower_intensity=2.,upper_intensity=3.,
            mu=1.,depths=9,solver='piecewise_linear_1d')


@pytest.mark.parametrize('temperature',[4500.,6000.,9000.])
@pytest.mark.parametrize('initial_weights',[[1.,1.,1.],[100.,1.,.01],[.01,1.,100.]])
@pytest.mark.parametrize('rate_scale',[.01,1.,100.])
def test_actual_collision_operators_and_se_recover_independent_lte(temperature,initial_weights,rate_scale):
    from pipeline.bride_numerical_validation import collision_lte_control
    d=collision_lte_control(temperature=temperature,initial_weights=initial_weights,
                            rate_scale=rate_scale)
    # Forward error depends on conditioning of the conservation-augmented
    # linear system; 9000 K is deliberately an ill-conditioned ion limit.
    # This is the standard eps*kappa arithmetic bound, not a fitted tolerance.
    bound=64*np.finfo(float).eps*d['condition_inf']*np.max(d['expected_fraction'])
    np.testing.assert_allclose(d['actual_fraction'],np.broadcast_to(
        d['expected_fraction'][:,None],d['actual_fraction'].shape),rtol=0,atol=bound)
    assert np.max(np.abs(d['conservation'])) < TOL
    assert d['normwise_backward_error'] < TOL
    assert np.min(d['actual_fraction']) > 0


def test_directional_rate_defect_is_detected_by_independent_lte_control(monkeypatch):
    from lightweaver.collisional_rates import CE
    from pipeline.bride_numerical_validation import collision_lte_control
    original=CE.compute_rates
    def wrong(self,atmos,eqpops,matrix):
        original(self,atmos,eqpops,matrix)
        matrix[self.i,self.j,:] *= 10.
    monkeypatch.setattr(CE,'compute_rates',wrong)
    d=collision_lte_control(temperature=6000.,initial_weights=[1.,1.,1.])
    # The wrong rates still conserve particles and solve to tiny residuals;
    # only an independent LTE population reference discriminates the error.
    assert np.max(np.abs(d['conservation'])) < TOL
    assert np.max(d['normalized_rate_residual']) < TOL
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(d['actual_fraction'],np.broadcast_to(
            d['expected_fraction'][:,None],d['actual_fraction'].shape),rtol=1e-11,atol=0)


@pytest.mark.parametrize('initial_weights',[[1.,1.,1.],[100.,1.,.01],[.01,1.,100.]])
@pytest.mark.parametrize('lambda_iterate,ng',[(False,False),(True,False),(True,True)])
def test_coupled_line_transfer_and_se_recover_thermal_equilibrium(initial_weights,lambda_iterate,ng):
    from pipeline.bride_numerical_validation import collision_lte_control
    # A 1e-9 population accuracy target is fixed for this iterative control;
    # iteration change/J thresholds are stricter (1e-11), not abundance targets.
    d=collision_lte_control(temperature=6000.,initial_weights=initial_weights,radiative=True,
                            lambda_iterate=lambda_iterate,ng=ng)
    np.testing.assert_allclose(d['actual_fraction'],np.broadcast_to(
        d['expected_fraction'][:,None],d['actual_fraction'].shape),rtol=1e-9,atol=0)
    assert np.max(d['normalized_rate_residual']) < 1e-9
    assert np.max(np.abs(d['conservation'])) < TOL
    assert d['history'][-1]['population_change'] < 1e-11
