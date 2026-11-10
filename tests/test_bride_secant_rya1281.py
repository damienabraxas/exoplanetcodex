import numpy as np
from pipeline.bride_secant import InverseSecant,log_coordinates,populations_from_coordinates


def test_coordinates_recover_populations_with_particle_constraint_and_no_floor():
    n=np.array([[.6,.1],[.3,.2],[.1,.7]])*np.array([5.,9.])
    x=log_coordinates(n);actual=populations_from_coordinates(x,n.sum(axis=0))
    np.testing.assert_allclose(actual,n,rtol=1e-14)
    np.testing.assert_allclose(actual.sum(axis=0),[5.,9.],rtol=1e-14)
    assert np.all(actual>0)


def test_inverse_update_satisfies_independent_secant_equation_for_coupled_matrix():
    matrix=np.array([[.002,.001],[.001,1.]])
    h=InverseSecant(memory=4);s=np.array([.3,.2]);y=matrix@s
    assert h.update(s,y)['updated']
    np.testing.assert_allclose(h.apply(y),s,rtol=1e-12,atol=1e-12)
    # Also check the transpose action used by the update against the explicit H.
    explicit=np.eye(2)+np.outer(h.u[0],h.v[0])
    np.testing.assert_allclose(h.apply(s,transpose=True),explicit.T@s,rtol=1e-14)
    result=h.update(s,np.zeros(2));assert not result['updated']and result['restarted']
    np.testing.assert_array_equal(h.apply(s),s)
