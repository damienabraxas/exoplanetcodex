import numpy as np
import pytest
from pipeline.bride_column_checks import integrated_column,candidate_column


def test_linear_density_integral_exact_in_si_with_outward_height():
    z=np.array([10.,8.,3.,0.]);rho=2.+(10.-z)*3.
    np.testing.assert_allclose(integrated_column(z,rho,7.),7+2*(10-z)+1.5*(10-z)**2,rtol=1e-15)


def test_column_node_identity_and_convex_monotonicity():
    a=np.array([[1.,3.,8.],[2.,6.,16.]])
    for scale in ('linear','log'):
        np.testing.assert_allclose(candidate_column(a,[1.,0.],scale=scale),a[0],rtol=1e-15)
        result=candidate_column(a,[.5,.5],scale=scale)
        assert np.all(np.diff(result)>0)
        assert np.all((result>=a[0])&(result<=a[1]))
    np.testing.assert_allclose(candidate_column(a,[.5,.5],scale='log'),a[0]*np.sqrt(2),rtol=1e-15)


@pytest.mark.parametrize('weights',[[-.1,1.1],[.1,.1],[np.nan,1.]])
def test_no_extrapolation_or_invalid_weights(weights):
    with pytest.raises(ValueError):candidate_column([[1,2],[2,4]],weights,scale='linear')


def test_geometry_and_density_validation():
    for z,rho,top in [([0,1],[1,1],1),([1,0],[-1,1],1),([1,0],[1,1],0)]:
        with pytest.raises(ValueError):integrated_column(z,rho,top)
