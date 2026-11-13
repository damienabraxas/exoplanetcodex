"""Discrete native ghost formula vs an independent optical-depth formal integral."""
import numpy as np
import pytest
pytest.importorskip('lightweaver')
from pipeline.native_ghost_boundary import linear_ghost_intensity

def test_linear_ghost_matches_formal_integral_in_thin_and_thick_limits():
    from scipy.integrate import quad
    tau=np.array([0.,1e-8,1e-4,.1,1.,10.,1e4]);result=linear_ghost_intensity(tau,np.where(tau>0,2*tau,0),2.,np.array([1.]))[:,0]
    expected=[]
    for t in tau:
        # Integrate measured-backwards distance from the node: S=2*(1-u/t).
        expected.append(0 if t==0 else quad(lambda x:2*(1-x)*np.exp(-t*x)*t,0,1,epsabs=1e-13)[0]if t<100 else 2*(1-1/t))
    np.testing.assert_allclose(result,expected,rtol=2e-12,atol=1e-15)
    assert result[-1]<2 and result[-1]>1.99

def test_invalid_ghost_input_refuses_without_clipping():
    for chi,eta,step,mu in [([-1.],[1.],1.,[1.]),([0.],[1.],1.,[1.]),([1.],[1.],0.,[1.]),([1.],[1.],1.,[0.])]:
        with pytest.raises(ValueError):linear_ghost_intensity(chi,eta,step,mu)
