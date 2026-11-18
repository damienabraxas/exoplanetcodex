import numpy as np
import pytest
from pipeline.absolute_rate_integrals import bb_rates,bf_rates,trapezoid_measure


def test_bb_preserves_angular_normalization_and_spontaneous_limit():
    wave=np.array([499.,500.,501.]);phi=np.ones((3,2));w=1.0001
    up,down=bb_rates(wave,phi,np.zeros_like(phi),lambda0_nm=500.,Aji=3.,Bij=4.,Bji=2.,angular_sum=w)
    np.testing.assert_array_equal(up,0)
    expected=3*np.sum(trapezoid_measure(wave)*500/wave)/2
    np.testing.assert_allclose(down,expected,rtol=1e-15)
    J=np.full_like(phi,7*w)
    up,down=bb_rates(wave,phi,J,lambda0_nm=500.,Aji=3.,Bij=4.,Bji=2.,angular_sum=w)
    np.testing.assert_allclose(up,4*7*expected/3)
    np.testing.assert_allclose(down,(2*7+3)*expected/3)


def test_bf_planck_reciprocity_at_actual_nonunit_angular_sum():
    wave=np.array([90.,100.,120.]);T=np.array([4000.,6000.]);h=6.6260755e-34;k=1.380658e-23;c=2.99792458e8;w=1.0001
    nu=c/(wave[:,None]*1e-9);J=w*2*h*nu**3/c**2/np.expm1(h*nu/(k*T))
    up,down=bf_rates(wave,np.full(3,1e-22),J,nstar_ratio=np.array([2.,3.]),T_K=T,h=h,k=k,c=c,angular_sum=w)
    np.testing.assert_allclose(down/up,[2.,3.],rtol=2e-15)
    with pytest.raises(ValueError):trapezoid_measure([100.,90.])
    with pytest.raises(ValueError):bb_rates(wave,np.zeros_like(J),J,lambda0_nm=100.,Aji=1.,Bij=1.,Bji=1.,angular_sum=w)
