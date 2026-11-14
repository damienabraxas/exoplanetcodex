import numpy as np
import pytest
from pipeline.fixed_population_thermal import thermal_correction


def test_corrected_opacity_emission_satisfy_kirchhoff_and_nlte_departures():
    wave=np.array([130.,844.,1300.]);T=np.array([4000.,6000.]);shape=(3,2);h=6.626176e-34;k=1.380662e-23;c=2.99792458e8
    lower=np.array([4.,6.]);upper=np.array([2.,3.]);ratio=lower/upper;V=np.ones(shape);nu=c/(wave[:,None]*1e-9);g=ratio[None,:]*np.exp(-h*nu/(k*T[None,:]))*1.001;bare=np.broadcast_to(2*h*nu**3/c**2,shape)*.999
    dc,de,new_g,new_bare=thermal_correction(V,g,bare,ratio,upper,wave,T,h,k,c)
    chi=lower[None,:]*V-upper[None,:]*g*V+dc;eta=upper[None,:]*g*bare*V+de
    B=new_bare/np.expm1(h*c/(wave[:,None]*1e-9*k*T[None,:]))
    np.testing.assert_allclose(eta/chi,B,rtol=1e-13)
    # Off-LTE populations preserve the prescribed departure-dependent formula.
    dc,de,ng,nb=thermal_correction(V,g,bare,ratio,upper*1.2,wave,T,h,k,c)
    np.testing.assert_allclose(lower[None,:]*V-upper[None,:]*1.2*g*V+dc,lower[None,:]*V-upper[None,:]*1.2*ng*V)
    with pytest.raises(ValueError):thermal_correction(V,g,bare,ratio,upper,wave,T*0,h,k,c)
