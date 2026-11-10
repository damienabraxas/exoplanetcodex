import numpy as np
import pytest


def test_lte_voigt_integrated_strength_kirchhoff_and_zero_opacity():
    pytest.importorskip('lightweaver')
    from scipy.integrate import quad
    import lightweaver.constants as c
    from pipeline.bride_lte_blend import lte_line_coefficients
    T=np.array([6000.]);ion=np.array([1e15]);xi=np.array([1000.]);nu0=c.CLight/630e-9
    kwargs=dict(lambda0_nm=630.,mass_amu=58.69,excitation_eV=0.,loggf=0.,gamma_s_inv=1e8)
    def integrand(x):
        nu=nu0+x*1e9
        chi,eta=lte_line_coefficients(np.array([c.CLight/nu/1e-9]),T,ion,xi,**kwargs)
        # Remove frequency-dependent stimulated-emission term for the sum rule.
        return chi.item()/(-np.expm1(-c.HPlanck*nu/(c.KBoltzmann*T.item())))*1e9
    integral=quad(integrand,-1e4,1e4,epsabs=0,epsrel=1e-7,points=[-10,0,10],limit=300)[0]
    expected=c.QElectron**2/(4*c.Epsilon0*c.MElectron*c.CLight)*ion.item()
    # Finite +/-1e13 Hz integral omits bounded Lorentzian tails ~Gamma/(2*pi^2*extent).
    assert abs(integral/expected-1)<1e-6
    chi,eta=lte_line_coefficients(np.array([630.]),T,ion,xi,**kwargs)
    B=2*c.HPlanck*nu0**3/c.CLight**2/np.expm1(c.HPlanck*nu0/(c.KBoltzmann*T))
    np.testing.assert_allclose(eta/chi,B[None,:],rtol=2e-15)
    zero=lte_line_coefficients(np.array([630.]),T,ion*0,xi,**kwargs)
    np.testing.assert_array_equal(zero,np.zeros_like(zero))
