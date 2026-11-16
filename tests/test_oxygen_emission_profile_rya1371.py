import numpy as np
import pytest
pytest.importorskip('lightweaver')
from pipeline.oxygen_emission_profile import diffusion_intensity,profile_planck,EmissionProfileBoundary


def test_diffusion_gradient_and_angle_are_retained():
    result=diffusion_intensity([[2.,4.]],[[1.,3.]],2.,[.5,1.])
    np.testing.assert_array_equal(result,[[4.25,4.5]])
    np.testing.assert_array_equal(diffusion_intensity([[4.,4.]],[[1.,3.]],2.,[.5,1.]),[[4.,4.]])


@pytest.mark.parametrize('B,chi,spacing,mu',[([[2.,4.]],[[0.,0.]],2.,[1.]),([[4.,1.]],[[1.,1.]],.01,[1.]),([[2.,4.]],[[1.,1.]],2.,[0.]),([[2.,float('nan')]],[[1.,1.]],2.,[1.])])
def test_invalid_diffusion_refuses_without_clipping(B,chi,spacing,mu):
    with pytest.raises(ValueError):diffusion_intensity(B,chi,spacing,mu)


def test_held_profile_matches_independent_planck_at_test_grid():
    import lightweaver.constants as c
    wave=np.array([90.,777.,2000.]);T=np.array([4000.,6000.]);nu=c.CLight/(wave[:,None]*1e-9)
    expected=2*c.HPlanck*nu**3/c.CLight**2/np.expm1(c.HPlanck*nu/(c.KBoltzmann*T[None,:]))
    np.testing.assert_allclose(profile_planck(wave,T,'held-LW',None),expected,rtol=2e-14,atol=0)
    with pytest.raises(ValueError):profile_planck(wave,T,'unknown',None)


def test_fixed_diagnostic_gate_precedes_native_calls(monkeypatch):
    from pipeline.bride_hop_transfer import make_context
    monkeypatch.setattr('pipeline.bride_hop_transfer.compute_packet',lambda *a,**k:pytest.fail('Native call before refusal'))
    with pytest.raises(ValueError,match='fixed-population'):
        make_context({},None,None,None,None,extra_vac_nm=[],emission_profile_control=dict(profile='held-LW',native_source=None))


def test_callable_diffusion_reproduces_compiled_boundary_on_disposable_column():
    import lightweaver as lw
    from lightweaver.LwCompiled import BackgroundProvider
    from lightweaver import rh_atoms
    # Passive-only continuum: no atomic-data packet, native engine or solar fit.
    wave=np.array([100.,500.,777.,2000.]);height=np.array([200.,100.,0.]);T=np.array([4500.,5000.,6000.])
    chi=np.array([[.01,.02,.04],[.001,.002,.004],[.005,.01,.02],[.03,.04,.05]])
    class Provider(BackgroundProvider):
        def compute_background(self,atmos,out_chi,out_eta,out_sca):
            out_chi[:]=chi;out_eta[:]=chi*profile_planck(wave,T,'held-LW',None);out_sca[:]=0
    def build(bc):
        atmos=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=height.copy(),temperature=T.copy(),ne=np.full(3,1e16),nHTot=np.full(3,1e20),vturb=np.zeros(3),vlos=np.zeros(3),lowerBc=bc)
        atmos.quadrature(3)
        rad=lw.RadiativeSet([rh_atoms.H_6_atom()]);rad.set_detailed_static('H');eq=rad.compute_eq_pops(atmos,mols=lw.MolecularTable([]));spect=rad.compute_wavelength_grid(extraWavelengths=wave)
        # A passive radset still adds background-grid nodes; supply the explicit
        # no-transition grid used by this disposable continuum fixture.
        spect.radSet.set_passive('H');spect.wavelength=wave.copy()
        ctx=lw.Context(atmos,spect,eq,Nthreads=1,formalSolver='piecewise_linear_1d',backgroundProvider=Provider);ctx.depthData.fill=True
        return ctx
    expected=build(lw.ThermalisedRadiation());expected.formal_sol_gamma_matrices(lambdaIterate=True)
    boundary=EmissionProfileBoundary('held-LW',None);actual=build(boundary)
    actual.formal_sol_gamma_matrices(lambdaIterate=True);boundary.update_from_context(actual)
    actual.formal_sol_gamma_matrices(lambdaIterate=True);assert boundary.update_from_context(actual)==0
    np.testing.assert_allclose(actual.spect.J,expected.spect.J,rtol=3e-14,atol=0)
    np.testing.assert_allclose(actual.spect.I,expected.spect.I,rtol=3e-14,atol=0)
