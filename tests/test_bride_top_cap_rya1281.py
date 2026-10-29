"""Independent slab limits and actual backend checks for the finite-column cap."""
import copy
from types import SimpleNamespace as S

import numpy as np
import pytest

pytest.importorskip('lightweaver')
from pipeline.bride_top_cap import cap_intensity,FiniteColumnCap


def test_cap_zero_thin_thick_pure_emission_and_angular_limits():
    chi=np.array([0.,1e-12,1.,1e4]);source=7.;eta=source*chi
    mu=np.array([.2,.7,1.])
    np.testing.assert_array_equal(cap_intensity(chi,eta,0.,mu),0.)
    result=cap_intensity(chi,eta,2.,mu)
    np.testing.assert_allclose(result[1],eta[1]*2/mu,rtol=6e-12)
    np.testing.assert_allclose(result[-1],source,rtol=1e-15)
    np.testing.assert_allclose(result[2],source*(1-np.exp(-2/mu)),rtol=1e-15)
    np.testing.assert_allclose(cap_intensity([0.],[3.],2.,mu),[6/mu])
    # Equal source and incoming intensity remains unchanged through any slab.
    transmitted=source*np.exp(-chi[:,None]*2/mu)+result
    np.testing.assert_allclose(transmitted,source,rtol=1e-15)
    for bad in (-1.,np.nan):
        with pytest.raises(ValueError):cap_intensity([bad],[1.],1.,mu)
    with pytest.raises(ValueError):cap_intensity([1.],[1.],1.,[0.])


def test_snapshot_includes_scattering_and_is_deepcopy_safe():
    cap=FiniteColumnCap(2.)
    wave=np.array([400.,500.,600.]);chi=np.ones((3,2,2,4))*.5
    eta=np.ones_like(chi)*3.
    ctx=S(spect=S(wavelength=wave),depthData=S(chi=chi,eta=eta),
          background=S(sca=np.ones((3,4))*.2),atmos=S(pyAtmos=S(muz=np.array([.2,.7]))))
    previous=np.array([2.,4.,6.])
    assert cap.update_from_context(ctx,previous)==1.
    np.testing.assert_array_equal(cap.eta,3+.2*previous)
    assert cap.update_from_context(ctx,previous)==0.
    clone=copy.deepcopy(cap);clone.muz=np.array([1.])
    expected=cap_intensity([.5],[3.8],2.,[1.])
    np.testing.assert_allclose(clone.compute_bc(None,S(wavelength=np.array([500.])))[:,:,0],expected)
    with pytest.raises(ValueError,match='extrapolate'):
        clone.compute_bc(None,S(wavelength=np.array([399.])))


def test_actual_backend_isothermal_cap_matches_independent_slab():
    import lightweaver as lw
    from lightweaver.LwCompiled import BackgroundProvider

    class Bottom(lw.BoundaryCondition):
        def compute_bc(self,atmos,spect):
            return np.full((len(spect.wavelength),len(self.muz),1),9.)

    class Coefficients(BackgroundProvider):
        def compute_background(self,atmos,chi,eta,sca):
            chi[:]=2.;eta[:]=14.;sca[:]=0.

    cap=FiniteColumnCap(3.)
    atmos=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,
        depthScale=np.linspace(1.,0.,7),temperature=np.full(7,6000.),ne=np.full(7,1e16),
        nHTot=np.full(7,1e20),vturb=np.zeros(7),vlos=np.zeros(7),logG=2.,
        lowerBc=Bottom(),upperBc=cap)
    mus=np.array([.2,.7]);atmos.rays(muz=mus,upOnly=False)
    spect=S(wavelength=np.array([500.,600.]),radSet=S(activeAtoms=[],detailedAtoms=[]))
    ctx=lw.Context(atmos,spect,{},backgroundProvider=Coefficients,Nthreads=1,
                   formalSolver='piecewise_besser_1d')
    cap=ctx.atmos.pyAtmos.zUpperBc
    ctx.depthData.fill=True
    for _ in range(3):
        previous=np.array(ctx.spect.J)[:,0].copy();ctx.formal_sol_gamma_matrices(lambdaIterate=True)
        cap.update_from_context(ctx,previous)
    intensity=np.array(ctx.depthData.I)
    # At grid top, inward ray crossed cap length3 with chi2.
    expected_in=7*(1-np.exp(-6/mus))
    np.testing.assert_allclose(intensity[:, :, 0, 0],np.broadcast_to(expected_in,(2,2)),rtol=2e-14)
    # Outward grid-top ray crossed internal slab length1 with independent bottom9.
    expected_out=7+2*np.exp(-2/mus)
    np.testing.assert_allclose(intensity[:, :, 1, 0],np.broadcast_to(expected_out,(2,2)),rtol=2e-14)


def test_native_cap_constructor_and_diagnostic_ray_grid():
    import json
    from config.constants import codex_path,require_codex_path
    from pipeline.bride_reference import read_marcs_structure,make_reference_context
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_top_cap import emergent_rays
    if not codex_path('experiments.bride_oxygen_atom').exists():
        pytest.skip('Held benchmark not installed')
    atom_path=require_codex_path('experiments.bride_oxygen_atom')
    ref=json.loads((atom_path.parent/'reference-5750-4.5-A8.70-v1.json').read_text())
    structure=read_marcs_structure(require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name'])
    source=read_multi_atom(atom_path)
    ctx=make_reference_context(structure,source,rays=3,qscale=1,finite_column_top_cap=True)
    cap=ctx.atmos.pyAtmos.zUpperBc
    assert structure['column_mass_kg_m2'][0]==pytest.approx(.09156734,rel=1e-14)
    assert cap.length_m==pytest.approx(.09156734/1.0854e-6,rel=1e-14)
    for _ in range(3):
        old=np.array(ctx.spect.J)[:,0].copy();ctx.formal_sol_gamma_matrices(lambdaIterate=True)
        cap.update_from_context(ctx,old)
    before=cap.chi.copy()
    grid,outer=emergent_rays(ctx,np.linspace(777.2,777.3,7),np.array([.2,.7,1.]))
    assert grid.shape==outer.shape==(7,3)
    assert np.isfinite(outer).all() and np.all(outer>0)
    np.testing.assert_array_equal(cap.chi,before)
