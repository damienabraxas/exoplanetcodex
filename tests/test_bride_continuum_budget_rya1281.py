"""Actual LTE backend: exact diagnostic grid, Kirchhoff and component closure."""
import numpy as np
import pytest

pytest.importorskip('lightweaver')


def native_context(extra):
    from config.constants import require_codex_path
    from pipeline.bride_reference import make_reference_context,read_marcs_structure
    from pipeline.multi_atom import read_multi_atom
    folder=require_codex_path('experiments.bride_marcs_opacity')
    if not (folder/'manifest.json').exists():pytest.skip('Held native fixture not installed')
    import json
    manifest=json.loads((folder/'manifest.json').read_text())
    structure=read_marcs_structure(require_codex_path('grids.marcs_standard')/manifest['atmosphere']['name'])
    atom=read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))
    return make_reference_context(structure,atom,rays=2,qscale=1,lte=True,
        continuum_only=True,extra_wavelengths_vac_nm=extra)


def test_exact_grid_lte_kirchhoff_and_positive_component_remainder():
    from pipeline.bride_continuum_budget import continuum_budget
    from lightweaver import constants as c
    wave=np.array([130.237,135.394,276.976,500.072,777.995,1604.48])
    ctx=native_context(wave);before=np.array(ctx.detailedAtoms[0].n).copy()
    ctx.depthData.fill=True;ctx.formal_sol_gamma_matrices(lambdaIterate=True)
    np.testing.assert_array_equal(ctx.detailedAtoms[0].n,before)
    grid=np.asarray(ctx.spect.wavelength);index=np.searchsorted(grid,wave)
    np.testing.assert_array_equal(grid[index],wave)
    chi=np.asarray(ctx.depthData.chi)[:,0,0,:]
    eta=np.asarray(ctx.depthData.eta)[:,0,0,:]
    sca=np.asarray(ctx.background.sca);absorption=chi-sca
    temperature=ctx.atmos.pyAtmos.temperature
    planck=2*c.HPlanck*c.CLight/(grid[:,None]*1e-9)**3/np.expm1(
        c.HC/(grid[:,None]*1e-9*c.KBoltzmann*temperature[None,:]))
    np.testing.assert_allclose(eta[index],(absorption*planck)[index],rtol=2e-11,atol=0)
    bf,scattering=continuum_budget(ctx)
    oxygen=chi-np.asarray(ctx.background.chi)
    assert np.all(oxygen>=-1e-12*chi)
    assert np.all(absorption-oxygen-sum(bf.values())>=-1e-10*absorption)
    assert np.all(sca-sum(scattering.values())>=-1e-10*sca)
    # Ground-state red-wing scattering has a physical long-wave decline,
    # separate from wavelength-independent electron scattering.
    electron=scattering['electron_Thomson']
    np.testing.assert_array_equal(electron[index[0]],electron[index[-1]])
    hydrogen=scattering['H(Z=1)_resonant_Rayleigh']
    assert hydrogen[index[0],0]>hydrogen[index[-1],0]>0


@pytest.mark.parametrize('wave',[[0.],[np.nan],[[500.]],[]])
def test_extra_wavelengths_fail_closed(wave):
    with pytest.raises(ValueError,match='finite positive vacuum nm'):
        native_context(wave)
