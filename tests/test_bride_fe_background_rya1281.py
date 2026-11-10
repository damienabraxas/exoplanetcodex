"""Actual backend tests for source-derived Fe replacement without double counting."""
import copy
from functools import partial

import numpy as np
import pytest

pytest.importorskip('lightweaver')


def inputs():
    import json
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.multi_atom import read_multi_atom
    folder=require_codex_path('experiments.bride_fe_packet')
    if not (folder/'manifest.json').exists():pytest.skip('Pinned Fe packet absent')
    m=json.loads((folder/'manifest.json').read_text())
    s=read_marcs_structure(require_codex_path('grids.marcs_standard')/m['atmosphere']['name'])
    a=read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))
    return folder,s,a


@pytest.mark.parametrize('mode',['raw_specific','local_neutral'])
def test_actual_backend_replacement_emission_scattering_and_replay(mode):
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_fe_background import fe_factory,load_packet
    from pipeline.bride_continuum_budget import continuum_budget
    from lightweaver.atomic_table import PeriodicTable
    from lightweaver import constants as c
    path,s,a=inputs();extra=np.array([130.237,208.1,276.976,302.,777.3])
    base=make_reference_context(s,a,rays=2,qscale=1,extra_wavelengths_vac_nm=extra)
    modified=make_reference_context(s,a,rays=2,qscale=1,extra_wavelengths_vac_nm=extra,
        fe_background=fe_factory(path,mode))
    wave=np.asarray(base.spect.wavelength);packet=load_packet(path)
    mask=(wave>=packet['wavelength_vac_nm'][0])&(wave<=packet['wavelength_vac_nm'][-1])
    old=continuum_budget(base)[0]['Fe(Z=26)']
    field='opacity_m2_kg'if mode=='raw_specific'else 'opacity_net_per_neutral_Fe_m2'
    expected=np.column_stack([np.interp(wave[mask],packet['wavelength_vac_nm'],packet[field][:,d])for d in range(len(s['temperature_K']))])
    if mode=='raw_specific':expected*=s['rho_kg_m3'][None,:]
    else:
        fe=base.kwargs['eqPops'].atomicPops[PeriodicTable[26]]
        neutral=np.array([l.stage==0 for l in fe.model.levels]);expected*=fe.n[neutral].sum(axis=0)[None,:]
    change=np.asarray(modified.background.chi)-np.asarray(base.background.chi)
    np.testing.assert_allclose(change[mask],expected-old[mask],rtol=2e-12,atol=1e-30)
    np.testing.assert_array_equal(change[~mask],0.)
    np.testing.assert_array_equal(modified.background.sca,base.background.sca)
    planck=2*c.HPlanck*c.CLight/(wave[:,None]*1e-9)**3/np.expm1(c.HC/(wave[:,None]*1e-9*c.KBoltzmann*s['temperature_K'][None,:]))
    np.testing.assert_allclose(np.asarray(modified.background.eta)-np.asarray(base.background.eta),change*planck,rtol=2e-11,atol=1e-35)
    before=np.array(modified.background.chi).copy();modified.background.update_background(modified.atmos)
    np.testing.assert_array_equal(modified.background.chi,before)
    cloned=copy.deepcopy(modified)
    np.testing.assert_array_equal(cloned.background.chi,modified.background.chi)


def test_source_equal_existing_fe_is_identity_and_changed_node_rejected():
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_fe_background import FeIBackground,fe_factory,load_packet
    from pipeline.bride_continuum_budget import continuum_budget
    path,s,a=inputs();wave=np.array([208.1,276.976,302.])
    base=make_reference_context(s,a,rays=2,qscale=1,extra_wavelengths_vac_nm=wave)
    indices=np.searchsorted(np.asarray(base.spect.wavelength),wave)
    packet=load_packet(path);packet['wavelength_vac_nm']=wave
    packet['opacity_m2_kg']=continuum_budget(base)[0]['Fe(Z=26)'][indices]/s['rho_kg_m3'][None,:]
    same=make_reference_context(s,a,rays=2,qscale=1,extra_wavelengths_vac_nm=wave,
        fe_background=partial(FeIBackground,packet=packet,normalization='raw_specific'))
    np.testing.assert_allclose(np.asarray(same.background.chi)[indices],np.asarray(base.background.chi)[indices],rtol=2e-15,atol=0)
    shifted=copy.deepcopy(s);shifted['temperature_K'][0]+=10.
    with pytest.raises(ValueError,match='changed T/ne/nH'):
        make_reference_context(shifted,a,rays=2,qscale=1,fe_background=fe_factory(path,'local_neutral'))
