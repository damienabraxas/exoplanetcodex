"""Actual native lifetime and failure replay; no accuracy from gauge agreement."""
import hashlib
import json
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_lifetime_summary import load_parts,scf_state,decay,mchf_stop,summarize
from rya1354_lifetime_mchf import validate_hf_seed
from rya1354_two_method_summary import verify_artifacts
from rya1354_pack_native import pack
from rya1354_wave_probe import decode,difference


@pytest.fixture(scope='module')
def aux():
    return load_parts(ROOT/'data/audit/rya1354/lifetime_aux_v4b')[0]


@pytest.fixture(scope='module')
def ladder():
    return load_parts(ROOT/'data/audit/rya1354/lifetime_ladder_v4')[0]


def test_every_balanced_state_meets_native_criterion_before_promotion(ladder):
    root='runs/lifetime_ladder_v3'
    receipt=json.loads(ladder[root+'/receipt.json'])
    assert len(receipt['qualified_states'])==24
    for n in (0,3,4,5):
        for state in ('ground','quintet','triplet'):
            path=f'{root}/n{n}_{state}'
            for suffix in ('','_verify'):
                assert scf_state(ladder,path+suffix,state)['all_reported_functions_below_ACC']
            assert ladder[path+'_verify/rwfn.inp']==ladder[path+'/state.w']
            if n:
                previous=0 if n==3 else n-1
                assert ladder[path+'/previous.w']==ladder[f'{root}/n{previous}_{state}_verify/state.w']
        d=decay(ladder,f'{root}/n{n}_')
        assert len(d['channels'])==13
        for g in ('B','C'):
            assert d['model'][g]['total_A_s1']*d['model'][g]['model_one_photon_lifetime_ns']==pytest.approx(1e9)


def test_invalid_hf_and_unconverged_independent_state_are_rejected(aux):
    bad='runs/lifetime_mchf_v1/reference_3P/'
    with pytest.raises(ValueError):
        validate_hf_seed(aux[bad+'HF.stdout'].decode()+aux[bad+'HF.stderr'].decode(),'3P')
    good='runs/lifetime_mchf_v3/reference_3P/'
    validate_hf_seed(aux[good+'HF.stdout'].decode()+aux[good+'HF.stderr'].decode(),'3P')
    with pytest.raises(ValueError,match='stop missing'):
        mchf_stop(aux['runs/lifetime_mchf_v2/n5_3S/mchf.stdout'].decode())


def test_archive_cache_omission_is_explicit_and_cannot_hide_a_wavefunction(tmp_path):
    source=tmp_path/'source';source.mkdir()
    output=tmp_path/'bundle'
    (source/'state.w').write_bytes(b'wavefunction')
    (source/'rci.res').write_bytes(b'regenerable restart matrix')
    records={p.name:dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in source.iterdir()}
    (source/'receipt.json').write_text(json.dumps(dict(files=records)))
    pack(tmp_path,[Path('source')],output,omit_rci_restart_cache=True)
    files,manifest=load_parts(output)
    assert 'source/rci.res' not in files
    assert manifest['external_restart_caches']['source/rci.res']['sha256']==records['rci.res']['sha256']
    with pytest.raises(ValueError):
        verify_artifacts(files)
    missing=dict(files);del missing['source/state.w']
    forged=dict(manifest['external_restart_caches']);forged['source/state.w']=records['state.w']
    with pytest.raises(ValueError,match='identity differs'):
        verify_artifacts(missing,external_restart_caches=forged)


def test_complete_lifetime_study_replays_without_calibrated_uncertainty(tmp_path):
    out=tmp_path/'summary.json'
    summarize(out)
    assert out.read_bytes()==(ROOT/'data/audit/rya1354/source_curation_v2/lifetime_summary_v4.json').read_bytes()
    d=json.loads(out.read_text())
    assert [r['n'] for r in d['balanced_layers']]==[0,3,4,5]
    assert set(d['core_controls'])=={'3','4','5'}
    assert not d['adoption_allowed']
    assert d['uncertainty'] is None and d['calibrated_sigma'] is None and d['covariance'] is None


def test_saved_wave_diagnostic_refuses_truncation_and_grid_mixing(ladder):
    data=ladder['runs/lifetime_ladder_v3/n3_triplet_verify/state.w']
    parsed=decode(data)
    assert len(parsed)==9 and (1,-1) in parsed and (3,-1) in parsed
    assert all(r['max_saved_radial_change']==0 for r in difference(parsed,parsed))
    with pytest.raises(ValueError):
        decode(data[:-1])
    changed={key:dict(value) for key,value in parsed.items()}
    radius=list(changed[(1,-1)]['radius']);radius[1]*=1.001
    changed[(1,-1)]['radius']=tuple(radius)
    with pytest.raises(ValueError,match='grid changed'):
        difference(parsed,changed)
