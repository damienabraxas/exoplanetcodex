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


def test_core_relaxation_refuses_wrong_layer_and_modified_initializer(aux, tmp_path):
    from rya1354_lifetime_relax import validate_core_layer, variational_orbitals
    root='runs/lifetime_core_n4_v1'
    receipt=json.loads(aux[root+'/receipt.json'])
    (tmp_path/'receipt.json').write_bytes(aux[root+'/receipt.json'])
    for state in ('ground','quintet','triplet'):
        for name in ('state.c','state.w','isodata','state.lsj.lbl'):
            relative=f'CV_2s_DC_{state}/{name}'
            path=tmp_path/relative;path.parent.mkdir(exist_ok=True)
            path.write_bytes(aux[root+'/'+relative])
    assert validate_core_layer(tmp_path,4)==hashlib.sha256((tmp_path/'receipt.json').read_bytes()).hexdigest()
    assert variational_orbitals(3)=='2*,3*'
    assert variational_orbitals(4)=='2*,3*,4*'
    assert variational_orbitals(5)=='2*,3*,4*,5*'
    with pytest.raises(ValueError,match='layer/status'):
        validate_core_layer(tmp_path,3)
    (tmp_path/'CV_2s_DC_triplet/state.w').write_bytes(b'wrong wavefunction')
    with pytest.raises(ValueError,match='differs from native receipt'):
        validate_core_layer(tmp_path,4)
    with pytest.raises(ValueError):
        variational_orbitals(6)


def test_mr_inventory_checks_actual_native_electron_parity_and_duplicates(aux):
    from rya1354_mr_inventory import inventory, references
    root='runs/lifetime_core_n4_v1'
    odd=aux[root+'/CV_2s_DC_triplet/state.c'].decode()
    even=aux[root+'/CV_2s_DC_ground/state.c'].decode()
    assert inventory(odd,'odd',True)['CSFs']==2795
    assert inventory(even,'even',True)['closed_1s']
    with pytest.raises(ValueError,match='parity'):
        inventory(odd,'even',True)
    with pytest.raises(ValueError,match='occupation'):
        inventory(odd.replace('1s ( 2)','1s ( 1)',1),'odd',True)
    with pytest.raises(ValueError,match='Duplicate'):
        inventory(odd+'\n'.join(odd.splitlines()[5:8])+'\n','odd',True)
    assert len(references('odd',True,True))==4
    assert len(references('even',False,False))==1


def test_native_mr_inventory_superset_and_only_completed_subcases_are_reported():
    from rya1354_mr_summary import summarize_inventory, summarize as mr_summary
    files,_=load_parts(ROOT/'data/audit/rya1354/mr_inventory_v5')
    result=summarize_inventory(files)
    assert len(result['cases'])==8
    assert result['expansions']['even_CV']==dict(SR_CSFs=3459,MR_CSFs=29185,additional_CSFs=25726,SR_subset_MR=True)
    assert result['expansions']['odd_CV']['MR_CSFs']==25382
    report=mr_summary()
    assert report['fixed_orbital_MR'] is None and report['relaxed_n4'] is None
    completed=report['completed_MR_subcases']
    assert set(completed['cases'])=={'DC','VV_MR_DC'}
    assert completed['cases']['VV_MR_DC']['model']['B']['model_one_photon_lifetime_ns']==pytest.approx(2.4595659765366853)
    assert report==json.loads((ROOT/'data/audit/rya1354/source_curation_v2/mr_summary_v5.json').read_text())
    from rya1354_mr_summary import fixed_cases
    partial,_=load_parts(ROOT/'data/audit/rya1354/mr_partial_v5')
    partial['runs/lifetime_mr_n4_v1/DC_triplet/state.w']+=b'wrong wavefunction'
    with pytest.raises(ValueError,match='identity differs'):
        fixed_cases(partial,[dict(name='DC',core=False)],'runs/lifetime_mr_n4_v1')
    assert not report['adoption_allowed'] and report['calibrated_sigma'] is None
    partial,_=load_parts(ROOT/'data/audit/rya1354/mr_partial_v5')
    partial['runs/lifetime_mr_n4_v1/VV_MR_DC_triplet/state.c']=partial['runs/lifetime_mr_n4_v1/DC_triplet/state.c']
    with pytest.raises(ValueError,match='independently audited generator'):
        fixed_cases(partial,[dict(name='DC',core=False),dict(name='VV_MR_DC',core=False)],'runs/lifetime_mr_n4_v1')
    changed=dict(files)
    path='runs/mr_inventory_n4_v1/receipt.json'
    receipt=json.loads(changed[path]);receipt['cases']['odd_CV_MR']['CSFs']-=1
    changed[path]=json.dumps(receipt).encode()
    with pytest.raises(ValueError,match='differs from actual'):
        summarize_inventory(changed)


def test_pending_collector_keeps_failed_campaigns_unaccepted(tmp_path):
    from rya1354_collect_pending import collect, terminal_receipt
    name='lifetime_relax_n4_v1'
    source=tmp_path/'runs'/name;source.mkdir(parents=True)
    path=source/'receipt.json'
    path.write_text('{')
    assert terminal_receipt(path) is None
    path.write_text(json.dumps(dict(status='FAILED_RETAINED',error='strict SCF cap',files={})))
    out=tmp_path/'evidence';out.mkdir()
    result=collect(tmp_path,out)
    assert result[name]['native_status']=='FAILED_RETAINED'
    assert result[name]['scientific_admission']=='HOLD' and not result[name]['adoption_allowed']
    assert result['lifetime_mr_n4_v1']['collection']=='WAITING_FOR_TERMINAL_RECEIPT'
    assert collect(tmp_path,out)==result
    path.write_text(json.dumps(dict(status='CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED',files={})))
    with pytest.raises(ValueError,match='changed after collection'):
        collect(tmp_path,out)
