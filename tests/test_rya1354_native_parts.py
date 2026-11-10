"""Lossless oversized native file packing, corruption refusal and legacy replay."""
import hashlib
import json
import os
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_pack_native import pack
from rya1354_native_parts import load_parts,reassemble


def test_incompressible_single_file_is_fragmented_and_native_receipt_replays(tmp_path):
    source=tmp_path/'run';source.mkdir()
    data=os.urandom(32000);(source/'large.+2T').write_bytes(data)
    record=dict(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
    (source/'receipt.json').write_text(json.dumps(dict(files={'large.+2T':record})))
    out=tmp_path/'bundle';pack(tmp_path,[Path('run')],out,limit=4096)
    files,manifest=load_parts(out)
    assert manifest['schema'].endswith('.v2')
    assert files['run/large.+2T']==data and (source/'large.+2T').read_bytes()==data
    assert all((out/name).stat().st_size<=4096 for name in manifest['parts'])
    assert len(manifest['fragments']['run/large.+2T'])>1
    damaged=dict(manifest);damaged['fragments']=json.loads(json.dumps(manifest['fragments']))
    damaged['fragments']['run/large.+2T'].reverse()
    (out/'manifest.json').write_text(json.dumps(damaged))
    with pytest.raises(ValueError,match='noncontiguous'):
        load_parts(out)


def test_fragment_reassembly_refuses_missing_duplicate_and_corrupt_chunks():
    data=b'abc';r=dict(member='.rya1354_fragments/key/0',offset=0,bytes=3,sha256=hashlib.sha256(data).hexdigest())
    with pytest.raises(ValueError,match='Missing'):
        reassemble({}, {'native':[r]})
    with pytest.raises(ValueError,match='bytes/hash'):
        reassemble({r['member']:b'bad'},{'native':[r]})
    with pytest.raises(ValueError,match='Duplicate'):
        reassemble({r['member']:data},{'native':[r,r]})
    with pytest.raises(ValueError,match='Unclaimed'):
        reassemble({r['member']:data},{})


def test_new_reader_preserves_legacy_evidence():
    files,manifest=load_parts(ROOT/'data/audit/rya1354/mr_inventory_v5')
    assert manifest['schema'].endswith('.v1')
    assert json.loads(files['runs/mr_inventory_n4_v1/receipt.json'])['status']=='INVENTORY_COMPLETED_NO_ATOMIC_PREDICTION'


def test_collection_failure_replaces_stale_waiting_heartbeat(tmp_path,monkeypatch):
    import rya1354_collect_pending as collector
    out=tmp_path/'evidence';out.mkdir()
    (out/'status.json').write_text('{"stale":"WAITING"}')
    source=tmp_path/'runs'/collector.CAMPAIGNS[0];source.mkdir(parents=True)
    (source/'receipt.json').write_text(json.dumps(dict(status='FAILED_RETAINED')))
    def fail(*args,**kwargs):
        raise ValueError('forced packing error')
    monkeypatch.setattr(collector,'pack',fail)
    with pytest.raises(ValueError,match='forced packing error'):
        collector.observe(tmp_path,out)
    status=json.loads((out/'status.json').read_text())
    assert status['collection_error']=='ValueError: forced packing error'
    assert status['campaigns'][collector.CAMPAIGNS[0]]['native_status']=='FAILED_RETAINED'
    assert all(s['collection']=='COLLECTION_FAILED_REVIEW_REQUIRED' for s in status['campaigns'].values())
    assert status['adoption_allowed'] is False


def test_two_s_method3_changes_only_native_method_selection():
    from rya1354_lifetime_relax import scf_controls
    from rya1354_lifetime_ladder import STATES
    for state in STATES:
        for verify in (False,True):
            old=scf_controls(state,4,1500,.5,verify)
            changed=scf_controls(state,4,1500,.5,verify,two_s_method3=True)
            marker='1500\ny\nn\n'
            before,after=old.split(marker)
            assert after.startswith('n\n')
            assert changed==(old if verify else before+marker+'y\n\n\n2s\n\n'+after[2:])
            assert '439\ny\n5.9604644775390625e-08\n' in changed


def test_default_controls_match_preserved_failed_native_input():
    from rya1354_lifetime_relax import scf_controls
    files,_=load_parts(ROOT/'data/audit/rya1354/lifetime_relax_failed_v6')
    assert files['runs/lifetime_relax_n4_v1/ground/rmcdhf.stdin'].decode()==scf_controls('ground',4,1500,.5,False)


def test_explicit_campaign_watcher_cannot_escape_task_runs(tmp_path):
    from rya1354_collect_pending import collect
    name='lifetime_relax_n4_method3_v2'
    result=collect(tmp_path,tmp_path/'evidence',[name])
    assert set(result)=={name}
    assert result[name]['collection']=='WAITING_FOR_TERMINAL_RECEIPT'
    for invalid in ('../other','/tmp/other','.','..'):
        with pytest.raises(ValueError,match='task-owned'):
            collect(tmp_path,tmp_path/'evidence',[invalid])


def test_full_native_mr_terminal_audit_and_failed_relaxation_replay():
    from rya1354_terminal_summary import summarize
    base=ROOT/'data/audit/rya1354'
    report=summarize(base/'lifetime_mr_complete_v6',base/'lifetime_relax_failed_v6')
    assert report==json.loads((base/'source_curation_v2/terminal_summary_v6.json').read_text())
    cases=report['fixed_orbital_MR']['cases']
    assert len(cases)==4 and all(len(c['channels'])==13 for c in cases.values())
    assert cases['CV_2s_MR_Breit']['model']['B']['model_one_photon_lifetime_ns']==pytest.approx(2.29598414558)
    assert cases['CV_2s_MR_Breit']['model']['C']['model_one_photon_lifetime_ns']==pytest.approx(1.73268651030)
    assert report['failed_relaxation']['accepted_lifetime'] is None
    assert not report['adoption_allowed'] and report['uncertainty'] is report['covariance'] is None


def test_terminal_audit_refuses_misdeclared_native_hamiltonian(monkeypatch):
    import rya1354_terminal_summary as audit
    base=ROOT/'data/audit/rya1354'
    files,manifest=load_parts(base/'lifetime_mr_complete_v6')
    files['runs/lifetime_mr_n4_v1/CV_2s_MR_Breit_ground/rci.stdin']=files['runs/lifetime_mr_n4_v1/CV_2s_MR_DC_ground/rci.stdin']
    monkeypatch.setattr(audit,'load_parts',lambda directory:(files,manifest))
    with pytest.raises(ValueError,match='Hamiltonian/root input'):
        audit.summarize(base/'lifetime_mr_complete_v6',base/'lifetime_relax_failed_v6')
