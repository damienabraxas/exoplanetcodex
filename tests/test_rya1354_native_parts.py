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


def test_node_checked_radial_budget_changes_only_native_inner_limit():
    from rya1354_lifetime_relax import scf_controls
    from rya1354_lifetime_ladder import STATES
    for state in STATES:
        for verify in (False,True):
            original=scf_controls(state,4,1500,.5,verify)
            expanded=scf_controls(state,4,1500,.5,verify,radial_iteration_floor=24)
            assert original.endswith('n\nn\n1\n')
            assert expanded==original[:-len('n\nn\n1\n')]+'y\n24\nn\n1\n'
    with pytest.raises(ValueError,match='positive'):
        scf_controls('ground',4,1500,.5,False,radial_iteration_floor=0)


def test_method3_divergence_is_retained_and_cannot_be_relabeled_completed(monkeypatch):
    import rya1354_method3_failure as audit
    bundle=ROOT/'data/audit/rya1354/lifetime_method3_failed_v7'
    report=audit.summarize(bundle)
    assert report==json.loads((ROOT/'data/audit/rya1354/source_curation_v2/method3_failure_v7.json').read_text())
    assert report['diagnostics']['iterations']==307
    assert report['diagnostics']['max_reported_self_consistency']==47.8
    assert report['accepted_lifetime'] is None and not report['adoption_allowed']
    files,manifest=load_parts(bundle)
    key='runs/lifetime_relax_n4_method3_v2/receipt.json'
    receipt=json.loads(files[key]);receipt['status']='CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED'
    files[key]=json.dumps(receipt).encode()
    monkeypatch.setattr(audit,'load_parts',lambda directory:(files,manifest))
    with pytest.raises(ValueError,match='failure receipt'):
        audit.summarize(bundle)


def test_native_relativistic_component_selection_and_explicit_control():
    from rya1354_lifetime_relax import scf_controls
    from rya1354_lifetime_ladder import STATES
    files,_=load_parts(ROOT/'data/audit/rya1354/rsl_component_v8')
    prefix='validation/rsl_component_v8/'
    receipt=json.loads(files[prefix+'receipt.json'])
    assert receipt['status']=='NATIVE_COMPONENT_SELECTION_VERIFIED_NO_ATOMIC_PREDICTION'
    assert receipt['atomic_prediction'] is False
    assert files[prefix+'component.stdout'].decode().splitlines()==['2 1 3','3 1 2 3','2 2 3']
    assert receipt['selections']==[[2,1,3],[3,1,2,3],[2,2,3]]
    for state in STATES:
        for verify in (False,True):
            old=scf_controls(state,4,1500,.5,verify,radial_iteration_floor=24)
            new=scf_controls(state,4,1500,.5,verify,radial_iteration_floor=24,both_2p_spectroscopic=True)
            assert new==old.replace('\n2s,2p', '\n2s,2p-,2p',1)


def test_history_does_not_confuse_inner_failures_or_partial_cycles_with_convergence():
    from rya1354_scf_history import history
    text=('Iteration number 1\n2s 1.0 1 2.0 1.0D-3 0 0 1 2 0 1\n'
          'Weighted average energy of these levels = -75.0\nGeneralised occupation numbers:\n'
          'Iteration number 2\n2s 1.0 2 2.0 9.0D-2 0 0 1 2 0 3\n'
          'Iteration number: 6, limit: 6\nConvergence not obtained\n')
    rows=history(text)
    assert [r['iteration'] for r in rows]==[1,2]
    assert rows[0]['complete_printed_cycle'] and not rows[1]['complete_printed_cycle']
    assert rows[1]['orbitals']['2s']==dict(residual=.09,method=2,nodes=3)
    with pytest.raises(ValueError,match='nonmonotonic'):
        history(text+'Iteration number 2\n')


def test_native_running_snapshots_preserve_actual_component_methods_without_admission():
    from rya1354_scf_history import history
    files,_=load_parts(ROOT/'data/audit/rya1354/scf_logs_v8')
    for name,method in (('nodebudget',3),('pcomponents',1)):
        prefix=f'evidence/{name}_log_snapshot_v8/'
        capture=json.loads(files[prefix+'capture.json'])
        assert capture['accepted_lifetime'] is None and not capture['adoption_allowed']
        for path,record in capture['files'].items():
            assert hashlib.sha256(files[prefix+path]).hexdigest()==record['sha256']
            assert len(files[prefix+path])==record['bytes']
        assert history(files[prefix+'rmcdhf.stdout'].decode())==capture['history']
        complete=[r for r in capture['history'] if r['complete_printed_cycle']]
        assert complete and all(r['orbitals']['2p-']['method']==method for r in complete)
        assert all(r['max_reported_residual']>5.960464477539063e-8 for r in complete)


def test_terminal_node_failures_replay_and_expose_bracket_guard():
    from rya1354_node_failure_summary import summarize
    report=summarize(ROOT/'data/audit/rya1354')
    assert report==json.loads((ROOT/'data/audit/rya1354/source_curation_v2/node_failures_v9.json').read_text())
    assert report['cases']['nodebudget']['first_193_cycles_exactly_reproduce_default_budget']
    assert report['cases']['pcomponents']['last_printed_outer_iteration']==206
    for case in report['cases'].values():
        final=case['failed_attempts'][-1]
        assert final['iterations']<final['limit']==24
        assert final['width_guard_triggered'] and final['counted_nodes']==3 and final['required_nodes']==1
        assert case['accepted_lifetime'] is None


def test_staged_preconditioning_changes_only_varied_scope_and_preserves_frozen_records(tmp_path):
    import struct
    from rya1354_lifetime_relax import scf_controls
    from rya1354_lifetime_staged import require_frozen
    for state in ('ground','quintet','triplet'):
        for verify in (False,True):
            full=scf_controls(state,4,1500,.5,verify,radial_iteration_floor=24,both_2p_spectroscopic=True)
            pre=scf_controls(state,4,1500,.5,verify,radial_iteration_floor=24,both_2p_spectroscopic=True,freeze_2s=True)
            assert pre==full.replace('\n2*,3*,4*\n','\n2p*,3*,4*\n',1)
    def record(data):
        marker=struct.pack('<i',len(data));return marker+data+marker
    def wave(two_s_energy):
        data=record(b'G92RWF')
        for n,e in ((1,-20.),(2,two_s_energy)):
            data+=record(struct.pack('<iidi',n,-1,e,2))
            data+=record(struct.pack('<5d',1.,0.,1.,0.,.01))
            data+=record(struct.pack('<2d',0.,.1))
        return data
    a=tmp_path/'a.w';b=tmp_path/'b.w';a.write_bytes(wave(-2.));b.write_bytes(wave(-2.))
    require_frozen(a,b,True)
    b.write_bytes(wave(-3.))
    with pytest.raises(ValueError,match='fixed radial'):
        require_frozen(a,b,True)
    require_frozen(a,b,False)  #2s is explicitly released;1s remains fixed.


def test_failed_preconditioner_cannot_release_core_or_complete_phase(tmp_path,monkeypatch):
    from types import SimpleNamespace
    import rya1354_lifetime_staged as staged
    called=[]
    def fail(args,phase,state,seed,receipt):
        called.append((phase,state))
        if state=='quintet':raise ValueError('unqualified preconditioner')
        return tmp_path/state
    monkeypatch.setattr(staged,'optimize',fail)
    receipt=dict(phases=['fixed_2s','released_2s'])
    with pytest.raises(ValueError,match='unqualified preconditioner'):
        staged.optimize_stages(SimpleNamespace(core=tmp_path,out=tmp_path),receipt)
    assert called==[('fixed_2s','ground'),('fixed_2s','quintet')]
    assert 'last_completed_phase' not in receipt


def test_live_snapshot_reads_explicit_staged_directory(tmp_path,monkeypatch):
    import rya1354_scf_history as audit
    monkeypatch.setattr(audit.platform,'node',lambda:'sirius')
    run=tmp_path/'run';out=tmp_path/'snapshot'
    stage=run/'fixed_2s_ground';stage.mkdir(parents=True)
    for name in ('rmcdhf.stdin','rmcdhf.stderr'):(stage/name).write_text('staged')
    (stage/'rmcdhf.stdout').write_text('Iteration number 7\n')
    audit.snapshot(run,out,'fixed_2s_ground')
    capture=json.loads((out/'capture.json').read_text())
    assert capture['source_stage']=='fixed_2s_ground' and capture['latest_printed_outer_iteration']==7
    assert capture['accepted_lifetime'] is None and not capture['adoption_allowed']
    with pytest.raises(ValueError,match='stage directory'):
        audit.snapshot(run,tmp_path/'bad','../other')


def test_native_preconditioner_snapshot_uses_intended_fixed_core_scope():
    from rya1354_lifetime_relax import scf_controls
    from rya1354_scf_history import history
    files,_=load_parts(ROOT/'data/audit/rya1354/staged_logs_v9')
    prefix='evidence/staged_log_snapshot_v9/'
    capture=json.loads(files[prefix+'capture.json'])
    assert capture['source_stage']=='fixed_2s_ground' and capture['accepted_lifetime'] is None
    assert files[prefix+'rmcdhf.stdin'].decode()==scf_controls('ground',4,1500,.5,False,
        radial_iteration_floor=24,both_2p_spectroscopic=True,freeze_2s=True)
    rows=history(files[prefix+'rmcdhf.stdout'].decode())
    assert rows==capture['history']
    complete=[r for r in rows if r['complete_printed_cycle']]
    assert complete
    for row in complete:
        assert len(row['orbitals'])==14 and '1s' not in row['orbitals'] and '2s' not in row['orbitals']
        assert row['orbitals']['2p-']['method']==row['orbitals']['2p']['method']==1
    for name,record in capture['files'].items():
        assert hashlib.sha256(files[prefix+name]).hexdigest()==record['sha256']
