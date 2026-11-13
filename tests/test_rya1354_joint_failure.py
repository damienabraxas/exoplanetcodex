"""Retain timeout as failure; replay matched-CI controls against native evidence."""
import json
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts


def test_joint_timeout_retains_all_cycles_without_admitting_a_state():
    from rya1354_joint_failure import summarize
    base=ROOT/'data/audit/rya1354'
    result=summarize(base)
    assert result==json.loads((base/'source_curation_v2/joint_failure_v16.json').read_text())
    assert result['completed_cycle_count']==459
    assert result['last_printed_cycle']['iteration']==460
    assert result['best_complete_cycle']['max_reported_residual']>result['original_ACC']
    assert result['exact_fixed_n_le_3_records_in_retained_wave']==9
    assert result['qualified_states']==0 and result['accepted_lifetime'] is None


def test_timed_out_joint_campaign_cannot_be_relabelled_complete(monkeypatch):
    import rya1354_joint_failure as audit
    original=audit.load_parts
    def changed(path):
        files,manifest=original(path)
        if path.name=='joint_failed_v16':
            key=audit.RUN+'receipt.json';r=json.loads(files[key]);r['status']='BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED'
            files[key]=json.dumps(r).encode()
        return files,manifest
    monkeypatch.setattr(audit,'load_parts',changed)
    with pytest.raises(ValueError,match='failure disposition'):
        audit.summarize(ROOT/'data/audit/rya1354')


def test_new_ci_controls_match_retained_native_hamiltonian_and_root_inputs():
    from rya1354_qualified_ci import controls
    files,_=load_parts(ROOT/'data/audit/rya1354/lifetime_core_n5_v4')
    for state in ('ground','quintet','triplet'):
        for case in ('DC','Breit'):
            key='runs/lifetime_core_n5_v1/CV_2s_'+case+'_'+state+'/rci.stdin'
            assert controls(state,case)==files[key].decode()


def test_ci_startup_checkpoint_cannot_be_reported_as_a_completed_model():
    from rya1354_qualified_ci_audit import summarize
    base=ROOT/'data/audit/rya1354'
    with pytest.raises(ValueError,match='No terminal qualified CI receipt'):
        summarize(base,base/'qualified_ci_initial_v16')
