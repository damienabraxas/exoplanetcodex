"""Native n5 checkpoint replay and refusal to admit incomplete evidence."""
import json
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_extension_progress import summarize,LOG,RUN
from rya1354_extension_audit import stage


def test_native_checkpoint_loading_and_initialization_are_verified_not_optimization():
    base=ROOT/'data/audit/rya1354'
    report=summarize(base)
    assert report==json.loads((base/'source_curation_v2/n5_progress_v14.json').read_text())
    assert report['completed_cycle_count']==24
    assert report['last_complete_cycle']['max_reported_residual']>report['original_ACC']
    assert report['loading_control_verified'] and report['baseline_orbitals']==25
    assert not report['optimization_qualified'] and report['accepted_lifetime'] is None


def test_live_native_log_cannot_be_admitted_as_closed_stage():
    files,_=load_parts(ROOT/'data/audit/rya1354/n5_ground_snapshot_v14')
    for name in ('rmcdhf.stdout','rmcdhf.stderr'):
        files[RUN+'outer_n5_ground/'+name]=files[LOG+name]
    with pytest.raises(ValueError,match='Native final orbital residual table absent'):
        stage(files,{}, {},'ground',False)


def test_checkpoint_rejects_replaced_native_loading_input(monkeypatch):
    import rya1354_extension_progress as audit
    original=audit.load_parts
    def changed(path):
        files,manifest=original(path)
        if path.name=='n5_ground_snapshot_v14':
            files[RUN+'outer_n5_ground/fixed_load/rwfn.inp']=b'changed'
        return files,manifest
    monkeypatch.setattr(audit,'load_parts',changed)
    with pytest.raises(ValueError,match='loading inputs differ'):
        audit.summarize(ROOT/'data/audit/rya1354')


def test_full_n5_native_model_replays_six_states_and_thirteen_channels():
    from rya1354_n5_terminal import report
    base=ROOT/'data/audit/rya1354'
    result=report(base)
    assert result==json.loads((base/'source_curation_v2/n5_terminal_v15.json').read_text())
    assert len(result['n5']['states'])==6 and len(result['n5']['decay']['channels'])==13
    assert result['n5']['decay']['model']['B']['model_one_photon_lifetime_ns']==pytest.approx(2.28448374775)
    assert result['n5']['decay']['model']['C']['model_one_photon_lifetime_ns']==pytest.approx(1.71320681549)
    assert result['accepted_lifetime'] is result['uncertainty'] is None
    assert not result['n5']['fully_released_qualified'] and not result['adoption_allowed']


@pytest.mark.parametrize('replacement',['loading_wave','transition_wave'])
def test_n5_terminal_audit_refuses_substituted_native_waves(monkeypatch,replacement):
    import rya1354_extension_audit as audit
    original=audit.load_parts
    def changed(path):
        files,manifest=original(path)
        if path.name=='n5_complete_v15':
            if replacement=='loading_wave':
                files[RUN+'outer_n5_ground/fixed_load/rwfn.out']=files[RUN+'outer_n5_quintet/fixed_load/rwfn.out']
            else:
                files[RUN+'to_ground/upper.w']=files[RUN+'outer_n5_ground_verify/state.w']
        return files,manifest
    monkeypatch.setattr(audit,'load_parts',changed)
    with pytest.raises(ValueError,match='frozen radial records differ|transition input differs'):
        audit.summarize(ROOT/'data/audit/rya1354',ROOT/'data/audit/rya1354/n5_complete_v15')
