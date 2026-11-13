"""Pinned FAC stopping markers, thresholds and deliberately incomplete logs."""
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_controls import diagnostics,input_text


@pytest.fixture(scope='module')
def native_log():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_oxygen_complete_v17')
    return files['runs/fac_oxygen_levels_v17/sfac.stdout'].decode()


def test_native_stopping_replay_requires_marker_and_scaled_conditions(native_log):
    result=diagnostics(native_log,3.0,2048)
    assert result['iterations']==18 and result['last']['iteration']==17
    assert result['printed_rule_verified']
    assert not diagnostics('\n'.join(native_log.splitlines()[:-1]),3.0,2048)['printed_rule_verified']
    capped='\n'.join(native_log.splitlines()[:2])
    assert not diagnostics(capped,3.0,1)['printed_rule_verified']


def test_wrong_declared_threshold_is_refused(native_log):
    with pytest.raises(ValueError,match='threshold multiplier differs'):
        diagnostics(native_log,.3,2048)


def test_input_keeps_level_only_scope_and_explicit_controls():
    text=input_text(.03,256)
    assert 'SetOptimizeTolerance(0.03)' in text and 'SetOptimizeMaxIter(256)' in text
    assert "BasisTable('o.basis',0)" in text and "BasisTable('o_export',10)" in text
    assert 'TransitionTable' not in text and 'SetBreit(0)' in text


def test_native_capped_and_coherent_campaigns_replay_without_accuracy_admission():
    import json
    from rya1354_fac_control_audit import summarize
    base=ROOT/'data/audit/rya1354'
    report=summarize(base)
    assert report==json.loads((base/'source_curation_v2/fac_controls_v18.json').read_text())
    assert not report['cases']['cap_canary']['diagnostics']['printed_rule_verified']
    assert report['cases']['tight10']['diagnostics']['iterations']==257
    assert report['cases']['coherent10']['diagnostics']['iterations']==19
    assert report['cases']['coherent100']['diagnostics']['iterations']==20
    assert report['cases']['coherent100']['mixing_export']['levels']==15
    assert not report['physical_LS_terms_admitted'] and report['uncertainty'] is None


def test_real_zero_exit_capped_native_run_is_still_refused():
    import json
    from rya1354_fac_controls import diagnostics
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_controls_failed_v18')
    root='runs/fac_numerical_controls_v18/'
    receipt=json.loads(files[root+'receipt.json'])
    assert receipt['cases']['cap_canary']['returncode']==0
    assert root+'cap_canary/o.lev' in files
    assert not diagnostics(files[root+'cap_canary/sfac.stdout'].decode(),3.0,1)['printed_rule_verified']


def test_mixing_export_corruption_is_refused():
    from rya1354_fac_control_audit import mixing
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_coherent_complete_v18')
    raw=files['runs/fac_coherent_controls_v18/coherent100/o_export_000.cm']
    assert mixing(raw)['electrons']==8
    with pytest.raises(ValueError,match='record'):
        mixing(raw[:-1])


def test_changed_inner_tolerance_input_is_refused(monkeypatch):
    import rya1354_fac_control_audit as audit
    original=audit.load_parts
    def changed(path):
        files,manifest=original(path)
        if path.name=='fac_coherent_complete_v18':
            key=audit.COHERENT+'coherent100/input.sf'
            files[key]=files[key].replace(b"'orbital:eneabserr'",b"'orbital:enerelerr'")
        return files,manifest
    monkeypatch.setattr(audit,'load_parts',changed)
    with pytest.raises(ValueError,match='declared numerical controls differ'):
        audit.summarize(ROOT/'data/audit/rya1354')
