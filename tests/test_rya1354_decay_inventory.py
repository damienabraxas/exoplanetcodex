"""Missing or duplicate weak channels must not yield a complete model lifetime."""
import json
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_lifetime_validation import evaluate,model_lifetimes
from rya1354_evaluated_4d import evaluate as evaluate_4d
from rya1354_two_method_summary import summarize
from rya1354_prospective_screen import screen
from rya1354_solver_smoke import run
from rya1354_rci_summary import summarize as rci_summary
from rya1354_forbidden_summary import summarize as forbidden_summary,parse as parse_forbidden


def test_fixed_wavefunctions_replay_and_experiment_remains_nonblind(tmp_path):
    out=tmp_path/'lifetime.json'
    evaluate(out)
    held=ROOT/'data/audit/rya1354/source_curation_v2/lifetime_validation.json'
    assert out.read_bytes()==held.read_bytes()
    d=json.loads(out.read_text())
    assert len(d['channels'])==13
    assert d['experiment']['value_ns']==1.79
    assert d['chi_square'] is None and d['calibrated_sigma'] is None
    assert not d['adoption_allowed']
    for g in ('B','C'):
        assert sum(r[f'model_branching_fraction_{g}'] for r in d['channels'])==pytest.approx(1)


def test_missing_zero_or_weak_channel_and_duplicate_refuse_lifetime():
    d=json.loads((ROOT/'data/audit/rya1354/source_curation_v2/lifetime_validation.json').read_text())
    rows=d['channels']
    zero=next(r for r in rows if r['computed_zero_in_this_model'])
    with pytest.raises(ValueError,match='decay channel'):
        model_lifetimes([r for r in rows if r is not zero])
    with pytest.raises(ValueError,match='decay channel'):
        model_lifetimes(rows+[rows[0].copy()])


def test_all_physical_4d_components_replay_without_forced_grid_aliases(tmp_path):
    out=tmp_path/'4d.json'
    evaluate_4d(out)
    assert out.read_bytes()==(ROOT/'data/audit/rya1354/source_curation_v2/evaluated_4d_validation.json').read_bytes()
    d=json.loads(out.read_text())
    assert len({r['physical_reference_id'] for r in d['rows']})==9
    assert all(r['experimental_accuracy_score'] is None for r in d['rows'])
    assert all(r['grid_alias_status'].startswith('NOT_ASSIGNED') for r in d['rows'])


def test_two_method_capsules_replay_and_reject_first_transition_attempt(tmp_path):
    out=tmp_path/'methods.json'
    summarize(out)
    held=ROOT/'data/audit/rya1354/source_curation_v2/two_method_ladder.json'
    assert out.read_bytes()==held.read_bytes()
    d=json.loads(out.read_text())
    assert [r['n'] for r in d['rows']]==[3,4,5]
    assert all(r['same_occupation_multiset'] for r in d['basis_checks'])
    assert d['native_control']['byte_identical']
    assert 'no .ls' in d['rejected_attempt']['reason']
    assert d['uncertainty'] is None


@pytest.mark.parametrize('message,marker',[
    ('The order of 3s and 3p is inconsistent','INCONSISTENT_ORBITAL_ORDER'),
    ('Maximum iterations in SCF Exceeded.','MAXIMUM ITERATIONS IN SCF EXCEEDED')])
def test_native_zero_exit_orbital_order_refusal_remains_failure(tmp_path,message,marker):
    executable=tmp_path/'native_stub'
    executable.write_text(f'#!/bin/sh\nprintf "{message}\\n"\nexit 0\n')
    executable.chmod(0o700)
    receipt=run(executable,'',tmp_path)
    assert receipt['returncode']==0
    assert marker in receipt['semantic_errors']


def test_no_active_feed_is_not_called_never_processed_or_blind(tmp_path):
    out=tmp_path/'prospective.json'
    screen(out)
    held=ROOT/'data/audit/rya1354/source_curation_v2/prospective_screen.json'
    assert out.read_bytes()==held.read_bytes()
    d=json.loads(out.read_text())
    assert len(d['rows'])==28
    assert {r['target'] for r in d['selected']}=={'Li','P','S'}
    assert all(not r['blind'] and not r['primary_experiments_verified'] for r in d['selected'])
    assert all(r['active_feed_products']>0 for r in d['rows'] if r['target'] in ('Na','Mg','Fe II'))


def test_rci_warning_control_and_sensitivity_replay(tmp_path):
    out=tmp_path/'rci.json'
    rci_summary(out)
    assert out.read_bytes()==(ROOT/'data/audit/rya1354/source_curation_v2/rci_sensitivities.json').read_bytes()
    d=json.loads(out.read_text())
    assert set(d['cases'])=={'DC','Breit_static_limit','VP_self_energy','CV_2s_one_hole_DC'}
    assert d['uncertainty'] is None
    assert 'byte-identical' in d['numerical_warning']['positive_control']


def test_native_air_column_is_rejected_when_it_disagrees_with_ssot():
    import numpy as np
    from pipeline.wavelength_util import vac_to_air
    d=json.loads((ROOT/'data/audit/rya1354/source_curation_v2/two_method_ladder.json').read_text())
    for row in d['rows']:
        native=row['mchf']
        qualified=float(vac_to_air(np.array([native['native_vacuum_A']]))[0])
        assert abs(qualified-native['printed_air_A_rejected'])>1
        assert 'not used' in native['air_warning']


def test_forbidden_full_operator_inventory_replays_and_deduplicates(tmp_path):
    out=tmp_path/'forbidden.json'
    forbidden_summary(out)
    assert out.read_bytes()==(ROOT/'data/audit/rya1354/source_curation_v2/forbidden_convergence.json').read_bytes()
    d=json.loads(out.read_text())
    assert [r['n'] for r in d['layers']]==[3,4,5]
    assert not d['ladder_qualified_for_convergence']
    assert d['layers'][1]['numerical_status']=='ITERATION_LIMIT_HIT_UNQUALIFIED'
    assert d['layers'][2]['numerical_status']=='FROZEN_SEED_FROM_UNCONVERGED_N4_UNQUALIFIED'
    for layer in d['layers']:
        assert len(layer['channels'])==16
        assert {r['operator'] for r in layer['channels']}=={'M1','E2','M3','E4'}
        assert len(layer['reference_comparisons'])==4
    assert d['uncertainty'] is None


def test_forbidden_missing_or_disagreeing_duplicate_refuses_summary():
    import tarfile
    import re
    with tarfile.open(ROOT/'data/audit/rya1354/forbidden_ladder_v1.tar.gz') as archive:
        text=archive.extractfile('runs/forbidden_ladder_v1/n3_multipoles/lower.upper.t').read().decode()
    first=next(line for line in text.splitlines() if re.match(r'\s*f1\s+\d',line))
    with pytest.raises(ValueError):
        parse_forbidden(text.replace(first+'\n','',1))
    # Alter one directional partner only: unequal duplicate rates must not be averaged.
    tokens=first.split()
    changed=first.replace(tokens[10],'9.87654D+01',1)
    with pytest.raises(ValueError,match='Duplicate directions disagree'):
        parse_forbidden(text.replace(first,changed,1))


def test_forbidden_resumed_native_stop_does_not_claim_orbital_convergence(tmp_path):
    out=tmp_path/'resumed.json'
    forbidden_summary(out,qualified=True)
    assert out.read_bytes()==(ROOT/'data/audit/rya1354/source_curation_v2/forbidden_convergence_v2.json').read_bytes()
    d=json.loads(out.read_text())
    assert d['native_stopping_sequence_completed']
    assert not d['full_orbital_residual_convergence_validated']
    assert not d['ladder_qualified_for_convergence']
    assert [r['native_SCF_diagnostics']['iterations'] for r in d['layers']]==[9,44,129]
    for layer in d['layers']:
        diagnostics=layer['native_SCF_diagnostics']
        assert not diagnostics['iteration_limit_hit']
        assert not diagnostics['all_reported_functions_below_ACC']
        assert diagnostics['max_reported_self_consistency']>diagnostics['native_function_ACC']
        assert len(layer['channels'])==16
    assert not d['adoption_allowed']
