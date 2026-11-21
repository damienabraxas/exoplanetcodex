from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_radial_ls_audit import ls_blocks,review,N3,EXTENDED


def test_angular_math_guard_precedes_any_data_or_numpy_operation(monkeypatch):
    monkeypatch.setattr('platform.node',lambda:'Mac')
    with pytest.raises(RuntimeError,match='Sirius'):review({}, {})


def test_six_shell_metadata_accepts_closed_4s_and_rejects_undeclared_shell():
    text='1s(2) 2s(2) 3s(2) 4s(2)\n1S0 1S0 1S0 1S0 1S0 1S0 1S0\n*\n'
    assert len(ls_blocks(text,EXTENDED)[0])==1
    with pytest.raises(ValueError,match='Undeclared'):ls_blocks(text,N3)
    with pytest.raises(ValueError,match='Undeclared'):ls_blocks(text.replace('4s','5s'),EXTENDED)


def test_fast_canary_requires_exact_case_and_corruption_coverage():
    from rya1354_fac_radial_recouple import qualified_canary,CASES
    r=dict(vector_control_location=dict(field='actual_first_vector_coefficient',original='0.25000000',mutated='-0.25000000'),schema='rya1354.fac-ls-fast-canary.v39',status='EXACT_LEGACY_EQUIVALENCE_AND_NEGATIVE_CONTROLS_PASS',cases={n:dict(all_nondiagnostic_report_fields_exact=True,full_LS_coefficients=329198) for n in CASES},negative_controls={'d4_seniority':'Unsupported NK d-shell seniority','vector_sign':'Expanded LS coefficient/sign replay differs','matrix_sign':'Expanded LS coefficient/sign replay differs','explicit_stdout_ID':'Explicit stdout/binary ID order differs'},forced_proposal_recovery=dict(pass_exact=True))
    qualified_canary(r)
    r['negative_controls'].pop('matrix_sign')
    with pytest.raises(ValueError,match='unqualified'):qualified_canary(r)


def test_actual_vector_mutation_leaves_all_energy_identity_lines_unchanged():
    from rya1354_native_parts import load_parts
    from rya1354_fac_ls_fast_canary_v39 import flip_actual_vector_coefficient
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_full_angular_complete_v29')
    raw=files['runs/fac_full_recoupling_v29/ground_coherent100/state.lsj.j'].decode()
    changed,location=flip_actual_vector_coefficient(raw)
    a=raw.splitlines(keepends=True);b=changed.splitlines(keepends=True)
    assert len(a)==len(b)
    assert [i for i,(x,y) in enumerate(zip(a,b)) if x!=y]==[location['line']]
    import re
    assert re.fullmatch(r'\s*(?:[-+]?\d+\.\d{8}\s*)+',a[location['line']])
    assert float(location['original'])==-float(location['mutated']) and abs(float(location['original']))>.01


def test_completed_first_6600_level_angular_evidence_has_corrected_controls_and_native_ids():
    from rya1354_fac_radial_angular_metadata import summarize
    report=summarize(ROOT/'data/audit/rya1354')
    assert report['levels']==6600 and report['full_LS_coefficients_reconstructed']==4354454
    assert report['label_ID_overflows']==5600 and report['targets']=={'3S':[6],'5S':[5]}
    assert report['canary']['negative_controls']['vector_sign']=='Expanded LS coefficient/sign replay differs'
    assert report['prior_invalid_vector_control_rejection']=='Missing LS vector identity'
    assert not report['adoption_allowed'] and report['accepted_lifetime'] is None
