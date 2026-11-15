from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_radial_ci_v36 import input_text,CASES
from rya1354_fac_radial_ci_v35 import input_text as qualified
from rya1354_fac_radial_preflight import input_text as radial


def test_qualified_first_input_is_unchanged():
    assert input_text(10)==qualified()


def test_tighter_input_changes_only_declared_inner_controls():
    assert input_text(100).removeprefix(radial(100))==input_text(10).removeprefix(radial(10))
    assert input_text(100).startswith(radial(100))
    assert CASES==(('ground_coherent100',100),('equal_coherent10',10),('equal_coherent100',100))


def test_all_four_actual_cases_replay_with_original_failures_and_physical_invariants():
    from rya1354_fac_radial_ci_remaining_audit import summarize
    report=summarize(ROOT/'data/audit/rya1354')
    assert report['retained_v33_status']==report['retained_v34_status']=='FAILED_RETAINED'
    assert report['scientific_admission']=='HOLD' and not report['adoption_allowed']
    assert set(report['cases'])=={'ground_coherent10','ground_coherent100','equal_coherent10','equal_coherent100'}
    counts={'ground_coherent10':4118356,'ground_coherent100':4118342,'equal_coherent10':4078220,'equal_coherent100':4078216}
    for name,case in report['cases'].items():
        assert case['post_CI_orbitals_unchanged'] and case['potential']['protected_header_grid_potential_byte_identical']
        assert case['transport']['levels']==case['transport']['CSFs']==6600
        assert case['transport']['signed_printed_coefficients_verified']==counts[name]
        assert not case['transport']['LS_mapping_admitted']
        assert all(v['groups']==114 and v['max_abs_delta']==v['max_abs_shift']==0 for v in case['configuration_state'].values())
