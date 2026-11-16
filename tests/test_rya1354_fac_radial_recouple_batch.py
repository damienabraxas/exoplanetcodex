from pathlib import Path
from types import SimpleNamespace
import json
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_radial_recouple_batch import command,CASES


def test_remaining_case_commands_preserve_qualified_child_and_exact_cli():
    a=SimpleNamespace(**{k:Path('/native')/k for k in ('grasp','first','identity','remaining','old_converted','canary','out')})
    assert CASES==('ground_coherent100','equal_coherent10','equal_coherent100')
    for case in CASES:
        cmd=command(a,case)
        assert Path(cmd[1]).name=='rya1354_fac_radial_recouple.py'
        assert '--old_converted' in cmd and '--old-converted' not in cmd
        assert cmd[cmd.index('--case')+1]==case
        assert cmd[cmd.index('--out')+1]==str(a.out/case)
    with pytest.raises(ValueError,match='Undeclared'):command(a,'ground_coherent10')


def test_all_four_actual_ls_cases_have_complete_signed_replay_and_separate_target_identity():
    from rya1354_fac_radial_angular_remaining_metadata import summarize
    report=summarize(ROOT/'data/audit/rya1354')
    assert report==json.loads((ROOT/'data/audit/rya1354/source_curation_v2/fac_radial_angular_v40.json').read_text())
    assert report['scientific_admission']=='HOLD' and not report['adoption_allowed']
    assert set(report['cases'])=={'ground_coherent10',*CASES}
    assert sum(v['full_LS_coefficients_reconstructed'] for v in report['cases'].values())==17417816
    for case in report['cases'].values():
        assert case['levels']==6600 and case['label_ID_overflows']==5600
        assert all(len(ids)==1 for ids in case['targets'].values())
        assert case['target_metadata']['3S']['two_J']==2 and case['target_metadata']['5S']['two_J']==4
        assert all(v['parity']=='-' for v in case['target_metadata'].values())
