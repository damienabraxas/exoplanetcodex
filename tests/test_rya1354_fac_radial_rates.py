from pathlib import Path
import json
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_radial_rates import rate_input,source_path
from rya1354_fac_radial_ci_v36 import input_text


def mapping():
    return dict(targets={'3S':[6],'5S':[5]},levels={
        '0':dict(energy_au=-10.,two_J=4,parity='+'),
        '5':dict(energy_au=-9.,two_J=4,parity='-'),
        '6':dict(energy_au=-8.,two_J=2,parity='-'),
        '7':dict(energy_au=-7.,two_J=2,parity='+')})


def test_rate_extension_preserves_radial_protocol_and_complete_dynamic_operators():
    for scale in (10,100):
        for gauge in (1,2):
            text,requests,lower=rate_input(scale,gauge,mapping())
            assert text.startswith(input_text(scale))
            assert text.count('Structure(')==1
            assert f'SetTransitionOptions({gauge},0)\nSetTransitionCut(0,0)' in text
            assert lower==[0,5]
            assert {(r['lower'],r['operator']) for r in requests}=={
                (0,'E1'),(0,'M2'),(0,'E3'),(5,'M1'),(5,'E2'),(5,'M3')}
    # A newly lower state must be included; IDs/counts are not copied from n3.
    m=mapping();m['levels']['7']['energy_au']=-8.5
    _,requests,lower=rate_input(10,2,m)
    assert lower==[0,5,7]
    assert {r['operator'] for r in requests if r['lower']==7}=={'E1','M2'}


def test_rate_controls_and_first_native_source_are_explicit():
    for scale,gauge in [(1,2),(10,3)]:
        with pytest.raises(ValueError,match='Undeclared'):rate_input(scale,gauge,mapping())
    root=Path('/native')
    assert source_path(root,'ground_coherent10')==root/'runs/fac_radial_ci_first_v33/ground_coherent10'
    assert source_path(root,'equal_coherent100')==root/'runs/fac_radial_ci_remaining_v36/equal_coherent100'


def test_radial_review_rejects_rate_input_with_changed_physics(tmp_path):
    from rya1354_fac_radial_ci_v36 import review
    text=rate_input(10,2,mapping())[0].replace('SetBreit(0)','SetBreit(1)')
    (tmp_path/'input.sf').write_text(text)
    with pytest.raises(ValueError,match='Declared radial CI input differs'):
        review(tmp_path,{},10,declared_input=text)


def test_portable_basis_reader_keeps_decorated_native_section_boundary():
    from rya1354_fac_radial_rates_audit import basis_keys
    data=b'========Basis Table========\n0 0 0 0 2 3 x y z q\n============Mixing Coefficients===================\nnot a basis row\n'
    assert basis_keys(data)=={(2,3)}
    with pytest.raises(ValueError,match='boundary missing'):basis_keys(data.split(b'Mixing Coefficients')[0])


def test_actual_all_eight_radial_rates_have_complete_independent_accounting():
    from rya1354_fac_radial_rates_audit import summarize
    report=summarize(ROOT/'data/audit/rya1354')
    assert report==json.loads((ROOT/'data/audit/rya1354/source_curation_v2/fac_radial_rates_v41.json').read_text())
    assert len(report['cases'])==8 and report['scientific_admission']=='HOLD'
    assert not report['adoption_allowed'] and report['accepted_lifetime'] is None
    assert sum(case['audit']['channels'] for case in report['cases'].values())==104
    for case in report['cases'].values():
        audit=case['audit']
        assert audit['printed_rows_verified']+len(audit['writer_suppressed'])==audit['channels']
        assert case['lower_states']==list(range(6))
        assert audit['exploratory_model_total_A_s1']>0
        for cfg in case['configuration_state'].values():
            assert cfg['groups']==114 and cfg['mode']==-1
            assert cfg['max_abs_delta']==cfg['max_abs_shift']==0
