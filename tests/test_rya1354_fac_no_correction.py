"""Explicit disabling, read-only diagnostics and model-policy isolation."""
from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_no_correction import input_text,config_state
from rya1354_fac_space import input_text as ground_input
from rya1354_fac_potential import input_text as equal_input


@pytest.mark.parametrize('policy',['ground','equal'])
def test_disabled_changes_only_explicit_switches_and_read_only_diagnostics(policy):
    for mult,scale in ((.3,10),(.03,100)):
        text=input_text(policy,mult,scale)
        stripped=text.replace("SetConfigEnergyMode(-1)\nSetOption('radial:config_energy',-1)\n",'').replace("ConfigDiagnostics('configs.before.tsv')\n",'').replace("ConfigDiagnostics('configs.after.tsv')\n",'')
        expected=ground_input(2,mult,scale) if policy=='ground' else equal_input(mult,scale)
        assert stripped==expected and text.count('Config(')==15
        assert 'ConfigEnergy(0)' in text and 'ConfigEnergy(1)' in text


def test_actual_configuration_values_required_not_absent_print_lines():
    zero='# FAC_CONFIG_V1 mode=-1 groups=1\n# group config nele energy delta shift name\n0 0 8 -74 0 0 ground\n'
    r=config_state(zero,True);assert r['max_abs_delta']==r['max_abs_shift']==0
    for bad in (zero.replace('-74 0 0','-74 1e-9 0'),zero.replace('-74 0 0','-74 0 1e-9'),zero.replace('mode=-1','mode=0'),zero.split('0 0 8')[0]):
        with pytest.raises(ValueError):config_state(bad,True)
    enabled=zero.replace('mode=-1','mode=0').replace('-74 0 0','-74 -0.1 0')
    assert config_state(enabled,False)['max_abs_delta']==.1


def test_build_helper_refuses_mac_before_touching_output(tmp_path):
    import platform
    from rya1354_fac_config_probe import build
    if platform.node().lower().split('.')[0]=='sirius':pytest.skip('Mac policy guard')
    with pytest.raises(RuntimeError,match='Sirius'):build(tmp_path/'source',tmp_path/'probe',tmp_path/'evidence')
    assert not (tmp_path/'probe').exists()


def test_actual_all_four_disabled_cases_and_enabled_positive_control():
    from rya1354_fac_no_correction_audit import summarize
    r=summarize(ROOT/'data/audit/rya1354')
    assert len(r['cases'])==4 and r['accepted_lifetime'] is None
    assert r['enabled_auxiliary_optimization_convergence_qualified'] is False
    for state in r['enabled_minimal_diagnostic_canary'].values():assert state['mode']==0 and state['max_abs_delta']>0
    for case in r['cases'].values():
        assert case['transport']['levels']==788
        for s in case['configuration_state'].values():assert s['mode']==-1 and s['max_abs_delta']==s['max_abs_shift']==0
        for rate in case['rates'].values():
            assert rate['audit']['channels']==13
            assert rate['audit']['printed_rows_verified']+len(rate['audit']['writer_suppressed'])==13
            for s in rate['configuration_state'].values():assert s['max_abs_delta']==s['max_abs_shift']==0


def test_missing_final_configuration_rejected_against_basis_inventory():
    from rya1354_fac_no_correction_audit import verify_config_inventory
    text='# FAC_CONFIG_V1 mode=-1 groups=1\n0 0 8 -74 0 0 ground\n'
    with pytest.raises(ValueError,match='inventory'):verify_config_inventory(text,{(0,0),(0,1)})
