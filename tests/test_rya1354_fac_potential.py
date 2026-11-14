"""Fixed-space optimization policy and channel guard tests."""
from pathlib import Path
import sys
import json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_potential import input_text,rate_input
from rya1354_fac_space import input_text as ground_input
from rya1354_native_parts import load_parts


def test_only_declared_optimization_line_changes():
    for mult,scale in ((.3,10),(.03,100)):
        original=ground_input(2,mult,scale);text=input_text(mult,scale)
        assert text.replace("OptimizeRadial(['ground','odd3s'],[0.5,0.5])","OptimizeRadial(['ground'])")==original
        assert text.count('Config(')==15 and text.count('Structure(')==1


def test_rate_inventory_uses_mapping_and_balanced_prefix():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_expanded_angular_complete_v24')
    mapping=json.loads(files['evidence/fac_expanded_replay_v24/report.json'])['cases']['n3_r2_coherent100']
    for gauge in (1,2):
        text,requests,lower=rate_input(.03,100,gauge,mapping)
        assert text.startswith(input_text(.03,100)) and len(requests)==13 and lower==list(range(6))
        assert f'SetTransitionOptions({gauge},0)' in text


def test_actual_fixed_space_full_control_replay():
    from rya1354_fac_potential_audit import summarize
    r=summarize(ROOT/'data/audit/rya1354')
    assert r['scientific_admission']=='HOLD' and r['accepted_lifetime'] is None
    assert r['weight_default_discrepancy_preserved']
    for row in r['cases'].values():
        assert row['transport']['levels']==788 and row['transport']['declared_configurations_present']==15
        for rate in row['rates'].values():assert rate['channels']==13 and rate['printed_rows_verified']==12
    assert abs(r['cases']['coherent100']['rates']['B']['exploratory_model_lifetime_ns']-r['cases']['coherent10']['rates']['B']['exploratory_model_lifetime_ns'])<1e-5


def test_explicit_weights_bypass_documented_default_discrepancy():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_potential_complete_v26')
    root='runs/fac_potential_v26/source_documents/'
    radial=files[root+'faclib/radial.c'].decode();config=files[root+'faclib/config.c'].decode();manual=files[root+'doc/manual.tex'].decode()
    assert 'static int _acfg_wmode = 2;' in radial
    assert 'acfg->weight[i] += fabs(cfg->sweight);' in config
    assert 'acfg->weight[i] = weight[i];' in config
    assert 'configuration group is given an equal weight' in manual
    for mult,scale in ((.3,10),(.03,100)):
        assert "OptimizeRadial(['ground','odd3s'],[0.5,0.5])" in input_text(mult,scale)


def test_unlisted_temporary_packing_candidate_rejected(tmp_path):
    import pytest
    from rya1354_fac_potential_audit import verify_directory
    (tmp_path/'manifest.json').write_text(json.dumps({'parts':{'part001.tar.gz':{}}}))
    (tmp_path/'part001.tar.gz').write_bytes(b'part')
    verify_directory(tmp_path)
    (tmp_path/'.candidate_000.tar.gz').write_bytes(b'unlisted temporary packing candidate')
    with pytest.raises(ValueError,match='unlisted'):verify_directory(tmp_path)
