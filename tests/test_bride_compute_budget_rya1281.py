import copy,json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]


def test_only_declared_cost_cap_changes_and_original_is_preserved():
    from pipeline.bride_compute_budget import extended_controls
    from pipeline.bride_thermal_sweep import effective_controls
    plan=json.loads((ROOT/'config/experiments/rya1281_mac_thermal_slice.json').read_text())
    declaration=json.loads((ROOT/'config/experiments/rya1281_ng_extended.json').read_text())
    original=effective_controls(plan['controls'],True);saved=copy.deepcopy(original)
    result=extended_controls(original,declaration,'thermal_g4338','record')
    assert original==saved and result['seconds_per_case']==1200
    assert {k for k in original if original[k]!=result[k]}=={'seconds_per_case'}
    assert result['maxiter']==1500 and result['rtol']==1e-6
    for key,value in [('rtol',.01),('maxiter',3000)]:
        bad=dict(declaration,**{key:value})
        with pytest.raises(ValueError):extended_controls(original,bad,'thermal_g4338','record')
    for seconds in (True,600,0,1200.0):
        bad=dict(declaration,seconds_per_case=seconds)
        with pytest.raises(ValueError):extended_controls(original,bad,'thermal_g4338','record')
    with pytest.raises(ValueError):extended_controls(original,declaration,'thermal_T5572','record')
    with pytest.raises(ValueError):extended_controls(original,declaration,'thermal_g4338','off')


def test_comparison_requires_explicit_cost_extension_and_retains_accuracy_guards():
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_compute_budget import extended_controls
    from pipeline.bride_thermal_sweep import effective_controls
    from scripts.analyze_bride_thermal_refinement import verify_pair
    plan=json.loads((ROOT/'config/experiments/rya1281_mac_thermal_slice.json').read_text())
    path=ROOT/'config/experiments/rya1281_ng_extended.json';declaration=json.loads(path.read_text())
    base=json.loads((ROOT/'data/results/rya1281/thermal_slice_v1/runs/thermal_g4338/report.json').read_text())
    ref=copy.deepcopy(base);ref.update(numeric_refinement=True,ng_control='record',
        generator=fingerprint(ROOT/'scripts/run_bride_thermal_sweep.py'),compute_extension=fingerprint(path),
        budget_adapter=fingerprint(ROOT/'pipeline/bride_compute_budget.py'),
        initialization=dict(mode='matched_physical_baseline',populations=base['coupled']['populations'],
            helper=fingerprint(ROOT/'pipeline/bride_thermal_initial.py')))
    ref['coupled']['controls']=extended_controls(effective_controls(plan['controls'],True),declaration,'thermal_g4338','record')
    with pytest.raises(ValueError,match='explicitly'):verify_pair(base,ref,plan,ng_control='record')
    verify_pair(base,ref,plan,ng_control='record',extended_budget=True)
    for key,value in [('seconds_per_case',1500),('rtol',1e-5),('maxiter',3000)]:
        bad=copy.deepcopy(ref);bad['coupled']['controls'][key]=value
        with pytest.raises(ValueError,match='budgets'):verify_pair(base,bad,plan,ng_control='record',extended_budget=True)
