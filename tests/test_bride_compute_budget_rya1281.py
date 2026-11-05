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
