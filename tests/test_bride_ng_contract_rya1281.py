import copy,json
from pathlib import Path
import pytest
from scripts.run_bride_solar_matched_mali import validate_solver_controls


def test_explicit_ng_control_cannot_relax_physics_thresholds_or_compute_budget():
    root=Path(__file__).resolve().parents[1]
    warm=json.loads((root/'config/experiments/rya1281_solar_matched_columns.json').read_text())['controls']
    plan=json.loads((root/'config/experiments/rya1281_solar_matched_transition_no_ng.json').read_text())
    validate_solver_controls(plan,warm)
    for key,value in [('rtol',1e-5),('seconds_per_case',601),('cap_radiation_rtol',1e-7),('rays',3)]:
        bad=copy.deepcopy(plan);bad['controls'][key]=value
        with pytest.raises(ValueError):validate_solver_controls(bad,warm)
    bad=copy.deepcopy(plan);bad['numerical_algorithm']['method']='unknown'
    with pytest.raises(ValueError):validate_solver_controls(bad,warm)
