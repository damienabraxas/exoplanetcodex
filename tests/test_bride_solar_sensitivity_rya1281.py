"""Check physical field isolation and source immutability for the declared sweep."""
import copy,json
from pathlib import Path
import numpy as np
import pytest


def test_held_sensitivity_changes_xi_without_changing_thermal_density_or_geometry():
    pytest.importorskip('lightweaver')
    from config.constants import require_codex_path
    from pipeline.bride_solar_matched import prepare_columns,column_state
    plan=json.loads((Path(__file__).resolve().parents[1]/'config/experiments/rya1281_solar_sensitivity.json').read_text());before=copy.deepcopy(plan)
    columns,_=prepare_columns(plan,require_codex_path('grids.marcs_standard'))
    solar=column_state(columns['solar'])
    for case in plan['cases']:
        if not case['name'].startswith('xi_'):continue
        state=column_state(columns[case['name']])
        for key in ['height_m','temperature_K','electron_pressure_Pa','ne_m3','nH_m3']:
            np.testing.assert_array_equal(getattr(state,key),getattr(solar,key))
        np.testing.assert_array_equal(state.vturb_m_s,np.full(56,case['xi_km_s']*1000))
    assert plan==before
    assert columns['T_5762']['temperature_K'][0]<columns['solar']['temperature_K'][0]<columns['T_5782']['temperature_K'][0]
    bad=copy.deepcopy(plan);bad['cases'][-1]['xi_km_s']=-.1
    with pytest.raises(ValueError):prepare_columns(bad,require_codex_path('grids.marcs_standard'))
