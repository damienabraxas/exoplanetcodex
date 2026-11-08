import copy,json
from pathlib import Path
import numpy as np
import pytest


def test_departure_pattern_maps_to_local_lte_without_changing_sources():
 from pipeline.bride_solar_matched import departure_seed
 p=np.array([[2.,4.],[3.,6.]]);s=np.array([[1.,2.],[4.,8.]]);t=np.array([[2.,1.],[7.,3.]]);before=[a.copy()for a in [p,s,t]];out=departure_seed(p,s,t,[20.,30.]);np.testing.assert_allclose(out.sum(axis=0),[20.,30.]);actual=out/t;expected=p/s;np.testing.assert_allclose(actual/actual.sum(axis=0),expected/expected.sum(axis=0))
 for x,old in zip([p,s,t],before):np.testing.assert_array_equal(x,old)
 with pytest.raises(ValueError):departure_seed(p,s*0,t,[20.,30.])


def test_held_four_columns_keep_explicit_owners_and_all_native_depths():
 pytest.importorskip('lightweaver')
 from config.constants import require_codex_path,get_star_params
 from pipeline.bride_solar_matched import prepare_columns,column_state
 plan=json.loads((Path(__file__).resolve().parents[1]/'config/experiments/rya1281_solar_matched_columns.json').read_text());columns,receipt=prepare_columns(plan,require_codex_path('grids.marcs_standard'))
 solar=get_star_params('solar');assert plan['cases'][-1]['teff']==solar['teff']and plan['cases'][-1]['logg']==solar['logg']
 for name,s in columns.items():
  before=copy.deepcopy(s);state=column_state(s);assert len(state.nH_m3)==56 and state.case=='model_mass_H'and np.all(np.array(state.vturb_m_s)==1000)
  assert 'native'not in state.ne_owner and 'network'not in state.nH_owner
  for key in ['height_m','temperature_K','electron_pressure_Pa','rho_kg_m3']:np.testing.assert_array_equal(s[key],before[key])
  if name!='baseline':assert 'column_mass_kg_m2'not in s
 bad=copy.deepcopy(columns['baseline']);bad['rho_kg_m3'][0]=np.nan
 with pytest.raises(ValueError):column_state(bad)
