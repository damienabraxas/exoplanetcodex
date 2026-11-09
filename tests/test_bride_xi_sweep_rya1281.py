from tests.bride_test_requirements import requires_lightweaver
import json,copy
from pathlib import Path
import pytest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def inputs():
    return [json.loads((ROOT/'config/experiments'/name).read_text())for name in ('rya1281_mac_xi_slice.json','rya1281_mac_solar_bound_control.json')]


def test_axis_contract_and_extension_without_changing_fixed_fields():
    from pipeline.bride_xi_sweep import validate_plan
    plan,ref=inputs();validate_plan(plan,ref)
    extra=copy.deepcopy(plan['cases'][0]);extra.update(name='wider_xi050',xi=.5,bound_point_coefficients='wider_xi050')
    plan['cases'].append(extra);validate_plan(plan,ref)
    for key,value in (('teff',5782),('abundance',8.8),('drop_bottom',1)):
        bad=copy.deepcopy(plan);bad['cases'][0][key]=value
        with pytest.raises(ValueError,match='fixed inputs'):validate_plan(bad,ref)
    bad=copy.deepcopy(plan);bad['cases'][1]['xi']=.9
    with pytest.raises(ValueError,match='duplicate'):validate_plan(bad,ref)


def test_producer_xi_input_boundary_and_refusals():
    from pipeline.bride_xi_sweep import source_input_xi
    source="'OTHER:' 'value'\n'XIFIX:' 'T'\n1.000000000\n"
    assert float(source_input_xi(source,.9).splitlines()[-1])==.9
    for bad in (source.replace("'T'","'F'"),source+"trailer\n",source+source,source.replace('1.000000000','2.0')):
        with pytest.raises(ValueError):source_input_xi(bad,1.1)


@requires_lightweaver
def test_reconstructed_xi_requires_fresh_binding():
    from pipeline.bride_xi_sweep import reconstruct
    from pipeline.bride_opacity_binding import validate_binding,local_binding
    plan,_=inputs();low,ls,lc,_=reconstruct(plan['cases'][0]);mid,ms,mc,_=reconstruct(plan['cases'][1])
    for key in ('temperature_K','electron_pressure_Pa','rho_kg_m3','height_m'):
        np.testing.assert_array_equal(low[key],mid[key])
    np.testing.assert_array_equal(ls['column_mass_kg_m2'],ms['column_mass_kg_m2'])
    assert low['vturb_m_s']==900 and mid['vturb_m_s']==1000
    with pytest.raises(ValueError):validate_binding(local_binding(mid,8.7),low,8.7)


@requires_lightweaver
def test_fresh_coefficient_context_uses_xi_and_refuses_live_change():
    from config.constants import require_codex_path
    from pipeline.bride_xi_sweep import reconstruct
    from pipeline.bride_solar_control import bound_fe_factory
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_top_cap import FiniteColumnCap
    from pipeline.multi_atom import read_multi_atom
    from dataclasses import replace
    plan,_=inputs();case=plan['cases'][0];base,s,contract,_=reconstruct(case)
    atom=replace(read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')),abundance=case['abundance'])
    provider,_=bound_fe_factory(require_codex_path('experiments.bride_xi_coefficients')/case['name'],base,s,atom.abundance,column_contract=contract)
    ctx=make_reference_context(s,atom,rays=2,qscale=1,fe_background=provider,experimental_column_cap=FiniteColumnCap(s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]))
    np.testing.assert_array_equal(ctx.atmos.pyAtmos.vturb,np.full(56,900.))
    ctx.atmos.pyAtmos.vturb[0]=np.nextafter(900.,np.inf)
    with pytest.raises(ValueError,match='geometry/velocity'):ctx.background.update_background(ctx.atmos)
