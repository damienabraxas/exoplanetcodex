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


def test_reconstructed_xi_requires_fresh_binding():
    from pipeline.bride_xi_sweep import reconstruct
    from pipeline.bride_opacity_binding import validate_binding,local_binding
    plan,_=inputs();low,ls,lc,_=reconstruct(plan['cases'][0]);mid,ms,mc,_=reconstruct(plan['cases'][1])
    for key in ('temperature_K','electron_pressure_Pa','rho_kg_m3','height_m'):
        np.testing.assert_array_equal(low[key],mid[key])
    np.testing.assert_array_equal(ls['column_mass_kg_m2'],ms['column_mass_kg_m2'])
    assert low['vturb_m_s']==900 and mid['vturb_m_s']==1000
    with pytest.raises(ValueError):validate_binding(local_binding(mid,8.7),low,8.7)
