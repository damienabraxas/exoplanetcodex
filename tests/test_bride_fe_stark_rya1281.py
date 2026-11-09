from tests.bride_test_requirements import requires_lightweaver
import json,copy
import numpy as np
import pytest


def row():return dict(line_id='controlled',wavelength_air_A=7777.86,excitation_eV=4.,ionization_eV=7.902,log_gamma_stark_10000=-5.)


@requires_lightweaver
def test_native_stark_cases_and_refusals():
    from pipeline.bride_fe_stark import stark_gamma
    from config.constants import require_codex_path
    r=json.loads((require_codex_path('experiments.bride_mac_grid')/'rya1281-fe-stark-source-v2/report.json').read_text())
    rows={v['line_id']:v for v in r['rows']}
    assert len(r['cases'])==738
    for case in r['cases']:
        got=stark_gamma(np.array([case['temperature_K']]),np.array([case['ne_m3']]),rows[case['line_id']],case['mode'])[0]
        expected=case['native_gamma_s_inv']
        if expected:assert abs(got/expected-1)<128*np.finfo(np.float32).eps
        else:assert got==0
    for mode in ('provided','cowley'):
        np.testing.assert_array_equal(stark_gamma(np.array([4000.,10000.]),np.zeros(2),row(),mode),0)
    for mode in ('bad',):
        with pytest.raises(ValueError):stark_gamma(np.array([5000.]),np.ones(1),row(),mode)
    bad=row();bad['log_gamma_stark_10000']=0.
    with pytest.raises(ValueError):stark_gamma(np.array([5000.]),np.ones(1),bad,'provided')
    bad=row();bad['excitation_eV']=10.
    with pytest.raises(ValueError):stark_gamma(np.array([5000.]),np.ones(1),bad,'cowley')


@requires_lightweaver
def test_actual_stark_none_identity_and_changed_ne_refusal():
    from tests.test_bride_solar_control_rya1281 import inputs
    from tests.test_bride_trace_fe_lines_rya1281 import line
    from tests.test_bride_fe_pressure_rya1281 import recipe
    from pipeline.bride_solar_control import bound_fe_factory
    from pipeline.bride_fe_pressure import pressure_fe_factory
    from pipeline.bride_fe_stark import stark_fe_factory
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_top_cap import FiniteColumnCap
    folder,base,s,atom,contract=inputs('fine_T5772');provider,_=bound_fe_factory(folder,base,s,atom.abundance,column_contract=contract)
    def ctx(factory):return make_reference_context(s,atom,rays=2,qscale=1,fe_background=factory,
        extra_wavelengths_vac_nm=np.array([777.5,778.,778.5]),experimental_column_cap=FiniteColumnCap(s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]))
    original=ctx(pressure_fe_factory(provider,[line()],(777.,779.),[recipe()]))
    zero=ctx(stark_fe_factory(provider,[line()],(777.,779.),[recipe()],[row()],'none'))
    for field in ('chi','eta','sca'):np.testing.assert_array_equal(getattr(original.background,field),getattr(zero.background,field))
    given=ctx(stark_fe_factory(provider,[line()],(777.,779.),[recipe()],[row()],'provided'))
    np.testing.assert_array_equal(given.background.sca,original.background.sca)
    clone=copy.deepcopy(given);np.testing.assert_array_equal(clone.background.chi,given.background.chi)
    ne=given.atmos.pyAtmos.ne;ne[0]=np.nextafter(ne[0],np.inf)
    with pytest.raises(ValueError):given.background.update_background(given.atmos)


def test_h2_full_domain_refuses_uncovered_and_invalid_temperatures():
    from scripts.audit_bride_h2_domain import require_full_domain,domain_mask
    require_full_domain(np.array([1000.,9000.]),[1000.,9000.])
    with pytest.raises(ValueError,match='complete column'):require_full_domain(np.array([5000.,np.nextafter(9000.,np.inf)]),[1000.,9000.])
    with pytest.raises(ValueError):domain_mask(np.array([np.nan]),[1000.,9000.])
