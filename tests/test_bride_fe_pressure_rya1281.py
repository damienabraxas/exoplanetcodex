from tests.bride_test_requirements import requires_lightweaver
import copy,json
import numpy as np
import pytest


def recipe():return dict(line_id='controlled',recipe='W',sigma_au=0.,alpha=0.,log_gamma6_10000=-7.5)


@requires_lightweaver
def test_pressure_units_density_and_temperature_scaling_and_native_oracle():
    from pipeline.bride_fe_pressure import pressure_gamma
    from config.constants import require_codex_path
    T=np.array([4000.,10000.]);H=np.array([1e23,1e23]);He=H*.1;H2=H*.001
    expected=10**(-7.5)*(T/10000.)**.3*(H+.41336*He+.85*H2)*1e-6
    np.testing.assert_allclose(pressure_gamma(T,H,He,H2,recipe()),expected,rtol=2e-15)
    np.testing.assert_array_equal(pressure_gamma(T,H*0,He*0,H2*0,recipe()),0)
    for density in [H*-1,np.full(2,np.nan)]:
        with pytest.raises(ValueError):pressure_gamma(T,density,He,H2,recipe())
    folder=require_codex_path('experiments.bride_mac_grid')/'rya1281-fe-pressure-source-v2'
    source=json.loads((folder/'report.json').read_text())
    for r in source['recipes']:
        actual=pressure_gamma(np.array([5772.]),np.array([1e23]),np.array([1e22]),np.array([1e20]),r)[0]
        assert abs(actual/r['native_reference_gamma_s_inv']-1)<128*np.finfo(np.float32).eps


@requires_lightweaver
def test_depth_width_kernel_equals_independent_single_depths_and_refuses_shape():
    from pipeline.bride_lte_blend import lte_line_coefficients
    wave=np.linspace(777.9,778.1,21);T=np.array([4000.,6000.]);ions=np.array([1e15,2e15]);xi=np.array([1000.,1000.])
    kwargs=dict(lambda0_nm=778.,mass_amu=55.845,excitation_eV=4.,loggf=-2.)
    widths=np.array([1e8,1e10]);result=lte_line_coefficients(wave,T,ions,xi,gamma_s_inv=widths,**kwargs)
    for j in range(2):
        one=lte_line_coefficients(wave,T[j:j+1],ions[j:j+1],xi[j:j+1],gamma_s_inv=widths[j],**kwargs)
        for array,single in zip(result,one):np.testing.assert_array_equal(array[:,j],single[:,0])
    for width in [np.array([1e8]),np.array([[1e8,1e8]])]:
        with pytest.raises(ValueError,match='Width'):lte_line_coefficients(wave,T,ions,xi,gamma_s_inv=width,**kwargs)


@requires_lightweaver
def test_actual_pressure_context_cache_and_live_neutral_refusal():
    from tests.test_bride_solar_control_rya1281 import inputs
    from tests.test_bride_trace_fe_lines_rya1281 import line,context
    from pipeline.bride_solar_control import bound_fe_factory
    from pipeline.bride_fe_pressure import pressure_fe_factory
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_top_cap import FiniteColumnCap
    from lightweaver.atomic_table import PeriodicTable
    folder,base,s,atom,contract=inputs('fine_T5772');provider,_=bound_fe_factory(folder,base,s,atom.abundance,column_contract=contract)
    factory=pressure_fe_factory(provider,[line()],(777.,779.),[recipe()])
    pressure=make_reference_context(s,atom,rays=2,qscale=1,fe_background=factory,
        extra_wavelengths_vac_nm=np.array([777.5,778.,778.5]),experimental_column_cap=FiniteColumnCap(s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]))
    natural=context([line()]);np.testing.assert_array_equal(pressure.background.sca,natural.background.sca)
    center=np.argmin(abs(np.asarray(pressure.spect.wavelength)-778.))
    assert np.any(np.asarray(pressure.background.chi)[center]<np.asarray(natural.background.chi)[center])
    clone=copy.deepcopy(pressure);np.testing.assert_array_equal(clone.background.chi,pressure.background.chi)
    h=pressure.kwargs['eqPops'].atomicPops[PeriodicTable['H']].n
    h[0,0]*=1.001
    with pytest.raises(ValueError,match='neutral colliders'):pressure.background.update_background(pressure.atmos)
