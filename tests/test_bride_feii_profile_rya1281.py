import numpy as np
import pytest


def test_fixed_contrast_preserves_scale_and_symmetric_linear_background():
    from pipeline.bride_feii_profile import feature_metrics
    w=np.linspace(-1,1,1001);f=1+.001*w
    kw=dict(core_half_A=.1,shoulder_inner_A=.35,shoulder_outer_A=.55)
    assert abs(feature_metrics(w,f,0.,**kw)['contrast'])<1e-14
    f=f-.01*np.exp(-(w/.08)**2)
    a=feature_metrics(w,f,0.,**kw);b=feature_metrics(w,f*1.03,0.,**kw)
    assert a['contrast']>0
    assert b['contrast']==pytest.approx(a['contrast'],abs=1e-14)
    assert b['relative_core_area_mA']==pytest.approx(a['relative_core_area_mA'],abs=1e-12)


def test_feature_truncation_refuses_instead_of_moving_shoulders():
    from pipeline.bride_feii_profile import feature_metrics
    with pytest.raises(ValueError,match='Truncated'):
        feature_metrics(np.linspace(-.4,.4,101),np.ones(101),0.,core_half_A=.1,shoulder_inner_A=.35,shoulder_outer_A=.55)


def test_ionic_stark_modes_zero_density_and_source_temperature_scaling():
    from pipeline.bride_feii_profile import feii_stark_gamma
    row=dict(ionization_eV=16.188,excitation_eV=3.199,wavelength_air_A=7777.099,log_gamma_stark_10000=-6.530)
    for mode in ('provided','cowley'):
        assert feii_stark_gamma(5772.,0.,row,mode)==0
    assert feii_stark_gamma(10000.,1e20,row,'provided')==pytest.approx(10**-6.530*1e14)
    assert feii_stark_gamma(10000.,1e20,row,'cowley')==pytest.approx(feii_stark_gamma(4000.,1e20,row,'cowley'))
    with pytest.raises(ValueError):feii_stark_gamma(5772.,1e20,dict(row,ionization_eV=4.),'cowley')


def test_native_recipe_requires_unambiguous_delimited_record():
    from pipeline.bride_feii_profile import native_recipe_record
    good=' FEII_RECIPE W 0.000000000000000E+00 0.000000000000000E+00\n'
    assert native_recipe_record(good)==('W',0.,0.)
    for text in (' FEII_RECIPEW 0.0 0.0\n',good+good,' FEII_RECIPE U 0.0 0.0\n'):
        with pytest.raises(ValueError):native_recipe_record(text)
