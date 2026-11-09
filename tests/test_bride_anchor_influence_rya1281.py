import numpy as np
import pytest


def test_accounting_reconstructs_independent_fit_difference_and_keeps_cancellation():
    from pipeline.bride_anchor_influence import fixed_support_blocks
    from pipeline.bride_continuum_audit import support_geometry
    x=np.linspace(-5,-4,101);g=support_geometry(x,0.,np.arange(len(x)))
    a=1+.001*np.sin(x*30);b=a+.002+.003*x
    r=fixed_support_blocks(x,a,b,g['centre_weights'],.1)
    expected=np.polyfit(x,b,1)[-1]-np.polyfit(x,a,1)[-1]
    assert r['centre_delta']==pytest.approx(expected,abs=2e-14)
    assert r['centre_delta']==pytest.approx(.002,abs=2e-14)
    assert r['absolute_contribution_sum']>abs(r['centre_delta'])
    assert sum(v['count'] for v in r['blocks'])==len(x)
    assert abs(r['closure'])<1e-15


def test_constant_and_slope_controls_do_not_create_false_centre_offsets():
    from pipeline.bride_anchor_influence import fixed_support_blocks
    from pipeline.bride_continuum_audit import support_geometry
    x=np.linspace(-5,2,141);h=support_geometry(x,0.,np.arange(len(x)))['centre_weights']
    one=np.ones_like(x)
    assert fixed_support_blocks(x,one,one+.001,h,.1)['centre_delta']==pytest.approx(.001)
    assert abs(fixed_support_blocks(x,one,one+.003*x,h,.1)['centre_delta'])<1e-15


def test_invalid_grouping_and_support_refuse():
    from pipeline.bride_anchor_influence import fixed_support_blocks
    x=np.arange(10.);one=np.ones(10)
    for width in (0.,-1.,np.nan):
        with pytest.raises(ValueError):fixed_support_blocks(x,one,one,one,width)
    with pytest.raises(ValueError):fixed_support_blocks(x[::-1],one,one,one,.1)
