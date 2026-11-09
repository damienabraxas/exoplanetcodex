import numpy as np
import pytest


def test_gaussian_quadrature_and_explicit_no_deconvolution():
    from pipeline.bride_resolution_audit import extra_gaussian_sigma_pixels
    w=np.linspace(7765,7780,1501);centre=7772.;target=500000.
    shipped=extra_gaussian_sigma_pixels(w,centre,target)
    for native in (700000.,1000000.):
        adjusted=extra_gaussian_sigma_pixels(w,centre,target,native_R=native)
        assert adjusted/shipped==pytest.approx(np.sqrt(1-(target/native)**2))
    for native in (0.,target,400000.,np.nan):
        with pytest.raises(ValueError):extra_gaussian_sigma_pixels(w,centre,target,native_R=native)
    assert extra_gaussian_sigma_pixels(w,centre,target,no_extra=True)==0


def test_shipped_filter_parity_and_no_extra_input_immutability():
    from scipy.ndimage import gaussian_filter1d
    from pipeline.bride_resolution_audit import broaden_reference
    w=np.linspace(7765,7780,1501);f=1-.1*np.exp(-((w-7772)/.03)**2);old=f.copy()
    actual,meta=broaden_reference(w,f,7772.,500000.)
    expected=gaussian_filter1d(f,7772./500000./2.3548/np.median(np.diff(w)))
    np.testing.assert_array_equal(actual,expected)
    unchanged,_=broaden_reference(w,f,7772.,500000.,no_extra=True)
    np.testing.assert_array_equal(unchanged,old)
    assert not np.shares_memory(unchanged,f)
    np.testing.assert_array_equal(f,old)
    assert meta['boundary_mode']=='reflect' and meta['truncate']==4.
