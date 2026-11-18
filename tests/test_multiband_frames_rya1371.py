import numpy as np
import pytest
from pipeline.spectral_frames import supported_air_axis


def test_air_conversion_never_receives_uv_and_missing_air_is_explicit(monkeypatch):
    received=[]
    def convert(w):received.extend(w);return w-1.
    monkeypatch.setattr('pipeline.wavelength_util.vac_to_air',convert)
    air,mask=supported_air_axis([1302.,2000.,6300.,53000.])
    np.testing.assert_array_equal(mask,[False,False,True,True]);assert np.isnan(air[:2]).all()
    np.testing.assert_array_equal(air[2:],[6299.,52999.]);assert received==[6300.,53000.]


@pytest.mark.parametrize('axis',[[],[2.,1.],[0.,1.],[1.,float('nan')]])
def test_invalid_vacuum_axes_refuse(axis):
    with pytest.raises(ValueError):supported_air_axis(axis)
