import numpy as np
import pytest

from pipeline.telluric.arl_gdas import (NX, NY, ArlHeader, _nearest_grid,
                                        interpolate_at_height, unpack)


def _pack(values, exponent, first):
    scale = 2.0 ** (7 - exponent)
    out = np.empty((NY, NX), dtype=np.uint8)
    first_column = first
    for j in range(NY):
        previous = first_column
        for i in range(NX):
            byte = int((values[j, i] - previous) * scale + 127.5)
            out[j, i] = byte
            previous += (byte - 127) / scale
            if i == 0:
                first_column = previous
    return out.tobytes()


def test_unpack_is_inverse_of_noaa_row_predictor():
    x = np.arange(NX)[None, :]
    y = np.arange(NY)[:, None]
    values = 250.0 + 0.5 * x + 0.5 * y
    h = ArlHeader(2023, 8, 2, 15, 1, "TEMP", 6, 0.0, values[0, 0])
    decoded = unpack(_pack(values, h.exponent, h.first_value), h)
    assert np.max(np.abs(decoded - values)) == 0.0


def test_la_silla_nearest_gdas_cell_is_29s_71w():
    i, j, lat, lon = _nearest_grid(-29.26, -70.73)
    assert (i, j, lat, lon) == (289, 61, -29.0, -71.0)


def test_site_interpolation_uses_log_pressure():
    profile = {"levels": [
        {"height_m": 2000, "press_hpa": 800, "temp_k": 294, "relhum_pct": 9},
        {"height_m": 2600, "press_hpa": 750, "temp_k": 292, "relhum_pct": 11},
    ]}
    got = interpolate_at_height(profile, 2300)
    assert got["pressure_hpa"] == pytest.approx((800 * 750) ** 0.5)
    assert got["temperature_c"] == pytest.approx(20.0 - 273.15 + 273.0)
    assert got["relative_humidity_pct"] == 10.0
