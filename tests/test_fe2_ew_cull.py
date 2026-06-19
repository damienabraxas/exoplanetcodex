"""RYA-352 — Fe II EW-quality cull gate (_apply_fe2_ew_quality_cull)."""
import numpy as np
import pandas as pd

from pipeline.abundances_derive import _apply_fe2_ew_quality_cull
from config.constants import PIPELINE

CEIL = PIPELINE['vmic_ew_ceiling_mA']
ERRMAX = PIPELINE['fe2_ew_err_frac_max']


def _row(wl, ew, err, blend=False, el='Fe', ion='II'):
    return {'element': el, 'ion': ion, 'wavelength_air_A': wl,
            'ew_mA': ew, 'ew_err_mA': err, 'blend_flag': blend}


def test_clean_fe2_lines_kept():
    # 6084/6456: clean EW lines a88ef0f wrongly dropped "because RYA-347" → must survive.
    df = pd.DataFrame([_row(6084.102, 28.5, 6.57), _row(6456.380, 85.4, 6.87)])
    out = _apply_fe2_ew_quality_cull(df)
    assert set(out['wavelength_air_A']) == {6084.102, 6456.380}


def test_saturated_line_culled():
    df = pd.DataFrame([_row(5234.623, CEIL + 20, 22.0)])  # EW above the COG ceiling
    assert len(_apply_fe2_ew_quality_cull(df)) == 0


def test_high_error_line_culled():
    df = pd.DataFrame([_row(5256.932, 10.8, 10.8 * (ERRMAX + 0.3))])  # err/EW > max
    assert len(_apply_fe2_ew_quality_cull(df)) == 0


def test_blend_line_culled():
    df = pd.DataFrame([_row(6149.258, 46.5, 7.2, blend=True)])
    assert len(_apply_fe2_ew_quality_cull(df)) == 0


def test_non_fe2_untouched():
    # A saturated Fe I / Ba II line must NOT be culled by the Fe II-only gate.
    df = pd.DataFrame([_row(5000.0, CEIL + 50, 5.0, el='Fe', ion='I'),
                       _row(5853.7, CEIL + 50, 5.0, el='Ba', ion='II')])
    assert len(_apply_fe2_ew_quality_cull(df)) == 2


def test_missing_err_column_no_hierr_cull():
    # Without ew_err_mA the HIERR criterion can't fire (no silent cull on missing data).
    df = pd.DataFrame([{'element': 'Fe', 'ion': 'II', 'wavelength_air_A': 6084.1,
                        'ew_mA': 28.5, 'blend_flag': False}])
    assert len(_apply_fe2_ew_quality_cull(df)) == 1
