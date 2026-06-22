"""
tests/test_nlte_out_of_hull_rya409.py
=====================================
RYA-409 Part A — the loud out-of-hull NLTE guard. Out-of-hull (Teff,logg,[Fe/H])
must NOT silently fall back to 1D-LTE: it flags loud + distinct (NLTE_OUT_OF_HULL),
warns, and (strict) raises. This is the 55-Cnc metal-rich clamp made loud.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from pipeline.nlte_corrections import apply_element_nlte_corrections, element_grid_in_bounds

# 55 Cnc A: cool, metal-rich — [Fe/H] +0.31 just past the +0.30 MPIA grid ceiling.
CNC = {'teff_K': 5172, 'logg': 4.43, 'feh': 0.31}
SOLAR = {'teff_K': 5772, 'logg': 4.44, 'feh': 0.0}
# the 7 registered MPIA/INSPECT elements whose grid ceiling is +0.30 (clamp at 55 Cnc)
CLAMPED = ['Ca', 'Ti', 'Cr', 'Na', 'Mg', 'Mn', 'Si']


def _dfs(element, ion='I', wave=6000.0):
    ab = pd.DataFrame([{'element': element, 'ion': ion, 'A_X': 5.0}])
    pl = pd.DataFrame([{'element': element, 'ion': ion,
                        'wavelength_air_A': wave, 'a_1dlte': 5.0}])
    return ab, pl


def test_out_of_hull_flags_loud_not_silent():
    ab, pl = _dfs('Ca')
    with pytest.warns(RuntimeWarning, match="OUTSIDE the NLTE grid hull"):
        out = apply_element_nlte_corrections(ab, CNC, per_line_df=pl, elements=['Ca'])
    row = out.iloc[0]
    assert row['nlte_flag'] == 'NLTE_OUT_OF_HULL'          # distinct, not generic 'NLTE_unavailable'
    assert row['A_X_nlte'] == pytest.approx(5.0)           # LTE retained — but WITH the flag
    assert not np.isfinite(row['delta_nlte_mean'])         # no delta silently applied


def test_strict_raises_out_of_hull():
    ab, pl = _dfs('Ca')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with pytest.raises(ValueError, match="out.of.hull|OUTSIDE the NLTE grid hull"):
            apply_element_nlte_corrections(ab, CNC, per_line_df=pl, elements=['Ca'], strict=True)


def test_all_seven_clamp_loud_at_55cnc():
    """A mock 55-Cnc run: all 7 ceiling-+0.30 elements report NLTE_OUT_OF_HULL loudly."""
    flagged = []
    for el in CLAMPED:
        ion = 'II' if el == 'Ba' else 'I'                  # (Ba not clamped; all 7 here are neutral)
        ab, pl = _dfs(el)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out = apply_element_nlte_corrections(ab, CNC, per_line_df=pl, elements=[el])
        if out.iloc[0]['nlte_flag'] == 'NLTE_OUT_OF_HULL':
            flagged.append(el)
    assert set(flagged) == set(CLAMPED)                    # all 7 loudly clamped, none silent


def test_in_hull_does_not_trip_the_guard():
    # At solar, Ca is in-hull -> the guard must NOT fire (flag is anything but OUT_OF_HULL).
    assert element_grid_in_bounds('Ca', **{'teff': 5772, 'logg': 4.44, 'feh': 0.0})
    ab, pl = _dfs('Ca')
    out = apply_element_nlte_corrections(ab, SOLAR, per_line_df=pl, elements=['Ca'])
    assert out.iloc[0]['nlte_flag'] != 'NLTE_OUT_OF_HULL'
