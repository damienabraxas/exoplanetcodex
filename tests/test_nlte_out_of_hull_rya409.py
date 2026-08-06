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
# RYA-409 found 7 MPIA/INSPECT elements clamped at +0.30. RYA-410 re-sourced Na/Mg/Si
# onto Amarsi-2020 PySME grids ([Fe/H] -> +0.6), CLOSING their clamp; Ca/Mn stayed on
# MPIA (PySME-Amarsi cross-check STOPPED — Ca model diff, Mn HFS), and Cr is not in
# Amarsi-2020 — so those three still clamp at 55 Cnc.
#
# RYA-545 REMOVED Ti from this list: production Ti moved onto Ti_Mallinson2024_PySME.csv,
# whose hull spans 2500-8000 K and reaches [Fe/H] +0.31, retiring the MPIA clamp (and the
# RYA-505 6500 K Ti wall with it). Ti clamping at 55 Cnc is now the *old* behaviour —
# asserting it would re-encode a closed defect.
STILL_CLAMPED = ['Ca', 'Cr', 'Mn']
RESOURCED = ['Na', 'Mg', 'Si']            # RYA-410: Amarsi grid now covers +0.31


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


def test_remaining_four_clamp_loud_at_55cnc():
    """A mock 55-Cnc run: the four still-MPIA elements (Ca/Ti/Cr/Mn) report
    NLTE_OUT_OF_HULL loudly — none silently drops to LTE (RYA-409 guard)."""
    flagged = []
    for el in STILL_CLAMPED:
        ab, pl = _dfs(el)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out = apply_element_nlte_corrections(ab, CNC, per_line_df=pl, elements=[el])
        if out.iloc[0]['nlte_flag'] == 'NLTE_OUT_OF_HULL':
            flagged.append(el)
    assert set(flagged) == set(STILL_CLAMPED)              # all four loudly clamped, none silent


def test_resourced_elements_no_longer_clamp_at_55cnc():
    """RYA-410: Na/Mg/Si moved to Amarsi grids covering [Fe/H]+0.31, so the 55-Cnc clamp
    is CLOSED — they are in-hull and must NOT trip the out-of-hull guard."""
    for el in RESOURCED:
        assert element_grid_in_bounds(el, teff=CNC['teff_K'], logg=CNC['logg'], feh=CNC['feh']), \
            f"{el} should be in-hull at 55 Cnc after RYA-410 re-source"
        ab, pl = _dfs(el)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out = apply_element_nlte_corrections(ab, CNC, per_line_df=pl, elements=[el])
        assert out.iloc[0]['nlte_flag'] != 'NLTE_OUT_OF_HULL', \
            f"{el} clamp should be closed by RYA-410, not OUT_OF_HULL"


def test_in_hull_does_not_trip_the_guard():
    # At solar, Ca is in-hull -> the guard must NOT fire (flag is anything but OUT_OF_HULL).
    assert element_grid_in_bounds('Ca', **{'teff': 5772, 'logg': 4.44, 'feh': 0.0})
    ab, pl = _dfs('Ca')
    out = apply_element_nlte_corrections(ab, SOLAR, per_line_df=pl, elements=['Ca'])
    assert out.iloc[0]['nlte_flag'] != 'NLTE_OUT_OF_HULL'
