"""
tests/test_threed_corrections_rya399.py
=======================================
RYA-399 — the metal 3D abundance-correction leg.

Pins the contract: the registry/grid load with the cited solar values, the
correction composes ON TOP of 1D-NLTE (no double-count) and touches only the
registered (element, ion), unregistered elements pass through, and — the headline
— the published 3D corrections do NOT close the RYA-398 graded solar residual
(validate-don't-tune: they are independent literature, demonstrably not fit to it).
"""
import numpy as np
import pandas as pd
import pytest

from pipeline import threed_corrections as td
from config.constants import THREED_CORRECTION_ELEMENTS


def test_registry_and_grid_values():
    assert set(THREED_CORRECTION_ELEMENTS) == {'Si', 'Ti', 'Cr'}
    assert td.solar_threed_delta('Si') == pytest.approx(-0.01)
    assert td.solar_threed_delta('Ti') == pytest.approx(+0.06)
    assert td.solar_threed_delta('Cr') == pytest.approx(+0.03)
    assert np.isnan(td.solar_threed_delta('Fe'))          # unregistered → NaN


def test_all_corrections_small_and_within_ceiling():
    grid = td._load_solar3d_grid()
    for el, node in grid.items():
        assert abs(node['delta_3d']) <= td._THREED_MAG_CEILING_DEX
        assert abs(node['delta_3d']) <= 0.1                # the finding: 3D is small here


def test_apply_composes_on_top_of_nlte():
    df = pd.DataFrame({
        'element': ['Si', 'Ti', 'Cr', 'Fe', 'Cr'],
        'ion':     ['I',  'I',  'I',  'I',  'II'],   # Cr II must be untouched (ion filter)
        'A_X':      [7.9, 5.5, 6.0, 7.46, 6.0],
        'A_X_nlte': [7.888, 5.471, 6.022, 7.50, 6.0],
    })
    out = td.apply_threed_corrections(df)
    si = out[(out.element == 'Si') & (out.ion == 'I')].iloc[0]
    assert si['A_X_3dnlte'] == pytest.approx(7.888 - 0.01, abs=1e-3)
    ti = out[(out.element == 'Ti')].iloc[0]
    assert ti['A_X_3dnlte'] == pytest.approx(5.471 + 0.06, abs=1e-3)
    # Fe is unregistered → no 3D row written
    fe = out[out.element == 'Fe'].iloc[0]
    assert fe['threed_flag'] == 'no_3D' and not np.isfinite(fe['delta_3d'])
    # Cr II (wrong ion) untouched
    cr2 = out[(out.element == 'Cr') & (out.ion == 'II')].iloc[0]
    assert cr2['threed_flag'] == 'no_3D'


def test_apply_falls_back_to_A_X_without_nlte():
    df = pd.DataFrame({'element': ['Ti'], 'ion': ['I'], 'A_X': [5.5]})
    out = td.apply_threed_corrections(df)
    assert out.iloc[0]['A_X_3dnlte'] == pytest.approx(5.5 + 0.06, abs=1e-3)


def test_unavailable_when_no_base_abundance():
    df = pd.DataFrame({'element': ['Cr'], 'ion': ['I'], 'A_X': [np.nan]})
    out = td.apply_threed_corrections(df)
    assert out.iloc[0]['threed_flag'] == '3D_unavailable'


def test_3D_does_not_close_residual():
    """The headline: 3D leaves the +0.4 dex graded residual essentially intact and
    for Ti/Cr makes it WORSE (positive correction) — so 3D is not the lever."""
    rep = td.residual_after_3d()
    assert set(rep['element']) == {'Si', 'Ti', 'Cr'}
    assert (rep['verdict'] == '3D_NOT_THE_LEVER').all()
    assert (~rep['closes_residual']).all()
    # Ti/Cr: positive 3D → residual magnitude grows (the wrong way)
    for el in ('Ti', 'Cr'):
        r = rep[rep.element == el].iloc[0]
        assert abs(r['resid_3Dnlte']) > abs(r['resid_1Dnlte'])


def test_validate_dont_tune_values_are_literature_not_fit():
    """The deltas are the published solar corrections; none was chosen to minimise
    our residual (if they had been, resid_3Dnlte would be ~0 — it is ~+0.4)."""
    rep = td.residual_after_3d()
    assert (rep['resid_3Dnlte'].abs() > 0.30).all()
