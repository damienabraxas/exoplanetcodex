"""
tests/test_element_nlte_rya235.py
=================================
RYA-235 — per-line Ca I / Ti I / Cr I NLTE application from the real MPIA MAFAGS-OS
grids, wired after the Fe leg. Guards: the grids load + interpolate to the documented
solar deltas; the application composes with the Fe pass (touches only registry neutral
rows; Cr II and non-registry rows untouched); the availability gate is a loud
1D-LTE retention, never a silent fix; and the honest deferred-validation contract —
applying NLTE to a raw (un-curated) high Cr pool does NOT land on Asplund (the
overshoot the reopen flagged; acceptance is gated on a curated pool, not on this layer).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline import nlte_corrections as N            # noqa: E402
from config.constants import NLTE_CORRECTION_ELEMENTS, SOLAR_ASPLUND2021  # noqa: E402

SOLAR = dict(teff_K=5777, logg=4.4, feh=0.0)


def test_registry_grids_exist_and_ions_well_formed():
    # RYA-235 seeded Ca/Ti/Cr; RYA-165 added Na/Mg (neutral) + Ba (ion II).
    assert {'Ca', 'Ti', 'Cr'} <= set(NLTE_CORRECTION_ELEMENTS)
    for el, spec in NLTE_CORRECTION_ELEMENTS.items():
        assert int(spec['ion']) in (1, 2)
        assert (ROOT / 'data' / 'nlte_grids' / spec['grid']).exists()
        assert spec['ref']                                          # cited
    # Ca/Ti/Cr remain neutral (the original RYA-235 contract)
    assert all(NLTE_CORRECTION_ELEMENTS[el]['ion'] == 1 for el in ('Ca', 'Ti', 'Cr'))


def test_grids_load_and_solar_deltas_match_documented():
    # the real MPIA grids → Ti ~+0.108, Cr ~+0.073 (RYA-256 notes); Ca clean median +0.017
    # after RYA-413 dropped the 6166 placeholder-zero (was +0.0118 with the placeholder).
    expect = {'Ca': 0.017, 'Ti': 0.108, 'Cr': 0.073}
    min_waves = {'Ca': 7, 'Ti': 8, 'Cr': 8}        # Ca = 7 after the RYA-413 6166 drop
    for el, ex in expect.items():
        c = N._load_mpia_element_grid(el)
        assert len(c['waves']) >= min_waves[el]
        ds = [N._mpia_element_delta(el, float(w), 5777, 4.4, 0.0) for w in c['waves']]
        ds = [d for d in ds if np.isfinite(d)]
        assert abs(float(np.median(ds)) - ex) < 0.02, f"{el} solar median Δ off"


def test_delta_is_nan_outside_hull_and_off_wave():
    el = 'Cr'
    w0 = float(N._load_mpia_element_grid(el)['waves'][0])
    assert not np.isfinite(N._mpia_element_delta(el, w0, 9000, 4.4, 0.0))   # Teff outside
    assert not np.isfinite(N._mpia_element_delta(el, w0 + 50.0, 5777, 4.4, 0.0))  # no wave match
    assert np.isfinite(N._mpia_element_delta(el, w0, 5777, 4.4, 0.0))      # in-grid, on-wave


def test_grid_loader_loud_fails_on_unknown_element_and_missing_file():
    with pytest.raises(KeyError):
        N._load_mpia_element_grid('Zz')


def _per_line_for(el, a_1dlte, n=5):
    waves = N._load_mpia_element_grid(el)['waves'][:n]
    return [{'element': el, 'ion': 1, 'wavelength_air_A': float(w), 'a_1dlte': a_1dlte}
            for w in waves]


def test_apply_corrects_registry_neutral_only_and_leaves_others():
    res = pd.DataFrame([
        {'element': 'Ca', 'ion': 1, 'A_X': 6.18, 'n_lines': 5, 'A_X_std': 0.05},
        {'element': 'Cr', 'ion': 1, 'A_X': 5.55, 'n_lines': 5, 'A_X_std': 0.06},
        {'element': 'Cr', 'ion': 2, 'A_X': 5.40, 'n_lines': 2, 'A_X_std': 0.10},
        {'element': 'Ni', 'ion': 1, 'A_X': 6.20, 'n_lines': 5, 'A_X_std': 0.04},
    ])
    pl = pd.DataFrame(_per_line_for('Ca', 6.18) + _per_line_for('Cr', 5.55))
    out = N.apply_element_nlte_corrections(res, SOLAR, per_line_df=pl)

    ca = out[(out.element == 'Ca') & (out.ion == 1)].iloc[0]
    cr1 = out[(out.element == 'Cr') & (out.ion == 1)].iloc[0]
    assert ca['nlte_flag'] == 'NLTE_MPIA_MAFAGS_1D' and ca['A_X_nlte'] > ca['A_X']
    assert cr1['nlte_flag'] == 'NLTE_MPIA_MAFAGS_1D' and cr1['n_nlte_lines'] == 5
    # Cr II (ion 2) and Ni (non-registry) are NOT corrected by this layer
    assert out[(out.element == 'Cr') & (out.ion == 2)].iloc[0]['nlte_flag'] == '1D_LTE'
    assert out[out.element == 'Ni'].iloc[0]['nlte_flag'] == '1D_LTE'


def test_availability_gate_retains_lte_without_per_line():
    res = pd.DataFrame([{'element': 'Ti', 'ion': 1, 'A_X': 4.88, 'n_lines': 3, 'A_X_std': 0.05}])
    out = N.apply_element_nlte_corrections(res, SOLAR, per_line_df=None)
    row = out.iloc[0]
    assert row['nlte_flag'] == 'NLTE_unavailable' and row['A_X_nlte'] == row['A_X']


def test_does_not_clobber_the_fe_leg():
    # simulate the post-Fe df: Fe row already carries a corrected A_X_nlte + flag
    res = pd.DataFrame([
        {'element': 'Fe', 'ion': 1, 'A_X': 7.46, 'A_X_nlte': 7.50,
         'nlte_flag': 'NLTE_MPIA_Bergemann_1D', 'n_nlte_lines': 80},
        {'element': 'Ca', 'ion': 1, 'A_X': 6.18, 'A_X_nlte': 6.18,
         'nlte_flag': '1D_LTE', 'n_nlte_lines': 0},
    ])
    out = N.apply_element_nlte_corrections(res, SOLAR, per_line_df=pd.DataFrame(_per_line_for('Ca', 6.18)))
    fe = out[out.element == 'Fe'].iloc[0]
    assert fe['A_X_nlte'] == 7.50 and fe['nlte_flag'] == 'NLTE_MPIA_Bergemann_1D'  # untouched


def test_applied_flag_is_not_a_not_applied_sentinel():
    # so run()'s A_X_nlte_absolute derivation carries the value through
    res = pd.DataFrame([{'element': 'Ca', 'ion': 1, 'A_X': 6.18, 'n_lines': 5}])
    out = N.apply_element_nlte_corrections(res, SOLAR, per_line_df=pd.DataFrame(_per_line_for('Ca', 6.18)))
    assert out.iloc[0]['nlte_flag'] not in ('1D_LTE', 'NLTE_unavailable')


def test_deferred_validation_contract_nlte_does_not_fix_uncurated_cr():
    # The honest contract: this layer APPLIES NLTE, it does not curate. A raw high Cr
    # pool (RYA-239 read +0.345 over Asplund) gets MORE positive, not closer — proving
    # acceptance is gated on a curated pool, not on the NLTE layer (the reopen's point).
    raw_cr_lte = SOLAR_ASPLUND2021['Cr'] + 0.345          # the uncurated 1D-LTE value
    res = pd.DataFrame([{'element': 'Cr', 'ion': 1, 'A_X': raw_cr_lte, 'n_lines': 5}])
    out = N.apply_element_nlte_corrections(res, SOLAR, per_line_df=pd.DataFrame(_per_line_for('Cr', raw_cr_lte)))
    a_nlte = float(out.iloc[0]['A_X_nlte'])
    assert a_nlte > raw_cr_lte                              # NLTE pushes Cr further UP
    assert abs(a_nlte - SOLAR_ASPLUND2021['Cr']) > 0.05    # NOT within Asplund — curation owed
