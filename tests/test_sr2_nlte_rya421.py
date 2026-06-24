"""
tests/test_sr2_nlte_rya421.py
=============================
RYA-421 — Sr II resonance-doublet NLTE wiring. Closes the RYA-401 ion-mismatch gap at
the grid level: the Sr II 4077/4215 delta grid (Bergemann 2012a via INSPECT) is fetched,
registered, and near-LTE at our metallicities, with a LOUD [Fe/H]-ceiling flag past +0.0
and the 4215 Fe I-blend cool-star discipline encoded. (The Sr II *measurement* is the
remaining Step-1 piece — not asserted here.)
"""
import pandas as pd

from config.constants import NLTE_CORRECTION_ELEMENTS, SR2_LINES, SR2_COOLSTAR_TEFF_K
from pipeline.nlte_corrections import (_load_mpia_element_grid, _mpia_element_delta,
                                       element_grid_in_bounds, detect_placeholder_zero_lines)


def test_sr_registered_as_sr_ii_delta_grid():
    sr = NLTE_CORRECTION_ELEMENTS['Sr']
    assert sr['ion'] == 2 and sr['grid'] == 'Sr_Bergemann2012_INSPECT.csv'
    assert 'INSPECT' in sr['flag']
    assert 'Mashonkina' in sr['ref'] and 'Bergemann' in sr['ref']   # source-flag pattern


def test_grid_loads_clean_no_placeholder():
    df = pd.read_csv('data/nlte_grids/Sr_Bergemann2012_INSPECT.csv')
    # both resonance lines present, varying (not a placeholder-zero column)
    assert set(round(w, 1) for w in df.wave_A.unique()) >= {4077.7, 4215.5}
    assert detect_placeholder_zero_lines(df) == []
    c = _load_mpia_element_grid('Sr')
    assert len(c['waves']) == 2


def test_solar_anchor_is_small_near_lte():
    # the sanity gate: solar Sr II NLTE correction is SMALL (~0.01), not the metal-poor regime
    for w in (4077.709, 4215.519):
        d = _mpia_element_delta('Sr', w, 5772, 4.44, 0.0)
        assert abs(d) < 0.02, f"solar Sr II {w} delta {d} too large — wrong grid orientation?"


def test_feh_ceiling_flag_metal_rich_out_of_hull():
    # grid [Fe/H] ceiling is +0.0; the metal-rich targets fall out-of-hull (loud, RYA-409)
    assert element_grid_in_bounds('Sr', 5772, 4.44, 0.0) is True           # Sun in range
    assert element_grid_in_bounds('Sr', 5172, 4.43, 0.32) is False         # 55 Cnc above ceiling
    assert element_grid_in_bounds('Sr', 6400, 4.27, 0.25) is False         # tau Boo above ceiling


def test_4215_coolstar_blend_rule_encoded():
    assert SR2_LINES['primary'] == [4077.709]          # clean line, all targets
    assert SR2_LINES['crosscheck'] == [4215.519]       # Fe I-blended -> cool-star caution
    assert SR2_COOLSTAR_TEFF_K == 5300.0               # below this 4215 is cross-check/dropped
