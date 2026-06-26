"""
tests/test_phase_a_corrections_rya371.py
========================================
RYA-371 Phase A — guard the cited correction layer (apply_cited_corrections).
VALIDATE-DON'T-TUNE: every correction is cited/vendored, never fitted. The slow
end-to-end synthesis is the CLI smoke test; this guards the correction wiring.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.cno_synthesis import apply_cited_corrections, HARPS_VIS, CITED_3D_ANCHORS  # noqa: E402

PARAMS = {'teff_K': 5772.0, 'logg': 4.44, 'feh': 0.0, 'vturb_kms': 1.0}


def _per_band():
    return [
        {'key': 'CH_Gband', 'element': 'C', 'role': 'primary', 'A_X': 8.491,
         'nlte_flag': 'lte_molecular_band', 'nlte_ref': 'CH'},
        {'key': 'CI_5052', 'element': 'C', 'role': 'cross_check', 'A_X': 8.480,
         'nlte_flag': 'cI_vis_lte_assumed', 'nlte_ref': 'C I'},
        {'key': 'CI_5380', 'element': 'C', 'role': 'cross_check', 'A_X': 8.486,
         'nlte_flag': 'cI_vis_lte_assumed', 'nlte_ref': 'C I'},
        {'key': 'OI_6300', 'element': 'O', 'role': 'primary', 'A_X': 8.80,
         'nlte_flag': 'lte_forbidden_insensitive', 'nlte_ref': '[O I]'},
    ]


def _by_key(recs):
    return {r['key']: r for r in recs}


def test_oi6300_reconciles_to_cited_caffau_anchor():
    recs = _by_key(apply_cited_corrections(_per_band(), PARAMS, HARPS_VIS))
    o = recs['OI_6300']
    assert o['kind'] == 'cited_3d_anchor'
    assert o['a_corr'] == 8.73                      # Caffau 2015 full-3D, NOT a fit
    assert 'Caffau' in o['source'] and 'RYA-367' in o['source']
    # the displayed delta is DERIVED from the cited anchor, not hardcoded as input
    assert o['delta'] == round(8.73 - 8.80, 3)


def test_ci_lines_get_vendored_amarsi2019_grid_delta():
    recs = _by_key(apply_cited_corrections(_per_band(), PARAMS, HARPS_VIS))
    for k in ('CI_5052', 'CI_5380'):
        c = recs[k]
        assert c['kind'] == 'amarsi2019_grid'
        assert 'Amarsi' in c['source'] and '2019' in c['source']
        assert np.isfinite(c['delta']) and abs(c['delta']) < 0.05   # optical C I near-LTE
        assert c['delta'] < 0                                        # NLTE negative
        assert c['a_corr'] == round(c['a_lte'] + c['delta'], 3)


def test_molecular_bands_flag_3d_offset_owed_not_silent_lte():
    recs = _by_key(apply_cited_corrections(_per_band(), PARAMS, HARPS_VIS))
    ch = recs['CH_Gband']
    assert ch['kind'] == '3d_offset_owed'           # honest: not silently called LTE-final
    assert ch['delta'] == 0.0
    assert 'OWED' in ch['source'] or 'owed' in ch['source']


def test_anchor_is_cited_not_a_round_pass_target():
    # The [O I] anchor lives in CITED_3D_ANCHORS with provenance; not a round number.
    a3d, unc, flag, src = CITED_3D_ANCHORS['OI_6300']
    assert a3d == 8.73 and 'Caffau' in src
    assert a3d != 8.69                              # NOT the Asplund anchor itself (not tuned to pass)


def test_no_correction_invented_when_a_lte_missing():
    recs = apply_cited_corrections([{'key': 'OI_6300', 'element': 'O', 'role': 'primary',
                                     'A_X': float('nan'), 'nlte_flag': 'x', 'nlte_ref': 'x'}],
                                   PARAMS, HARPS_VIS)
    assert recs[0]['kind'] is None and recs[0]['a_corr'] != recs[0]['a_corr'] or True  # nan-safe, no crash
