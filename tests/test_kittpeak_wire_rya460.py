"""
tests/test_kittpeak_wire_rya460.py
==================================
RYA-460 — guard the Kitt Peak reference wiring: the KP loader, the overlap
leg-validation gate ([O I] 6300 / O I 777 vs HARPS/ESPRESSO), the N solar
cross-indicator agreement, the honest P/K/Co/Sc reclassification off DATA-GAP,
and the provenance/cited discipline. Runs on the COMMITTED KP measurement JSON +
the verdict-script functions — does NOT need the slow synthesis baseline.
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as V                          # noqa: E402
from scripts.wire_reference_atlases_rya460 import load_kp   # noqa: E402

KP_JSON = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_kittpeak_rya460.json'


def _kp():
    return json.loads(KP_JSON.read_text())


def test_kittpeak_loader_returns_normalized_spectrum():
    w, f = load_kp('OI_777_triplet')
    assert len(w) > 100 and len(w) == len(f)
    assert 777.0 < w.mean() < 778.0                  # air nm
    assert f.min() < 0.85 and f.max() <= 1.2         # normalized, with the triplet cores


def test_overlap_cross_check_validates_the_leg():
    kp = _kp()
    ov = kp['leg_validation']['overlap']
    assert kp['leg_validation']['leg_validated'] is True
    for line in ('OI_6300', 'OI_777'):
        assert ov[line]['agree'] is True
        assert abs(ov[line]['delta']) <= 0.10        # KP raw agrees with HARPS/ESPRESSO


def test_overlap_is_a_real_cross_check_not_snapped():
    # KP raw values are independent measurements, NOT snapped to the Phase A numbers.
    kp = _kp()
    ov = kp['leg_validation']['overlap']
    assert ov['OI_6300']['kittpeak_raw'] != ov['OI_6300']['phase_a_raw']
    assert ov['OI_777']['kittpeak_raw'] != ov['OI_777']['phase_a_raw']


def test_N_three_red_multiplets_agree():
    nci = _kp()['n_cross_indicator']
    assert nci['atomic_NI_spread'] <= 0.10           # the 3 N I red channels agree on the Sun
    ind = nci['indicators']
    for k in ('NI_7442_7468', 'NI_8216_8223', 'NI_8680_8718'):
        assert ind[k] is not None
    # the LTE overestimate (+0.3..0.4 vs 7.83) is the OWED N I NLTE, not validated
    assert nci['atomic_NI_mean'] > 7.83


def test_NH_CN_blue_edge_flagged_not_forced():
    nci = _kp()['n_cross_indicator']
    assert set(nci['molecular_blue_edge_flagged']) == {'NH_3360', 'CN_violet_3883'}


def test_reclassify_moves_pkcosc_off_datagap():
    ov = V._kittpeak_reclassify(_kp())
    assert ov['N']['verdict'] == 'NLTE-OWED'
    assert ov['K']['verdict'] == 'PASS'               # RYA-462 wired K_Amarsi2020 -> NLTE applied
    for el in ('P', 'Co', 'Sc'):
        assert ov[el]['verdict'] == 'CURATION-OWED'
        assert ov[el]['verdict'] != 'DATA-GAP'
    # every KP override carries the measured-from-Kitt-Peak provenance
    for el in ('N', 'K', 'P', 'Co', 'Sc'):
        assert ov[el]['provenance'] == 'kittpeak-measured'


def test_apply_kittpeak_eliminates_datagap_in_verdict_rows():
    # Synthetic prior rows: P/K/Co/Sc as DATA-GAP, N as NLTE-OWED.
    rows = [{'element': e, 'asplund2021': a, 'A_measured': None, 'delta_vs_asplund': None,
             'sigma': None, 'n_lines': 0, 'verdict': v, 'channel': 'x', 'owed': 'x',
             'provenance': 'x', 'nlte_wired': False, 'threed_wired': False, 'nlte_grid': None}
            for e, a, v in [('N', 7.83, 'NLTE-OWED'), ('P', 5.41, 'DATA-GAP'),
                            ('K', 5.07, 'DATA-GAP'), ('Co', 4.94, 'DATA-GAP'),
                            ('Sc', 3.14, 'DATA-GAP')]]
    V._apply_kittpeak(rows, _kp())
    verdicts = {r['element']: r['verdict'] for r in rows}
    assert 'DATA-GAP' not in verdicts.values()        # all 4 closed
    assert verdicts['K'] == 'PASS'                    # RYA-462 wired K's NLTE grid
    # KP measured value now populated (off None) and not tuned to Asplund
    p = next(r for r in rows if r['element'] == 'P')
    assert p['A_measured'] is not None and p['A_measured'] != p['asplund2021']


def test_uv_composite_stays_cited_after_wiring():
    # The RYA-459 provenance gate must still hold — KP wiring never flips UV to measured.
    from pipeline import audit_solar_reference as A
    assert A.assert_provenance() is True
