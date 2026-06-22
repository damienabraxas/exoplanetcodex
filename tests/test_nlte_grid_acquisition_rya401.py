"""
tests/test_nlte_grid_acquisition_rya401.py
==========================================
RYA-401 — Beast round-2 NLTE-grid vendoring (Al/K/S/Cu/V/Sr). The honest outcome: 0 of 6
vendorable, all documented gaps (no fabrication). These guards pin that result so it
can't silently regress: the 6 are out of GET-GRID, none is registered as NLTE-applied,
each carries its `rya401:` probe/gap/path provenance, S's indicator decision is recorded,
the Sr ion-mismatch is captured, K cross-refs RYA-380, and the RYA-400 audit still PASSes.
"""
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import audit_physics_regime_rya400 as A            # noqa: E402
from config.constants import NLTE_CORRECTION_ELEMENTS  # noqa: E402

MAP = yaml.safe_load((ROOT / 'config' / 'physics_regime_rya400.yaml').read_text())['elements']
SIX = ['Al', 'K', 'S', 'Cu', 'V', 'Sr']


def test_no_get_grid_verdicts_remain():
    # round-2 flipped every GET-GRID element out; nothing should still claim GET-GRID
    assert not [el for el, s in MAP.items() if s.get('verdict') == 'GET-GRID']


def test_six_flipped_to_documented_gaps():
    for el in SIX:
        assert MAP[el]['verdict'] in ('HARD-carry-forward', 'GET-DATA'), el
        r = MAP[el].get('rya401')
        assert r and r.get('probe') and r.get('gap_reason') and r.get('path'), \
            f"{el} missing rya401 probe/gap/path provenance"


def test_none_of_the_six_is_registered_no_false_done():
    for el in SIX:
        assert el not in NLTE_CORRECTION_ELEMENTS, f"{el} must not be registered (no grid vendored)"


def test_s_indicator_decision_is_optical():
    s = MAP['S']
    assert 'indicator_decision' in s
    assert 'optical' in s['indicator_decision'].lower()
    assert '6743' in s['indicator_decision'] or '6757' in s['indicator_decision']


def test_sr_ion_mismatch_recorded():
    sr = MAP['Sr']
    assert sr['ion'] == 'II' and sr['verdict'] == 'GET-DATA'
    assert 'mismatch' in sr['rya401']['gap_reason'].lower()
    assert 'Sr I' in sr['rya401']['gap_reason']        # our line is Sr I, grid is Sr II


def test_k_cross_refs_rya380_not_resolved_here():
    assert MAP['K'].get('xref') == 'RYA-380'
    assert 'RYA-380' in MAP['K']['rya401']['path']


def test_cu_and_v_are_genuine_grid_gaps():
    assert MAP['Cu']['grid']['status'] == 'GAP'
    assert MAP['V']['grid']['status'] == 'GAP'


def test_rya400_audit_still_passes():
    assert A.run() == 0                                 # map ↔ code still consistent
