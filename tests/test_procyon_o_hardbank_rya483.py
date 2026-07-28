"""
tests/test_procyon_o_hardbank_rya483.py
=======================================
RYA-483 — guard the Procyon O hard-bank result (reads the committed JSON, does NOT
re-run the slow synthesis). The findings that must not silently regress:
  * solar control unmoved (the RYA-478 guard — method must not move the Sun);
  * the cross-instrument zero-point is the LARGE continuum-lever value (~0.18), the
    finding that corrects RYA-478's 0.04-0.08 estimate — not silently shrunk;
  * [O I] 6300 is terminal (solar control fails) -> the corroboration leg is unavailable;
  * the verdict is PROVISIONAL / not hard-banked, with the residual named (not a forced bank).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

J = ROOT / 'data' / 'results' / 'procyon_O_hardbank_rya483.json'


def _d():
    return json.loads(J.read_text())


def test_solar_control_unmoved():
    d = _d()
    # KPNO solar O I 777 (1D-NLTE) reproduces our Sun 8.735 to ~1 mdex — method doesn't move it
    assert d['solar_control_unmoved'] is True
    assert abs(d['part_a_zeropoint']['solar']['a_nlte_local'] - d['our_solar_O']) < 0.02


def test_primary_is_oi777_differential_vs_our_sun():
    d = _d()
    assert d['primary_indicator'] == 'OI_777_1D_NLTE'
    assert d['our_solar_O'] == 8.735                         # our measured Sun, never Asplund
    assert abs(d['procyon_OH_vs_our_sun'] - (d['procyon_O_AO_nlte'] - d['our_solar_O'])) < 1e-6


def test_zeropoint_is_the_large_continuum_lever_finding():
    # the headline correction to RYA-478: the cross-instrument zero-point is ~0.18 dex,
    # dominated by the UVES continuum-localization lever (solar -0.05 vs procyon -0.23).
    zp = _d()['part_a_zeropoint']
    assert zp['continuum_diff_zp'] < -0.12                    # big differential, not 0.04-0.08
    assert zp['solar']['shift'] > -0.10                       # Sun barely moves
    assert zp['procyon']['shift'] < -0.15                     # UVES swings a lot
    assert _d()['banked_sigma'] >= 0.15                       # σ carries the zero-point, not hidden


def test_oi6300_leg_is_terminal_unavailable():
    b = _d()['part_b_oi6300_leg']
    assert b['solar_terminal'] is True                        # solar [O I] 6300 control failed
    assert b['leg_available'] is False                        # cannot corroborate O I 777
    # the solar control is far off our Sun (the terminal signature, RYA-447/448/455)
    assert abs(b['solar_kpno_oi6300_lte'] - _d()['our_solar_O']) > 0.30


def test_verdict_is_provisional_not_hard_banked_with_named_blockers():
    d = _d()
    assert d['hard_bankable'] is False
    assert 'PROVISIONAL' in d['verdict'] or 'NOT_HARD' in d['verdict']
    # both residuals explicitly named (not a silent bank)
    blockers = ' '.join(d['blockers']).lower()
    assert 'continuum' in blockers and '6300' in blockers
