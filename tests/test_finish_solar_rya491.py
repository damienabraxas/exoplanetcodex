"""
tests/test_finish_solar_rya491.py
=================================
RYA-491 finish-solar — guard the committed result (reads the JSON, no re-run).
  * the IR O I 844 measurement exists and corroborates 777 with a formation-depth offset;
  * O I 926 is excised as unusable (weak + flag-limited) — NOT averaged in;
  * O I 777 stays the bankable primary; both RT legs recorded;
  * C holds on the new reference; S still high on the correct line (gf-floor finding); N flagged.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

J = ROOT / 'data' / 'results' / 'finish_solar_rya491.json'


def _d():
    return json.loads(J.read_text())


def test_frame_declared_and_validate_dont_tune():
    d = _d()
    assert d['validate_dont_tune'] is True
    assert 'Salami' in d['frame_declared'] and 'Birch' in d['frame_declared']   # RYA-481


def test_ir_oi844_measured_and_corroborates_777_with_formation_offset():
    d = _d()
    o844 = next(r for r in d['ir_O'] if r['indicator'] == 'OI_844')
    assert o844['reliable'] is True and o844['flag_frac'] < 0.05
    # 844 sits below 777 (forms deeper) — a real offset, the formation-depth finding
    assert d['O_844_formation_offset'] < 0 and abs(d['O_844_formation_offset']) < 0.30
    assert d['O_primary_777_3d'] == 8.736


def test_oi926_excised_as_unusable_not_averaged():
    d = _d()
    o926 = next(r for r in d['ir_O'] if r['indicator'] == 'OI_926')
    assert o926['reliable'] is False                 # weak + flag-limited -> excluded
    assert d['O_926_unusable'] is True
    # the usable-indicator spread excludes the railed 926 (so it is small, a real finding)
    assert d['O_reliable_spread'] < 0.4


def test_C_holds_on_new_reference():
    cc = _d()['C_confirm']
    assert cc['c_5052_3d'] is not None
    assert abs(cc['c_5052_3d'] - cc['banked']) < 0.10    # ~8.45 vs banked 8.491


def test_S_still_high_on_correct_line_gf_finding():
    d = _d()
    su = d['S_fix']
    # measured the CORRECT line (6743.5), not the discarded 6748.68
    assert any(abs(p['line'] - 6743.53) < 0.1 for p in su['per_line'])
    assert su['still_high'] is True                      # +0.40 over Asplund -> gf lever
    assert 'Costa Silva' in d['S_gf_caveat'] and 'gf' in d['S_gf_caveat']


def test_N_flagged_and_crires_blocked():
    d = _d()
    assert 'RYA-369' in d['N_status'] and 'DEFERRED' in d['N_status']
    assert 'STAGGER' in d['blocked'] and 'CRIRES' in d['blocked']
