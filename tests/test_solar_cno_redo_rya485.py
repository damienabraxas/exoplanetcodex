"""
tests/test_solar_cno_redo_rya485.py
===================================
RYA-485/486 solar CNO redo — guard the committed result (reads the JSON, does NOT re-run
the synthesis). The findings that must not silently regress:
  * the O I 777 four-reference swap-test agrees (continuum NOT reference-driven) -> the
    RYA-483 Procyon σ is intrinsic, and the solar O denominator is re-confirmed at 8.736;
  * both RT legs recorded (3D for the Sun, 1D for warm-star regime-matching);
  * S measured by synthesis (supersedes the EW-pool 7.753), reported as a finding not tuned.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

J = ROOT / 'data' / 'results' / 'solar_cno_redo_rya485.json'
EXTRACTS = ROOT / 'data' / 'solar_reference_swaptest'


def _d():
    return json.loads(J.read_text())


def test_validate_dont_tune_flag_and_asplund_is_comparison_only():
    d = _d()
    assert d['validate_dont_tune'] is True
    # the Asplund values are recorded as comparison, the banked 371 as the re-check baseline
    assert d['asplund_comparison_only']['O'] == 8.69 and d['banked_371']['O'] == 8.735


def test_oi777_swap_test_agrees_continuum_not_reference_driven():
    s = _d()['oi777_swap_test']
    assert len(s['per_reference']) == 4                  # Kurucz + Reiners + Baker + Wallace
    names = {r['reference'] for r in s['per_reference']}
    assert names == {'kurucz_kpno', 'reiners_iag', 'baker_iag_tel', 'wallace_kpno'}
    # all four agree to within ~0.05 dex -> not reference-driven
    assert s['spread'] < 0.05
    assert 'NOT reference-driven' in s['verdict']


def test_solar_O_denominator_reconfirmed_and_both_legs():
    d = _d()
    s = d['oi777_swap_test']
    # combined 3D-NLTE matches the 371-banked 8.735 to a few mdex (re-confirmed, not moved)
    assert abs(s['combined_a_nlte_3d'] - d['banked_371']['O']) <= 0.01
    # both RT legs recorded; 1D denominator is the warm-star regime-match (higher than 3D)
    rm = d['regime_matched_O']
    assert rm['solar_A_O_3D'] == s['combined_a_nlte_3d']
    assert rm['solar_A_O_1D'] > rm['solar_A_O_3D']       # 1D leg ~+0.02 above 3D
    assert abs(rm['solar_A_O_1D'] - rm['solar_A_O_3D']) < 0.05


def test_sulphur_measured_by_synthesis_supersedes_ew_and_is_a_finding():
    su = _d()['sulphur']
    assert su['n'] == 2 and su['a_nlte'] is not None
    # synthesis value, distinct from (lower than) the EW-pool 7.753 banked
    assert su['a_nlte'] < _d()['banked_371']['S']
    # sits above Asplund 7.12 -> a finding, surfaced not snapped
    assert su['a_nlte'] - su['asplund'] > 0.15
    assert su['sigma'] > 0.0                              # real line-to-line scatter


def test_swaptest_reference_extracts_present():
    for fn in ('reiners_oi777_vac.txt', 'baker_oi777_vac.txt', 'wallace_kpno_oi777_vac.txt'):
        assert (EXTRACTS / fn).exists()
