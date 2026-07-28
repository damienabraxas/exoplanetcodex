"""
tests/test_procyon_cno_redo_rya493.py
=====================================
RYA-493 — guard the Procyon CNO redo (reads the committed JSON, no re-run).
  * O re-banked regime-matched (777 1D-NLTE - solar 1D 8.755), confirms provisional +0.066;
  * C holds + the 0.253 spread re-stated arm-matched (VIS C I separated from atomic-vs-molecular);
  * N additive leg measured with the warm-Teff NLTE finding (Procyon delta materially > solar);
  * metals untouched (surgical CNO).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

J = ROOT / 'data' / 'results' / 'procyon_cno_redo_rya493.json'


def _d():
    return json.loads(J.read_text())


def test_O_rebanked_regime_matched_confirms_provisional():
    o = _d()['O']
    # [O/H] = 777 1D-NLTE (8.82) - solar 1D denominator (8.755)
    assert abs(o['oi777_a_1d'] - (o['oh_rebanked'] + 8.755)) < 0.002
    assert _d()['solar_denominators']['O_1D'] == 8.755          # regime-matched, new-solar
    # confirms the provisional +0.066 (tiny delta)
    assert abs(o['delta_from_provisional']) < 0.02
    assert o['continuum_sigma'] == 0.186                        # carried, not re-litigated


def test_C_holds_and_spread_restated_arm_matched():
    c = _d()['C']
    assert c['holds'] is True
    # the spread is split, NOT collapsed: VIS C I (Bruntt-comparable) vs atomic-vs-molecular (ours)
    assert c['vis_ci_spread'] > 0.1                             # ~0.197, the real C I internal spread
    assert c['atomic_vs_molecular_spread'] < 0.10               # ~0.038, atomic/molecular agree
    assert c['full_spread'] != c['vis_ci_spread']               # not the same number


def test_N_additive_warm_teff_nlte_confirmed():
    n = _d()['N']
    assert n['a_nlte'] is not None and n['a_lte'] is not None
    # the strategy premise: Procyon (6554 K) N I NLTE is materially LARGER than solar (~-0.015)
    assert n['nlte_delta_procyon'] is not None
    assert abs(n['nlte_delta_procyon']) > abs(n['nlte_solar']) * 2   # ~-0.059 vs -0.015
    assert n['warm_teff_larger'] is True
    # honestly flagged: solar N not yet banked
    assert n['solar_n_bank_pending'] is True


def test_metals_untouched_surgical_cno():
    assert _d()['metals_untouched'] is True
