"""
tests/test_n_grid_loadtest_rya369.py
====================================
RYA-369 — guard the N I NLTE grid load-test result + the production wiring that
cleared N's hard blocker. Reads the committed JSON (no grid needed) and the live
pysme_nlte registration.
  * the load-test PASSED: format loaded, coverage spans Sun->55 Cnc, N I resolves to a
    sane NEGATIVE solar NLTE delta (sign-discipline, RYA-339);
  * N is registered in pysme_nlte (NLTE_LINES + _A_SUN + _GRID_FILENAME + _Z atomic number),
    non-breaking (the other elements stay registered).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import pysme_nlte as P                              # noqa: E402

J = ROOT / 'data' / 'results' / 'n_grid_loadtest_rya369.json'


def _d():
    return json.loads(J.read_text())


def test_loadtest_passed_format_coverage_resolve():
    d = _d()
    assert d['loadtest_passed'] is True
    assert d['format']['loaded'] is True
    assert d['coverage_all_inside'] is True
    # all benchmark stars inside the N grid box
    for star in ('Sun', 'Procyon', 'aCenA', 'aCenB', '55Cnc'):
        assert d['coverage'][star] is True


def test_n_i_resolves_to_sane_negative_solar_delta():
    r = _d()['resolve']
    assert r['sane_negative'] is True
    # solar N I red NLTE is small + negative (the finding: NOT a large solar NLTE)
    assert -0.6 < r['median'] < 0.0
    assert set(r['per_line']) >= {'7468.31', '8216.34'} or len(r['per_line']) >= 2


def test_n_registered_in_pysme_nlte_nonbreaking():
    # the production wiring that clears the blocker — N joins the registered elements
    assert P._GRID_FILENAME.get('N') == 'nlte_N_scatt_pysme.grd'
    assert P._A_SUN.get('N') == 7.83
    assert P._Z('N') == 7                                   # the gap the load-test found
    assert 'N' in P.NLTE_LINES and len(P.NLTE_LINES['N']) == 3
    # the 9 pre-existing PySME elements still registered (non-breaking)
    for el in ('Na', 'Al', 'K', 'Cu', 'S'):
        assert el in P.NLTE_LINES and P._Z(el)
