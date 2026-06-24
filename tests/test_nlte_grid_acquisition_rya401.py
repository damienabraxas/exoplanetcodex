"""
tests/test_nlte_grid_acquisition_rya401.py
==========================================
RYA-401 — Beast round-2 NLTE-grid vendoring (Al/K/S/Cu/V/Sr). Round-2's honest outcome was
0 of 6 vendorable via a scrape. RYA-412 + RYA-421 reconcile these guards to the current
adjudication: **Al and S were SUPERSEDED by RYA-402** (departure grids via PySME → registered,
LOCKED) and **Sr by RYA-421** (the ion mismatch was resolved at the grid level — the Sr II
4077/4215 NLTE grid is fetched + registered → LOCKED; measurement owed, recorded in rya401).
**K/Cu/V remain unregistered gaps** — K and Cu HAVE a grid (RYA-402 PySME) but are data-blocked
(GET-DATA: K telluric/RYA-380, Cu line-quality/RYA-395); V is a genuine no-public-grid gap
(HARD-carry-forward). No element still claims GET-GRID, none is falsely registered without a
grid, each carries `rya401:` provenance, and the RYA-400 audit still PASSes.
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
# RYA-412 + RYA-421 + RYA-428 adjudication of the round-2 six:
LOCKED_REGISTERED  = ['Al', 'S']    # PySME grid + measured prescribed-ion line → LOCKED (RYA-402)
PENDING_REGISTERED = ['Sr']         # grid registered (RYA-421) but Sr II measurement OWED →
                                    # GET-DATA-pending, NOT LOCKED (RYA-428 reverted premature handled)
REGISTERED         = LOCKED_REGISTERED + PENDING_REGISTERED
STILL_GAP          = ['K', 'Cu', 'V']     # unregistered: data-blocked or no public grid


def test_no_get_grid_verdicts_remain():
    # round-2 + RYA-402 flipped every GET-GRID element out; nothing should still claim GET-GRID
    assert not [el for el, s in MAP.items() if s.get('verdict') == 'GET-GRID']


def test_six_adjudicated_post_merge():
    # Al/S LOCKED (RYA-402: grid + measured line). Sr is grid-REGISTERED but GET-DATA-pending
    # (RYA-428: measurement owed -> not LOCKED). K/Cu/V are documented gaps. Each of the six
    # carries its rya401 provenance (probe + a resolution 'outcome' or a gap_reason/path).
    for el in LOCKED_REGISTERED:
        assert MAP[el]['verdict'] == 'LOCKED', f"{el} should be LOCKED"
    for el in PENDING_REGISTERED:
        assert MAP[el]['verdict'] == 'GET-DATA', f"{el} should be GET-DATA-pending (measurement owed)"
    for el in STILL_GAP:
        assert MAP[el]['verdict'] in ('HARD-carry-forward', 'GET-DATA'), el
    for el in SIX:
        r = MAP[el].get('rya401')
        assert r and r.get('probe'), f"{el} missing rya401 probe provenance"
        assert r.get('outcome') or (r.get('gap_reason') and r.get('path')), \
            f"{el} missing rya401 outcome (resolved) or gap_reason/path (genuine gap)"


def test_registration_matches_adjudication_no_false_done():
    # Al/S registered via RYA-402, Sr via RYA-421; the still-gap three must NOT be registered.
    for el in REGISTERED:
        assert el in NLTE_CORRECTION_ELEMENTS, f"{el} should be registered"
    for el in STILL_GAP:
        assert el not in NLTE_CORRECTION_ELEMENTS, f"{el} must not be registered (unresolved gap)"


def test_s_indicator_decision_is_optical():
    s = MAP['S']
    assert 'indicator_decision' in s
    assert 'optical' in s['indicator_decision'].lower()
    assert '6743' in s['indicator_decision'] or '6757' in s['indicator_decision']


def test_sr_ion_mismatch_grid_resolved_measurement_owed():
    # RYA-421 resolved the Sr ion mismatch at the GRID level (Sr II 4077/4215 NLTE grid
    # fetched + registered). RYA-428: that is NOT "handled" — the Sr II MEASUREMENT is owed,
    # so the verdict stays GET-DATA-pending, not LOCKED (registration != measurement).
    sr = MAP['Sr']
    assert sr['ion'] == 'II' and sr['verdict'] == 'GET-DATA'
    assert 'Sr' in NLTE_CORRECTION_ELEMENTS and NLTE_CORRECTION_ELEMENTS['Sr']['ion'] == 2
    assert 'RYA-421' in sr['rya401']['outcome']
    assert 'OWED' in sr['rya401']['outcome'] or 'owed' in sr['rya401']['outcome'].lower()
    assert 'measurement' in sr['rya401']['outcome'].lower()   # Step-1 measurement owed, recorded


def test_k_cross_refs_rya380_not_resolved_here():
    # K's grid is now HAVE-VALIDATED (RYA-402) but it stays GET-DATA, data-blocked on RYA-380
    # (telluric/measurement). The RYA-380 cross-ref moved into the rya401 'outcome'.
    assert MAP['K'].get('xref') == 'RYA-380'
    assert 'RYA-380' in MAP['K']['rya401']['outcome']
    assert MAP['K']['verdict'] == 'GET-DATA'


def test_v_genuine_grid_gap_cu_data_blocked():
    # V is the genuine no-public-grid gap (HARD-carry-forward). Cu's grid is now HAVE
    # (RYA-402 PySME) but unregistered — the blocker moved from the grid to line quality
    # (GET-DATA → RYA-395), so Cu no longer claims a grid gap.
    assert MAP['V']['grid']['status'] == 'GAP'
    assert MAP['V']['verdict'] == 'HARD-carry-forward'
    assert MAP['Cu']['grid']['status'] == 'HAVE'
    assert MAP['Cu']['verdict'] == 'GET-DATA' and 'Cu' not in NLTE_CORRECTION_ELEMENTS


def test_rya400_audit_still_passes():
    assert A.run() == 0                                 # map ↔ code still consistent
