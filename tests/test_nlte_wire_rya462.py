"""
tests/test_nlte_wire_rya462.py
==============================
RYA-462 Part B/D — wiring the K NLTE grid into the production path + verdict.

Pins: K is registered in NLTE_CORRECTION_ELEMENTS; the applied solar delta is the
VENDORED grid value read through the existing interpolation subsystem (validate-don't-
tune, not fitted to Asplund); the verdict reclassifies K NLTE-OWED -> PASS off the
cited correction; the PASS is the grid value (not snapped to Asplund); an unregistered
element is never silently corrected; and the Cr canary holds — Cr's small NLTE delta
cannot reach Asplund, so it stays at its curation floor (a Cr PASS would be the RYA-398
red flag), and the K wiring never touches Cr.

Runs on the COMMITTED Kitt Peak JSON + grid CSV + verdict JSON — no slow synthesis.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as V              # noqa: E402
from config.constants import NLTE_CORRECTION_ELEMENTS  # noqa: E402
from pipeline import nlte_corrections as N       # noqa: E402

KP_JSON = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_kittpeak_rya460.json'
VERDICT_JSON = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'


def _kp():
    return json.loads(KP_JSON.read_text())


# ── registration ──────────────────────────────────────────────────────────────
def test_K_registered_in_nlte_elements():
    assert 'K' in NLTE_CORRECTION_ELEMENTS
    spec = NLTE_CORRECTION_ELEMENTS['K']
    assert spec['grid'] == 'K_Amarsi2020_PySME.csv'
    assert int(spec['ion']) == 1
    assert (ROOT / 'data' / 'nlte_grids' / spec['grid']).exists()


# ── the delta is the vendored grid value (validate-don't-tune) ─────────────────
def test_K_solar_delta_is_vendored_not_fitted():
    N._mpia_element_cache.clear()
    d = N._mpia_element_delta('K', 7698.964, 5772.0, 4.44, 0.0)
    assert np.isfinite(d)
    assert -0.40 < d < -0.20                 # severe negative K I resonance NLTE
    assert abs(d - (-0.3121)) < 1e-3         # the exact vendored node value


def test_apply_nlte_grid_delta_helper_uses_subsystem():
    a_nlte, d, flag = V._apply_nlte_grid_delta('K', 7698.964, 5.411)
    assert np.isfinite(d) and d < 0
    assert abs(a_nlte - 5.099) < 0.01        # 5.411 + (-0.312)
    assert 'NLTE' in flag and 'unavailable' not in flag.lower()


def test_unregistered_element_never_silently_corrected():
    a_nlte, d, flag = V._apply_nlte_grid_delta('Xx', 5000.0, 1.0)
    assert a_nlte == 1.0 and not np.isfinite(d)
    assert 'unavailable' in flag.lower()


# ── the verdict reclassifies K to PASS off the cited grid ──────────────────────
def test_K_verdict_now_pass():
    ov = V._kittpeak_reclassify(_kp())
    assert ov['K']['verdict'] == 'PASS'
    assert abs(ov['K']['A_measured'] - 5.099) < 0.01
    assert abs(ov['K']['A_measured'] - 5.07) <= V.TOL_PASS   # reconciles within TOL


def test_K_pass_is_grid_value_not_snapped_to_asplund():
    ov = V._kittpeak_reclassify(_kp())
    # PASS because the cited NLTE delta lands it in band — NOT because we snapped to 5.07.
    assert ov['K']['A_measured'] != 5.07
    assert ov['K']['provenance'] == 'kittpeak-measured'


# ── Cr canary ──────────────────────────────────────────────────────────────────
def test_cr_nlte_delta_cannot_reach_asplund():
    # Cr's NLTE grid solar delta is SMALL (<=~0.13) across every grid line — NLTE alone
    # cannot pull the +0.40 raw Cr offset to Asplund. The residual is a CURATION floor.
    N._mpia_element_cache.clear()
    c = N._load_mpia_element_grid('Cr')
    deltas = [N._mpia_element_delta('Cr', w, 5772.0, 4.44, 0.0) for w in c['waves']]
    deltas = [d for d in deltas if np.isfinite(d)]
    assert deltas and max(abs(d) for d in deltas) < 0.20


def test_K_wiring_does_not_touch_cr():
    ov = V._kittpeak_reclassify(_kp())
    assert 'Cr' not in ov                     # KP / K wiring never reclassifies Cr


# ── committed verdict artifact: the re-emitted counts + invariants ─────────────
def test_committed_verdict_counts_and_invariants():
    d = json.loads(VERDICT_JSON.read_text())
    by = {r['element']: r for r in d['verdicts']}
    assert by['K']['verdict'] == 'PASS'
    assert by['Cr']['verdict'] == 'CURATION-OWED'      # canary: held at floor, not PASS
    assert abs(by['Cr']['delta_vs_asplund'] - 0.402) < 0.01
    for el in ('O', 'C', 'Fe'):
        assert by[el]['verdict'] == 'PASS'             # unchanged
    # RYA-476: Mn certifies PASS via HFS synthesis + the live triplet-exact Amarsi
    # NLTE δ (A(Mn) 5.470, +0.050 vs Asplund — inside TOL). Earned by the measurement,
    # not a tune (the δ is the small line-exact value, not the MPIA high-EP +0.108).
    assert by['Mn']['verdict'] == 'PASS'
    assert abs(by['Mn']['A_measured'] - 5.470) < 0.01
    assert by['N']['verdict'] == 'NLTE-OWED'           # the one remaining NLTE-OWED
    counts = d['summary']['counts']
    assert counts == {'PASS': 5, 'NLTE-OWED': 1, 'CURATION-OWED': 20}
