"""
tests/test_co_synthesis_rya564.py
=================================
RYA-564 regression guards for the cobalt decision.

Two invariants, one per half of the ticket:

  A. The blue-edge Co I 3845 value (A=6.128, +1.188; Kitt Peak SNR~24, chi2r~3100) is
     DEMOTED to diagnostic-only and can NEVER be the reported A(Co) — not when the red-line
     synthesis is missing, not when it runs, not when it fails to clear reliability. This is
     the ticket's CRITICAL condition and the whole reason the artifact must not reach the
     RYA-527 gold-v3 freeze.

  B. The red-line HFS synthesis governs: a reliable measurement lands as the value (PASS
     inside TOL_PASS, CURATION-OWED outside), and an unreliable one lands as NO VALUE —
     never as a fallback to 3845.

Plus provenance guards on the committed measurement artifact (gf single source, NLTE
departures actually engaged, reliability floor honoured).
"""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as pc  # noqa: E402

CO_JSON = ROOT / 'data' / 'results' / 'co_synthesis_rya564.json'


def _kp_payload(a_1dlte=6.128):
    """Minimal Kitt Peak payload carrying the blue-edge Co I 3845 measurement."""
    return {
        'leg_validation': {'leg_validated': True, 'overlap': []},
        'n_cross_indicator': {'atomic_NI_mean': None, 'atomic_NI_spread': None},
        'measurements': [{
            'key': 'CoI_3845', 'element': 'Co', 'a_1dlte': a_1dlte, 'cont_snr': 24.0,
            'red_chi2': 3126.9, 'blue_edge': True, 'asplund2021': 4.94,
            'delta_vs_asplund': round(a_1dlte - 4.94, 3), 'status': 'ok',
        }],
    }


def _co_row(a=6.128):
    return [{'element': 'Co', 'asplund2021': 4.94, 'A_measured': a,
             'delta_vs_asplund': None if a is None else round(a - 4.94, 3), 'sigma': None,
             'n_lines': 1, 'channel': 'x', 'verdict': 'CURATION-OWED', 'owed': '',
             'provenance': 'kittpeak-measured'}]


def _synth(a_co, n=5, scatter=0.069):
    return {'_summary': {
        'A_Co': a_co, 'n_reliable': n, 'scatter': scatter,
        'median_nlte_delta': 0.0951, 'nlte_anchor_consistent': True,
        'reliable_lines': {'5352.040': 4.875, '6771.034': 5.008} if a_co else {},
        'primary_only': {'median': 4.924, 'n': 4, 'scatter': 0.077, 'mean': 4.926, 'lines': {}},
        'iag_crosscheck': {'median': 4.939, 'n': 5, 'scatter': 0.061, 'mean': 4.93, 'lines': {}},
        'reason': None if a_co else 'no red Co I line cleared the reliability floor',
    }}


# ── A. the 3845 artifact is demoted, unconditionally ────────────────────────────────
def test_kittpeak_co_emits_no_value():
    """The KP channel still records 3845, but stops reporting it as A(Co)."""
    ov = pc._kittpeak_reclassify(_kp_payload())['Co']
    assert ov.get('A_measured') is None, "3845 must not be handed back as a value"
    assert ov['A_measured_blank'] is True
    assert 'DIAGNOSTIC-ONLY' in ov['channel']
    assert ov['n_lines'] == 0


def test_kittpeak_apply_blanks_a_standing_co_value():
    """The demotion CLEARS a value already on the row — a demoted artifact must not
    survive just because something upstream had written it."""
    rows = _co_row(6.128)
    pc._apply_kittpeak(rows, _kp_payload())
    assert rows[0]['A_measured'] is None
    assert rows[0]['delta_vs_asplund'] is None


def test_3845_never_survives_even_with_no_synthesis():
    """Full order-of-operations with the synthesis artifact ABSENT: Co ends with no value."""
    rows = _co_row(6.128)
    pc._apply_kittpeak(rows, _kp_payload())
    pc._apply_co_synthesis(rows, None)          # synthesis hasn't run
    assert rows[0]['A_measured'] is None
    assert rows[0]['verdict'] == 'CURATION-OWED'


@pytest.mark.parametrize('a_co', [None, 4.965])
def test_3845_value_never_reappears(a_co):
    """Whatever the synthesis says, 6.128 is never the reported value."""
    rows = _co_row(6.128)
    pc._apply_kittpeak(rows, _kp_payload())
    pc._apply_co_synthesis(rows, _synth(a_co))
    assert rows[0]['A_measured'] != pytest.approx(6.128)


# ── B. the red-line synthesis governs ───────────────────────────────────────────────
def test_reliable_red_line_value_lands_and_passes():
    rows = _co_row(6.128)
    pc._apply_kittpeak(rows, _kp_payload())
    pc._apply_co_synthesis(rows, _synth(4.965))
    r = rows[0]
    assert r['A_measured'] == pytest.approx(4.965)
    assert r['delta_vs_asplund'] == pytest.approx(0.025, abs=1e-3)
    assert abs(r['delta_vs_asplund']) <= pc.TOL_PASS
    assert r['verdict'] == 'PASS'
    assert r['n_lines'] == 5
    assert 'RYA-564' in r['provenance']


def test_offset_red_line_value_is_owed_not_passed():
    """Validate-don't-tune: a synthesised Co that sits outside tol is a finding, not a PASS."""
    rows = _co_row(6.128)
    pc._apply_kittpeak(rows, _kp_payload())
    pc._apply_co_synthesis(rows, _synth(5.30))
    assert rows[0]['verdict'] == 'CURATION-OWED'
    assert rows[0]['A_measured'] == pytest.approx(5.30)


def test_no_reliable_line_means_no_value_not_a_fallback():
    rows = _co_row(6.128)
    pc._apply_kittpeak(rows, _kp_payload())
    ov = pc._apply_co_synthesis(rows, _synth(None))['Co']
    assert rows[0]['A_measured'] is None
    assert rows[0]['verdict'] == 'CURATION-OWED'
    assert 'not a fallback' in ov['owed']


# ── provenance guards on the committed measurement artifact ─────────────────────────
@pytest.mark.skipif(not CO_JSON.exists(), reason="RYA-564 measurement artifact not present")
def test_measurement_artifact_provenance():
    d = json.loads(CO_JSON.read_text())
    meta, summ = d['_meta'], d['_summary']
    assert meta['nlte_engine'].startswith('Engine-B TS-Gerber')
    assert meta['reliable_dEW_dA_floor'] == 40.0
    lines = {k: v for k, v in d.items() if not k.startswith('_')}
    assert lines, "no per-line records"
    for w, rec in lines.items():
        # NLTE departures must have ENGAGED for every measured line — never silent LTE.
        assert rec['nlte']['engaged'] is True, f"{w}: departures not engaged"
        # gf carried through from the canonical single source, and the HFS components
        # reconcile with it (verbatim, or explicitly rescaled onto it — never silently off).
        assert rec['hfs']['canonical_loggf'] == rec['loggf']
        if not rec['hfs']['rescaled_to_canonical']:
            assert abs(rec['hfs']['ges_minus_canonical']) <= 0.005, f"{w}: off-SSOT gf"
    # every line that counts toward the value really did clear the floor
    for w in summ.get('reliable_lines', {}):
        h = lines[w]['harps']
        assert h['reliable'] is True
        assert h['dEW_dA_mA_dex'] >= meta['reliable_dEW_dA_floor']
        assert h['railed'] is False


@pytest.mark.skipif(not CO_JSON.exists(), reason="RYA-564 measurement artifact not present")
def test_gate_lines_are_diagnostic_not_the_value():
    """RYA-534's 5000/5013 gated the NLTE atom but are sub-mA — they must not be carrying
    the abundance. Recorded, excluded by the reliability floor, and labelled diagnostic."""
    d = json.loads(CO_JSON.read_text())
    for w in ('5000.875', '5013.324'):
        assert d[w]['role'] == 'diagnostic'
        assert d[w]['harps']['reliable'] is False
        assert w not in (d['_summary'].get('reliable_lines') or {})


@pytest.mark.skipif(not CO_JSON.exists(), reason="RYA-564 measurement artifact not present")
def test_nlte_delta_reconciles_with_the_rya534_anchor():
    """The per-line deltas this run READ from the grid must reproduce the RYA-534 gate."""
    s = json.loads(CO_JSON.read_text())['_summary']
    assert s['nlte_anchor_consistent'] is True
    assert abs(s['median_nlte_delta'] - 0.10) <= 0.12
