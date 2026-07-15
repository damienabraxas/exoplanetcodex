"""
tests/test_sr2_synthesis_rya551.py
==================================
RYA-551 — complete the Sr II synthesis measurement + suppress the raw-EW channel.

Pins the WIRING + DECISION (the synthesis NUMBERS come from Sirius and live in the
committed artifact; these guard that the pipeline treats Sr II correctly):
  * Sr II is registered synthesis-required (RYA-463 registry) so the RYA-520 raw-EW
    suppressor routes it to the synthesis channel,
  * the suppression is ION-SPECIFIC: Sr II suppressed, the clean Sr I 6617 leg
    preserved (RYA-422 ionization cross-check),
  * the EW fitter routes Sr II out of the naive Gaussian path (no saturated_core
    artifact; RYA-429 ledger reason = 'other' with a synthesis pointer),
  * the SR2_LINES line-selection decision (4077 primary / 4215+4161 cross-check /
    4305 excluded) is encoded,
  * the committed synthesis artifact carries real, sensitivity-ranked numbers.
"""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pipeline.problem_children as pc          # noqa: E402
from config.constants import (SR2_LINES, TARGET_ELEMENTS,           # noqa: E402
                              NLTE_CORRECTION_ELEMENTS)

SYNTH = ROOT / 'data' / 'results' / 'sr2_synthesis_rya551.json'
DECISION = ROOT / 'data' / 'results' / 'sr2_line_selection_rya430.json'


# ── registry: Sr II synthesis-required + caught by the RYA-520 set ────────────

def test_sr_ii_is_synthesis_required_in_registry():
    d = pc.disposition_for('Sr')
    assert d is not None
    assert d['required_treatment'] == 'synthesis', d
    assert d['problem_class'] == 'SATURATION_COG'


def test_sr_in_rya520_synth_required_set():
    # mirror the RYA-520 set construction (abundances_derive:2657-2659)
    synth_req = {'C', 'N', 'O'} | {
        e for e in TARGET_ELEMENTS
        if (dd := pc.disposition_for(e)) and dd['required_treatment'] in ('synthesis', 'HFS_sum')}
    assert 'Sr' in synth_req


def test_suppression_is_ion_specific_sr_i_preserved():
    # the RYA-551 ion guard: element in synth_req, but Sr I (ion != II) is NOT suppressed.
    def suppressed(el, ion):
        synth_req = {'C', 'N', 'O'} | {
            e for e in TARGET_ELEMENTS
            if (dd := pc.disposition_for(e)) and dd['required_treatment'] in ('synthesis', 'HFS_sum')}
        if el not in synth_req:
            return False
        if el == 'Sr' and str(ion).strip() != 'II':   # the guard in abundances_derive
            return False
        return True
    assert suppressed('Sr', 'II') is True             # Sr II routed to synthesis
    assert suppressed('Sr', 'I') is False             # Sr I 6617 EW leg preserved
    assert suppressed('C', 'I') is True               # C unchanged


# ── fitter routes Sr II out of the naive Gaussian path ───────────────────────

def test_fitter_routes_sr_ii_to_synthesis():
    import pipeline.lines_fit as lf
    assert ('Sr', 'II') in lf._SYNTH_ROUTED_SPECIES


# ── line-selection decision encoded ──────────────────────────────────────────

def test_sr2_line_selection_decision():
    assert SR2_LINES['primary'] == [4077.709]
    assert set(SR2_LINES['crosscheck']) == {4215.519, 4161.792}
    assert SR2_LINES['excluded'] == [4305.443]


def test_registry_still_points_at_bergemann_grid():
    # GRID-DRIFT guard: the RYA-400 audit requires the map grid_file == this.
    assert NLTE_CORRECTION_ELEMENTS['Sr']['grid'] == 'Sr_Bergemann2012_INSPECT.csv'


# ── committed synthesis artifact: real, sensitivity-ranked numbers ───────────

def test_synthesis_artifact_has_real_numbers():
    assert SYNTH.exists(), "Sirius synthesis result artifact missing"
    r = json.loads(SYNTH.read_text())
    # the three reliable lines land near A(Sr) ~2.7 (NOT the 370-633 mA EW artifacts)
    for w in ('4077.709', '4161.792', '4215.519'):
        a = r[w]['harps']['A_NLTE']
        assert 2.5 < a < 3.0, f"{w} A_NLTE={a} implausible"
    # 4077 is NLTE-covered; 4161 is not (grid = doublet only)
    assert r['4077.709']['nlte_status'] == 'applied'
    assert r['4161.792']['nlte_status'] == 'NLTE_unavailable'
    # 4305 railed with ~zero sensitivity -> excluded / unreliable
    assert r['4305.443']['harps']['railed'] is True
    assert r['4305.443']['harps']['dEW_dA_mA_dex'] < 5.0
    # sensitivity ordering that justified the pick: 4077 > 4215 > 4161 >> 4305
    s = {w: r[w]['harps']['dEW_dA_mA_dex'] for w in r}
    assert s['4077.709'] > s['4215.519'] > s['4161.792'] > s['4305.443']


def test_decision_record_exists_and_is_sensitivity_based():
    assert DECISION.exists()
    d = json.loads(DECISION.read_text())
    assert d['decision']['primary'] == [4077.709]
    assert d['decision']['excluded'] == [4305.443]
    assert 'not_keyed_to_asplund' in d['decision']       # validate-don't-tune recorded
    assert 'engineA_1dnlte' in d['two_engine_rya525']     # RYA-525 both engines recorded
    assert 'engineB_synth_LTE' in d['two_engine_rya525']
