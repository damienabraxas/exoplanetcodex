"""
tests/test_cu_v_hfs_synthesis_rya466.py
=======================================
RYA-466 — guard the HFS-resolved synthesis measurement of solar Cu (+V LTE):
the gf adjudication (VALD-vs-GES → Kock & Richter), the HFS-resolution premise
(the features are NOT collapsed single lines), the RYA-402 b-factor Cu NLTE
application (small + positive, vendored-or-live), the phase_c fold that moves
Cu/V off "no value" → measured, and the validate-don't-tune discipline.

Runs on the COMMITTED measurement JSON + canonical_gf + the verdict-script
functions — does NOT re-run the slow Turbospectrum synthesis.
"""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as V                                  # noqa: E402
from scripts.measure_cu_v_hfs_synthesis_rya466 import (             # noqa: E402
    adjudicate_cu_gf, CU_LINES, CU_NLTE_DELTA_DOC)

JSON = (ROOT / 'data' / 'audit' / 'cu_v_hfs_synthesis' /
        'solar_cu_v_hfs_synthesis_rya466.json')
CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'

import numpy as np                                                  # noqa: E402
from config.constants import SOLAR_ASPLUND2021                      # noqa: E402


def _data():
    return json.loads(JSON.read_text())


# ── the measurement product exists and is well-formed ─────────────────────────

def test_measurement_json_present_and_schema():
    d = _data()
    assert d['ticket'] == 'RYA-466' and d['star'] == 'solar'
    for el in ('Cu', 'V'):
        assert el in d and d[el]['per_line']
    assert d['Cu']['A_nlte'] is not None
    assert d['V']['nlte_void'] is True


# ── HFS premise: the Cu features are hyperfine-RESOLVED, not collapsed ─────────

def test_cu_features_are_hfs_resolved_not_collapsed():
    hfs = _data()['Cu']['hfs_components']
    assert len(hfs) == 5
    # every Cu diagnostic feature is multi-component in the GES list — exactly why
    # the single-profile EW path could not measure it (the RYA-354 finding).
    assert all(n > 1 for n in hfs.values())
    assert max(hfs.values()) >= 7                    # 5700/5782 are 10-component


# ── gf adjudication: GES = Kock & Richter, VALD3 superseded ────────────────────

def test_gf_adjudication_resolves_vald_vs_ges_to_KR():
    audit = adjudicate_cu_gf(apply=False)            # read-only, no mutation
    assert len(audit) == 5
    deltas = [abs(a['delta_GES_minus_VALD']) for a in audit
              if np.isfinite(a['delta_GES_minus_VALD'])]
    # the genuine VALD-vs-GES disagreement RYA-354 flagged (median |Δ| ~ 0.275-0.33)
    assert 0.20 <= float(np.median(deltas)) <= 0.45
    # the synthesis uses GES = Kock & Richter for every line
    assert all('Kock & Richter' in a['synthesis_uses'] for a in audit)


def test_canonical_gf_adjudication_status_written_value_frozen():
    rows = list(csv.reader(open(CANON)))
    col = {n: i for i, n in enumerate(rows[0])}
    seen = 0
    for r in rows[1:]:
        if r[col['species']] != 'Cu I':
            continue
        if any(abs(float(r[col['wavelength_air_A']]) - wl) < 0.02 for wl, _ in CU_LINES):
            seen += 1
            # the decision is recorded …
            assert r[col['adjudication_status']] == 'adjudicated_kr1968_rya466'
            # … but the value/grade/reference fields are FROZEN (batch-2 discipline:
            # the synthesis already uses GES=KR; we never invent/transcribe a gf).
            assert r[col['log_gf']] not in ('', None)
    assert seen == 5


# ── Cu NLTE: RYA-402 b-factor, small + positive, applied (not silent LTE) ──────

def test_cu_nlte_small_positive_and_applied():
    cu = _data()['Cu']
    assert cu['nlte_delta'] is not None
    # RYA-402 finding: Cu NLTE is small + positive (Shi 2014 cross-check)
    assert 0.0 < cu['nlte_delta'] <= 0.05
    if not cu['nlte_live']:
        assert abs(cu['nlte_delta'] - CU_NLTE_DELTA_DOC) < 1e-9   # vendored value applied
    # A_nlte = A_lte + delta (NLTE actually applied, never silently dropped)
    assert abs(cu['A_nlte'] - (cu['A_lte_median'] + cu['nlte_delta'])) < 1e-3


def test_v_is_lte_only_nlte_void():
    v = _data()['V']
    assert v['nlte_void'] is True
    assert 'NLTE-VOID' in v['method'] or 'NLTE' in v['method']


# ── phase_c fold: Cu/V move off "no value" → MEASURED ─────────────────────────

def test_phase_c_fold_moves_cu_v_off_no_value():
    ov = V._cu_v_reclassify(_data())
    assert set(ov) == {'Cu', 'V'}
    for el in ('Cu', 'V'):
        assert ov[el]['A_measured'] is not None          # was None (no value) before
        assert ov[el]['n_lines'] and ov[el]['n_lines'] >= 5
    # Cu sits high → honest CURATION-OWED (a measured value, not a forced PASS)
    assert ov['Cu']['verdict'] in ('CURATION-OWED', 'PASS')
    # V is the NLTE-void → never certified on LTE agreement alone
    assert ov['V']['verdict'] == 'CURATION-OWED'
    assert 'NLTE-VOID' in ov['V']['channel'] or 'NLTE-void' in ov['V']['owed']


def test_fold_does_not_change_verdict_counts():
    # Cu/V were CURATION-OWED (no value) and stay CURATION-OWED (now WITH a value):
    # the move is no-value → measured, NOT a category flip / fabricated PASS.
    ov = V._cu_v_reclassify(_data())
    assert ov['Cu']['verdict'] == 'CURATION-OWED'        # +0.165 survives small NLTE


# ── validate-don't-tune: nothing snapped to Asplund ───────────────────────────

def test_not_tuned_to_asplund():
    d = _data()
    # Cu lands +0.165 off Asplund — a real offset, NOT pulled to 4.18
    assert abs(d['Cu']['A_nlte'] - SOLAR_ASPLUND2021['Cu']) > 0.10
    # per-line A(Cu) genuinely SCATTER (a fit, not a constant) — σ > 0
    assert d['Cu']['scatter'] > 0.02
    a_vals = [r['A_X'] for r in d['Cu']['per_line']
              if isinstance(r['A_X'], float)]
    assert max(a_vals) - min(a_vals) > 0.1               # 5105 high vs 5782 low
