"""
tests/test_mn_hfs_synthesis_rya473.py
=====================================
RYA-473 — guard the HFS-resolved synthesis measurement of solar Mn on the Den
Hartog e6S→z6P triplet (6013/16/21). The REAL unlock RYA-468 pointed at: the gf is
graded (Den Hartog, MED) but the EW path SAT-culls these HFS-split lines, so Mn
stayed no-value; synthesis on the GES HFS-resolved line list measures it.

Guards: the measurement product + schema; the HFS premise (6 components/feature,
NOT a collapsed single line); the gf provenance (GES gf-sum == the RYA-468 Den
Hartog canonical value — cited, not invented/tuned); the NLTE application (the MPIA/
Bergemann grid δ, applied + LOUDLY flagged with the RYA-411 high-EP-vs-triplet
caveat, never silent LTE); the phase_c fold that moves Mn off "no value" → measured
without flipping the verdict counts; and validate-don't-tune (a real fit that
scatters, no value snapped to Asplund).

Runs on the COMMITTED measurement JSON + the verdict-script functions — does NOT
re-run the slow Turbospectrum synthesis.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as V                                  # noqa: E402
from scripts.measure_mn_hfs_synthesis_rya473 import MN_LINES, MN_DH_GF  # noqa: E402
from config.constants import SOLAR_ASPLUND2021                      # noqa: E402

JSON = (ROOT / 'data' / 'audit' / 'mn_hfs_synthesis' /
        'solar_mn_hfs_synthesis_rya473.json')


def _data():
    return json.loads(JSON.read_text())


# ── the measurement product exists and is well-formed ─────────────────────────

def test_measurement_json_present_and_schema():
    d = _data()
    assert d['ticket'] == 'RYA-473' and d['star'] == 'solar'
    mn = d['Mn']
    assert mn['A_lte_median'] is not None
    assert len(mn['per_line']) == 3
    assert {round(float(p['line']), 2) for p in mn['per_line']} == {6013.51, 6016.67, 6021.82}


# ── HFS premise: the triplet is hyperfine-RESOLVED, gf = Den Hartog (cited) ────

def test_triplet_is_hfs_resolved_with_den_hartog_gf():
    hfs = _data()['Mn']['hfs_components']
    assert len(hfs) == 3
    for wl, _ in MN_LINES:
        d = hfs[str(round(wl, 3))]
        # multi-component = exactly why the single-profile EW path can't measure it
        assert d['n_components'] == 6
        # the GES gf-sum IS the RYA-468 Den Hartog canonical value (cited, not invented)
        assert d['gf_matches_dh'] is True
        assert abs(d['gf_sum_loggf'] - MN_DH_GF[wl]) < 0.02


def test_gf_and_hfs_provenance_is_den_hartog_not_invented():
    mn = _data()['Mn']
    assert 'Den Hartog' in mn['gf_provenance'] and 'RYA-468' in mn['gf_provenance']
    assert 'Den Hartog' in mn['hfs_provenance'] and 'Table 4' in mn['hfs_provenance']
    # the EW-saturation rationale (why synthesis, not EW) is recorded
    assert 'SAT' in mn['why_synthesis'] and 'RYA-468' in mn['why_synthesis']


# ── NLTE: MPIA/Bergemann grid δ, applied + flagged (never silent LTE) ──────────

def test_mn_nlte_applied_and_flagged():
    mn = _data()['Mn']
    assert mn['nlte_delta'] is not None
    # Mn NLTE is positive (Bergemann); the MPIA grid solar δ ~ +0.1
    assert 0.05 < mn['nlte_delta'] <= 0.15
    # A_nlte = A_lte + delta (actually applied, never silently dropped)
    assert abs(mn['A_nlte'] - (mn['A_lte_median'] + mn['nlte_delta'])) < 1e-3
    # the RYA-411 caveat is carried LOUD: the triplet is not a grid node, the
    # vendored MPIA δ is the high-EP reference value (the line-exact δ would differ)
    src = mn['nlte_source']
    assert 'RYA-411' in src and 'high-EP' in src
    assert mn['nlte_live'] is False              # Amarsi .grd offline in this worktree


# ── phase_c fold: Mn moves off "no value" → MEASURED, counts unchanged ─────────

def test_phase_c_fold_moves_mn_off_no_value():
    ov = V._mn_reclassify(_data())
    assert set(ov) == {'Mn'}
    assert ov['Mn']['A_measured'] is not None     # was None (no value) before
    assert ov['Mn']['n_lines'] == 3
    assert 'HFS synthesis' in ov['Mn']['channel'] and 'Den Hartog' in ov['Mn']['channel']
    # vendored MPIA δ + the RYA-411 caveat → honest CURATION-OWED, NOT a forced PASS
    assert ov['Mn']['verdict'] == 'CURATION-OWED'
    assert 'do NOT certify PASS' in ov['Mn']['owed']


def test_fold_keeps_verdict_counts_4_1_21_0():
    # Mn was CURATION-OWED (no value) and stays CURATION-OWED (now WITH a value):
    # no-value → measured, NOT a category flip / fabricated PASS.
    ab, ew, phase_a = V._load()
    rows = V.build_verdicts(ab, ew, phase_a)
    V._apply_kittpeak(rows, V._load_kittpeak())
    V._apply_cu_v_synthesis(rows, V._load_cu_v_synthesis())
    V._apply_mn_synthesis(rows, _data())
    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
    assert counts['PASS'] == 4
    assert counts['NLTE-OWED'] == 1
    assert counts['CURATION-OWED'] == 21
    assert counts.get('DATA-GAP', 0) == 0
    # Mn carries the measured value now
    mn = next(r for r in rows if r['element'] == 'Mn')
    assert mn['A_measured'] is not None and mn['n_lines'] == 3
    assert 'synthesis' in mn['provenance']


# ── validate-don't-tune: a real fit that scatters, nothing snapped to Asplund ──

def test_not_tuned_to_asplund():
    mn = _data()['Mn']
    # per-line A(Mn) genuinely scatters (a fit, not a constant) — σ > 0
    assert mn['scatter'] > 0.02
    a_vals = [r['A_X'] for r in mn['per_line'] if isinstance(r['A_X'], float)]
    assert max(a_vals) - min(a_vals) > 0.1        # 6016 low vs 6013/6021 high
    # the reported value is the fitted median, not the Asplund anchor
    assert mn['A_lte_median'] != SOLAR_ASPLUND2021['Mn']
    # the +offset survives the NLTE correction (not pulled to 5.42)
    assert mn['A_nlte'] - SOLAR_ASPLUND2021['Mn'] > V.TOL_PASS
