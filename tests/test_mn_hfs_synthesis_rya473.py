"""
tests/test_mn_hfs_synthesis_rya473.py
=====================================
RYA-473 — guard the HFS-resolved synthesis measurement of solar Mn on the Den
Hartog e6S→z6P triplet (6013/16/21). The REAL unlock RYA-468 pointed at: the gf is
graded (Den Hartog, MED) but the EW path SAT-culls these HFS-split lines, so Mn
stayed no-value; synthesis on the GES HFS-resolved line list measures it.

RYA-476 update: the Amarsi GALAH Mn grid is now staged, so the NLTE δ is the LIVE
triplet-exact HFS-resolved value (+0.024), replacing the vendored MPIA high-EP
+0.108. RYA-411 predicted this (the low-EP strong-HFS triplet desaturates to a
smaller δ than the grid's high-EP nodes) → A(Mn)_NLTE 5.470 (+0.050 vs Asplund),
inside TOL → Mn certifies **PASS**. Verdict 4/1/21/0 → 5/1/20/0.

Guards: the measurement product + schema; the HFS premise (6 components/feature,
NOT a collapsed single line); the gf provenance (GES gf-sum == the RYA-468 Den
Hartog canonical value — cited, not invented/tuned); the LIVE triplet-exact NLTE δ
(small + positive, the line-exact Amarsi value, never silent LTE); the phase_c fold
that moves Mn no-value → PASS; and validate-don't-tune (a real fit that scatters;
the small physical δ — not a tuned one — lands it in range).

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

def test_mn_nlte_applied_live_triplet_exact():
    mn = _data()['Mn']
    assert mn['nlte_delta'] is not None
    # RYA-476: the Amarsi grid is staged → LIVE triplet-exact δ, small + positive,
    # well below the MPIA high-EP +0.108 (RYA-411: low-EP HFS triplet desaturates).
    assert mn['nlte_live'] is True
    assert 0.0 < mn['nlte_delta'] < 0.06
    # A_nlte = A_lte + delta (actually applied, never silently dropped)
    assert abs(mn['A_nlte'] - (mn['A_lte_median'] + mn['nlte_delta'])) < 1e-3
    # provenance names the live Amarsi triplet-exact path (not the MPIA high-EP value)
    src = mn['nlte_source']
    assert 'Amarsi' in src and 'triplet' in src and 'live' in src


# ── phase_c fold: Mn moves no-value → PASS (live triplet-exact NLTE) ───────────

def test_phase_c_fold_moves_mn_to_pass():
    ov = V._mn_reclassify(_data())
    assert set(ov) == {'Mn'}
    assert ov['Mn']['A_measured'] is not None     # was None (no value) before
    assert ov['Mn']['n_lines'] == 3
    assert 'HFS synthesis' in ov['Mn']['channel'] and 'Den Hartog' in ov['Mn']['channel']
    # RYA-476: live triplet-exact NLTE δ brings A(Mn) within TOL → certified PASS
    # (the vendored-MPIA caveat that withheld PASS in RYA-473 no longer applies).
    assert ov['Mn']['verdict'] == 'PASS'
    assert 'closes the Mn gap' in ov['Mn']['owed']


def test_fold_moves_verdict_counts_to_5_0_21_0():
    # RYA-476: Mn was CURATION-OWED (no value), now MEASURED + reconciled via the
    # live triplet-exact NLTE → PASS. The category flip is earned by the measurement,
    # NOT fabricated: it rests on the line-exact Amarsi δ (+0.024), not a tune.
    # RYA-556: the N I NLTE grid is now wired in the KP channel, so the NLTE-OWED
    # bucket is empty (N -> CURATION-OWED); counts are 5 / 0 / 21 / 0.
    ab, ew, phase_a = V._load()
    rows = V.build_verdicts(ab, ew, phase_a)
    V._apply_kittpeak(rows, V._load_kittpeak())
    V._apply_cu_v_synthesis(rows, V._load_cu_v_synthesis())
    V._apply_mn_synthesis(rows, _data())
    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
    assert counts['PASS'] == 5
    assert counts.get('NLTE-OWED', 0) == 0                   # RYA-556: N I NLTE grid wired
    assert counts['CURATION-OWED'] == 21                     # +N (was NLTE-OWED)
    assert counts.get('DATA-GAP', 0) == 0
    # Mn carries the measured value + PASS now
    mn = next(r for r in rows if r['element'] == 'Mn')
    assert mn['verdict'] == 'PASS'
    assert mn['A_measured'] is not None and mn['n_lines'] == 3
    assert 'synthesis' in mn['provenance']


# ── validate-don't-tune: a real fit that scatters; small PHYSICAL δ, not a tune ─

def test_not_tuned_to_asplund():
    mn = _data()['Mn']
    # per-line A(Mn) genuinely scatters (a fit, not a constant) — σ > 0
    assert mn['scatter'] > 0.02
    a_vals = [r['A_X'] for r in mn['per_line'] if isinstance(r['A_X'], float)]
    assert max(a_vals) - min(a_vals) > 0.1        # 6016 low vs 6013/6021 high
    # the reported value is the fitted median, not the Asplund anchor
    assert mn['A_lte_median'] != SOLAR_ASPLUND2021['Mn']
    # NOT tuned: the A(Mn)_LTE measurement itself sits +0.026 off Asplund (a real
    # offset), and the LIVE NLTE δ is the small line-exact Amarsi value (+0.024) —
    # the PASS is earned by physics, not by snapping A(Mn) to 5.42.
    assert mn['A_lte_median'] - SOLAR_ASPLUND2021['Mn'] != 0
    assert mn['nlte_delta'] == round(mn['A_nlte'] - mn['A_lte_median'], 3)
    assert mn['nlte_delta'] < 0.06                # line-exact, not pulled to fit
