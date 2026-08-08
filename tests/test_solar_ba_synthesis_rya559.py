"""
tests/test_solar_ba_synthesis_rya559.py
=======================================
RYA-559 — wire the Ba II 5853 synthesis + Korotin2015 NLTE route so Ba emits a value.

Ba was owed-NO-value: the EW path SAT-culls the strong Ba II 5853 (REW -4.90, at the
knee, +HFS/isotope). The fix folds an HFS-resolved LTE synthesis (Turbospectrum bsyn,
inverting the observed solar EW) + the production Engine-A Korotin2015 1D-NLTE delta
into the phase_c verdict — moving Ba off blank while staying CURATION-OWED (the +0.14
is the blend-inflated pool EW, a gf/blend floor -> RYA-161/162, NOT tuned).

CI-safe (no synthesis): the measurement lives in scripts/ba5853_synth.py (Sirius); this
tests the committed result JSON + the phase_c fold. Mirrors test_mn_hfs_synthesis_rya473.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as V                                  # noqa: E402
from config.constants import SOLAR_ASPLUND2021                      # noqa: E402

JSON = ROOT / 'data' / 'results' / 'solar_ba_synthesis_rya559.json'


from tests.gold_scale_blocker import (  # noqa: E402  RYA-681/669, parameterised RYA-674
    verdict_gold_version, xfail_if_regeneration_blocked)

def _data():
    return json.loads(JSON.read_text())


# ── the measurement product exists and is well-formed ─────────────────────────

def test_measurement_json_present_and_schema():
    d = _data()
    assert d['element'] == 'Ba' and d['ion'] == 'II' and d['line_A'] == 5853.668
    assert d['A_nlte'] is not None and d['A_lte'] is not None
    # HFS-resolved synthesis, total gf validated to the canonical Ba II 5853 value.
    assert d['hfs_components'] >= 20
    assert abs(d['total_loggf'] - (-1.0)) < 0.02           # Davidson 1992 / NIST -1.00
    # Engine-A production Korotin delta, read at the solar node (validate-don't-tune).
    assert d['engineA_in_bounds'] is True
    assert -0.05 < d['engineA_korotin_delta'] < 0.0        # solar Ba II 5853 small-negative
    # A_nlte = A_lte + Korotin delta (the delta is APPLIED, not a separate guess).
    assert abs(d['A_nlte'] - (d['A_lte'] + d['engineA_korotin_delta'])) < 0.005


# ── the fold lands the value, off no-value, CURATION-OWED ──────────────────────

def test_ba_reclassify_lands_value_curation_owed():
    ov = V._ba_reclassify(_data())
    assert set(ov) == {'Ba'}
    assert ov['Ba']['A_measured'] is not None              # was None (no value) before
    assert abs(ov['Ba']['A_measured'] - _data()['A_nlte']) < 1e-6
    assert ov['Ba']['verdict'] == 'CURATION-OWED'          # sits high -> owed, NOT PASS
    assert 'synthesis' in ov['Ba']['channel'] and '5853' in ov['Ba']['channel']
    assert 'Korotin' in ov['Ba']['channel']


def test_owed_text_is_honest_not_tuned():
    ov = V._ba_reclassify(_data())
    owed = ov['Ba']['owed']
    # validate-don't-tune: the blend caveat + clean-EW cross-check must be surfaced,
    # and the value must NOT be shrunk toward Asplund.
    assert 'blend' in owed.lower()
    assert 'do NOT tune' in owed
    # Engine-B corroboration recorded (departure-level; 5853 GES-blocked).
    assert '4554' in owed or 'Engine-B' in owed


def test_ba_stays_off_pass_high():
    # The measured Ba sits ABOVE Asplund (blend-inflated EW). A Ba PASS / Ba≈0 delta
    # would mean tuning crept in — the value must carry the raw positive offset.
    d = _data()
    ov = V._ba_reclassify(d)
    delta = ov['Ba']['A_measured'] - SOLAR_ASPLUND2021['Ba']
    assert delta > 0.10                                    # off-anchor, curation-owed
    assert ov['Ba']['verdict'] != 'PASS'


# ── the fold keeps the verdict counts at 5 / 0 / 21 (Ba was already CURATION-OWED) ─

def test_fold_keeps_counts_5_0_21():
    # Ba was CURATION-OWED with NO value; the fold gives it a value but the category is
    # unchanged, so the 5 / 0 / 21 counts (RYA-556 N-wired baseline) are preserved.
    # RYA-674: regenerate against the gold version the COMMITTED artifact declares it
    # was built from — not a hardcoded CURRENT. Same named input, guard at full strength.
    gold = verdict_gold_version()
    xfail_if_regeneration_blocked(gold)
    ab, ew, phase_a, _gold = V._load(gold)
    rows = V.build_verdicts(ab, ew, phase_a)
    # the full main() fold chain — the 5th PASS is Mn (its synthesis fold, RYA-476).
    V._apply_kittpeak(rows, V._load_kittpeak())
    V._apply_cu_v_synthesis(rows, V._load_cu_v_synthesis())
    V._apply_mn_synthesis(rows, V._load_mn_synthesis())
    V._apply_s_synthesis(rows, V._load_s_synthesis())
    V._apply_ba_synthesis(rows, _data())
    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
    assert counts['PASS'] == 5
    assert counts.get('NLTE-OWED', 0) == 0
    assert counts['CURATION-OWED'] == 21
    assert counts.get('DATA-GAP', 0) == 0
    ba = next(r for r in rows if r['element'] == 'Ba')
    assert ba['verdict'] == 'CURATION-OWED'
    assert ba['A_measured'] is not None                    # off blank — the RYA-559 deliverable
    assert ba['delta_vs_asplund'] == round(ba['A_measured'] - ba['asplund2021'], 3)
