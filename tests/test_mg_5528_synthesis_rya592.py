"""RYA-592 — Mg I 5528 second-line measurement: provenance + verdict invariants.

These lock the things that must not drift silently:
  * the Step-0 gf provenance in the SINGLE source (value, grade, reference);
  * the NIST accuracy ladder used to tier a gf (the '+' tiers were missing);
  * the committed measurement artefact's shape and its OWN internal consistency;
  * the disposition logic — a DISCORDANT second line must NOT promote Mg.

Nothing here re-runs the Sirius synthesis (grids stay on Sirius, RYA-567); they read
the committed artefact, which is the point of committing it.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RESULT = REPO / 'data' / 'results' / 'mg_5528_synthesis_rya592.json'
CANON = REPO / 'data' / 'linelists' / 'canonical_gf.csv'

TARGET = '5528.405'
REFLINE = '5711.088'


@pytest.fixture(scope='module')
def res():
    if not RESULT.exists():
        pytest.skip(f"{RESULT.name} not committed")
    return json.loads(RESULT.read_text())


# ── Step 0: gf provenance in the single source ───────────────────────────────

def _canon_row(wave):
    import csv
    with CANON.open() as f:
        for row in csv.DictReader(f):
            if row['species'] == 'Mg I' and abs(float(row['wavelength_air_A']) - wave) < 0.01:
                return row
    return None


@pytest.mark.parametrize('wave,loggf,grade', [(5528.405, -0.498, 'B+'),
                                              (5711.088, -1.724, 'B')])
def test_canonical_gf_value_and_grade(wave, loggf, grade):
    """The values re-derive from NIST ASD (log gf = log10(g_i * f_ik)); the GRADES are
    what ASD actually reports. Both lines carried grade 'A' before RYA-592 — the source
    says B+ / B. If this fails, the single source drifted off its cited primary."""
    row = _canon_row(wave)
    assert row is not None, f"Mg I {wave} missing from canonical_gf.csv"
    assert abs(float(row['log_gf']) - loggf) < 1e-6
    assert row['nist_grade'] == grade
    assert 'NIST ASD' in row['loggf_reference']
    # the cross-check source must stay recorded alongside, so the gf-scale systematic
    # is never silently forgotten
    assert 'Pehlivan Rhodin' in row['loggf_reference']


def test_both_lines_clear_the_grade_policy():
    """B+ (<=7%) and B (<=10%) both clear PIPELINE['min_nist_grade'] = 'B', so the
    grade correction must not have cost either line its usability."""
    from pipeline.curate_nonfe_pools import NIST_GRADE_HIGH
    for wave in (5528.405, 5711.088):
        assert _canon_row(wave)['nist_grade'] in NIST_GRADE_HIGH


def test_nist_grade_ladder_is_monotonic():
    """B+ is <=7%, B is <=10% — B+ must never tier BELOW B. The shipped set omitted the
    '+' tiers, which silently demoted a B+ line out of HIGH (an inverted ladder)."""
    from pipeline import curate_nonfe_pools as cur
    for g in ('AAA', 'AA', 'A+', 'A', 'B+', 'B'):
        assert cur._gf_tier(g, 'NIST ASD') == 'HIGH', g
    for g in ('D+', 'D', 'E'):
        assert cur._gf_tier(g, 'NIST ASD') == 'CULL', g
    for g in ('C+', 'C'):
        assert cur._gf_tier(g, 'NIST ASD') == 'MED', g


# ── the measurement artefact ─────────────────────────────────────────────────

def test_target_line_is_the_second_line(res):
    assert res['_meta']['target_line'] == 5528.405
    assert res['_meta']['concordance_reference_line'] == 5711.088


def test_both_engine_controls_reproduced(res):
    """validate-don't-tune: each engine re-derived its OWN committed 5711 delta before
    its new 5528 delta was used. A drifting control invalidates the 5528 number."""
    for eng in ('engine_a', 'engine_b'):
        e = res['_meta'][eng]
        assert abs(e['control_rederived'] - e['control_committed']) <= e['control_tol'], eng


def test_blend_actually_modelled_in_window(res):
    """The RYA-592 CRITICAL check: this is an in-window blend FIT, not an isolated-line
    EW inversion. If the blend were culled the component count would collapse."""
    assert res[TARGET]['n_blend_components_in_window'] > 50
    assert 'NOT an isolated-line EW inversion' in res['_meta']['method']


def test_fit_quality_and_no_railing(res):
    """A railed fit or a blown-up chi2 means the value is not measured — do not ship it."""
    for line in (TARGET, REFLINE):
        for arm in ('harps', 'iag'):
            a = res[line][arm]
            assert not a['railed'], f"{line}/{arm} railed"
            assert not a['gsig_railed'], f"{line}/{arm} gsig railed"
            assert a['red_chi2'] < 5.0, f"{line}/{arm} red_chi2 {a['red_chi2']}"


def test_fitted_velocity_matches_the_independent_measurement(res):
    """dv is a nuisance parameter, so it must agree with the abundance-BLIND centroid
    measurement — otherwise the fit is absorbing a profile mismatch into a shift."""
    for line in (TARGET, REFLINE):
        for arm in ('harps', 'iag'):
            a = res[line][arm]
            assert abs(a['dv_fitted_kms'] - a['dv_measured_kms']) < 0.35, f"{line}/{arm}"


def test_target_cleared_the_reliability_floor(res):
    t = res[TARGET]['harps']
    assert t['dEW_dA_mA_dex'] >= res['_meta']['reliable_dEW_dA_floor']
    assert t['reliable'] is True


def test_two_arms_agree_on_each_line(res):
    """HARPS and the IAG FTS atlas are independent solar arms; a per-line disagreement
    between them would mean the ARM, not the line, is driving the result."""
    for line in (TARGET, REFLINE):
        d = abs(res[line]['harps']['A_LTE'] - res[line]['iag']['A_LTE'])
        assert d < 0.05, f"{line} arm-to-arm {d:.3f}"


# ── the verdict: a discordant second line must NOT promote Mg ────────────────

def test_lines_are_discordant_and_mg_is_not_promoted(res):
    """The measured outcome: 5528 and 5711 disagree by ~0.21-0.23 dex, far outside the
    0.10 band. The second line therefore CONTRADICTS rather than confirms 5711, and Mg
    must stay owed. This test is the guard against a future silent promotion on this
    evidence — if the discordance is ever resolved, this expectation changes DELIBERATELY."""
    v = res['_verdict']
    assert v['concordance_worst_abs_dex'] > v['concordance_band']
    assert v['lines_concordant'] is False
    assert v['disposition'].startswith('CURATION-OWED')


def test_no_real_cross_engine_dce_from_a_saturated_line(res):
    """RYA-561 gate 3 (STRICT) needs a REAL cross-engine dCE. 5528 is far above the
    ratified saturation knee, so the Engine-A (EW->A) route is inadmissible and Mg
    stays single-engine. Gate 3 must read FAIL, never satisfied by atom-delta."""
    v = res['_verdict']
    assert v['engineA_ew_route_admissible'] is False
    assert v['real_cross_engine_dCE_available'] is False
    assert res[TARGET]['harps']['obs_ew_mA'] > res['_meta']['saturation_knee_mA']


def test_atom_delta_agreement_is_labelled_diagnostic_only(res):
    """The Engine-A/Engine-B NLTE atoms agree closely on 5528 — that agreement is
    explicitly NOT gate-3 evidence (RYA-561 rejected the atom-delta fallback)."""
    assert 'DIAGNOSTIC ONLY' in res['_meta']['gate3_note']
    assert res[TARGET]['nlte_delta_atom_spread'] < 0.02


# ── the RYA-592 evidence against the MERGED RYA-561 guard ────────────────────

def test_rya592_evidence_still_holds_mg_under_the_ratified_guard(res):
    """End-to-end tie-in: feed this ticket's actual measurement to the ratified
    three-gate rule and confirm Mg is still HELD. RYA-561's guard names RYA-592 as
    "Mg's path to PASS"; this is the honest answer — the second line does not open it.

    Gate 3 fails for a sourced reason: 5528 is saturated far past the ratified knee, so
    the Engine-A EW route is inadmissible and dCE stays None."""
    from pipeline.engine_selection import evaluate_floor_promotion
    from config.constants import SOLAR_ASPLUND2021
    rec = evaluate_floor_promotion(
        'Mg',
        a_value=res['_verdict']['target_A_NLTE_engineB'],
        reference_value=SOLAR_ASPLUND2021['Mg'],
        cross_engine_delta=None,          # no Engine-A on a saturated line
        species='Mg I')
    assert rec.promoted is False
    assert rec.gate1_atom_validated is True      # the 534 atom IS validated
    assert rec.gate3_cross_engine is False       # ...but there is still no real dCE
    assert 'gate3' in rec.reason


def test_a_concordant_second_line_would_still_not_self_promote(res):
    """Guard against the tempting shortcut: even if the two lines HAD agreed, the
    >=2-concordant-line variant of gate 3 is explicitly NOT ratified (RYA-561 flagged it
    for separate ratification). A concordant line count must never be passed to the
    guard as a cross-engine delta."""
    from pipeline.engine_selection import evaluate_floor_promotion
    from config.constants import SOLAR_ASPLUND2021
    rec = evaluate_floor_promotion('Mg', a_value=SOLAR_ASPLUND2021['Mg'],
                                   reference_value=SOLAR_ASPLUND2021['Mg'],
                                   cross_engine_delta=None, species='Mg I')
    assert rec.promoted is False, "a perfect value with no dCE must still be held"
