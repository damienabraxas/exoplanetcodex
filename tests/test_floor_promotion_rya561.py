"""
tests/test_floor_promotion_rya561.py
====================================
RYA-561 — the ratified three-gate floor-promotion rule (CURATION-OWED -> PASS).

Ryan's ruling (2026-07-27, RYA-561 comment) — STRICT gate 3. A two-engine-FLOOR-
governed element earns PASS iff ALL THREE hold:

  gate 1  its NLTE atom is RYA-534 anchor-validated (grid PASS on record)
  gate 2  |A(X) - Asplund2021| <= 0.10   (== the phase_c TOL_PASS, RYA-371)
  gate 3  a REAL cross-engine delta exists AND |dCE| <= 0.10

The firewall these tests exist to pin: **a MISSING dCE FAILS gate 3 and may never be
substituted by the atom delta.** The atom delta reproducing the published anchor IS
gate 1, so reusing it as a gate-3 fallback is gate 1 under a second name — and
substituting it in order to make Mg promote is validate-don't-tune territory
(RYA-161). Mg's path to PASS is the real Mg I 5528 second line (RYA-592).

On the committed two-engine record the rule promotes **Ca and only Ca**. Mg is a
ratified HOLD at CURATION-OWED-with-value 7.614 (atom validated, +0.064 inside the
band, but Engine-A is suppressed for Mg by design -> no dCE exists).
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

from pipeline import engine_selection as es

ROOT = Path(__file__).resolve().parent.parent
ARTIFACT = ROOT / 'data' / 'audit' / 'rya527_reemit' / 'proposed_gold_v3_diff.json'

# The committed two-engine record (data/audit/rya527_two_engine/solar_two_engine_records.json)
# — element, A(X), Asplund-2021, dCE, expected promotion.
CA = dict(element='Ca', a_value=6.372, reference_value=6.30, cross_engine_delta=0.016)
MG = dict(element='Mg', a_value=7.614, reference_value=7.55, cross_engine_delta=None)
NI = dict(element='Ni', a_value=6.253, reference_value=6.20, cross_engine_delta=-1.253)
SI = dict(element='Si', a_value=7.639, reference_value=7.51, cross_engine_delta=-0.756)
TI = dict(element='Ti', a_value=4.880, reference_value=4.97, cross_engine_delta=-1.326)
SR = dict(element='Sr', a_value=2.769, reference_value=2.83, cross_engine_delta=None)


# ── gate 1: the RYA-534 atom-validation source ─────────────────────────────────

def test_gate1_reads_the_534_provenance_not_a_hardcoded_list():
    """Single source of truth = the md5-pinned RYA-534 grid provenance gate verdict."""
    ok, citation = es.nlte_atom_validation('Mg')
    assert ok is True
    assert 'prov.json' in citation and 'PASS' in citation


def test_gate1_fails_on_a_check_verdict():
    """Ti's Engine-B atom is a CHECK (not PASS) — RYA-548 xfail. Gate 1 must fail."""
    ok, citation = es.nlte_atom_validation('Ti')
    assert ok is False
    assert 'CHECK' in citation


def test_gate1_fails_when_no_grid_is_on_record():
    """An element with no Engine-B Gerber grid (S) fails gate 1 honestly — never a
    silent pass, and never a raise (an absent grid is a real, documented state)."""
    ok, citation = es.nlte_atom_validation('S')
    assert ok is False
    assert 'no RYA-534' in citation


# ── the ratified thresholds ────────────────────────────────────────────────────

def test_declared_thresholds():
    assert es.FLOOR_PROMOTION['tol_pass_dex'] == 0.10       # == phase_c TOL_PASS (RYA-371)
    assert es.FLOOR_PROMOTION['cross_engine_dex'] == 0.10   # == RYA-525 cross_engine_mix_gate


# ── the rule, element by element ───────────────────────────────────────────────

def test_ca_promotes_on_all_three_gates():
    p = es.evaluate_floor_promotion(**CA, species='Ca I')
    assert (p.gate1_atom_validated, p.gate2_within_tol, p.gate3_cross_engine) == (True, True, True)
    assert p.promoted is True
    # the promotion must carry the 1D-vs-3D scale caveat (Ryan, RYA-561)
    assert es.FLOOR_PROMOTION_SCALE_CAVEAT in p.reason


def test_mg_is_held_by_gate3_alone():
    """Mg passes gates 1+2 — atom validated, +0.064 inside the band — and is held
    ONLY because the record has no cross-engine delta. This is the ratified hold."""
    p = es.evaluate_floor_promotion(**MG, species='Mg I')
    assert (p.gate1_atom_validated, p.gate2_within_tol) == (True, True)
    assert p.gate3_cross_engine is False
    assert p.promoted is False
    assert p.cross_engine_delta is None


def test_missing_dce_never_passes_gate3_even_with_a_perfect_atom():
    """THE FIREWALL (RYA-561/161): a missing dCE fails gate 3 no matter how well the
    NLTE atom reproduces the anchor. Gate 1 may not stand in for gate 3."""
    p = es.evaluate_floor_promotion(element='Mg', a_value=7.55, reference_value=7.55,
                                    cross_engine_delta=None, species='Mg I')
    assert p.gate1_atom_validated is True      # atom validated
    assert p.gate2_within_tol is True          # dead-on the reference
    assert p.gate3_cross_engine is False       # ...and still not promoted
    assert p.promoted is False


@pytest.mark.parametrize('rec,species', [(NI, 'Ni I'), (SI, 'Si I'), (TI, 'Ti I'), (SR, 'Sr II')])
def test_the_guard_holds_every_other_candidate(rec, species):
    """Ni (dCE -1.253, INSIDE tol at +0.053), Si (+0.129), Ti (atom CHECK) and
    Sr II (no dCE, RYA-428/551 near-UV systematic) must all stay owed."""
    assert es.evaluate_floor_promotion(**rec, species=species).promoted is False


def test_ni_is_held_by_dce_despite_being_within_tol():
    """Ni is the guard working: +0.053 clears gate 2, the -1.253 engine split kills it."""
    p = es.evaluate_floor_promotion(**NI, species='Ni I')
    assert p.gate2_within_tol is True
    assert p.gate3_cross_engine is False


def test_upper_limit_and_excluded_species_are_vetoed():
    """A ratified disposition outranks the gates: Li (UPPER_LIMIT, RYA-563) and
    Cr II (ratified exclusion, RYA-558/240) can never be promoted."""
    li = es.evaluate_floor_promotion(element='Li', a_value=1.05, reference_value=1.05,
                                     cross_engine_delta=0.0, species='Li I')
    assert li.promoted is False and 'UPPER_LIMIT' in li.reason
    cr2 = es.evaluate_floor_promotion(element='Cr', a_value=5.62, reference_value=5.62,
                                      cross_engine_delta=0.0, species='Cr II')
    assert cr2.promoted is False and 'excluded' in cr2.reason


# ── end-to-end: the re-emitted verdict ────────────────────────────────────────

@pytest.fixture(scope='module')
def payload():
    subprocess.run([sys.executable, 'scripts/rya527_reemit_verdict.py'], cwd=ROOT, check=True)
    return json.loads(ARTIFACT.read_text())


def test_reemit_promotes_ca_and_only_ca(payload):
    assert payload['floor_promotion_rule']['promoted'] == ['Ca']


def test_reemit_ca_is_pass_with_a_sourced_citation(payload):
    ca = next(r for r in payload['diff_table'] if r['element'] == 'Ca')
    assert ca['verdict'] == 'PASS'
    assert ca['v3_proposed'] == 6.372
    assert ca['floor_promotion']['applied'] is True
    assert ca['floor_promotion']['n_lines'] == 4        # 4 line-winners, sourced (RYA-561 check)
    assert 'RYA-561 PASS' in ca['source'] and '3D' in ca['source']


def test_reemit_mg_stays_curation_owed_with_its_value(payload):
    """Mg keeps the honest value 7.614 but NOT the PASS tier."""
    mg = next(r for r in payload['diff_table'] if r['element'] == 'Mg')
    assert mg['verdict'] == 'CURATION-OWED'
    assert mg['v3_proposed'] == 7.614
    assert mg['floor_promotion']['applied'] is False
    assert mg['floor_promotion']['cross_engine_delta'] is None


def test_reemit_leaves_the_cr_canary_untouched(payload):
    """The validate-don't-tune Cr canary tracks the phase_c graded-EW-pool Cr I
    (+0.402, RYA-398). Cr is a RATIFIED channel (RYA-558) so the promotion rule
    cannot reach it — the tuning detector must not drift."""
    cr = next(r for r in payload['diff_table'] if r['element'] == 'Cr')
    assert cr['v3_proposed'] == 6.022
    assert abs(cr['v3_delta_vs_asplund'] - 0.402) < 0.01
    assert cr['verdict'] == 'CURATION-OWED'
    assert cr.get('floor_promotion') is None           # never even evaluated


def test_reemit_holds_ni_si_na_ti_sr(payload):
    by = {r['element']: r for r in payload['diff_table']}
    for el in ('Ni', 'Si', 'Na', 'Ti', 'Sr'):
        assert by[el]['verdict'] == 'CURATION-OWED', el
        assert not (by[el].get('floor_promotion') or {}).get('applied'), el


def test_reemit_counts_move_by_exactly_one(payload):
    """Before: 5 PASS / 1 NLTE-OWED / 20 CURATION-OWED. The rule moves exactly Ca."""
    assert payload['counts'] == {'PASS': 6, 'NLTE-OWED': 1, 'CURATION-OWED': 19}
