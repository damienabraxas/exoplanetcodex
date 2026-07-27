"""
tests/test_floor_promotion_rya561.py
====================================
RYA-561 — the ratified three-gate floor-promotion rule (durable helper).

Ryan's ruling (2026-07-27, RYA-561) — STRICT gate 3. A two-engine-FLOOR-governed
element earns PASS (CURATION-OWED -> PASS) iff ALL THREE hold:

  gate 1  its NLTE atom is RYA-534 anchor-validated (grid PASS on record)
  gate 2  |A(X) - Asplund2021| <= 0.10   (== the phase_c TOL_PASS, RYA-371)
  gate 3  a REAL cross-engine delta exists AND |dCE| <= 0.10

The firewall these pin: **a MISSING dCE FAILS gate 3 and may never be substituted by
the atom delta.** The atom delta reproducing the published anchor IS gate 1, so reusing
it as a gate-3 fallback is gate 1 under a second name — and substituting it in order to
make Mg promote is validate-don't-tune territory (RYA-161). Mg's path to PASS is the
real Mg I 5528 second line (RYA-592), not a rule relaxation.

These pin the LAW of the durable helpers — `nlte_atom_validation` (single source of
truth = the md5-pinned RYA-534 grid provenance, no hardcoded element list) and
`evaluate_floor_promotion`. The consumer-side enforcement in the RYA-527 re-emit ladder
(Ca promoted, Mg held at 7.614, Cr canary untouched, counts 6/1/19) is exercised
separately on the RYA-527 stack, where that script lives.
"""
import pytest

from pipeline import engine_selection as es

# The committed two-engine record (RYA-527 stack:
# data/audit/rya527_two_engine/solar_two_engine_records.json) — A(X), Asplund-2021, dCE.
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
    ONLY because the record has no cross-engine delta (Engine-A is suppressed for Mg
    by design: RYA-520 synth-required + the SAT-culled b-triplet). Ratified hold."""
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
    """Ni (dCE -1.253, INSIDE tol at +0.053), Si (+0.129), Ti (atom CHECK, and -0.090
    is INSIDE tol) and Sr II (no dCE, and -0.061 is INSIDE tol — the ratified
    RYA-428/551 near-UV hold) must all stay owed."""
    assert es.evaluate_floor_promotion(**rec, species=species).promoted is False


def test_ni_is_held_by_dce_despite_being_within_tol():
    """Ni is the guard working: +0.053 clears gate 2, the -1.253 engine split kills it."""
    p = es.evaluate_floor_promotion(**NI, species='Ni I')
    assert p.gate2_within_tol is True
    assert p.gate3_cross_engine is False


def test_sr_would_have_been_promoted_by_the_rejected_fallback():
    """Sr II is why STRICT is right: it clears gates 1+2 with dCE None, so an
    atom-delta fallback would have promoted it — silently reversing the ratified
    RYA-428/551 hold (~0.05-0.1 dex near-UV systematic)."""
    p = es.evaluate_floor_promotion(**SR, species='Sr II')
    assert (p.gate1_atom_validated, p.gate2_within_tol) == (True, True)
    assert p.gate3_cross_engine is False and p.promoted is False


def test_upper_limit_and_excluded_species_are_vetoed():
    """A ratified disposition outranks the gates: Li (UPPER_LIMIT, RYA-563) and
    Cr II (ratified exclusion, RYA-558/240) can never be promoted."""
    li = es.evaluate_floor_promotion(element='Li', a_value=1.05, reference_value=1.05,
                                     cross_engine_delta=0.0, species='Li I')
    assert li.promoted is False and 'UPPER_LIMIT' in li.reason
    cr2 = es.evaluate_floor_promotion(element='Cr', a_value=5.62, reference_value=5.62,
                                      cross_engine_delta=0.0, species='Cr II')
    assert cr2.promoted is False and 'excluded' in cr2.reason


def test_promotion_record_is_auditable():
    """Every decision serialises with its gate values + declared thresholds, so a
    verdict artifact can show WHY an element promoted or was held."""
    d = es.evaluate_floor_promotion(**CA, species='Ca I').as_dict()
    assert d['promoted'] is True
    assert d['delta_vs_reference'] == pytest.approx(0.072, abs=1e-4)
    assert d['cross_engine_delta'] == pytest.approx(0.016, abs=1e-4)
    assert d['thresholds'] == dict(es.FLOOR_PROMOTION)
    assert 'prov.json' in d['atom_citation']
