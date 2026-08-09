"""
tests/test_ratified_constraints.py — RYA-674
============================================
THE THREE KNOWN VIOLATIONS, EXERCISED AS VIOLATIONS.

Each of the three ratified constraints in `pipeline/ratified_constraints.py` was
violated for real, in a committed artifact, by a module that did not know the
constraint existed (RYA-669). This file drives each of those three exact rows through
the emission-time gate and asserts it now raises — the failing direction exercised, not
argued for in a comment.

It also drives the LIVE artifacts through the gate and asserts they are clean, because
a guard that fires on everything is not a guard.
"""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import corrections_registry as cr                     # noqa: E402
from pipeline.ratified_constraints import (                          # noqa: E402
    RATIFIED_CONSTRAINTS, ConstraintType, RatifiedConstraintViolation, RowKind,
    assert_ratified_constraints_satisfied, check_ratified_constraints, constraint_by_id)
from pipeline.solar_scale_provenance import (                        # noqa: E402
    CORRECTIONS_APPLIED_COLUMN, FE_1D3D_CORRECTION_ID, SCALE_3D_NLTE,
    SCALE_STATE_COLUMN, encode_corrections_applied, method_scale_label)

EXPECTED_IDS = {'Li_6707_veto_1_409', 'Cr_II_species_exclusion',
                'Fe_1D_3D_correction_required_on_solar_report'}

TWO_ENGINE = ROOT / 'data' / 'audit' / 'rya527_two_engine' / 'solar_two_engine_records.json'
VERDICT = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'
DISPOSITION = ROOT / 'data' / 'audit' / 'element_disposition_rya663.json'
PROPOSED_V3 = ROOT / 'data' / 'audit' / 'rya527_reemit' / 'proposed_gold_v3_diff.json'


# ── the registry itself ──────────────────────────────────────────────────────

def test_all_three_known_constraints_are_registered():
    ids = {c.constraint_id for c in RATIFIED_CONSTRAINTS}
    assert EXPECTED_IDS <= ids, f"missing: {EXPECTED_IDS - ids}"


def test_every_constraint_cites_a_ratifying_ticket():
    """A constraint with no ratifying RYA-# is an opinion, not a ratification."""
    for c in RATIFIED_CONSTRAINTS:
        assert 'RYA-' in c.provenance_ticket, c.constraint_id
        assert c.constraint_type in vars(ConstraintType).values(), c.constraint_id
        assert c.scope and c.check_fn is not None


def test_constraints_are_active_and_lookup_works():
    for cid in EXPECTED_IDS:
        c = constraint_by_id(cid)
        assert c.active, f"{cid} is revoked by {c.revoked_by} — a revocation needs a ticket"
    with pytest.raises(KeyError):
        constraint_by_id('not_a_constraint')


# ── failure case 1: Li 1.409 (RYA-563 / RYA-669) ─────────────────────────────

LI_FLOOR_RECORD = {'element': 'Li', 'ion': 'I', 'reported': 1.409, 'engineA': 0.727,
                   'engineB': 1.409, 'selected_engines': ['engineB_synth'],
                   'mean_cross_engine_delta': 0.682}


def test_li_1_409_as_a_floor_record_is_caught():
    """The exact row in the committed July two-engine artifact."""
    with pytest.raises(RatifiedConstraintViolation, match=r'Li_6707_veto_1_409'):
        assert_ratified_constraints_satisfied(
            [LI_FLOOR_RECORD], 'test: two-engine emitter', RowKind.SPECIES_RECORD)


def test_li_1_409_adopted_as_the_reported_element_value_is_caught():
    """The RYA-669 proposed-gold row: `1.409 NEW ... two-engine synthesis floor`."""
    row = {'element': 'Li', 'ion': 'I', 'A_measured': 1.409,
           'source': 'two-engine synthesis floor (I, B_synth)', 'verdict': 'CURATION-OWED'}
    with pytest.raises(RatifiedConstraintViolation, match=r'Li_6707_veto_1_409'):
        assert_ratified_constraints_satisfied([row], 'test: verdict re-emit ladder')


def test_the_ratified_demotion_is_what_makes_the_floor_record_legal():
    """RYA-563's own remedy: the floor value becomes a DIAGNOSTIC, not a deletion."""
    demoted = dict(LI_FLOOR_RECORD, reported=None, diagnostic_only=True,
                   diagnostic_value=1.409,
                   diagnostic_reason='upper_limit disposition (RYA-563)')
    assert_ratified_constraints_satisfied(
        [demoted], 'test: two-engine emitter', RowKind.SPECIES_RECORD)


def test_the_phase_c_upper_limit_itself_is_NOT_a_violation():
    """Li 0.727 carried as an upper limit is the ratified treatment, not the leak. A
    guard that forbade it would forbid the correct answer."""
    row = {'element': 'Li', 'ion': 'I', 'A_measured': 0.727, 'verdict': 'CURATION-OWED',
           'channel': 'EW: Li I 6707 (single line, UPPER LIMIT)',
           'provenance': 'harps-measured (EW pool)'}
    assert_ratified_constraints_satisfied([row], 'test: phase_c')


def test_an_upper_limit_element_may_never_be_emitted_as_PASS():
    row = {'element': 'Li', 'ion': 'I', 'A_measured': 0.727, 'verdict': 'PASS',
           'channel': 'EW: Li I 6707'}
    with pytest.raises(RatifiedConstraintViolation, match=r'never a PASS point value'):
        assert_ratified_constraints_satisfied([row], 'test: phase_c')


def test_upper_limit_membership_is_registry_sourced_not_hardcoded():
    """RYA-563's permanent rule: single source = problem_children.csv, never `{'Li'}`."""
    src = (ROOT / 'pipeline' / 'ratified_constraints.py').read_text(encoding='utf-8')
    assert 'is_upper_limit_disposition' in src
    # the check body must not test element identity against a literal
    body = src.split('def _check_upper_limit_veto')[1].split('\ndef ')[0]
    assert "== 'Li'" not in body and '"Li"' not in body
    fe_ok = {'element': 'Fe', 'ion': 'I', 'A_measured': 7.466,
             'source': 'two-engine synthesis floor'}
    # Fe has no upper_limit disposition, so the same floor-sourced shape is fine for it
    # (its own constraint is checked separately).
    assert not any('Li_6707_veto' in v for v in
                   check_ratified_constraints([fe_ok], 'test'))


# ── failure case 2: Cr II 5.676 (RYA-240/558 / RYA-669) ──────────────────────

CR_II_FLOOR_RECORD = {'element': 'Cr', 'ion': 'II', 'reported': 5.676, 'engineA': 8.354,
                      'engineB': 5.676, 'selected_engines': ['engineB_synth']}


def test_cr_ii_5_676_as_a_floor_record_is_caught():
    with pytest.raises(RatifiedConstraintViolation, match=r'Cr_II_species_exclusion'):
        assert_ratified_constraints_satisfied(
            [CR_II_FLOOR_RECORD], 'test: two-engine emitter', RowKind.SPECIES_RECORD)


def test_cr_ii_proposed_as_the_reported_element_value_is_caught():
    row = {'element': 'Cr', 'ion': 'II', 'A_X': 5.676, 'confidence': 'gold',
           'source': 'two-engine floor'}
    with pytest.raises(RatifiedConstraintViolation, match=r'RATIFIED-EXCLUDED'):
        assert_ratified_constraints_satisfied([row], 'test: gold builder')


def test_cr_i_is_the_ratified_reported_ion_and_passes():
    row = {'element': 'Cr', 'ion': 'I', 'A_X': 6.022, 'confidence': 'gf_floor'}
    assert_ratified_constraints_satisfied([row], 'test: gold builder')


def test_cr_ii_kept_as_a_marked_diagnostic_is_legal():
    """RYA-558 ratified keeping Cr II as a cross-engine DIAGNOSTIC. The guard must
    permit exactly that, or it would forbid the ratified treatment."""
    demoted = dict(CR_II_FLOOR_RECORD, reported=None, diagnostic_only=True,
                   diagnostic_value=5.676)
    assert_ratified_constraints_satisfied(
        [demoted], 'test: two-engine emitter', RowKind.SPECIES_RECORD)


def test_exclusion_scope_is_single_sourced_from_engine_selection():
    from pipeline.engine_selection import RATIFIED_EXCLUDED_SPECIES
    assert 'Cr II' in RATIFIED_EXCLUDED_SPECIES


# ── failure case 3: the Fe correction (RYA-553 / RYA-681 / RYA-669) ──────────

A_CORRECT, A_DOUBLED, A_UNCORRECTED = 7.466, 7.416, 7.516


def test_the_doubled_fe_anchor_is_caught():
    """RYA-669's 7.416 — on neither scale. FE_GATE cannot catch it (the correction
    magnitude equals the gate half-width); this does."""
    row = {'element': 'Fe', 'ion': 'I', 'A_measured': A_DOUBLED}
    with pytest.raises(RatifiedConstraintViolation, match=r'NEITHER recognised'):
        assert_ratified_constraints_satisfied([row], 'test: phase_c')


def test_gold_v3s_shape_a_3d_value_under_a_1d_label_is_caught():
    row = {'element': 'Fe', 'ion': 'I', 'A_X': A_CORRECT,
           'method_scale': '1D-NLTE (Fe I)'}
    with pytest.raises(RatifiedConstraintViolation, match=r'CONTRADICTS ITSELF'):
        assert_ratified_constraints_satisfied([row], 'test: gold builder')


def test_a_self_consistent_1d_anchor_is_still_a_violation_when_REPORTED():
    """gold v2's Fe row is honest but pre-correction. Frozen history is fine; EMITTING
    it as a reported solar anchor is the missing-correction case."""
    row = {'element': 'Fe', 'ion': 'I', 'A_measured': A_UNCORRECTED,
           'method_scale': '1D-NLTE (Fe I)'}
    with pytest.raises(RatifiedConstraintViolation, match=r'must carry the tabulated'):
        assert_ratified_constraints_satisfied([row], 'test: phase_c')


def test_the_v4_shape_passes():
    """What RYA-677 should freeze: the 3D value, the corrections_applied declaration,
    and a derived label."""
    row = {'element': 'Fe', 'ion': 'I', 'A_X': A_CORRECT,
           'method_scale': method_scale_label('Fe', SCALE_3D_NLTE),
           SCALE_STATE_COLUMN: SCALE_3D_NLTE,
           CORRECTIONS_APPLIED_COLUMN: encode_corrections_applied([FE_1D3D_CORRECTION_ID])}
    assert_ratified_constraints_satisfied([row], 'test: gold builder')


def test_a_species_record_is_not_a_claim_about_the_reported_anchor():
    """The floor's raw Fe I leg sits at 7.58 — on neither scale, and correctly so: it is
    a pre-report engine measurement, not the reported anchor. The Fe constraint applies
    to ELEMENT_VALUE rows only, or it would fire on every honest diagnostic."""
    rec = {'element': 'Fe', 'ion': 'I', 'reported': 7.58, 'engineA': 7.475}
    assert_ratified_constraints_satisfied(
        [rec], 'test: two-engine emitter', RowKind.SPECIES_RECORD)


# ── the corrections registry (RYA-674 §2A/§2B) ───────────────────────────────

def test_corrections_registry_loads_and_every_pointer_resolves():
    cr.assert_registry_consistent()
    assert '1D_3D_solar_Fe_Magic2013' in cr.correction_ids()


def test_the_registry_tabulates_no_magnitude_of_its_own():
    """THE anti-pattern this ticket exists to kill. The YAML may name where a number
    lives; it may not carry a copy of one. A tabulated magnitude / pre_range /
    post_range would be a second hand-maintained copy of the Magic-2013 offset."""
    entry = cr.correction('1D_3D_solar_Fe_Magic2013')
    for banned in ('magnitude', 'pre_range', 'post_range', 'value'):
        assert banned not in entry, (
            f"{banned!r} is tabulated in the corrections registry — it must be a "
            f"*_source pointer into config/constants.py (RYA-674)")
    from config.constants import CORRECTIONS_3D
    assert cr.magnitude('1D_3D_solar_Fe_Magic2013') == CORRECTIONS_3D['Fe_1D3D_solar_dex']


def test_the_scale_bands_are_derived_and_exactly_adjacent():
    """RYA-674 §2B's pre_range / post_range, computed rather than typed. The two bands
    must touch and not overlap, so no value is ever on both scales or in a gap."""
    bands = cr.scale_bands('1D_3D_solar_Fe_Magic2013')
    post_lo, post_hi = bands['3D-NLTE']
    pre_lo, pre_hi = bands['1D-NLTE']
    assert post_hi == pytest.approx(pre_lo)
    assert post_lo < A_CORRECT < post_hi
    assert pre_lo < A_UNCORRECTED < pre_hi
    assert not (post_lo <= A_DOUBLED <= post_hi or pre_lo <= A_DOUBLED <= pre_hi)


def test_value_check_is_required_for_the_fe_correction():
    """§2B's trigger, COMPUTED: |−0.05| >= 0.5 × FE_GATE's 0.05 half-width, so the gate
    alone cannot catch a doubled application and the value must be classified."""
    assert cr.value_check_required('1D_3D_solar_Fe_Magic2013') is True


def test_an_unknown_correction_identifier_raises():
    with pytest.raises(cr.CorrectionsRegistryError):
        cr.correction('not_a_correction')


# ── the live artifacts are CLEAN (a guard that fires on everything is no guard) ──

@pytest.mark.skipif(not VERDICT.exists(), reason='verdict artifact absent')
def test_the_committed_phase_c_verdict_satisfies_every_constraint():
    rows = json.loads(VERDICT.read_text(encoding='utf-8'))['verdicts']
    assert_ratified_constraints_satisfied(rows, 'committed phase_c verdict artifact')


@pytest.mark.skipif(not DISPOSITION.exists(), reason='disposition artifact absent')
def test_the_committed_disposition_report_satisfies_every_constraint():
    rows = json.loads(DISPOSITION.read_text(encoding='utf-8'))['dispositions']
    assert_ratified_constraints_satisfied(rows, 'committed disposition report')


@pytest.mark.skipif(not TWO_ENGINE.exists(), reason='two-engine artifact absent')
def test_the_committed_two_engine_artifact_still_carries_both_known_leaks():
    """The July artifact is the EVIDENCE, not a regression: it predates the guard and is
    read as a review artifact (RYA-663 marks it stale). Pinning both leaks here means
    the day it is regenerated, this test tells us whether the demotion actually ran."""
    records = json.loads(TWO_ENGINE.read_text(encoding='utf-8'))['records']
    found = check_ratified_constraints(
        records, 'committed two-engine artifact', RowKind.SPECIES_RECORD)
    ids = {cid for cid in EXPECTED_IDS if any(cid in v for v in found)}
    assert ids == {'Li_6707_veto_1_409', 'Cr_II_species_exclusion'}, found
    assert len(found) == 2, found


@pytest.mark.skipif(not PROPOSED_V3.exists(), reason='RYA-527 re-emit artifact absent')
def test_the_historical_proposed_gold_v3_table_is_caught_at_element_level():
    """THE artifact the ticket was written about. `data/audit/rya527_reemit/
    proposed_gold_v3_diff.json` proposed A(Li) 1.409 "from the two-engine synthesis
    floor" and A(Cr II) 5.676 as element-level gold values. Both are now caught by the
    gate that `scripts/rya527_reemit_verdict.py` calls before writing that table."""
    rows = json.loads(PROPOSED_V3.read_text(encoding='utf-8'))['diff_table']
    found = check_ratified_constraints(rows, 'RYA-527 proposed gold v3')
    ids = {cid for cid in EXPECTED_IDS if any(cid in v for v in found)}
    assert ids == {'Li_6707_veto_1_409', 'Cr_II_species_exclusion'}, found
    with pytest.raises(RatifiedConstraintViolation):
        assert_ratified_constraints_satisfied(rows, 'RYA-527 proposed gold v3')


# ── the gate's own contract ──────────────────────────────────────────────────

def test_a_row_the_gate_cannot_understand_is_never_silently_skipped():
    with pytest.raises(RatifiedConstraintViolation, match=r'cannot check a candidate row'):
        assert_ratified_constraints_satisfied(['not a row'], 'test')


def test_the_failure_message_names_module_constraint_ticket_and_row():
    with pytest.raises(RatifiedConstraintViolation) as exc:
        assert_ratified_constraints_satisfied(
            [CR_II_FLOOR_RECORD], 'my emitter', RowKind.SPECIES_RECORD)
    msg = str(exc.value)
    for needle in ('my emitter', 'Cr_II_species_exclusion', 'RYA-240', 'Cr II', '5.676'):
        assert needle in msg, needle


def test_every_emission_path_calls_the_gate():
    """The whole design: the check is a property of emitting, so a path that writes a
    per-element result and does not call it is the next RYA-669."""
    paths = ('scripts/phase_c_verdict_rya371.py',
             'scripts/build_solar_reference_v2_rya522.py',
             'pipeline/element_disposition.py',
             'scripts/rya527_two_engine_run.py',
             'scripts/rya527_reemit_verdict.py',
             'scripts/rya669_phase2_reemit.py')
    for rel in paths:
        src = (ROOT / rel).read_text(encoding='utf-8')
        assert 'assert_ratified_constraints_satisfied(' in src, (
            f"{rel} emits per-element results but never invokes the RYA-674 gate")


# ── RYA-691 / RYA-701: an unreliable value must not be emitted ────────────────
# The rule is RYA-679's (`not railed AND dEW_dA >= floor`). RYA-691 made the
# two-engine loader honour it at read time; Ryan ratified it as a RYA-674 constraint
# on 2026-08-08 so it also holds at EMISSION time, in any module.

from pipeline.reliability_contract import (                          # noqa: E402
    GATED_PREFIX, ReliabilityState, classify_reliability_basis, reliability_basis)


def _sr_record(basis):
    """Sr II 4077 as the floor emits it — the species RYA-679 nearly demoted, and the
    one RYA-691 found reading a `reliable` flag that was sitting right there unused."""
    return dict(element='Sr', ion='II', species='Sr II', reported=2.759,
                engineB_source='RYA-551 Sr II 4077.709 synth',
                engineB_reliability=basis)


def test_a_demoted_measurement_cannot_be_emitted():
    with pytest.raises(RatifiedConstraintViolation, match=r'reliable=False'):
        assert_ratified_constraints_satisfied(
            [dict(_sr_record(None), engineB_reliability=False)],
            'test emitter', RowKind.SPECIES_RECORD)


def test_a_gated_measurement_passes():
    assert_ratified_constraints_satisfied(
        [_sr_record(reliability_basis(True))], 'test emitter', RowKind.SPECIES_RECORD)


def test_an_ungated_measurement_passes_only_with_its_reason():
    """'No flag' is legitimate — the CNO cross-arm artifact genuinely has none. It is
    legitimate only when the record says which artifact and why (RYA-691 §3A)."""
    ok = reliability_basis(None, absent_reason='RYA-491/237 cross-arm artifact is not a '
                                               'profile fit: no dEW_dA / railed')
    assert_ratified_constraints_satisfied(
        [_sr_record(ok)], 'test emitter', RowKind.SPECIES_RECORD)
    with pytest.raises(RatifiedConstraintViolation, match=r'unreadable reliability basis'):
        assert_ratified_constraints_satisfied(
            [_sr_record('UNGATED')], 'test emitter', RowKind.SPECIES_RECORD)


def test_a_basis_outside_the_vocabulary_is_refused_not_assumed_good():
    """Assuming good is the entire RYA-691 defect."""
    with pytest.raises(RatifiedConstraintViolation, match=r'unreadable reliability basis'):
        assert_ratified_constraints_satisfied(
            [_sr_record('looked fine to me')], 'test emitter', RowKind.SPECIES_RECORD)


def test_a_row_declaring_nothing_is_out_of_scope():
    """Verdict and gold rows carry no reliability key. Failing them would be a schema
    change dressed as a constraint — and RYA-691 §3A forbids fabricating a uniform
    check over artifacts with genuinely different semantics."""
    assert_ratified_constraints_satisfied(
        [dict(element='Mn', ion='I', A_X=5.466)], 'test emitter')


def test_a_demoted_measurement_may_still_be_kept_as_a_diagnostic():
    assert_ratified_constraints_satisfied(
        [dict(_sr_record(None), engineB_reliability=False, diagnostic_only=True)],
        'test emitter', RowKind.SPECIES_RECORD)


def test_the_producer_cannot_write_a_basis_the_gate_would_refuse():
    """The round trip is the point: builder and classifier are one module, so the
    string the loader writes is by construction the string the gate can read."""
    assert classify_reliability_basis(reliability_basis(True)) == ReliabilityState.GATED
    assert classify_reliability_basis(
        reliability_basis(None, absent_reason='why')) == ReliabilityState.UNGATED
    with pytest.raises(ValueError, match=r'requires an absent_reason'):
        reliability_basis(None)
    with pytest.raises(ValueError, match=r'must not be emitted at all'):
        reliability_basis(False)


def test_the_loader_writes_its_basis_through_the_shared_vocabulary():
    """RYA-701's actual change: the two basis strings stopped being private to one
    script. If they drift back into f-strings there, the gate silently stops matching."""
    src = (ROOT / 'scripts' / 'rya527_two_engine_run.py').read_text(encoding='utf-8')
    assert 'from pipeline.reliability_contract import reliability_basis' in src
    assert f'"{GATED_PREFIX}' not in src and f"'{GATED_PREFIX}" not in src
