"""RYA-632 — results-ledger consistency guard.

Engine tests are fixture-based and need no live files. The final test pins the LIVE
ledger state and is a loud strict-xfail: the guard is PRE-DECLARED red on the
un-ratified Co reconciliation (RYA-564 left it to "a deliberate pass, not a drive-by")
plus the five same-class rows the RYA-632 audit surfaced. When those are ratified the
strict marker fires and forces this test to be un-marked deliberately.
"""

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline import provenance_honesty  # noqa: E402
from pipeline.ledger_consistency_guard import (  # noqa: E402
    ArtifactState,
    check_counts,
    check_element,
    collect_element_states,
    load_allowed_divergences,
    run_check,
)

ROOT = Path(__file__).resolve().parents[1]
# Taken from the module that OWNS the claim (RYA-596/653), not re-typed: if the
# emitters reword it, these fixtures follow and the coupling test below fails loudly
# rather than this suite quietly testing a string nothing writes any more.
ZERO_SURVIVOR = provenance_honesty.ZERO_SURVIVOR_CHANNEL


def _allowed(*entries):
    """entries: (element, artifactA, artifactB)."""
    out = {}
    for el, a, b in entries:
        out.setdefault(el, set()).add(frozenset((a, b)))
    return out


# ── C1 — the RYA-596 phantom: a zero-survivor claim next to graded survivors ─────────

def test_c1_zero_survivor_claim_beside_survivors_fails():
    """The exact shape RYA-596 found in the live verdict: phase_c asserts 'no line
    survives' on a row the tracker shows carrying 10 lines."""
    states = [
        ArtifactState("phase_c", "Ti", verdict="CURATION_OWED", raw_cause=ZERO_SURVIVOR),
        ArtifactState("tracker", "Ti", verdict="CURATION_OWED", n_lines=10),
    ]
    errs = check_element(states, allowed_divergences={})
    assert len(errs) == 1
    assert "claim no data while" in errs[0]
    assert "Ti" in errs[0]


def test_c1_passes_when_the_pair_is_annotated_and_ratified():
    states = [
        ArtifactState("phase_c", "Ti", verdict="CURATION_OWED", raw_cause=ZERO_SURVIVOR),
        ArtifactState("tracker", "Ti", verdict="CURATION_OWED", n_lines=10),
    ]
    allowed = _allowed(("Ti", "phase_c", "tracker"))
    assert check_element(states, allowed) == []


def test_c1_annotation_for_a_different_element_does_not_transfer():
    states = [
        ArtifactState("phase_c", "Ti", verdict="CURATION_OWED", raw_cause=ZERO_SURVIVOR),
        ArtifactState("tracker", "Ti", verdict="CURATION_OWED", n_lines=10),
    ]
    assert check_element(states, _allowed(("Ni", "phase_c", "tracker"))) != []


def test_c1_self_contradiction_inside_one_artifact_can_never_be_annotated():
    """One artifact claiming 'no line survives' on a row it also gives survivors for is
    not a cross-artifact divergence — no exceptions entry can license it."""
    states = [ArtifactState("gold", "Ba", n_lines=1, raw_cause=ZERO_SURVIVOR)]
    assert check_element(states, _allowed(("Ba", "gold", "gold"))) != []


def test_c1_clean_when_every_artifact_agrees_there_is_no_data():
    """Mg/Y/Zr/Eu are CURATION-OWED *and* verified zero-survivor (RYA-596) — a real,
    correct state that must not read as a contradiction."""
    states = [
        ArtifactState("phase_c", "Mg", verdict="CURATION_OWED", n_lines=0,
                      raw_cause=ZERO_SURVIVOR),
        ArtifactState("gold", "Mg", verdict="CURATION_OWED", n_lines=0, tier="owed",
                      raw_cause=ZERO_SURVIVOR),
        ArtifactState("tracker", "Mg", verdict="CURATION_OWED", tier="owed"),
    ]
    assert check_element(states, {}) == []


def test_c1_recognises_the_claim_in_the_emitters_own_words():
    """The guard must read the zero-survivor claim from pipeline.provenance_honesty
    (RYA-653), which is what the phase_c verdict and the gold builder actually assert
    with. If it kept a private copy of the wording, a reword in those emitters would
    silently blind C1 — the guard would stop seeing the very phantom it exists for.

    Both live spellings are covered: the gold builder writes the bare claim, the
    verdict channel writes it behind an "EW present; " prefix.
    """
    for cause in (provenance_honesty.ZERO_SURVIVOR_CLAIM,
                  provenance_honesty.ZERO_SURVIVOR_CHANNEL):
        states = [
            ArtifactState("gold", "Ba", verdict="CURATION_OWED", raw_cause=cause),
            ArtifactState("phase_c", "Ba", verdict="CURATION_OWED", n_lines=1),
        ]
        assert check_element(states, {}) != [], f"C1 went blind to {cause!r}"


# ── C2 — an owed tier freezes no value (RYA-522, Ryan 2026-07-05) ────────────────────

def test_c2_owed_tier_with_a_frozen_value_fails():
    states = [ArtifactState("gold", "Ti", verdict="NLTE_OWED", tier="nlte_owed",
                            frozen_value=5.471, n_lines=10)]
    errs = check_element(states, {})
    assert any("is owed but froze a value (5.471)" in e for e in errs)


def test_c2_gold_tier_may_freeze_a_value():
    states = [ArtifactState("gold", "K", verdict="PASS", tier="gold",
                            frozen_value=5.099, n_lines=1)]
    assert check_element(states, {}) == []


# ── C3 — a stale tracker's counts do not match the tally ────────────────────────────

def test_c3_stale_tracker_counts_fail():
    """Co measured (PASS) but the tracker never bumped — the live RYA-632 case."""
    errs = check_counts({"PASS": 5, "CURATION_OWED": 20},
                        {"PASS": 6, "CURATION_OWED": 20})
    assert errs == ["count mismatch [PASS]: tracker=5 tallied=6"]


def test_c3_matching_counts_are_clean():
    assert check_counts({"PASS": 6, "CURATION_OWED": 20},
                        {"PASS": 6, "CURATION_OWED": 20}) == []


def test_c3_a_bucket_missing_from_one_side_still_fails():
    assert check_counts({"NLTE_OWED": 1}, {}) == ["count mismatch [NLTE_OWED]: tracker=1 tallied=0"]


# ── C4 — cross-artifact verdict disagreement ────────────────────────────────────────

def test_c4_unannotated_verdict_disagreement_fails():
    states = [
        ArtifactState("physics_regime", "Co", verdict="GET_DATA"),
        ArtifactState("phase_c", "Co", verdict="PASS", n_lines=5),
    ]
    errs = check_element(states, {})
    assert any("verdict disagreement across artifacts" in e for e in errs)


def test_c4_annotated_and_ratified_disagreement_passes():
    states = [
        ArtifactState("physics_regime", "Sr", verdict="GET_DATA"),
        ArtifactState("phase_c", "Sr", verdict="PASS", n_lines=5),
    ]
    assert check_element(states, _allowed(("Sr", "physics_regime", "phase_c"))) == []


def test_c4_agreeing_pairs_need_no_annotation():
    """Only the pairs that actually disagree need a ratified entry — annotating an
    agreement would be meaningless paperwork."""
    states = [
        ArtifactState("physics_regime", "Sc", verdict="GET_DATA"),
        ArtifactState("phase_c", "Sc", verdict="CURATION_OWED", n_lines=1),
        ArtifactState("gold", "Sc", verdict="CURATION_OWED", n_lines=1),
    ]
    allowed = _allowed(("Sc", "physics_regime", "phase_c"),
                       ("Sc", "physics_regime", "gold"))
    assert check_element(states, allowed) == []


def test_c4_partial_annotation_still_fails():
    states = [
        ArtifactState("physics_regime", "Sc", verdict="GET_DATA"),
        ArtifactState("phase_c", "Sc", verdict="CURATION_OWED", n_lines=1),
        ArtifactState("gold", "Sc", verdict="CURATION_OWED", n_lines=1),
    ]
    allowed = _allowed(("Sc", "physics_regime", "phase_c"))
    assert check_element(states, allowed) != []


def test_no_states_is_clean():
    assert check_element([], {}) == []


# ── exceptions file loading ─────────────────────────────────────────────────────────

def test_load_allowed_divergences_rejects_an_entry_without_a_ratifying_ticket(tmp_path):
    p = tmp_path / "known_verdict_divergences.yaml"
    p.write_text(textwrap.dedent("""
        - element: Co
          between: [physics_regime, phase_c]
          reason: "looks fine to me"
    """), encoding="utf-8")
    with pytest.raises(ValueError, match="no ratified_by"):
        load_allowed_divergences(p)


def test_load_allowed_divergences_rejects_a_malformed_pair(tmp_path):
    p = tmp_path / "known_verdict_divergences.yaml"
    p.write_text(textwrap.dedent("""
        - element: Co
          between: [phase_c]
          reason: "half a pair"
          ratified_by: RYA-000
    """), encoding="utf-8")
    with pytest.raises(ValueError, match="two distinct artifacts"):
        load_allowed_divergences(p)


def test_load_allowed_divergences_reads_a_valid_entry(tmp_path):
    p = tmp_path / "known_verdict_divergences.yaml"
    p.write_text(textwrap.dedent("""
        - element: Sr
          between: [physics_regime, phase_c]
          reason: "LOCKED verifier is EW-pool-based; Sr is synthesis-sourced"
          ratified_by: RYA-551
    """), encoding="utf-8")
    assert load_allowed_divergences(p) == {"Sr": {frozenset(("physics_regime", "phase_c"))}}


def test_load_allowed_divergences_on_a_missing_file_is_empty(tmp_path):
    assert load_allowed_divergences(tmp_path / "nope.yaml") == {}


def test_the_committed_exceptions_file_seeds_only_the_ratified_axis():
    """RYA-632 seeded Sr alone; RYA-654 ratified Co/N/K/Cu after re-cutting the GET-DATA
    class on the EW-measurable-vs-synthesis-sourced axis. The two ABSENCES are the
    load-bearing part of this test:

      * Sc — measurement-not-trusted (its only value is the KP Sc II 4246 blue-edge HFS
        single line, held LOW_CONFIDENCE by phase_c itself). Nothing to diverge from, so
        annotating it would launder an un-done element. It stays red.
      * P  — EW-measurable and merely stale, so RYA-654 FIXED the physics_regime row
        instead. A fix must never leave an annotation behind.
    """
    allowed = load_allowed_divergences()
    assert set(allowed) == {"Sr", "Co", "N", "K", "Cu"}
    assert "Sc" not in allowed, "annotating Sc would launder an un-done element"
    assert "P" not in allowed, "P was reconciled by fixing the yaml, not by annotating"
    for element in ("Co", "N", "K", "Cu"):
        # every pair that actually disagrees needs its own entry (RYA-632 rule)
        assert frozenset(("physics_regime", "phase_c")) in allowed[element]
        assert frozenset(("physics_regime", "tracker")) in allowed[element]


def test_documented_owed_annotates_but_never_suppresses():
    """The RYA-654 message mechanism must not become a second exceptions file.

    _OWED_NOT_LAUNDERED explains a red; it may never remove one. If this ever inverts,
    the laundering the guard exists to stop has been re-introduced through the door
    marked 'documentation'.
    """
    from pipeline.ledger_consistency_guard import (_OWED_NOT_LAUNDERED,
                                                   _annotate_documented)
    assert "Sc" in _OWED_NOT_LAUNDERED
    errs = check_element(
        [ArtifactState("physics_regime", "Sc", verdict="GET_DATA"),
         ArtifactState("phase_c", "Sc", verdict="CURATION_OWED", n_lines=1)],
        load_allowed_divergences())
    assert errs, "Sc must still FAIL — documenting a red does not clear it"
    annotated = _annotate_documented(errs[0])
    assert "DOCUMENTED-OWED, still failing" in annotated
    assert errs[0] in annotated          # the original failure text is never replaced


# ── live adapter ────────────────────────────────────────────────────────────────────

def test_adapter_reads_all_four_live_artifacts():
    states_by_element, tracker_counts, tallied_counts = collect_element_states()
    assert len(states_by_element) == 26, "TARGET_ELEMENTS grain (Fe II arbiter row dropped)"
    for el, states in states_by_element.items():
        arts = {s.artifact for s in states}
        assert arts == {"tracker", "phase_c", "gold", "physics_regime"}, f"{el}: {arts}"
    assert sum(tracker_counts.values()) == 26
    assert sum(tallied_counts.values()) == 26


def test_live_gold_never_freezes_a_value_at_an_owed_tier():
    """C2 against the live ledger — the RYA-522 invariant, independently of the
    un-ratified divergences below."""
    states_by_element, _, _ = collect_element_states()
    offenders = [e for e, states in states_by_element.items()
                 if any("is owed but froze a value" in m
                        for m in check_element(states, load_allowed_divergences()))]
    assert offenders == []


@pytest.mark.xfail(strict=True, reason=(
    "PRE-DECLARED RED, now down to ONE documented-owed element (was six un-ratified "
    "divergences plus two count mismatches at RYA-632). RYA-654 cleared Co/N/K/Cu by "
    "ratification, P by fixing the stale yaml row, and both count mismatches by "
    "generating the tracker; RYA-665's gold v3 freeze cleared Ba, whose red was the "
    "phantom blank cause frozen into gold v2. What remains is red ON PURPOSE: Sc "
    "(measurement-not-trusted — cleared by measuring Sc, never by an exceptions entry). "
    "When Sc clears, this strict marker fires and must be removed deliberately."))
def test_live_ledger_is_consistent():
    assert run_check() == 0


def test_the_live_reds_are_exactly_the_documented_ones():
    """The stronger companion to the xfail above: not merely 'red', but red for EXACTLY
    the adjudicated reason. A new contradiction cannot hide behind the known one.

    The triple equality is the load-bearing part: the live offender set must equal the
    documented-owed table, so an element cannot be red without an explanation NOR carry a
    stale explanation without being red. RYA-665's v3 freeze cleared Ba on both sides at
    once — it left the offender set and its _OWED_NOT_LAUNDERED entry was removed."""
    from pipeline.ledger_consistency_guard import _OWED_NOT_LAUNDERED
    states_by_element, tracker_counts, tallied_counts = collect_element_states()
    allowed = load_allowed_divergences()
    offenders = {e for e in states_by_element
                 if check_element(states_by_element[e], allowed)}
    assert offenders == set(_OWED_NOT_LAUNDERED) == {"Sc"}
    # the RYA-632 count mismatches (PASS 5 vs 6, NLTE_OWED 1 vs 0) are gone for good:
    # the tracker is generated from the same phase_c the tally is taken over.
    assert check_counts(tracker_counts, tallied_counts) == []


def test_guard_is_runnable_as_a_module():
    """`python -m pipeline.ledger_consistency_guard` is the merge-gate entry point; it
    must exit non-zero and name the offending elements, not traceback."""
    r = subprocess.run([sys.executable, "-m", "pipeline.ledger_consistency_guard"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 1
    assert "RESULTS-LEDGER CONSISTENCY GUARD FAILED" in r.stderr
    assert "Traceback" not in r.stderr
    assert "Sc:" in r.stderr
    # a documented red still prints its explanation AND still counts as a failure
    assert "DOCUMENTED-OWED, still failing" in r.stderr
    assert "0 undocumented" in r.stderr
