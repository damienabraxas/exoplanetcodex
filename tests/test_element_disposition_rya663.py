"""RYA-663 — per-element disposition report (+ RYA-675 staleness-detector narrow).

The engine tests are fixture-based. The live tests pin the two facts that make this
report worth trusting: that it applies the RATIFIED gate rather than a private copy,
and that it distinguishes "no two-engine record exists" from "the delta failed".

RYA-675 adds a third: that `artifact_age_stale` and `cross_channel_disagreement` are
DISTINCT signals with distinct scopes and distinct remedies, and can never be collapsed
back into one "stale" verdict. Those tests live here rather than in a second
`tests/test_element_disposition.py` (the name RYA-675 §4 guessed at) because a module
with two test files grows two contradictory ideas of what it does — which is the exact
failure mode the narrow exists to end.
"""

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline import element_disposition as ed  # noqa: E402
from pipeline.engine_selection import FLOOR_PROMOTION  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def _row(element, verdict="CURATION-OWED", n_lines=1, channel="EW: 1 curated line(s)",
         a_measured=None, ref=6.0):
    return {"element": element, "verdict": verdict, "n_lines": n_lines,
            "channel": channel, "A_measured": a_measured, "asplund2021": ref}


# ── bucketing ────────────────────────────────────────────────────────────────
def test_bucket_pass_wins_over_everything():
    assert ed.classify_owed_reason(_row("Fe", verdict="PASS")) == ed.OwedReason.PASS


def test_bucket_zero_survivors_is_owed_blank():
    assert ed.classify_owed_reason(_row("Y", n_lines=0)) == ed.OwedReason.OWED_BLANK


def test_bucket_held_is_detected_from_the_channel_not_a_name_list():
    """A hardcoded element list would go stale the moment a measurement lands; the
    bucket has to come from the row's own channel string."""
    held = _row("Ca", channel="EW: 2 curated line(s); value HELD at gold tier 'owed' (RYA-522)")
    assert ed.classify_owed_reason(held) == ed.OwedReason.OWED_HELD


def test_bucket_synthesis_channel_is_measured_awaiting_freeze():
    syn = _row("Ba", channel="synthesis: Ba II 5853.668 HFS-resolved LTE COG")
    assert ed.classify_owed_reason(syn) == ed.OwedReason.MEASURED_UNFROZEN


def test_bucket_plain_ew_pool():
    assert ed.classify_owed_reason(_row("Si", channel="EW: 7 line(s)")) == ed.OwedReason.EW_POOL


def test_upper_limit_veto_beats_the_channel():
    """Li has a real EW line and a real value; the ratified UPPER_LIMIT disposition
    (RYA-563) must still take precedence over the channel-derived bucket."""
    li = _row("Li", n_lines=1, channel="EW: Li I 6707 (single line, UPPER LIMIT)")
    assert ed.classify_owed_reason(li) == ed.OwedReason.UPPER_LIMIT


# ── the STRICT gate-3 distinction (the reason this report exists) ─────────────
def test_gate3_unevaluable_is_not_gate3_failed():
    """'No two-engine record' points at the RYA-527 re-run; 'a record with no delta'
    points at the element. Collapsing them would hide which reds the re-run clears."""
    report = ed.build_report()
    states = {d["element"]: d["gate3_state"] for d in report["dispositions"]}
    assert ed.GATE3_UNEVALUABLE in states.values()
    assert ed.GATE3_FAILED in states.values()
    assert states != {}


def test_vetoed_elements_report_gate3_not_applicable():
    """A short-circuited element never reached gate 3 — reporting UNEVALUABLE there
    would claim we looked for a two-engine record and found none."""
    report = ed.build_report()
    for d in report["dispositions"]:
        if d["owed_reason"] in (ed.OwedReason.PASS, ed.OwedReason.UPPER_LIMIT):
            assert d["gate3_state"] == ed.GATE3_NA, d["element"]


# ── the report must not invent a gate ────────────────────────────────────────
def test_thresholds_come_from_the_ratified_module():
    """If this report ever grows its own thresholds, promotion stops being ratified.
    Pin them to engine_selection so a divergence is a test failure, not a silent drift."""
    assert ed.build_report()["thresholds"] == dict(FLOOR_PROMOTION)


def test_every_element_lands_in_exactly_one_bucket():
    report = ed.build_report()
    elements = [d["element"] for d in report["dispositions"]]
    assert len(elements) == len(set(elements))
    assert len(elements) == report["phase_c_summary"]["n_elements"]
    for d in report["dispositions"]:
        assert d["owed_reason"] in ed.OWED_REASON_PLAN


def test_promotions_are_exactly_the_promoted_rows():
    report = ed.build_report()
    assert report["can_flip_now"] == [d["element"] for d in report["dispositions"]
                                      if d["promoted"]]


# ── honesty about the inputs ─────────────────────────────────────────────────
def test_a_promotion_on_a_stale_input_is_marked_provisional():
    """The one flip available today rests on a cross-engine delta from a REVIEW
    artifact that demonstrably lags the live channel. Presenting that as settled is
    exactly the overclaim this report is supposed to prevent."""
    report = ed.build_report()
    if not report["gate3_provisional"]:
        pytest.skip("inputs are no longer stale — the caveat is correctly absent")
    for d in report["dispositions"]:
        if d["promoted"] and d["cross_engine_delta"] is not None:
            assert "PROVISIONAL" in d["promotion_reason"], d["element"]
            # RYA-675: a doubt with no exit is what left Ca stuck. Every provisional
            # stamp must name the remedy that clears it.
            assert d["remedy"] in (ed.REMEDY_REGENERATE, ed.REMEDY_ADJUDICATE,
                                   ed.REMEDY_REGENERATE_THEN_ADJUDICATE), d["element"]
            assert d["remedy"] in d["promotion_reason"], d["element"]


def test_cross_channel_disagreement_is_concrete_not_a_date_heuristic():
    """Each entry must name an element and both numbers — a claim that cannot be
    checked is not evidence."""
    for e in ed.build_report()["cross_channel_disagreement"]:
        assert {"element", "two_engine", "live_phase_c", "delta"} <= set(e)
        assert abs(e["two_engine"] - e["live_phase_c"]) > 0.001


# ── RYA-675: the two signals, and their two remedies ─────────────────────────
# The defect these pin: `detect_stale_inputs` read ANY cross-channel value difference as
# proof the artifact was old, so RYA-669's fresh in-run artifact was reported stale on
# six elements and Ca's gate-3 flag could not be cleared by the re-run that existed to
# clear it. The three cases below are age-only, disagreement-only, and both.

def _prov(path, when, sha="abc1234"):
    return {"path": path, "commit": sha, "committed": when}


NEWER = "2026-08-07T17:11:53-06:00"
OLDER = "2026-07-27T00:50:37-06:00"


def test_artifact_age_only():
    """A genuinely OLD artifact, with both channels agreeing on every value.

    Only the age signal may fire, and the remedy must be REGENERATE — the one remedy a
    re-run actually executes."""
    age = ed.detect_artifact_age_stale(
        _prov("two_engine.json", OLDER),
        [_prov("phase_c.json", NEWER), _prov("gold.csv", OLDER)])
    assert len(age) == 1, age
    assert age[0]["input"] == "phase_c.json"
    assert OLDER in age[0]["note"] and NEWER in age[0]["note"]

    rows = [_row("Ca", a_measured=6.324)]
    two_engine = {"Ca": {"reported": 6.324, "mean_cross_engine_delta": 0.016}}
    assert ed.detect_cross_channel_disagreement(rows, two_engine) == []

    c = ed.classify_disposition("Ca", artifact_age_stale=True, disagreeing_elements=[])
    assert c.signals == (ed.SIGNAL_ARTIFACT_AGE,)
    assert c.remedy == ed.REMEDY_REGENERATE
    assert c.provisional is True


def test_cross_channel_disagreement_only():
    """The RYA-669 shape: an artifact generated during the very run that reads it, whose
    values still disagree with the live channel. Nothing here is old, so the age signal
    must stay silent and the remedy must be ADJUDICATE — re-running changes nothing."""
    age = ed.detect_artifact_age_stale(
        _prov("two_engine.json", NEWER),
        [_prov("phase_c.json", OLDER), _prov("gold.csv", NEWER)])
    assert age == [], "a same-run artifact predates nothing — this is the RYA-675 defect"

    rows = [_row("Fe", a_measured=7.466), _row("Ca", a_measured=6.324)]
    two_engine = {"Fe": {"reported": 7.580}, "Ca": {"reported": 6.324}}
    dis = ed.detect_cross_channel_disagreement(rows, two_engine)
    assert [d["element"] for d in dis] == ["Fe"]
    assert dis[0]["delta"] == pytest.approx(0.114)

    c = ed.classify_disposition("Fe", artifact_age_stale=False,
                                disagreeing_elements=["Fe"])
    assert c.signals == (ed.SIGNAL_CROSS_CHANNEL,)
    assert c.remedy == ed.REMEDY_ADJUDICATE
    # Ca is untouched: a disagreement about Fe is evidence about Fe and nothing else.
    # The pre-RYA-675 blanket stamp is precisely what this asserts is gone.
    assert ed.classify_disposition("Ca", artifact_age_stale=False,
                                   disagreeing_elements=["Fe"]).remedy == ed.REMEDY_NONE


def test_both_signals_active():
    """Old AND disagreeing. Both signals must be reported separately — never merged into
    one 'stale' verdict — and the remedy is ordered: regenerate first, because the fresh
    emission may resolve the disagreement without any human decision."""
    age = ed.detect_artifact_age_stale(
        _prov("two_engine.json", OLDER),
        [_prov("phase_c.json", NEWER), _prov("gold.csv", NEWER)])
    assert len(age) == 2

    rows = [_row("Si", a_measured=7.888)]
    dis = ed.detect_cross_channel_disagreement(rows, {"Si": {"reported": 7.639}})
    assert [d["element"] for d in dis] == ["Si"]

    c = ed.classify_disposition("Si", artifact_age_stale=True,
                                disagreeing_elements=["Si"])
    assert c.signals == (ed.SIGNAL_ARTIFACT_AGE, ed.SIGNAL_CROSS_CHANNEL)
    assert c.remedy == ed.REMEDY_REGENERATE_THEN_ADJUDICATE
    assert c.provisional is True


def test_a_row_that_read_nothing_is_not_made_provisional_by_the_artifacts_age():
    """Age is a property of the FILE, so it can only bear on a row that read the file.
    A PASS row whose gates were short-circuited never consulted the artifact."""
    c = ed.classify_disposition("Mn", artifact_age_stale=True, disagreeing_elements=[],
                                reads_artifact=False)
    assert c.signals == () and c.remedy == ed.REMEDY_NONE and c.provisional is False


def test_a_disagreement_is_reported_even_where_the_row_read_nothing():
    """Fe's shape: PASS on the ratified channel, gates short-circuited, yet the two
    channels still differ. The reconciliation is owed to the freeze — so ADJUDICATE is
    reported, while the row itself is NOT provisional."""
    c = ed.classify_disposition("Fe", artifact_age_stale=True,
                                disagreeing_elements=["Fe"], reads_artifact=False)
    assert c.signals == (ed.SIGNAL_CROSS_CHANNEL,)
    assert c.remedy == ed.REMEDY_ADJUDICATE
    assert c.provisional is False


def test_the_two_signals_are_never_collapsed_into_one_key():
    """The whole point of RYA-675. A single combined 'stale' key would let the next
    consumer re-derive the conflation this ticket removed."""
    report = ed.build_report()
    assert "stale_input_evidence" not in report
    assert isinstance(report["artifact_age_stale"], bool)
    assert isinstance(report["cross_channel_disagreement"], list)
    # gate3_provisional is the GLOBAL taint, so it may only track the GLOBAL signal.
    assert report["gate3_provisional"] == report["artifact_age_stale"]
    assert (report["gate3_provisional_signal"] == ed.SIGNAL_ARTIFACT_AGE
            if report["artifact_age_stale"] else
            report["gate3_provisional_signal"] is None)


def test_every_remedy_names_a_plan_and_every_provisional_row_names_a_remedy():
    report = ed.build_report()
    for d in report["dispositions"]:
        assert d["remedy"] in ed.REMEDY_PLAN, d["element"]
        assert d["remedy_plan"] == ed.REMEDY_PLAN[d["remedy"]]
        assert bool(d["staleness_signals"]) == (d["remedy"] != ed.REMEDY_NONE), d["element"]
        if d["provisional"]:
            assert d["remedy"] != ed.REMEDY_NONE, d["element"]


def test_age_detection_uses_commit_times_not_mtimes():
    """A checkout or a `touch` rewrites mtimes and would make a stale artifact look
    fresh (RYA-659). The detector must read only the git provenance handed to it."""
    art = _prov("two_engine.json", OLDER)
    assert ed.detect_artifact_age_stale(art, [_prov("phase_c.json", OLDER)]) == []
    assert ed.detect_artifact_age_stale(art, [_prov("phase_c.json", NEWER)])


def test_ion_mismatch_is_reported_not_gated():
    """gold v1 holds Sr I while the registry ratifies Sr II. Gating across species
    would report a ~2 dex miss as physics; it is a species mismatch and must say so."""
    report = ed.build_report()
    sr = next(d for d in report["dispositions"] if d["element"] == "Sr")
    assert "ION MISMATCH" in sr["value_source"]
    assert sr["value"] is None


def test_reported_ion_prefers_the_registry_over_gold():
    """Gold's ion column reads Sr I / Ba I, the superseded raw-EW leg. The ratified
    registry lock must win, or gate 3 would read the wrong species' delta."""
    gold_ions = {"Sr": "I", "Ba": "I", "Fe": "I"}
    assert ed.reported_ion("Sr", gold_ions) == "II"
    assert ed.reported_ion("Ba", gold_ions) == "II"
    assert ed.reported_ion("Fe", gold_ions) == "I"      # not locked -> gold answers


def test_value_disagreements_surface_what_the_ledger_guard_cannot_see():
    """The RYA-632 guard compares verdicts and counts, never values: gold's frozen Fe
    against the live Fe is invisible to it. This list is that blind spot."""
    rows = {v["element"]: v for v in ed.build_report()["value_disagreements"]}
    assert "Fe" in rows, "the known Fe 7.516-vs-7.466 split must be surfaced"
    assert rows["Fe"]["spread"] > 0.1


# ── the generator contract ───────────────────────────────────────────────────
def test_committed_report_matches_a_fresh_run():
    """Same contract as the RYA-654 tracker generator: a report that has drifted from
    the ledger must not be able to sit in the repo."""
    r = subprocess.run([sys.executable, "scripts/gen_element_disposition.py", "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_markdown_renders_from_the_same_report_object():
    md = ed.render_markdown(ed.build_report())
    assert "do not hand-edit" in md
    assert "| El | verdict |" in md
