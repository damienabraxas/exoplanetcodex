"""RYA-663 — per-element disposition report.

The engine tests are fixture-based. The live tests pin the two facts that make this
report worth trusting: that it applies the RATIFIED gate rather than a private copy,
and that it distinguishes "no two-engine record exists" from "the delta failed".
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


def test_stale_input_evidence_is_concrete_not_a_date_heuristic():
    """Each entry must name an element and two numbers — a staleness claim that cannot
    be checked is not evidence."""
    for line in ed.build_report()["stale_input_evidence"]:
        assert "two-engine reports" in line and "live verdict channel measures" in line


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
