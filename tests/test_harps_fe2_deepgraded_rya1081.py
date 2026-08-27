"""
tests/test_harps_fe2_deepgraded_rya1081.py — RYA-1081
=====================================================
HARPS Fe II VIS DEEPGRADED reads 7.966 against KP's ~7.57. The suspicion was a few skewing
lines. It is not: the offset is pool-wide and arm-specific, on gf that are identical between
the arms. These tests pin the findings that would be easy to lose — above all that the one
excluded line is NOT what closes the gap, so nobody later reads the exclusion as the repair.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

VERDICT = ROOT / "data" / "results" / "rya1081" / "rya1081_verdict.json"
TABLE = ROOT / "data" / "results" / "rya1081" / "rya1081_perline_audit.csv"
PROBLEM = ROOT / "data" / "registry" / "problem_children.csv"


@pytest.fixture(scope="module")
def v():
    if not VERDICT.exists():
        pytest.skip("RYA-1081 verdict absent")
    return json.loads(VERDICT.read_text())


def test_the_offset_is_arm_specific_on_identical_gf(v):
    """THE DIAGNOSTIC. Same nine lines, same PRIMARY LAB gf, and HARPS sits above KP on
    eight of nine. A gf error moves both arms together; it cannot open a gap between them."""
    assert v["gf"]["all_primary_lab"] is True
    assert v["gf"]["so_gf_cannot_explain_the_arm_gap"] is True
    for arm in ("kp_molecfit", "kp_kurucz2005"):
        c = v["kp_crosscheck"][arm]
        assert c["median_delta"] > 0.2
        assert c["n_harps_above_kp"] == 8 and c["n"] == 9


def test_the_published_median_is_not_driven_by_the_high_puller(v):
    """🔴 The ticket's premise, refuted. Removing the worst line in the pool moves the
    published median by 0.035 dex — a median is robust to exactly the suspected failure."""
    r = v["robustness"]
    assert abs(r["shift"]) < 0.05, (
        "the outlier now moves the median materially — the pool-wide reading may no longer "
        "hold, so re-derive before repeating it")
    assert r["spread_excluding_outlier_dex"] > 1.0


def test_the_excluded_line_is_not_presented_as_the_fix(v):
    """🔴 THE MISREADING THIS EXISTS TO PREVENT. 4303.170 is dispositioned on a measured fit
    failure present on BOTH arms. It is not saturation, not a diagnosis, and — the part that
    matters — not what closes the HARPS-KP gap."""
    row = [ln for ln in PROBLEM.read_text().splitlines() if "4303.170" in ln]
    assert len(row) == 1, "the 4303.170 disposition is missing or duplicated"
    note = row[0]
    assert "BAD_FIT" in note and "exclude" in note
    assert "UNDIAGNOSED" in note, "the mechanism must be recorded as unknown, not implied"
    assert "-0.035" in note or "0.035" in note, (
        "the note must carry how little the exclusion moves the value")
    # and the arm gap must survive the exclusion, which is the substantive claim
    c = v["kp_crosscheck"]["kp_molecfit"]
    assert c["median_delta_excluding_outlier"] > 0.2


def test_the_saturation_question_is_unanswered_not_answered(v):
    """RYA-833. The ratified COG covers Fe I only, so 'not flagged' is not 'not saturated',
    and EW magnitude must not be substituted for it (a CRITICAL failure — and moot, since
    this route measures no EW at all)."""
    s = v["saturation_cog"]
    assert s["runnable"] is False
    assert s["n_covered"] == 0
    assert "UNANSWERED" in s["reason"]
    t = pd.read_csv(TABLE)
    assert t.ew_mA.isna().all() and t.rew.isna().all(), (
        "an EW appeared in this product — the saturation reasoning must be revisited")


def test_the_graded_reference_does_not_exist_and_is_not_faked(v):
    """The deep contribution cannot be isolated because there is no GRADED Fe II tier
    anywhere in the feed. Manufacturing one from a depth proxy is the CRITICAL conflation."""
    g = v["graded_reference_gap"]
    assert g["runnable"] is False
    assert g["fe2_tiers_present_in_feed"] == ["DEEPGRADED"]


def test_it_is_not_the_rya911_residual(v):
    """RYA-911 was the EW leg running LOW. This is the synthesis leg running HIGH, with no
    EW leg present at all — routing it to 911 would attach it to the wrong mechanism."""
    r = v["rya911"]
    assert r["is_the_rya911_residual"] is False
    assert r["any_measured_ew_in_this_product"] is False
    assert v["kp_crosscheck"]["kp_molecfit"]["median_delta"] > 0


def test_no_value_was_tuned_toward_kp(v):
    """🔴 CRITICAL (RYA-161 / RYA-523): find the bad line, never tune. The published value
    must still be the pool median of the unmodified per-line abundances."""
    t = pd.read_csv(TABLE)
    assert abs(float(t.abundance.median()) - v["published"]["value"]) < 1e-9
    assert v["published"]["value"] == pytest.approx(7.966, abs=1e-3)
    # exactly one line is dispositioned; the other eight stand as measured
    assert v["published"]["n_lines"] == 9


def test_the_fit_quality_correlation_was_not_overread(v):
    """A correction to this audit's own first pass, pinned so it cannot come back: the
    A-vs-red_chi2 correlation of +0.86 is driven entirely by the outlier. On the other
    eight it is ~+0.31, which does NOT support 'fit quality explains the offset'."""
    c = v["kp_crosscheck"]["kp_molecfit"]
    assert c["corr_A_vs_red_chi2_harps"] > 0.7
    assert c["corr_A_vs_red_chi2_harps_excluding_outlier"] < 0.5
