"""RYA-1117 — the RYA-1110 counts, re-validated under a source-derived tolerance.

RYA-1109 was bitten by `line_match`'s 0.005 Å default on a table printed in NANOMETRES.
This suite pins the answer to "did RYA-1110 have the same problem?" — it did not — and,
more usefully, pins the METHOD that settles it, because "I derived the tolerance" is a
claim and a plateau is a measurement.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import rya1117_tolerance_recheck as R          # noqa: E402
from pipeline import line_match                # noqa: E402


@pytest.fixture(scope="module")
def doc():
    return R.run()


# ── the premise ──────────────────────────────────────────────────────────────
def test_rya1110_never_used_the_fixed_default():
    """🔴 THE TICKET'S PREMISE, CHECKED FIRST. RYA-1117 was raised on the belief that
    RYA-1110 computed its counts with the 0.005 Å default. It did not — 0.015 Å, derived
    there from Jofré's 2-dp print precision. Pinned so the refutation cannot rot: if
    RYA-1110's tolerance is ever quietly reset to the default, this fails."""
    assert R.DERIVED_TOL_A == 0.015
    assert R.BEFORE_TOL_A == line_match.MATCH_TOL_A == 0.005
    assert R.DERIVED_TOL_A == pytest.approx(3 * R.BEFORE_TOL_A)


def test_the_precision_is_measured_not_read_off_a_format_string(doc):
    """Jofré prints ÅNGSTRÖMS to 2 dp where AGSS21 printed NANOMETRES to 2 dp — a 10x
    difference in resolution and the whole reason the two tickets diverge."""
    p = doc["measured_precision"]
    for k in ("jofre_table4_fe1", "jofre_table5_fe2", "vizier_ew_dat"):
        assert p[k]["decimals"] == {2: p[k]["n"]}, k
        assert p[k]["print_half_width_A"] == pytest.approx(0.005)
    # table6 is the coarse one, and it is NOT joined by tolerance at all
    assert doc["measured_precision"]["vizier_table6_dat"]["decimals"] == {1: 242}
    assert doc["measured_precision"]["vizier_table6_dat"]["print_half_width_A"] == pytest.approx(0.05)
    # the other side of every tolerance join is far finer, so Jofré dominates the budget
    assert p["heiter2021_geslines"]["print_half_width_A"] <= 1e-4
    assert p["canonical_gf"]["print_half_width_A"] <= 1e-4


# ── the plateau, and why the naive rule is wrong here ─────────────────────────
def test_every_published_count_is_its_clean_region_plateau(doc):
    assert R.check(doc) == []
    for key, published in R.PUBLISHED.items():
        pl = doc["plateaus"][key]
        assert pl["value"] == published, key
        assert pl["derived_tol_inside"], key
        assert not pl["still_climbing_at_clean_edge"], key


def test_the_naive_longest_run_would_have_picked_a_SATURATED_region(doc):
    """🔴 THE METHOD FINDING. "The longest run of equal values" is not sufficient on its
    own: past the point where the key stops discriminating, the count goes flat because the
    matcher is REFUSING ambiguous lines while the null admits coincidences. That flatness
    is saturation, not agreement — and here it would have given two counts the wrong value.
    """
    for key in ("canonical_gf_resolved", "gf_comparable"):
        nv = doc["plateaus"][key]["naive_longest_run"]
        assert not nv["region_is_clean"], key
        assert nv["value"] != R.PUBLISHED[key], (
            f"{key}: the naive rule happens to agree here, so this test no longer "
            f"demonstrates anything — re-derive it before trusting it")
    assert doc["plateaus"]["canonical_gf_resolved"]["naive_longest_run"]["value"] == 157
    assert doc["plateaus"]["gf_comparable"]["naive_longest_run"]["value"] == 119


def test_the_clean_region_is_bounded_by_ambiguity_and_by_the_null(doc):
    """The region has to END, or "clean" means nothing."""
    clean = [s for s in doc["sweep"] if R._is_clean(s)]
    dirty = [s for s in doc["sweep"] if not R._is_clean(s)]
    assert clean and dirty, "the sweep must span both sides of the clean edge"
    assert max(s["tol_A"] for s in clean) < min(s["tol_A"] for s in dirty)
    # ambiguity is what ends it first; the null leaks later
    first_dirty = min(dirty, key=lambda s: s["tol_A"])
    assert first_dirty["canonical_gf_ambiguous"] > 0
    assert all(n == 0 for s in clean for n in s["displaced_null"])


def test_the_widest_window_in_the_sweep_is_visibly_broken(doc):
    """A sweep that never breaks has not bracketed the answer."""
    worst = max(doc["sweep"], key=lambda s: s["tol_A"])
    assert worst["tol_A"] == 0.5
    assert max(worst["displaced_null"]) > 100, "the null must blow up at a silly window"


# ── the counts themselves ────────────────────────────────────────────────────
def test_the_before_after_table_shows_what_the_default_would_have_cost(doc):
    b, a = doc["before"], doc["after"]
    assert (b["canonical_gf_resolved"], a["canonical_gf_resolved"]) == (154, 159)
    assert (b["gf_comparable"], a["gf_comparable"]) == (117, 121)
    assert (b["gf_identical"], a["gf_identical"]) == (65, 67)
    assert (b["heiter2021_resolved"], a["heiter2021_resolved"]) == (152, 159)


def test_the_graded_overlap_counts_do_not_depend_on_the_tolerance_at_all(doc):
    """16 Fe I / 0 Fe II is the same at every clean window — the one result in RYA-1110
    that could not have been a tolerance artifact under any choice."""
    clean = [s for s in doc["sweep"] if R._is_clean(s)]
    assert {s["graded_LAB_overlap_FeI"] for s in clean} == {16}
    assert {s["graded_LAB_overlap_FeII"] for s in clean} == {0}
    # it survives even the broken windows, right out to 0.25 A
    upto = [s for s in doc["sweep"] if s["tol_A"] <= 0.25]
    assert {s["graded_LAB_overlap_FeI"] for s in upto} == {16}


def test_the_fe2_structural_zero_is_a_wavelength_fact(doc):
    """Not a tolerance artifact and not close to being one."""
    z = doc["structural_zero_FeII"]
    assert z["our_graded_FeII_VIS"]["hi_A"] < z["gbs_FeII_selected"]["lo_A"]
    assert z["closest_cross_pair_A"] == pytest.approx(409.53, abs=0.01)
    assert z["window_multiple_needed"] > 10_000
    assert z["gbs_FeII_selected"]["n"] == 9 and z["our_graded_FeII_VIS"]["n"] == 10


def test_the_recheck_changed_no_value_and_no_line():
    """RYA-1117 re-validates counts. It must not touch the artifact it is checking."""
    import subprocess
    r = subprocess.run(["git", "status", "--porcelain",
                        "data/linelists/reference_sets/gbs_solar_fe_rya1110.csv",
                        "data/linelists/canonical_gf.csv",
                        "data/reference/jofre2014_gbs",
                        "data/reference/heiter2021_ges"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.stdout.strip() == "", f"RYA-1117 modified data it was only meant to read:\n{r.stdout}"
