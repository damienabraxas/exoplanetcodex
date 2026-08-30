"""RYA-1104 -- the three read-only tests, held against the committed artifacts.

WHAT THIS TEST IS FOR. RYA-1104 is a diagnostic: it ships no value, so there is no product
whose number could regress. What CAN regress is the diagnostic's own reasoning, and every
one of its three verdicts is a claim about a committed artifact that a later commit could
silently invalidate:

  * Test 1 says the Amarsi leg's 1D-LTE base is the SYNTHESIS pool. If the leg is re-run
    against a different base, or the provenance string stops naming one, the ticket's
    conclusion changes and nothing else in the repo would notice.
  * Test 2 says the Kurucz-floored subset is EMPTY. canonical_gf moves (RYA-822 extended it
    blueward mid-flight), and a tier demotion on one of these 50 lines would make the split
    real without anyone re-reading this ticket.
  * Test 3 says the feed's 0.000 is a difference-of-medians artefact and not a wiring bug.

⚠️ AND THE VACUITY GUARD IS THE POINT OF `test_test3_sides_are_different_objects`. A paired
differential of a product against ITSELF is 0.000 on every line and reads exactly like a
real null (RYA-853, RYA-1080: twice in one session a guard's two sides converged and the
guard went silent). So the ⟨3D⟩ pair is asserted to be two DISTINCT objects before its
difference is believed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.rya1104_amarsi_base_diagnostic import (  # noqa: E402
    HOLDINGS, base_stem_from_provenance, amarsi_leg, diagnose_base, diagnose_gf_split, diagnose_nlte_wiring)

HOLDING = "kpno_molecfit"


@pytest.fixture(scope="module")
def t1():
    return diagnose_base(HOLDING)


@pytest.fixture(scope="module")
def t2():
    return diagnose_gf_split(HOLDING)


@pytest.fixture(scope="module")
def t3():
    return diagnose_nlte_wiring(HOLDING)


# ── TEST 1 ────────────────────────────────────────────────────────────────────

def test_amarsi_base_is_the_synthesis_pool(t1):
    """The ticket's premise, held as an assertion. The base is SYNTH, not EW."""
    assert t1["declared_base"]["is_synthesis"] is True
    assert t1["declared_base"]["handler"] == "SynthesisHandler"
    assert t1["declared_base"]["route_axis"] == "synth"
    assert "_SYNTH_GRADED_" in t1["declared_base"]["artifact"]


def test_base_is_verified_value_for_value_not_taken_from_the_string(t1):
    """A provenance line is a CLAIM; the base file is the referee (RYA-1042)."""
    v = t1["base_verified_line_by_line"]
    assert v["n_matched_in_base"] == v["n_pool"] == 50
    assert v["max_abs_residual_dex"] == 0.0


def test_the_mlp_takes_no_line_strength_input(t1):
    """Test 1 step 4's STOP branch does not fire: the network is EW-agnostic."""
    f = t1["mlp_input_features"]
    assert f["takes_ew_or_reduced_ew"] is False
    assert "ew" not in [p.lower() for p in f["signature"]]


def test_the_offset_is_not_a_route_difference(t1):
    """+0.065 = line-subset + correction. Neither term is a route."""
    d = t1["decomposition"]
    assert d["synth_base_full_pool"]["n"] == 67
    assert d["synth_base_amarsi_subset"]["n"] == d["post_amarsi"]["n"] == 50
    assert d["line_subset_step_dex"] + d["amarsi_correction_step_dex"] == pytest.approx(
        d["total_vs_full_synth_base_dex"], abs=1e-3)


def test_the_ew_product_is_a_different_population_and_is_not_differenced(t1):
    """An EW-route product exists for this holding -- over 13 lines, not these 50.

    It is reported BESIDE the synthesis numbers and never subtracted from them. That is
    the whole lesson of the ticket: two aggregates over different populations look
    comparable and are not.
    """
    ew = t1["ew_route_product"]
    assert ew["exists"] is True
    assert ew["handler"] == "ProfileFitHandler"
    assert ew["n"] != t1["decomposition"]["post_amarsi"]["n"]
    assert ew["per_line_file_committed"] is False


def test_the_network_swings_with_the_line_population(t1):
    """The offset's actual locus: same network, two pools, ~+0.045 dex apart.

    The control side is not ours -- it is the RYA-817 reactivation control reproducing
    Amarsi's own published solar row from Amarsi's own inputs, and it PASSES. That is what
    makes the swing attributable to the pool rather than to the reactivation.
    """
    n = t1["network_on_two_pools"]
    assert n["control_pool"]["control_passed"] is True
    assert n["control_pool"]["mean_correction_dex"] == pytest.approx(-0.0024, abs=5e-4)
    assert n["our_pool"]["median_correction_dex"] == pytest.approx(0.043, abs=5e-4)
    assert n["swing_dex"] > 0.04
    # the populations really are disjoint at the low end -- that is the mechanism's premise
    assert n["our_pool"]["elo_eV"]["min"] > n["control_pool"]["elo_eV"]["max"] / 2
    assert n["control_pool"]["elo_eV"]["min"] < 0.5 < n["our_pool"]["elo_eV"]["min"]


def test_elo_is_not_claimed_to_explain_the_whole_swing(t1):
    """⚠️ The diagnostic must NOT over-attribute. Elo transport covers a minority of it.

    A tolerance, or a mechanism, needs its own null (RYA-1042). This holds the honest
    fraction: if a later change makes Elo appear to explain everything, that is a finding
    to re-derive, not a number to quietly accept.
    """
    frac = t1["network_on_two_pools"]["elo_trend"]["fraction_of_swing_explained"]
    assert 0.0 < frac < 0.6


def test_the_axis_rail_is_closed_as_a_confound(t1):
    """The band run rails A(Fe;3N) at 7.5 and the control runs at 7.46 -- immaterial."""
    ax = t1["network_on_two_pools"]["axis_confound_closed"]
    assert ax["railed_at_ceiling"] is True
    assert ax["sensitivity_across_whole_axis_dex"] < 0.005


# ── TEST 2 ────────────────────────────────────────────────────────────────────

def test_every_line_resolves_uniquely_on_the_dual_key(t2):
    """λ+EP, gf_resolver's own tolerances, no argmin, no silent drop (RYA-1033/RYA-429)."""
    m = t2["match"]
    assert m["n_unmatched"] == 0 and m["n_ambiguous"] == 0
    assert m["n_unique_match"] == t2["n_pool"] == 50


def test_the_kurucz_floored_subset_is_empty(t2):
    """Test 2's prediction is refuted: the split does not exist to be made."""
    assert t2["n_kurucz_floored"] == 0
    assert t2["n_lab_tier"] == 50
    assert t2["splittable"] is False
    assert sorted(t2["by_tier"]) == ["LAB"]


def test_the_gf_axis_now_AGREES_with_the_products_own_budget(t2):
    """RYA-1104 found one product making two claims about its oscillator strengths;
    RYA-1106 resolved it in favour of the one that was measured.

    🔴 THIS ASSERTION IS INVERTED ON PURPOSE. It read `product_gf_axis == "kurucz"` while
    the contradiction stood, which is what a diagnostic ticket's test should say: it pinned
    the defect so it could not be lost. The defect is now fixed at its source
    (`treatment_axes.LEGACY["ENGINE-A-3DNLTE"]["gf"]`), so the same test becomes the
    regression guard for the fix -- the budget still measures a lab-gf pool, and the axis
    now says so too. Leaving it asserting "kurucz" would have pinned the bug in place.
    """
    assert t2["product_gf_axis"] == "lab"
    assert "rung 3" in t2["budget_gf_rung"]
    assert "GF-LAB" in t2["budget_gf_rung"]


# ── TEST 3 ────────────────────────────────────────────────────────────────────

def test_test3_sides_are_different_objects(t3):
    """⚠️ THE VACUITY CONTROL. Assert the two sides DIFFER before believing their difference.

    RYA-853's referee ingested itself and RYA-1080's guard diffed a ref against itself;
    both produced a clean, well-formed, plausible null. A ⟨3D⟩ pair that had silently
    written to one path (RYA-915) would do the same here.
    """
    assert t3["ran"] is True
    w = t3["witnesses"]
    assert w["per_line_values_differ"]["distinct"] is True
    assert w["rya915_one_path_overwrite"]["distinct_paths"] is True
    assert w["rya915_one_path_overwrite"]["distinct_bytes"] is True
    assert t3["paired_differential"]["n_nonzero"] > 0


def test_departures_are_engaged(t3):
    """Not a wiring bug. Four independent witnesses, one of them the CODE."""
    assert t3["departures_engaged"] is True
    c = t3["witnesses"]["code_path"]
    assert c["passes_nlte_deck_only_when_nlte"] is True
    # RYA-1049: `for_node(..., node=...)` silently returns the 1D deck. The shipped call
    # site passes the ⟨3D⟩ deck key FIRST-POSITIONALLY, which is the reachable way.
    assert c["deck_key_is_first_positional"] is True


def test_the_feed_zero_is_the_wrong_statistic_not_a_zero(t3):
    """difference-of-medians 0.000 vs a paired +0.032 -- and RYA-1051 is corroborated."""
    q = t3["paired_differential"]
    assert t3["published"]["difference_of_published_values"] == 0.0
    assert q["median"] == pytest.approx(0.032, abs=1e-3)
    assert q["mean"] == pytest.approx(0.0267, abs=1e-3)
    assert q["n_nonzero"] == 66 and q["n_paired"] == 67
    assert q["collision"] is True


def test_the_collision_is_a_coincidence_not_a_quantiser(t3):
    """If it were a rounding floor it would reproduce; a coincidence will not (RYA-1083)."""
    nq = t3["not_a_quantiser"]
    assert nq["min_gap_dex"] <= 0.001
    assert nq["n_unique_abundance_values"] > 100


def test_the_nlte_products_own_departure_columns_read_as_lte(t3):
    """⚠️ A REAL BOOKKEEPING DEFECT, held so it is not lost when this ticket closes.

    The abundances DO carry the departures -- that is what Test 3 establishes -- but
    `nlte_delta_dex` is 0.0 and `nlte_source` says "no departure applied" on the NLTE
    member. A later reader asking the artifact whether departures ran gets the wrong
    answer from the column built to tell them.

    Asserted as PRESENT rather than fixed: RYA-1104 ships nothing.
    """
    b = t3["bookkeeping_defect"]
    assert b["present"] is True
    assert b["contradicts_the_products_own_provenance"] is True


# ── the whole thing, on every holding ─────────────────────────────────────────

@pytest.mark.parametrize("holding", sorted(HOLDINGS))
def test_every_amarsi_leg_names_a_synthesis_base(holding):
    """Not one arm -- all four. A conclusion drawn on kpno_molecfit alone would be a
    property of that holding, and the ticket asks for a check holding for exactly that
    reason."""
    prod, _lines, _s = amarsi_leg(holding)
    _stem, base = base_stem_from_provenance(prod.iloc[0])
    assert "_SYNTH_" in base.name
    assert "_PROFILEFIT_" not in base.name


@pytest.mark.parametrize("holding", sorted(HOLDINGS))
def test_no_amarsi_pool_has_a_kurucz_floored_line(holding):
    t = diagnose_gf_split(holding)
    assert t["n_kurucz_floored"] == 0, t["by_tier"]
