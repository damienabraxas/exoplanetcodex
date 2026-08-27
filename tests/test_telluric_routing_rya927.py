"""RYA-927 — one telluric route contract for every catalog instrument."""
from pathlib import Path

import pandas as pd
import pytest

from pipeline.telluric.routing import (classify_via_policy, correction_needed_at,
                                       route_for, route_for_holding)
from pipeline.telluric_policy import TelluricCatalogContradiction

ROOT = Path(__file__).resolve().parents[1]


def test_every_catalog_instrument_has_an_EXPLICIT_DISPOSITION():
    """RYA-927's contract, kept and STRENGTHENED by RYA-1078.

    The original assertion was "every instrument returns one of six known engines". That
    stayed true for `delbouille_liege_intensity` only because `route_for` SETTLED its
    contradictory row (`telluric_required=yes` AND `telluric_basis=line_selection`) by
    preferring the requirement column -- so the test was passing ON the defect, certifying
    a route the catalog does not actually license.

    The contract RYA-927 wanted is not "always returns a route"; it is NEVER SILENTLY
    GUESSES. So every instrument must reach an explicit disposition, and there are now two
    legitimate kinds: a named route, or a named REFUSAL. What is still forbidden -- and is
    what this test exists to catch -- is a third outcome sneaking in, or a refusal that
    does not say what it is refusing.
    """
    df = pd.read_csv(ROOT / "data/catalog/instrument_catalog.csv")
    allowed = {"molecfit_gdas", "molecfit_gdas+clean_line_selection", "corrected_product",
               "clean_line_selection", "none", "audit_required"}
    routed, refused = [], []
    for instrument in df.instrument_id:
        try:
            route = route_for(instrument)
        except TelluricCatalogContradiction as exc:
            assert instrument in str(exc), "a refusal must name the instrument it refuses"
            assert "telluric_required" in str(exc) and "telluric_basis" in str(exc), (
                "a contradiction must name BOTH columns -- otherwise the reader cannot "
                "tell which two values disagree")
            refused.append(instrument)
            continue
        assert route.engine in allowed
        assert route.reason
        routed.append(instrument)

    assert refused == ["delbouille_liege_intensity"], (
        f"expected exactly the one known contradictory catalog row to be refused, got "
        f"{refused}. A NEW name here means a catalog row was edited into disagreement; a "
        f"MISSING name means route_for went back to settling a contradiction silently.")
    assert len(routed) == len(df) - 1


def test_harps_and_kitt_peak_are_clean_line_routes():
    assert route_for("harps").engine == "molecfit_gdas+clean_line_selection"
    assert route_for("kpno_solar_atlas").engine == "clean_line_selection"
    assert not correction_needed_at(6631.218, "harps")
    assert not correction_needed_at(6696.185, "harps")
    assert correction_needed_at(6872.0, "harps")


def test_correction_required_instruments_use_observation_aware_engine():
    for instrument in ("crires_plus", "nirps", "spirou", "apogee", "nirspec"):
        route = route_for(instrument)
        assert route.engine == "molecfit_gdas"
        assert route.correction_required is True
        assert route.analysis_ready_required is True
    assert correction_needed_at(6872.0, "crires_plus")


def test_unresolved_basis_is_a_stop_not_a_fallback():
    route = route_for("feros")
    assert route.engine == "audit_required"
    assert route.analysis_ready_required is True


def test_crires_holding_state_overrides_instrument_capability():
    corrected = route_for_holding("solar_crires_plus_y_rya794")
    assert corrected.engine == "corrected_product"
    assert corrected.correction_required is False
    raw = route_for_holding("solar_vesta_crires_plus_idp")
    assert raw.engine == "molecfit_gdas"
    assert raw.correction_required is True


# ── RYA-1078: the two modules read ambiguity from ONE place ──────────────────

def test_routing_delegates_the_reading_and_does_not_own_a_second_one():
    """The split closed. Not "the two lists agree" -- the same function, called twice."""
    from pipeline import telluric_policy as tp
    for value in ("mode_dependent", "unspecified", "", "correction_required",
                  "line_selection", "yes", "no", "not_applicable", "corrected", "wat"):
        assert classify_via_policy(value) == tp.classify(value)


def test_unspecified_resolves_the_SAME_WAY_in_both_modules():
    """The divergence RYA-1078 exists to close, asserted at the module boundary.

    Before: `unspecified` meant `audit_required` in routing and ran on clean-line
    selection in `telluric_policy.gate_holding`. One value, two answers, in the two
    modules that both claim to own the decision.
    """
    from pipeline import telluric_policy as tp
    assert tp.classify("unspecified") == tp.AMBIGUOUS
    assert route_for("espresso").engine == "audit_required"      # basis unspecified
    with pytest.raises(tp.TelluricStateUnknown):
        tp.gate_holding("alpha_cen_b_espresso", "espresso")


def test_route_for_REFUSES_the_contradiction_instead_of_preferring_a_column():
    """🔴 THE DEFECT. Before RYA-1078 this returned engine='molecfit_gdas' -- the
    requirement column silently winning over the basis column."""
    with pytest.raises(TelluricCatalogContradiction):
        route_for("delbouille_liege_intensity")


def test_the_refusal_did_not_cost_any_resolved_instrument_its_route():
    """Positive control: refusing must not be what route_for now does in general."""
    assert route_for("crires_plus").engine == "molecfit_gdas"
    assert route_for("kpno_solar_atlas").engine == "clean_line_selection"
    assert route_for("iag_fts_solar_atlas").engine == "corrected_product"
    assert route_for("hst_stis").engine == "none"
    assert route_for("harps").engine == "molecfit_gdas+clean_line_selection"
    assert route_for("uves").engine == "audit_required"
