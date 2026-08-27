"""RYA-1072: a declared AMBIGUITY on the telluric requirement never resolves to `no`.

THE DEFECT. `requires_correction` was `return v in ('yes','true','1','required')` -- an
allow-list for the TRUE case with everything else falling into FALSE. `mode_dependent`,
a value `instrument_catalog.csv` actually uses, therefore took the identical branch to a
declared `no`, and `gate_holding` served UVES on the clean-line-selection basis with no
mode ever consulted. UVES RED860 reaches 10427 A, inside the registered H2O complexes.

EVERY TEST HERE CARRIES ITS POSITIVE CONTROL, in the same test where that is possible.
A function that always raised, or always returned a constant, must fail -- otherwise
this file could pass while gating everything forever and prove nothing about the rule.
"""
import csv
from pathlib import Path

import pytest

from pipeline import telluric_policy as TP

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "data" / "catalog" / "instrument_catalog.csv"


def _catalog_rows():
    return list(csv.DictReader(CATALOG.open(encoding="utf-8")))


# ── the three-outcome reading ────────────────────────────────────────────────

@pytest.mark.parametrize("value", ["mode_dependent", "unspecified", "", "   ",
                                   "maybe", "tbd", None, "YES please"])
def test_an_unrecognised_value_is_AMBIGUOUS_and_never_NOT_REQUIRED(value):
    assert TP.correction_requirement_from_value(value) == TP.AMBIGUOUS
    assert TP.correction_requirement_from_value(value) != TP.NOT_REQUIRED


@pytest.mark.parametrize("value,expected", [
    ("yes", TP.REQUIRED), ("YES", TP.REQUIRED), (" Yes ", TP.REQUIRED),
    ("true", TP.REQUIRED), ("1", TP.REQUIRED), ("required", TP.REQUIRED),
    ("no", TP.NOT_REQUIRED), ("NO", TP.NOT_REQUIRED), ("false", TP.NOT_REQUIRED),
    ("0", TP.NOT_REQUIRED), ("not-required", TP.NOT_REQUIRED),
])
def test_the_resolved_values_still_read_exactly_as_before(value, expected):
    """The positive control for the parametrisation above: if AMBIGUOUS had swallowed
    everything, these would fail and the ambiguity tests would be vacuous."""
    assert TP.correction_requirement_from_value(value) == expected


def test_ambiguous_yes_and_no_in_ONE_test_so_a_constant_cannot_pass():
    """The ticket's explicit requirement. All three outcomes, one assertion block."""
    with pytest.raises(TP.TelluricStateUnknown):
        TP.requires_correction_from_value("mode_dependent")
    assert TP.requires_correction_from_value("yes") is True
    assert TP.requires_correction_from_value("no") is False


def test_the_ambiguous_refusal_names_the_value_and_the_fix():
    with pytest.raises(TP.TelluricStateUnknown) as exc:
        TP.requires_correction_from_value("mode_dependent")
    msg = str(exc.value)
    assert "mode_dependent" in msg
    assert "RYA-1072" in msg


# ── the recognised sets are CLOSED, and must stay closed ─────────────────────

def test_mode_dependent_is_in_neither_recognised_set():
    """🔴 Widening a set to swallow `mode_dependent` MASKS the bug -- it does not fix it.

    This is the regression guard the ticket asks for by name. If someone ever `fixes` a
    future refusal by adding the ambiguous value to one of these sets, this fails.
    """
    assert "mode_dependent" not in TP._REQUIRED_VALUES
    assert "mode_dependent" not in TP._NOT_REQUIRED_VALUES
    assert "unspecified" not in TP._REQUIRED_VALUES
    assert "unspecified" not in TP._NOT_REQUIRED_VALUES
    assert not (TP._REQUIRED_VALUES & TP._NOT_REQUIRED_VALUES), "the sets must be disjoint"


def test_every_catalog_value_reads_as_one_of_the_three_outcomes():
    """Nothing in the live catalog falls off the end of the reading."""
    seen = {TP.correction_requirement_from_value(r["telluric_required"])
            for r in _catalog_rows()}
    assert seen <= {TP.REQUIRED, TP.NOT_REQUIRED, TP.AMBIGUOUS}
    assert TP.AMBIGUOUS in seen, (
        "the catalog is expected to still carry at least one unresolved requirement "
        "(uves = mode_dependent); if it no longer does, this guard has gone vacuous")


def test_uves_is_the_ambiguous_row_and_the_IR_arms_are_not():
    """Measured against the live catalog, both directions in one test."""
    by_id = {r["instrument_id"]: r for r in _catalog_rows()}
    assert TP.correction_requirement("uves") == TP.AMBIGUOUS
    assert by_id["uves"]["telluric_required"] == "mode_dependent"
    for inst in ("crires_plus", "nirps", "spirou", "apogee"):
        assert TP.correction_requirement(inst) == TP.REQUIRED
    for inst in ("harps", "kpno_solar_atlas", "hst_stis"):
        assert TP.correction_requirement(inst) == TP.NOT_REQUIRED


# ── the GATE refuses, and does not fall through to clean-line selection ──────

def test_requires_correction_refuses_uves_and_still_answers_the_resolved_arms():
    with pytest.raises(TP.TelluricStateUnknown):
        TP.requires_correction("uves")
    assert TP.requires_correction("crires_plus") is True
    assert TP.requires_correction("harps") is False


def test_gate_holding_REFUSES_the_ambiguous_holding_rather_than_serving_it():
    """🔴 THE CRITICAL CASE. Before RYA-1072 this returned

        (True, "procyon_uves: uves is registered telluric_required=no on basis
                'mode_dependent', so tellurics are handled by per-line clean-line
                selection ...")

    -- a permission quoting a value the catalog does not contain.
    """
    with pytest.raises(TP.TelluricStateUnknown) as exc:
        TP.gate_holding("procyon_uves", "uves")
    assert "mode_dependent" in str(exc.value)


def test_gate_REFUSES_the_ambiguous_instrument():
    with pytest.raises(TP.TelluricStateUnknown):
        TP.gate("uves")


def test_the_gate_still_passes_and_still_refuses_every_resolved_holding():
    """Positive control for the two refusals above, over all four resolved outcomes."""
    assert TP.gate_holding("solar_harps", "harps")[0] is True            # line_selection
    assert TP.gate_holding("alpha_cen_a_stis", "hst_stis")[0] is True    # not_applicable
    assert TP.gate_holding("alpha_cen_b_nirps", "nirps")[0] is True      # applied
    assert TP.gate_holding("alpha_cen_a_crires_plus", "crires_plus")[0] is False
    assert TP.gate("crires_plus", analysis_ready=False)[0] is False
    assert TP.gate("crires_plus", analysis_ready=True)[0] is True
    assert TP.gate("harps")[0] is True


def test_no_gate_verdict_anywhere_claims_telluric_required_no_for_an_ambiguous_row():
    """The exact false sentence the defect produced must be unreachable.

    Sweeps every registered holding: a permissive verdict may only be returned for an
    instrument whose requirement actually resolves.
    """
    rows = list(csv.DictReader(
        (ROOT / "data" / "catalog" / "holdings_manifest_registry.csv").open(encoding="utf-8")))
    for r in rows:
        inst = r["instrument_id"]
        try:
            ok, why = TP.gate_holding(r["holding_id"], inst)
        except (TP.TelluricStateUnknown, TP.TelluricCatalogContradiction):
            continue
        if ok and "telluric_required=no" in why:
            assert TP.correction_requirement(inst) == TP.NOT_REQUIRED, (
                f"{r['holding_id']} was served on a 'telluric_required=no' basis while "
                f"{inst} reads {TP.correction_requirement(inst)} -- this is the RYA-1072 "
                f"defect exactly")


# ── ledger warning 2: the two catalog axes must not contradict each other ────

def test_the_contradictory_catalog_row_is_refused_not_silently_resolved():
    """`delbouille_liege_intensity` is telluric_required=yes AND basis=line_selection.

    The guard does NOT decide which value is stale -- that is Ryan's data call. It
    refuses so the contradiction cannot be settled by whichever column a consumer reads.
    """
    with pytest.raises(TP.TelluricCatalogContradiction) as exc:
        TP.reconcile_axes("delbouille_liege_intensity")
    assert "line_selection" in str(exc.value)
    with pytest.raises(TP.TelluricCatalogContradiction):
        TP.gate_holding("solar_delbouille_liege", "delbouille_liege_intensity")


def test_delbouille_is_the_ONLY_contradictory_row_in_the_catalog():
    """Positive control: reconcile_axes must not be a function that refuses everything."""
    bad = []
    for r in _catalog_rows():
        try:
            TP.reconcile_axes(r["instrument_id"])
        except TP.TelluricCatalogContradiction:
            bad.append(r["instrument_id"])
    assert bad == ["delbouille_liege_intensity"], (
        f"expected exactly the one known contradictory row, got {bad}")


def test_an_ambiguous_basis_asserts_nothing_and_so_cannot_contradict():
    """uves is mode_dependent on BOTH axes -- a consistent ambiguity, not a contradiction.

    The two failure modes need different fixes (determine a value vs. decide which of two
    is stale), so they must not collapse into one error type.
    """
    TP.reconcile_axes("uves")          # must not raise
    TP.reconcile_axes("espresso")      # telluric_required=no + basis unspecified
    assert TP.correction_requirement("uves") == TP.AMBIGUOUS


# ── the callers that must survive an unresolved instrument ───────────────────

def test_the_router_routes_the_ambiguity_to_audit_rather_than_crashing():
    from pipeline.telluric.routing import route_for
    assert route_for("uves").engine == "audit_required"
    assert route_for("crires_plus").engine == "molecfit_gdas"
    assert route_for("kpno_solar_atlas").engine == "clean_line_selection"
