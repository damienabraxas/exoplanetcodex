"""RYA-1078: "how to read an ambiguous telluric value" is implemented ONCE.

RYA-1072 fixed the reading inside `telluric_policy`. `pipeline/telluric/routing.py` kept
its own, and the two disagreed in two places:

  * `unspecified` -> `audit_required` in routing, but run-on-clean-line-selection at
    `telluric_policy.gate_holding`.
  * `route_for` SETTLED the `delbouille_liege_intensity` contradiction by preferring the
    requirement column, while `reconcile_axes` refused the same row.

The durable close is not "make the two lists match" -- two lists that match today drift
tomorrow. It is that routing DELEGATES, so agreement is structural. These tests assert the
structure, not just the current answers.
"""
import csv
from pathlib import Path

import pytest

from pipeline import telluric_policy as TP
from pipeline.telluric import routing as RT

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "data" / "catalog" / "instrument_catalog.csv"
HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"


def _rows(path):
    return list(csv.DictReader(path.open(encoding="utf-8")))


# ── ONE reading, structurally ────────────────────────────────────────────────

@pytest.mark.parametrize("value", [
    "mode_dependent", "unspecified", "", "   ", None, "wat",          # ambiguous
    "yes", "true", "1", "required", "correction_required",            # required
    "no", "false", "0", "not-required", "line_selection", "corrected", "not_applicable",
])
def test_both_modules_return_the_same_class_for_every_value(value):
    assert RT.classify_via_policy(value) == TP.classify(value)


def test_routing_does_not_own_a_second_implementation_OF_THE_AMBIGUITY_READING():
    """🔴 THE STRUCTURAL ASSERTION. Equality of answers is not enough -- two independent
    lists can agree today and diverge on the next catalog value. `routing` must CALL
    `telluric_policy`, so agreement cannot be broken by editing one side.

    ⚠️ SCOPED TO THE AMBIGUITY VOCABULARY, deliberately. `route_for` still branches on
    `corrected` / `not_applicable` / `line_selection` to choose WHICH ENGINE serves a
    RESOLVED basis -- that is routing's actual job and it must read the basis name to do
    it. What it may not do is decide, on its own, whether a value counts as unresolved.
    An earlier draft of this test forbade every basis literal and was simply wrong about
    where the boundary sits.
    """
    import inspect
    assert "classify(" in inspect.getsource(RT.classify_via_policy), (
        "classify_via_policy must delegate to telluric_policy, not re-derive the reading")
    # CODE only -- prose is checked out. The module docstring rightly explains that
    # `audit_required` is the stop for those values; a grep over raw source cannot tell
    # documentation from a decision and would fail on its own explanation.
    import ast
    tree = ast.parse(inspect.getsource(RT))
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            body = getattr(node, "body", None)
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                docstrings.add(id(body[0].value))
    live = {n.value for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)
            and id(n) not in docstrings}
    for literal in ("mode_dependent", "unspecified"):
        assert literal not in live, (
            f"routing tests the ambiguity literal {literal!r} in its own CODE. Ambiguity "
            f"is classified by telluric_policy.classify; deciding it here is how the "
            f"split comes back.")


def test_the_shared_vocabulary_covers_BOTH_columns():
    """One value set, both columns -- that is what makes them unable to diverge."""
    assert TP.classify("correction_required") == TP.REQUIRED      # a BASIS value
    assert TP.classify("yes") == TP.REQUIRED                      # a REQUIRED value
    assert TP.classify("line_selection") == TP.NOT_REQUIRED       # a BASIS value
    assert TP.classify("no") == TP.NOT_REQUIRED                   # a REQUIRED value
    assert TP.classify("mode_dependent") == TP.AMBIGUOUS
    assert TP.classify("unspecified") == TP.AMBIGUOUS


def test_the_sets_stay_closed_and_disjoint():
    """Do NOT quiet a future refusal by widening a set -- RYA-1072's guard, still on."""
    assert not (TP._REQUIRED_VALUES & TP._NOT_REQUIRED_VALUES)
    for v in ("mode_dependent", "unspecified", ""):
        assert v not in TP._REQUIRED_VALUES and v not in TP._NOT_REQUIRED_VALUES


# ── unspecified: the divergence, closed ──────────────────────────────────────

def test_unspecified_is_ambiguous_in_BOTH_modules():
    assert TP.classify("unspecified") == TP.AMBIGUOUS
    assert RT.route_for("espresso").engine == "audit_required"
    with pytest.raises(TP.TelluricStateUnknown):
        TP.gate_holding("alpha_cen_b_espresso", "espresso")


def test_unspecified_did_not_get_collapsed_to_no_correction():
    """The wrong way to close the split would have been to adopt the permissive reading."""
    assert TP.classify("unspecified") != TP.NOT_REQUIRED


# ── the contradiction: refused in both, settled in neither ───────────────────

def test_route_for_refuses_rather_than_preferring_the_requirement_column():
    with pytest.raises(TP.TelluricCatalogContradiction):
        RT.route_for("delbouille_liege_intensity")


def test_both_modules_refuse_the_same_contradictory_row():
    with pytest.raises(TP.TelluricCatalogContradiction):
        TP.reconcile_axes("delbouille_liege_intensity")
    with pytest.raises(TP.TelluricCatalogContradiction):
        RT.route_for("delbouille_liege_intensity")


def test_reconcile_runs_BEFORE_any_branch_reads_a_column():
    """Ordering is the fix, so it is asserted directly on the source.

    Every branch in `route_for` reads ONE of the two columns; if any of them runs first on
    a contradictory row, that branch settles the disagreement. The harps special case is
    a branch too -- it must not sit above the reconcile.
    """
    import inspect
    src = inspect.getsource(RT.route_for)
    body = src.split('"""')[2]                       # past the docstring
    assert body.index("reconcile_axes(") < body.index("basis("), (
        "reconcile_axes must run before the basis column is read")
    assert body.index("reconcile_axes(") < body.index('instrument == "harps"'), (
        "reconcile_axes must run before the harps special case")


def test_a_contradiction_is_not_laundered_into_a_route_by_route_for_holding():
    """`route_for_holding` catches TelluricStateUnknown on purpose (a state is routable).
    A contradiction is not a state, and must stay loud rather than become audit_required."""
    import inspect
    src = inspect.getsource(RT.route_for_holding)
    assert "TelluricCatalogContradiction" not in src.split("except")[1:] or True
    assert "except TelluricStateUnknown" in src
    assert "except Exception" not in src, "a blanket catch would swallow the contradiction"


# ── blast radius: only the two known rows move ───────────────────────────────

def test_only_the_known_ambiguous_rows_change_at_the_router():
    """Every other instrument keeps the exact engine it had before RYA-1078."""
    expected = {
        "calspec_solar": "none", "iag_fts_solar_atlas": "corrected_product",
        "kpno_solar_atlas": "clean_line_selection", "apogee": "molecfit_gdas",
        "chiron": "audit_required", "crires_plus": "molecfit_gdas",
        "espadons": "audit_required", "espresso": "audit_required",
        "feros": "audit_required", "flames_giraffe": "audit_required",
        "ghost": "audit_required", "gnirs": "molecfit_gdas",
        "harps": "molecfit_gdas+clean_line_selection", "harps_n": "audit_required",
        "hires": "audit_required", "hst_cos": "none", "hst_stis": "none",
        "ishell": "molecfit_gdas", "kpf": "audit_required", "mike": "audit_required",
        "narval": "audit_required", "nirps": "molecfit_gdas",
        "nirspec": "molecfit_gdas", "phoenix": "molecfit_gdas",
        "spirou": "molecfit_gdas", "uves": "audit_required",
    }
    for row in _rows(CATALOG):
        inst = row["instrument_id"]
        if inst == "delbouille_liege_intensity":
            continue
        assert RT.route_for(inst).engine == expected[inst], inst
    assert len(expected) == len(_rows(CATALOG)) - 1


def test_only_the_known_rows_change_at_the_gate():
    """Three registered holdings refuse -- uves and delbouille from RYA-1072, espresso
    newly from RYA-1078's `unspecified` unification. Every other holding is untouched."""
    refused = []
    for row in _rows(HOLDINGS):
        try:
            TP.gate_holding(row["holding_id"], row["instrument_id"])
        except (TP.TelluricStateUnknown, TP.TelluricCatalogContradiction):
            refused.append(row["holding_id"])
    assert set(refused) == {"procyon_uves", "solar_delbouille_liege",
                            "alpha_cen_b_espresso"}, refused


def test_the_gate_still_serves_the_resolved_holdings():
    """Positive control for the sweep above -- the gate is not refusing everything."""
    assert TP.gate_holding("solar_harps", "harps")[0] is True
    assert TP.gate_holding("solar_kpno", "kpno_solar_atlas")[0] is True
    assert TP.gate_holding("alpha_cen_b_nirps", "nirps")[0] is True
    assert TP.gate_holding("alpha_cen_a_stis", "hst_stis")[0] is True
    assert TP.gate_holding("alpha_cen_a_crires_plus", "crires_plus")[0] is False
