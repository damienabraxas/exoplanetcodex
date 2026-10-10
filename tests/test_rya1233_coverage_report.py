"""RYA-1233: the coverage table is the orchestrator's definition of DONE.

Two lookups it got wrong on Si, each of which reported finished work as NOT_RUN:
  - the run report keys a cell's pool by the FEED token (pool_tier: GRADED / SET-<NAME>),
    not the orchestrator's pool name ("codex");
  - a product the budget stage built and HELD (RYA-587 verdict) or skipped as empty is HELD
    with that reason, read from the stage's own report.json.
"""
import json

from pipeline import coverage_report, run_matrix as rm


def test_si_owes_nothing_unrun_and_budget_holds_carry_their_reason():
    cov = coverage_report.coverage("solar", "Si")
    assert cov["counts"]["NOT_RUN"] == 0, [r for r in cov["rows"] if r["outcome"] == "NOT_RUN"][:5]
    bud = rm.REPORT_DIR / "budget" / "solar_Si" / "report.json"
    held = [o for o in json.loads(bud.read_text())["products"]
            if (o.get("verdict") or o.get("skip")) and o.get("key_treatment") != "ENGINE-B"]
    assert held, "fixture: the committed Si budget holds some products"
    reasons = {r["reason"] for r in cov["rows"] if r["outcome"] == "HELD"}
    assert any(x.startswith("budget: ") for x in reasons)


def test_a_gerber_cell_without_a_deck_is_held_not_owed():
    """Gerber 2023 decks are per element (Fe, Al): a Si gerber cell is HELD with that reason."""
    cov = coverage_report.coverage("solar", "Si")
    g = [r for r in cov["rows"] if r["treatment"] == "ENGINE-B-NLTE"]
    assert g and all(r["outcome"] == "HELD" and "no Gerber deck for Si" in r["reason"] for r in g)


def test_done_needs_a_build_record_on_todays_inputs():
    """A published cell with no orchestrator build record, or one built on other inputs, is
    STALE and the element is not done."""
    cov = coverage_report.coverage("solar", "Si")
    for r in cov["rows"]:
        if r["outcome"] == "STALE":
            assert r["reason"] and r["A"] is not None
    assert cov["done"] == (cov["counts"]["NOT_RUN"] == 0 and cov["counts"]["STALE"] == 0)


def test_every_published_si_product_is_a_published_cell():
    cov = coverage_report.coverage("solar", "Si")
    assert cov["counts"]["PUBLISHED"] + cov["counts"]["STALE"] > 0
    for r in cov["rows"]:
        if r["outcome"] in ("PUBLISHED", "STALE"):
            assert r["A"] is not None and r["sigma"] is not None
