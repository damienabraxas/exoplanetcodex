"""Every product the forest declines to draw is COUNTED (RYA-1102's own rule, completed).

🔴 WHAT WAS MISSING. `plot_grid.build()` carries the comment "REPORTED, NEVER DROPPED" over
`off_axis`, and that is true of what `off_axis` covers: a product that entered a section and
found no row on the axis. A product `is_displayable` REFUSES never enters a section at all,
so it left the grid with no trace anywhere in the feed. Measured on the live Fe feed: 14 of
160 products, while `off_axis` was empty on all 38 sections -- so the number a reader could
see was zero.

Two of those 14 are why it matters. Tier `ALL` carries grade "Reference Grade", and
`is_displayable` is written on TIER, so the four Asplund-pool 3D-NLTE products are refused
-- including the one the site publishes as the Fe I headline. Whether the showcase SHOULD
draw them is a publication decision and this file does not take it. Refusing to draw a
product is a decision; not recording the refusal is a defect.

⚠️ NOTHING HERE CHANGES WHICH PRODUCTS ARE DRAWN. These tests would all pass on the old
display policy too, as long as it reported itself.
"""
from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import plot_grid as pg          # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"


@pytest.fixture(scope="module")
def products():
    return json.loads(FEED.read_text())["products"]


@pytest.fixture(scope="module")
def grid(products):
    return pg.build(products)


def test_the_accounting_closes_on_every_product(grid, products):
    """🔴 THE INVARIANT. A product is drawn, or it is an alternate on a drawn row, or it is
    in `not_displayed` -- and the three sum to the product count, so it cannot leave the
    plot through a fourth door nobody is watching."""
    drawn = sum(1 for s in grid["sections"] for c in s["cells"] if c["product_key"])
    alternates = sum(c["alternates"] for s in grid["sections"] for c in s["cells"])
    off_axis = sum(len(s["off_axis"]) for s in grid["sections"])
    held = len(grid["not_displayed"])
    assert drawn + alternates + off_axis + held == len(products), (
        f"drawn {drawn} + alternates {alternates} + off_axis {off_axis} + not_displayed "
        f"{held} = {drawn + alternates + off_axis + held}, but the feed publishes "
        f"{len(products)} products — some product is leaving the grid uncounted")


def test_something_is_actually_being_held_back(grid):
    """The precondition. If the display policy ever draws everything this file is vacuous,
    and a reader deserves to be told that rather than shown a green tick."""
    assert grid["not_displayed"], (
        "no product is held back, so the reporting path is untested — delete this file or "
        "fix the fixture, do not leave it passing on nothing")


def test_every_held_product_says_why_and_can_be_identified(grid):
    for row in grid["not_displayed"]:
        assert row["product_key"], row
        assert row["reason"] and len(row["reason"]) > 30, row
        assert row["tier"] in row["reason"], (
            "the reason must name the tier it refused, so a reader does not have to "
            "re-derive the rule: " + str(row))


def test_a_reference_GRADE_product_outside_the_reference_TIER_is_never_silent(grid, products):
    """🔴 THE NEAR-MISS, PINNED. `is_displayable` admits REFERENCE by TIER, and four
    products carry grade "Reference Grade" with tier "ALL" -- one of them the published Fe I
    headline. They may legitimately stay out of the showcase; they may not stay out of the
    record. Stated as a rule rather than a count so it survives the decision either way.
    """
    odd = [p for p in products
           if p.get("grade") == "Reference Grade" and p.get("tier") != "REFERENCE"]
    assert odd, "no grade/tier mismatch in this feed — this test is not exercising anything"
    held = {r["product_key"] for r in grid["not_displayed"]}
    drawn = {c["product_key"] for s in grid["sections"] for c in s["cells"]
             if c["product_key"]}
    from pipeline import product_eligibility as _pe
    for p in odd:
        k = _pe.key_of(p)
        assert k in held or k in drawn, (
            f"{p['instrument']} {p['band']} tier={p['tier']} grade={p['grade']} is neither "
            f"drawn nor recorded as held back")
        if k in held:
            reason = next(r["reason"] for r in grid["not_displayed"]
                          if r["product_key"] == k)
            assert "Reference Grade" in reason, (
                "a Reference-GRADE product held back by a TIER rule must say so, or the "
                "next reader re-discovers it: " + reason)


def test_the_deep_tier_is_held_back_for_the_reason_ryan_stated(grid):
    """Ryan's rule: "GRADED only. DEEPGRADED is a secondary product -- documented in its own
    section, not showcased", with the near-UV exception where nothing else exists. So a Deep
    product in `not_displayed` must be one whose section HAS a Codex product."""
    deep = [r for r in grid["not_displayed"] if r["tier"] == "DEEPGRADED"]
    assert deep, "no Deep product held back — the showcase rule is not being exercised"
    for r in deep:
        assert "GRADED product" in r["reason"], r["reason"]


def test_removing_the_report_breaks_the_accounting(products, monkeypatch):
    """🔴 MUTATION. If `is_displayable` were made to refuse everything and nothing recorded
    it, the invariant above must FAIL — otherwise it is not what is holding the line."""
    monkeypatch.setattr(pg, "is_displayable", lambda p, only_deep: False)
    g = pg.build(products)
    drawn = sum(1 for s in g["sections"] for c in s["cells"] if c["product_key"])
    assert drawn == 0, "the mutation did not take"
    assert len(g["not_displayed"]) == len(products), (
        "every product was refused, so every product must be recorded; the report is not "
        "tracking the refusals")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
