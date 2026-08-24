"""The published product store keeps its four promises — RYA-1034.

A real four-arm 3D-NLTE result (A(Fe I) VIS = 7.511, n=50) was produced on Sirius and
then could not be found: it lived in `~/out_g3d_9e651a/` on the only machine that can run
that leg, while every Mac copy still held the superseded ungraded 7.604 and
`data/results/band_products/` is gitignored. The value survived only in a Linear ticket.

These tests pin the properties that make `data/products/<star>/<El>.json` a source of
truth rather than another copy. They test the SHAPE, not the values -- the values are
supposed to change; the guarantees are not.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT / "data" / "products" / "solar"
KEY_FIELDS = ("element", "ion", "band", "instrument", "holding",
              "tier", "selector", "route", "treatment")


def _docs():
    return [(p, json.loads(p.read_text())) for p in sorted(STORE.glob("*.json"))]


def test_there_is_at_least_one_published_element():
    assert _docs(), (
        f"no element product published under {STORE}. This store is the source of "
        "truth; an empty one means results are living in ad-hoc directories again.")


@pytest.mark.parametrize("field", ["schema", "element", "version", "products",
                                   "archive", "quarantine"])
def test_every_document_carries_the_required_sections(field):
    for p, d in _docs():
        assert field in d, f"{p.name} has no `{field}`"


def test_keys_are_unique_among_current_products():
    """🔴 THE COLLISION THIS CAUGHT ON DAY ONE. The key first omitted route and selector,
    so IAG VIS GRADED ENGINE-A was 7.481 (n=4, profile fit) AND 7.484 (n=37, synthesis) --
    two different pools measured two different ways, one silently overwriting the other.
    RYA-1026 lists `1D-LTE` and `1D-LTE Synth` as separate rows, and RYA-984 makes two
    runs differing only in selector two products. The key must say so."""
    for p, d in _docs():
        seen = {}
        for row in d["products"]:
            k = "|".join(str(row.get(f) or "") for f in KEY_FIELDS)
            assert k not in seen, (
                f"{p.name}: duplicate current product for key {k}\n"
                f"    {seen.get(k)} vs {row.get('A')} — two products cannot share a key")
            seen[k] = row.get("A")


def test_every_product_states_where_it_came_from():
    """Provenance travels with the number, or the number is unverifiable. A Sirius-only
    result copied here without its origin path claims to have been produced here, which
    is the laundering RYA-772 names."""
    for p, d in _docs():
        for row in d["products"] + d["quarantine"]:
            pr = row.get("provenance") or {}
            for f in ("host", "path", "sha256"):
                assert pr.get(f), f"{p.name}: product {row.get('A')} has no provenance.{f}"
            assert len(pr["sha256"]) == 64, f"{p.name}: sha256 is not a sha256"


def test_nothing_is_withdrawn_or_replaced_without_a_stated_reason():
    """RYA-711: quarantined, never culled -- and a withdrawal with no reason is
    indistinguishable from a loss."""
    for p, d in _docs():
        for row in d["quarantine"]:
            assert row.get("quarantine_reason"), (
                f"{p.name}: a quarantined product carries no reason. It cannot be "
                "defended in the appendix, which is the only thing that makes a "
                "rejection legitimate (RYA-844).")
            assert row.get("quarantined_at")
        for row in d["archive"]:
            assert row.get("superseded_reason"), (
                f"{p.name}: an archived value carries no reason for being replaced.")


def test_a_published_value_is_a_number_with_a_line_count():
    """n is not decoration: it differs per arm for real reasons (IAG starts at 5001 A;
    the NLTE grid does not serve every line), and a value without it invites reading two
    different pools as the same measurement."""
    for p, d in _docs():
        for row in d["products"]:
            assert isinstance(row.get("A"), (int, float)), f"{p.name}: non-numeric A"
            assert row.get("n_lines") is not None, (
                f"{p.name}: product A={row['A']} publishes no n_lines")


def test_quarantined_products_are_not_also_current():
    """A withdrawn product that is still in `products` is not withdrawn at all -- that is
    exactly how the near-UV 8.529 came back after being 'quarantined' into a directory
    the collector still globbed."""
    for p, d in _docs():
        cur = {"|".join(str(r.get(f) or "") for f in KEY_FIELDS) for r in d["products"]}
        for row in d["quarantine"]:
            k = "|".join(str(row.get(f) or "") for f in KEY_FIELDS)
            assert k not in cur, f"{p.name}: {k} is quarantined AND current"
