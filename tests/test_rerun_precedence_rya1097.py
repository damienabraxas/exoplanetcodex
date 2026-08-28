"""
RYA-1097 — "new takes precedent UNLESS it is bad", and identity is what makes that decidable.

🔴 PRECEDENCE WAS RECENCY-ONLY. `publish_product` archived the existing record and installed
the new one whenever a value changed, with no reference to whether the new one is
publishable. So a bad re-run -- `sigma_syst` null, a pre-continuum-fix artifact, an
uninterpretable error bar -- would SUPERSEDE a good live product, and being newer is exactly
what makes it look authoritative. The gate now decides: a re-run that fails goes to
`quarantine[]` with its reasons and the existing record stays live.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline import product_eligibility as pe

ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data" / "products" / "solar" / "Fe.json"


@pytest.fixture(scope="module")
def doc():
    return json.loads(FEED.read_text())


def _live(doc):
    p = next((x for x in doc["products"] if x.get("sigma_syst") is not None), None)
    if p is None:
        pytest.skip("no eligible live product to build on")
    return p


# ── A. identity ──────────────────────────────────────────────────────────────────

def test_identity_is_the_nine_fields_INCLUDING_route_and_selector(doc):
    """RYA-1097's spec lists (element, ion, band, instrument, holding, tier, treatment,
    selector). `route` is kept as well, and that is not pedantry: IAG VIS GRADED ENGINE-A
    is 7.481 on n=4 by profile fit and 7.484 on n=37 by synthesis. A key without `route`
    collides two different measurements of two different pools, and the store's own
    overwrite guard caught exactly that on its first backfill."""
    for f in ("element", "ion", "band", "instrument", "holding", "tier", "selector",
              "treatment", "route"):
        assert f in pe.KEY_FIELDS, f
    # ⚠️ Perturb AWAY FROM whatever the sample already carries. A first version set
    # selector="DEEPGRADED" on a product whose selector was already DEEPGRADED and
    # asserted the key changed -- a test that could only pass by accident of which
    # product it picked.
    a = _live(doc)
    for field in ("route", "selector"):
        b = dict(a, **{field: f"NOT-{a.get(field)}"})
        assert pe.key_of(a) != pe.key_of(b), f"{field} must be part of the identity"


def test_every_live_product_is_unique_on_that_identity(doc):
    assert pe.duplicate_live_cells(doc) == {}


# ── D. gated re-run precedence ───────────────────────────────────────────────────

def test_a_BAD_rerun_does_not_supersede_a_good_product():
    """🔴 THE RULE, checked by structure. The publish path must consult the gate BEFORE
    replacing a live record, and must route a failing product into quarantine."""
    import ast
    src = (ROOT / "scripts" / "publish_product.py").read_text()
    tree = ast.parse(src)
    main = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    body = ast.get_source_segment(src, main)
    assert "pe.evaluate(" in body, (
        "publish must evaluate each pending row against the eligibility gate")
    # The gate call has to come BEFORE the supersede/replace, or it decides nothing.
    assert body.index("pe.evaluate(") < body.index('doc["products"][by_key[k]] = row'), (
        "the gate must be consulted before the live record is replaced")
    assert "quarantine_codes" in body and "quarantine_reason" in body, (
        "a refused re-run must be recorded with its reason codes, not dropped")


def test_a_refused_rerun_is_MOVED_not_discarded():
    """A re-run that fails the gate is still evidence: it is what the pipeline produced.
    Dropping it would leave no record that the attempt happened (RYA-833/711)."""
    import ast
    src = (ROOT / "scripts" / "publish_product.py").read_text()
    main = next(n for n in ast.walk(ast.parse(src))
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    body = ast.get_source_segment(src, main)
    assert 'doc.setdefault("quarantine", []).append(stamped)' in body


def test_a_GOOD_rerun_still_supersedes_and_archives(doc):
    """The other half. Gating precedence must not freeze the feed: an ELIGIBLE re-run has
    to replace its predecessor, and the predecessor has to be archived, not deleted."""
    import ast
    src = (ROOT / "scripts" / "publish_product.py").read_text()
    main = next(n for n in ast.walk(ast.parse(src))
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    body = ast.get_source_segment(src, main)
    assert 'doc.setdefault("archive", []).append(old)' in body
    assert 'doc["products"][by_key[k]] = row' in body
    # And it happened for real: RYA-1096 republished eight cells this way.
    archived = [p for p in doc.get("archive") or [] if p.get("superseded_reason")]
    assert archived, "nothing has ever been superseded -- the path is untested by use"


def test_the_gate_would_refuse_a_syst_incomplete_rerun(doc):
    """The specific bad re-run the rule is written for: a product with no systematic."""
    bad = dict(_live(doc), sigma_syst=None)
    assert "SYST_INCOMPLETE" in {r.code for r in pe.evaluate(bad)}
