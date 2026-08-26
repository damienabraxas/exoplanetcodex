"""RYA-1048 — the CONSISTENT tier is WIRED, and deliberately CLOSED.

An entire column of the product grid — CONSISTENT, every band — had no product path.
`gf_empirical` classified lines into `TIER_CONSISTENT` fine; `derive_band_products` had
zero references to it, so the third of the three sub-products (Graded / Deep Graded /
Consistent) could never be produced. These tests pin BOTH halves of the fix:

  * the path exists and is reachable (CLI choice, selector tag, selection branch), and
  * asking for it RAISES, because the saturation gate it depends on is undeclared.

🔴 THE RAISE IS THE FEATURE. RYA-1041 measured that the ensemble curve of growth for solar
Fe I has no resolvable knee, so the scalar `rew_max` is RETIRED and stage 2 gates each line
on its own measured dREW/dA. Until RYA-1043's per-line COG sets `min_d_rew_dA`, a CONSISTENT
run must refuse. A default — and above all Gray's -4.9, which RYA-1041 explicitly refused —
would not fail loudly on one line: it would silently admit or reject an ENTIRE GRID COLUMN
across every band, carrying the authority of a measurement nobody made (RYA-161).
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_empirical import (  # noqa: E402
    TIER_CONSISTENT, ThresholdNotDeclared, Thresholds,
)

SCRIPT = ROOT / "scripts" / "derive_band_products.py"


@pytest.fixture(scope="module")
def src() -> str:
    return SCRIPT.read_text(encoding="utf-8")


def test_the_cli_actually_offers_the_tier(src):
    """The defect was that you could not ASK for it. `--lines-tier` must accept it."""
    tree = ast.parse(src)
    found = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "add_argument"):
            continue
        if not (node.args and getattr(node.args[0], "value", "") == "--lines-tier"):
            continue
        for kw in node.keywords:
            if kw.arg == "choices":
                found = [e.value for e in kw.value.elts]
    assert found, "--lines-tier not found in derive_band_products"
    assert "consistent" in found, f"--lines-tier choices are {found}"


def test_the_product_gets_its_own_name(src):
    """Two tiers sharing a stem is how RYA-1006 destroyed an anchor. CONSISTENT products
    must be distinguishable on disk from GRADED and DEEPGRADED ones."""
    assert '"_CONSISTENT"' in src
    fn = src[src.index("def _selector_tag"):src.index("def conditioning_tag")]
    assert "consistent" in fn and "_CONSISTENT" in fn


def test_the_selection_branch_is_reachable(src):
    """A CLI choice that falls through to `select_lines` is the RYA-933 defect: the run
    asks for one pool and silently measures another."""
    assert "_cand_consistent(" in src
    assert 'getattr(a, "lines_tier", "all") == "consistent"' in src


def test_asking_for_it_RAISES_because_the_gate_is_undeclared():
    """🔴 THE HEART OF THE TICKET. Wired is not open."""
    from scripts.derive_band_products import _cand_consistent  # noqa: E402
    with pytest.raises(ThresholdNotDeclared) as e:
        _cand_consistent(None, lo_A=4200.0, hi_A=6910.0, species="Fe I")
    msg = str(e.value)
    assert "has not been declared" in msg
    # It must name a threshold, so the refusal points at the missing MEASUREMENT.
    assert any(t in msg for t in ("min_depth", "min_d_rew_dA", "max_red_chi2"))


def test_it_raises_BEFORE_touching_any_artifact():
    """The gate is checked before I/O, so the refusal is about the missing threshold and
    not an incidental file error that would mask it."""
    from scripts.derive_band_products import _cand_consistent  # noqa: E402
    with pytest.raises(ThresholdNotDeclared):
        _cand_consistent(None, lo_A=4200.0, hi_A=6910.0, species="Fe I",
                         ew_artifact="/nonexistent/path/that/would/fail/differently.csv")


def test_no_default_saturation_threshold_exists_anywhere():
    """⚠️ The trap this ticket exists to prevent. `rew_max` is RETIRED (RYA-1041): the
    ensemble COG has no resolvable knee, so no scalar can express saturation."""
    th = Thresholds()
    assert not hasattr(th, "rew_max"), (
        "rew_max is back. It was retired because the ensemble curve of growth has no "
        "resolvable knee — a scalar here invites Gray's -4.9 to return.")
    assert th.min_d_rew_dA is None, "the per-line saturation gate must stay undeclared"
    for name in ("min_depth", "max_red_chi2", "rew_min", "consistency_bound_dex"):
        assert getattr(th, name) is None, f"{name} acquired a default"


def test_the_borrowed_constant_is_not_in_the_source(src):
    """RYA-1041 refused Gray-16's -4.9 by name. It must not reappear as a LITERAL in the
    product path.

    ⚠️ Checked on the AST, not the text. A first version of this test grepped the source
    with `#`-comments stripped and failed on its own warning: the docstring above
    `_cand_consistent` NAMES -4.9 precisely so nobody reinstates it. Prose that warns
    against a constant is the opposite of code that uses one, and a text search cannot
    tell them apart.
    """
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            assert abs(float(node.value) - 4.9) > 1e-9, (
                f"the literal 4.9 appears in executable code at line {node.lineno} — "
                f"Gray-16's borrowed saturation constant, refused by RYA-1041")
        for attr in ("id", "attr", "arg"):
            if getattr(node, attr, None) == "rew_max":
                raise AssertionError(
                    f"`rew_max` is referenced in executable code at line "
                    f"{getattr(node, 'lineno', '?')} — the scalar was RETIRED because the "
                    f"ensemble curve of growth has no resolvable knee (RYA-1041)")


def test_the_tier_name_is_single_sourced():
    """The store, the grading tree and the product path must agree on the string."""
    assert TIER_CONSISTENT == "UNGRADED-CONSISTENT"
