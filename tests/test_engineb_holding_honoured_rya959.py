"""Engine-B must honour `--holding` — RYA-959.

RYA-933/934 added `--holding` so a caller can name ONE product instead of taking the
instrument's first covering candidate, because `solar_harps` and
`solar_harps_molecfit_corrected` are the one pair this repo must never collapse: they
differ precisely by whether tellurics were removed (RYA-911/927/931).

The synthesis route honoured it. **Engine-B did not.** A RYA-959 run launched with
`--holding solar_harps_molecfit_corrected` printed

    [holding] harps -> solar_harps (pre_normalised=True)
    [provenance] 1 holding(s) served this product: ['solar_harps'] — matches --instrument harps

and derived its whole Engine-B product on the UNCORRECTED sibling, while the 1D-LTE and
Engine-A legs — which read the EW artifact, and that artifact WAS measured on the
corrected holding — were correct. One product, two holdings, no complaint.

🔴 AND THE EXISTING PROVENANCE GUARD COULD NOT HAVE CAUGHT IT. It built its expectation
from the same `select_holding` call the lines were loaded with, so it proved the leg was
self-consistent, never that it honoured the request. A guard whose expected value is
derived from the thing it is checking cannot fail — that is the shape being pinned here.
"""
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DERIVER = ROOT / "scripts" / "derive_band_products.py"


def _calls_named(tree: ast.AST, func_name: str) -> list[ast.Call]:
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            f = node.func
            name = getattr(f, "id", None) or getattr(f, "attr", None)
            if name == func_name:
                out.append(node)
    return out


def _passes_holding(call: ast.Call) -> bool:
    return any(kw.arg == "holding" for kw in call.keywords)


def test_every_holding_selection_forwards_the_requested_holding():
    """Asserted over EVERY call site, not the one that was broken (RYA-870).

    Pinning the specific line that failed would let the next new call site reintroduce it
    silently — which is exactly how this survived RYA-933/934 in the first place.
    """
    tree = ast.parse(DERIVER.read_text(encoding="utf-8"))
    for fn in ("select_holding", "load_window_ex"):
        calls = _calls_named(tree, fn)
        assert calls, f"no {fn} calls found — has the deriver been restructured?"
        missing = [c.lineno for c in calls if not _passes_holding(c)]
        assert not missing, (
            f"{fn}() is called at {DERIVER.name} line(s) {missing} without forwarding "
            f"`holding=`. A leg that drops it silently falls back to the instrument's "
            f"first covering candidate, which for HARPS is the UNCORRECTED sibling of "
            f"the telluric-corrected product the caller asked for (RYA-931).")


def test_the_provenance_guard_compares_against_the_request_not_itself():
    """The guard must raise when a named holding is not the one that served the lines."""
    src = DERIVER.read_text(encoding="utf-8")
    assert "HOLDING IGNORED" in src, (
        "Engine-B must refuse to emit a product whose lines were measured on a holding "
        "other than the one `--holding` named. Without that check the leg can only ever "
        "verify it agreed with itself.")
    guard = src[src.index("HOLDING IGNORED") - 400: src.index("HOLDING IGNORED") + 200]
    assert "a.holding" in guard and "_b_holdings" in guard, (
        "the check must compare the REQUESTED holding against the holdings that actually "
        "served the lines")


def test_the_two_harps_holdings_are_distinct_products():
    """The premise. If these ever collapse to one id the guard above is pointless."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from measure_band_ew import _INSTRUMENT_HOLDINGS, PRE_NORMALISED

    ids = [h.holding_id for h in _INSTRUMENT_HOLDINGS["harps"]]
    assert "solar_harps" in ids and "solar_harps_molecfit_corrected" in ids
    assert len(set(ids)) == len(ids)
    # Both are pre-normalised, so `pre_normalised` cannot be used to tell them apart —
    # which is why the holding ID has to be carried explicitly rather than inferred.
    assert PRE_NORMALISED["solar_harps"] is True
    assert PRE_NORMALISED["solar_harps_molecfit_corrected"] is True
