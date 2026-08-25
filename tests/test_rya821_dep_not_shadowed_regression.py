"""The hoisted departure value must not be SHADOWED inside the fit closure (RYA-821).

`_solve` hoists Gerber departures once for decks that resolve a single A(X) (Fe), and
recomputes them per trial for decks with a real abundance axis (Al, 31 nodes). PR #383
wrote the per-trial result back to the SAME name the hoisted value lives in:

    _dep = None                      # outer
    ...
    def chi2(a_x):
        if _dep_per_abundance:
            _dep = _gn.for_node(...)  # <-- makes _dep a LOCAL of chi2 for the WHOLE body
        if _dep is not None:          # <-- therefore reads the local, unbound when the
                                      #     branch above did not run

Python decides local-vs-free at COMPILE time from the presence of an assignment, so the
branch not running does not fall back to the enclosing scope -- it raises
`UnboundLocalError`. Every element whose deck has no abundance axis therefore lost every
fit, INCLUDING plain 1D-LTE with no NLTE deck at all, while Al -- the one element that
takes the assigning branch -- kept working. Measured: 58/58 near-UV Fe I lines excluded
with "cannot access local variable '_dep'", n=0.

Pinned structurally rather than by running a synthesis, so it fails in milliseconds at
the fault instead of after an engine call. The invariant is the one that generalises:
a name the closure reads from the enclosing scope must not also be assigned in it.
"""
from __future__ import annotations

import ast
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "pipeline" / "abundances_derive.py"


def _closures(tree):
    """Every nested function, paired with the enclosing function that defines it."""
    for outer in ast.walk(tree):
        if not isinstance(outer, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for node in ast.iter_child_nodes(outer):
            for inner in ast.walk(node):
                if isinstance(inner, ast.FunctionDef) and inner is not outer:
                    yield outer, inner


def _assigned(fn) -> set[str]:
    out = set()
    for n in ast.walk(fn):
        if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store):
            out.add(n.id)
        elif isinstance(n, (ast.Global, ast.Nonlocal)):
            out -= set(n.names)          # explicitly rebound: not a shadow
    return out


def test_chi2_does_not_shadow_the_hoisted_departures():
    tree = ast.parse(SRC.read_text())
    checked = 0
    for outer, inner in _closures(tree):
        if inner.name != "chi2":
            continue
        checked += 1
        assigned = _assigned(inner)
        assert "_dep" not in assigned, (
            f"{inner.name}() (in {outer.name}) assigns `_dep`, which makes the hoisted "
            "outer `_dep` an unbound local on every path that does not take the "
            "assigning branch -- that is the RYA-821 regression, and it killed plain "
            "1D-LTE too. Write the per-trial value to a different name."
        )
    assert checked, "no chi2 closure found — the guard is vacuous; update the anchor"


def test_the_broken_pattern_really_does_raise():
    """Guard the guard: prove the shadowing actually raises, so the test above is not
    pinning a rule that Python does not enforce."""
    ns: dict = {}
    exec(
        "def make(dep_outer, per_abundance):\n"
        "    _dep = dep_outer\n"
        "    _per = per_abundance\n"
        "    def chi2(a):\n"
        "        if _per:\n"
        "            _dep = 'per-trial'\n"
        "        return _dep\n"
        "    return chi2\n",
        ns,
    )
    assert ns["make"](None, True)(1.0) == "per-trial"      # Al path: fine
    for outer in (None, "hoisted"):                        # Fe / plain-LTE paths
        try:
            ns["make"](outer, False)(1.0)
        except UnboundLocalError:
            continue
        raise AssertionError("shadowing did not raise — the invariant is not what we think")
