"""Every `--engine-b-deck` must be able to REACH its fit (RYA-1208).

🔴 RYA-1206 bound `_unlabelled` inside `if _nlte:` and read it at the bottom of the
Engine-B block, which EVERY deck reaches. So every non-NLTE Engine-B run died on
`UnboundLocalError: cannot access local variable '_unlabelled'` -- `ts-lte` (the DEFAULT,
and the production Engine B for all 27 species), `gerber-1d-lte` and `gerber-mean3d-lte`.

It survived a merge because the only runs after it either used an NLTE deck (which binds
the name) or passed `--skip-engine-b` (which skips the block). RYA-1203's NIR
re-derivation was the second kind, so that ticket's suite went green over a broken
default route.

These tests read the SOURCE rather than executing a synthesis: reaching the crash needs an
atlas, a deck and ~20 minutes, which is why no test had this. A name bound on one branch
and read on all of them is a static property, and static is where it should have been
caught.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "scripts" / "derive_band_products.py"


import functools


@functools.lru_cache(maxsize=1)
def _tree_for(text: str) -> ast.Module:
    return ast.parse(text)


def _synthesis_route() -> ast.FunctionDef:
    """`synthesis_route` from ONE cached parse, so node identity is comparable."""
    tree = _tree_for(SRC.read_text())
    return next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "synthesis_route")


def _engine_b_block() -> ast.If:
    """The `if not a.skip_engine_b:` block inside `synthesis_route`."""
    fn = _synthesis_route()
    for node in ast.walk(fn):
        if (isinstance(node, ast.If) and isinstance(node.test, ast.UnaryOp)
                and isinstance(node.test.op, ast.Not)
                and "skip_engine_b" in ast.dump(node.test)):
            return node
    raise AssertionError("could not locate the `if not a.skip_engine_b:` block")


def _assigned_here(stmt) -> set[str]:
    if isinstance(stmt, (ast.Assign, ast.AnnAssign)):
        tgts = stmt.targets if isinstance(stmt, ast.Assign) else [stmt.target]
        return {t.id for t in tgts if isinstance(t, ast.Name)}
    return set()


def _names_bound_unconditionally(block) -> set[str]:
    """Names bound on EVERY path through `block.body`.

    🔴 AN if/elif/else CHAIN THAT BINDS IN EVERY BRANCH IS UNCONDITIONAL. A first version
    looked only at top-level statements and flagged `eb_treatment`, which the deck chain
    assigns in all three branches -- a false positive that would have needed an exemption
    list, and an exemption list is how a guard quietly stops guarding.
    """
    out: set[str] = set()
    for stmt in block.body:
        out |= _assigned_here(stmt)
        if isinstance(stmt, ast.If):
            # collect the branches of the whole if/elif/.../else chain
            branches, node = [], stmt
            while True:
                branches.append(node.body)
                if len(node.orelse) == 1 and isinstance(node.orelse[0], ast.If):
                    node = node.orelse[0]
                    continue
                if node.orelse:
                    branches.append(node.orelse)
                else:
                    branches = []          # no else -> some path binds nothing
                break
            if branches:
                per = [set().union(*(_assigned_here(x) for x in b)) if b else set()
                       for b in branches]
                out |= set.intersection(*per) if per else set()
    return out


def _names_read_unconditionally(block: ast.If) -> set[str]:
    """Names LOADED at the block's own top level, outside any nested branch."""
    out: set[str] = set()
    for stmt in block.body:
        if isinstance(stmt, (ast.If, ast.Try, ast.For, ast.While, ast.With)):
            # a read inside a nested branch is guarded by that branch's own condition,
            # except for the `if <name>:` TEST itself, which always evaluates
            if isinstance(stmt, ast.If):
                for n in ast.walk(stmt.test):
                    if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load):
                        out.add(n.id)
            continue
        for n in ast.walk(stmt):
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load):
                out.add(n.id)
    return out


def test_every_name_read_on_all_engine_b_paths_is_bound_on_all_of_them():
    """The invariant, stated generally rather than pinned to `_unlabelled`: a local read
    unconditionally inside the Engine-B block must be bound unconditionally too."""
    block = _engine_b_block()
    bound = _names_bound_unconditionally(block)
    read = _names_read_unconditionally(block)
    # 🔴 SCOPED TO NAMES THIS BLOCK ACTUALLY BINDS. A first version flagged every read
    # name starting with "_" and caught `_fit_lines`, which is a closure defined earlier
    # in `synthesis_route` -- bound in the ENCLOSING scope, so reading it here is fine.
    # The real invariant is narrower and is the one that bites: a name the block binds on
    # SOME branch and reads on ALL of them.
    block_local = set()
    for stmt in block.body:
        for n in ast.walk(stmt):
            if isinstance(n, (ast.Assign, ast.AnnAssign)):
                tgts = n.targets if isinstance(n, ast.Assign) else [n.target]
                block_local |= {t.id for t in tgts if isinstance(t, ast.Name)}
    # Names bound EARLIER in `synthesis_route`, before this block, are already live when
    # the block runs -- `cand` is the candidate pool, built well above. Reading one here
    # is not a path hazard.
    #
    # 🔴 THE FUNCTION AND THE BLOCK MUST COME FROM THE **SAME PARSE**. A first version
    # re-parsed the source here and compared `stmt is block` against a node from the
    # other tree, so the identity check never matched, the loop never broke, and `before`
    # collected the WHOLE function -- including `_unlabelled` itself. The test then
    # passed on the very source it was written to catch. Verified by re-running it
    # against the pre-fix file, which is the only reason this was found.
    fn = _synthesis_route()
    before: set[str] = set()
    for stmt in fn.body:
        if stmt is block:
            break
        for n in ast.walk(stmt):
            before |= _assigned_here(n)
    else:                                    # pragma: no cover - block must be reachable
        raise AssertionError("the Engine-B block was not found in the same parse")
    missing = sorted(n for n in (read & block_local) - before if n not in bound)
    assert missing == [], (
        "read on every Engine-B path but bound on only some of them, so a deck that "
        f"skips the binding branch dies with UnboundLocalError: {missing}")


def test_unlabelled_specifically_is_bound_before_the_deck_branch():
    """The regression itself, pinned by name so the fix cannot be quietly reverted."""
    block = _engine_b_block()
    assert "_unlabelled" in _names_bound_unconditionally(block), (
        "_unlabelled must be initialised at the top of the Engine-B block: it is read "
        "after the fit on EVERY deck, and RYA-1206 bound it only under `if _nlte:`")


@pytest.mark.parametrize("deck", ["ts-lte", "gerber-1d-lte", "gerber-mean3d-lte"])
def test_the_non_nlte_decks_are_still_offered(deck):
    """A control: if these decks ever stop being选able the test above goes vacuous."""
    tree = ast.parse(SRC.read_text())
    choices = {c.value for n in ast.walk(tree)
               if isinstance(n, ast.Call)
               and getattr(n.func, "attr", "") == "add_argument"
               and any(getattr(k, "arg", "") == "choices" for k in n.keywords)
               for k in n.keywords if k.arg == "choices"
               for c in getattr(k.value, "elts", []) if isinstance(c, ast.Constant)}
    assert deck in choices, f"{deck} is no longer an --engine-b-deck choice"
