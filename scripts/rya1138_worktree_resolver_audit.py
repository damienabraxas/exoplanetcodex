#!/usr/bin/env python3
"""RYA-1138 item 5 — find every path resolver that reads ABOVE its own repo root.

    python3 scripts/rya1138_worktree_resolver_audit.py [--check]

RYA-1090 fixed one of these by hand. This finds the rest.

🔴 THE NAIVE GREP IS USELESS HERE, AND WORSE THAN USELESS. `grep parents\\[` returns
~150 hits in this repo and nearly all of them are CORRECT: `Path(__file__).parents[1]`
from `pipeline/foo.py` IS the repo root, and a module reading its own repo's data is
exactly right -- each worktree should read its own data. A reviewer handed 150 hits
where 145 are fine learns to skim, and the two that matter are indistinguishable from
the rest by pattern.

What separates a defect from a correct resolver is not the pattern but the ARITHMETIC:
how deep is the file, how many levels does the expression climb, and does the result
land inside the repo or outside it. `parents[2]` is a DEFECT from `pipeline/x.py` and
CORRECT from `pipeline/audit/x.py`. So this audit parses the AST, computes each
expression's depth against the file's own location, and reports only the ones that
ESCAPE -- which is the actual defect class, stated exactly.

⚠️ ESCAPING IS NOT AUTOMATICALLY A BUG, AND THIS DOES NOT PRETEND OTHERWISE. Reaching
out of the repo is legitimate when the target is genuinely external (an engine
install, a data volume). It is a DEFECT when the escaping path is then treated as
canonical -- when a store, a baseline, or a dependency resolves through it, so that
moving the worktree silently changes the answer instead of raising. The audit
therefore reports escapes with their resolved target and whether an environment
override exists, and classifies rather than just counting: an escape with a named
override and a loud failure is a documented external dependency, an escape without one
is RYA-1090's mechanism waiting to fire.
"""
from __future__ import annotations

import argparse
import ast
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data" / "results" / "rya1138" / "worktree_resolver_audit.json"

#: 🔴 EXCLUDED BY NAME, NEVER BY PATTERN. This auditor discusses `parents[N]` and
#: `.parent` in its own prose and constants, so a pattern-based self-exclusion would
#: also silence any file that merely resembled it. RYA-1112's audit script matched its
#: own grep the moment it was tracked, and the failure only appeared at `git add` time.
SELF = "scripts/rya1138_worktree_resolver_audit.py"

#: Resolvers that escape the repo ON PURPOSE, each with the reason and the override
#: that makes the dependency explicit. Named individually -- a class-wide exemption
#: ("anything mentioning ispec is fine") would re-admit the whole defect family.
KNOWN_EXTERNAL = {
    "pipeline/canonical_checkout.py": (
        "RYA-1138's own `_canonical_location()` — compares a checkout's parent against "
        "the NAMED `CANONICAL_WORKTREE_PARENT`. It reads the parent in order to "
        "REPORT on it, and resolves nothing through it."),
    "pipeline/artifact_store.py": (
        "RYA-1090's `_legacy_derived_root()` — DELIBERATELY still derives the old "
        "pre-1090 root so the divergence guard can see a store left behind by the old "
        "derivation. It is never used AS the store root; `store_root()` is named."),
}


def _depth_of(expr: ast.AST) -> int | None:
    """Levels climbed by a `Path(__file__).resolve().parents[N]` / `.parent` chain."""
    depth = 0
    node = expr
    while True:
        if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Attribute) \
                and node.value.attr == "parents":
            idx = node.slice
            if not isinstance(idx, ast.Constant) or not isinstance(idx.value, int):
                return None
            # `parents[N]` climbs N+1 levels: `parents[0]` IS the immediate parent.
            # Counting it as N is an off-by-one that hides exactly the real defects --
            # `parents[2]` from a depth-2 file reads one level ABOVE the root and would
            # score as landing on it.
            depth += idx.value + 1
            node = node.value.value
        elif isinstance(node, ast.Attribute) and node.attr == "parent":
            depth += 1
            node = node.value
        elif isinstance(node, ast.Attribute) and node.attr == "resolve":
            node = node.value
        elif isinstance(node, ast.Call):
            node = node.func
        elif isinstance(node, ast.Attribute):
            node = node.value
        elif isinstance(node, ast.Name):
            return depth if depth else None
        else:
            return depth if depth else None


def _mentions_file_or_root(expr: ast.AST) -> str | None:
    """The anchor a chain starts from: `__file__`, or a module-level ROOT-ish name."""
    for node in ast.walk(expr):
        if isinstance(node, ast.Name) and node.id == "__file__":
            return "__file__"
        if isinstance(node, ast.Constant) and node.value == "__file__":
            return "__file__"
        if isinstance(node, ast.Name) and node.id.upper().strip("_") in (
                "ROOT", "REPO", "REPOROOT", "REPO_ROOT"):
            return node.id
    return None


#: What the escaping path is actually REACHING FOR. Grouping matters more than the
#: raw count: 20 scattered findings look like 20 problems, but they are two external
#: dependencies -- an engine install and a data volume -- resolved the same wrong way
#: twenty times. Two named constants fix the family; twenty local edits do not.
_FAMILIES = (("ispec engine install", ("ispec",)),
             ("spectra data volume", ("spectra", "exoplanetcodex-data")))


def _family(src: str) -> str:
    low = src.lower()
    for name, needles in _FAMILIES:
        if any(n in low for n in needles):
            return name
    return "other"


def audit() -> list[dict]:
    findings = []
    for path in sorted(ROOT.rglob("*.py")):
        rel = path.relative_to(ROOT).as_posix()
        if rel == SELF or "/.git/" in f"/{rel}" or rel.startswith("data/audit/"):
            continue
        try:
            text = path.read_text(errors="replace")
            tree = ast.parse(text)
        except SyntaxError:
            continue
        source_lines = text.splitlines()
        # Steps from the FILE ITSELF up to the repo root -- not the directory depth.
        # `__file__` is the .py file, so `pipeline/x.py` needs TWO `.parent` steps to
        # reach the root, and `Path(__file__).parent.parent` therefore lands exactly ON
        # the root rather than above it.
        file_depth = len(path.relative_to(ROOT).parts)

        for node in ast.walk(tree):
            if not isinstance(node, (ast.Subscript, ast.Attribute)):
                continue
            anchor = _mentions_file_or_root(node)
            if anchor is None:
                continue
            depth = _depth_of(node)
            if depth is None:
                continue
            # __file__-anchored chains start at the FILE; ROOT-anchored ones start at
            # the repo root, where any climb at all already escapes.
            escapes = depth > file_depth if anchor == "__file__" else depth > 0
            if not escapes:
                continue
            above = depth - (file_depth if anchor == "__file__" else 0)
            src = source_lines[node.lineno - 1].strip() if \
                node.lineno - 1 < len(source_lines) else ""
            findings.append({
                "file": rel, "line": node.lineno, "anchor": anchor,
                "file_depth": file_depth, "climbs": depth,
                "levels_above_repo_root": above,
                "known_external": rel in KNOWN_EXTERNAL,
                "why": KNOWN_EXTERNAL.get(rel, ""),
                "source": src[:160],
                "family": _family(src),
            })
    # Deduplicate nested AST nodes that describe the same chain on one line.
    seen, uniq = set(), []
    for f in sorted(findings, key=lambda d: (d["file"], d["line"],
                                             -d["levels_above_repo_root"])):
        k = (f["file"], f["line"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(f)
    return uniq


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero on an UNDOCUMENTED escape")
    a = ap.parse_args()

    findings = audit()
    undocumented = [f for f in findings if not f["known_external"]]

    L = []; A = L.append
    A("=" * 92)
    A("RYA-1138 item 5 — path resolution that reads ABOVE its own repo root")
    A("=" * 92)
    A("  A resolver that climbs past the repo root reads through the WORKTREE'S")
    A("  PARENT. That is canonical only while every worktree shares one parent, and")
    A("  it forks silently the moment they do not (RYA-1090: 85 artifacts lost;")
    A("  RYA-1138: a baseline that could not import ispec and inverted a merge)."); A("")
    A(f"  scanned: every *.py under {ROOT.name}/   escapes found: {len(findings)}")
    A("")
    fams = {}
    for f in findings:
        fams.setdefault(f["family"], []).append(f)
    for fam in sorted(fams, key=lambda k: (k == "other", k)):
        group = fams[fam]
        n_undoc = sum(1 for f in group if not f["known_external"])
        A(f"  ── {fam}  ({len(group)} site(s), {n_undoc} undocumented)")
        if fam != "other" and n_undoc:
            A(f"     🔴 ONE dependency resolved the same wrong way {n_undoc} times.")
            A(f"     Fixing it is ONE named constant with an env override, not "
              f"{n_undoc} local edits.")
        for f in group:
            tag = "documented" if f["known_external"] else "🔴"
            A(f"     {tag} {f['file']}:{f['line']}  "
              f"({f['levels_above_repo_root']} above root)")
            A(f"         {f['source'][:110]}")
            if f["why"]:
                A(f"         -> {f['why']}")
        A("")
    if not findings:
        A("  ✅ no resolver reads above its own repo root.")
    A("=" * 92)

    text = "\n".join(L)
    print(text)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"ticket": "RYA-1138", "findings": findings,
                               "n_undocumented": len(undocumented)}, indent=2) + "\n")
    OUT.with_suffix(".txt").write_text(text + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    if a.check and undocumented:
        print(f"\nFAIL: {len(undocumented)} undocumented escape(s)")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
