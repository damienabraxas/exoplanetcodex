#!/usr/bin/env python3
"""
RYA-1037 — find every place a spectral line is matched on WAVELENGTH ALONE.
===========================================================================
Keying a line match on wavelength is the most-repeated defect in this pipeline. It has been
fixed one call site at a time — RYA-871, 1033, 1034, 818, 317, 703, 704 — and it keeps
coming back because nothing enforces the rule. This is the audit AND the enforcement: run
with `--check` it exits non-zero, which is what CI reads.

A wavelength is not a transition. Two levels of one species sit inside a tolerance window
routinely, so lambda alone selects a LINE ONLY BY LUCK. The dual key is lambda + excitation
potential, compared by tolerance on BOTH axes.

WHY AST AND NOT GREP. RYA-1033 shipped an AST guard precisely because grep missed ten
sites. `df.merge(other, on=WAVE)` where `WAVE = 'wavelength_air_A'` is invisible to a
pattern; so is a rounded key built two statements away from its use. Grep is kept only as a
cross-check, and the two counts are reported side by side.

WHAT IT FLAGS

  ROUNDED_WAVE_KEY   round()/np.round() applied to a wavelength-ish expression. RYA-1033:
                     Python round() and numpy round() disagree on ...x5, so a rounded
                     wavelength is not even a stable key against itself.
  WAVE_ONLY_MERGE    a pandas merge/join whose keys name a wavelength and no EP.
  WAVE_ONLY_TOL      an abs(a - b) < tol comparison on a wavelength with no EP term
                     anywhere in the enclosing function.
  WAVE_DICT_KEY      a dict/set keyed on a bare wavelength expression.
  ARGMIN_ON_WAVE     argmin/idxmin over a wavelength distance — "nearest wins". Two rows the
                     data cannot separate must not be separated by proximity; this is the
                     first-hit-in-file-order defect wearing a different hat.

WHAT IT DELIBERATELY DOES NOT FLAG, learned from the first pass (225 findings, mostly noise):
rounding a wavelength for DISPLAY or a reported coverage range — the rule is about a rounded
value USED AS A KEY, so the rounded result must reach a key position; band filters
(`w >= lo`, `w <= hi`), which select a window rather than identify a line; `abs(w - 6247.557)`
against a LITERAL, which names one known line rather than joining two tables; boolean-mask
subscripts (`out[np.isfinite(wave)]`); and type annotations.

Usage:
    python3 scripts/audit_line_keys_rya1037.py            # inventory
    python3 scripts/audit_line_keys_rya1037.py --check    # CI: exit 1 on any unwaived find
"""
from __future__ import annotations

import argparse
import ast
import collections
import json
import re
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "audit" / "rya1037"
SCAN_DIRS = ("pipeline", "scripts")

#: The canonical matcher itself is allowed to do all of this — it is the one place that may.
CANONICAL = "pipeline/line_match.py"

#: Findings that are KNOWN and OWNED, not fixed in this pass. Every entry needs a ticket and
#: a reason — a waiver without an owner is just a suppression, and suppression is how this
#: defect class survived seven tickets. Keyed (kind, file) WITH A COUNT: an ELEVENTH finding
#: in a file waived for ten still fails. Line numbers would be brittle (they shift on every
#: edit); a count survives edits and still refuses new occurrences.
WAIVERS = ROOT / "config" / "line_key_waivers.yaml"

#: Names that denote a wavelength. Deliberately broad: a false positive costs a waiver line,
#: a false negative costs another year of this defect.
#: `wl` alone is deliberately NOT here: CRIRES `wlen_id` is an ORDER id, not a wavelength,
#: and matching it produced 90 false positives on the first pass.
_WAVE = re.compile(r"(^|_)(wave|wavelength|lambda|lam_?air|wave_?air|air_A|wl_A|wl_air)"
                   r"($|_|[A-Z0-9])", re.I)
_EP = re.compile(r"(ep|e_?low|elo|excitation|chi|e_?i\b|lower_level)", re.I)


def _is_wave_name(name: str | None) -> bool:
    return bool(name) and bool(_WAVE.search(name))


def _expr_names(node: ast.AST) -> list[str]:
    out = []
    for n in ast.walk(node):
        if isinstance(n, ast.Name):
            out.append(n.id)
        elif isinstance(n, ast.Attribute):
            out.append(n.attr)
        elif isinstance(n, ast.Constant) and isinstance(n.value, str):
            out.append(n.value)
    return out


def _mentions_wave(node: ast.AST) -> bool:
    return any(_is_wave_name(n) for n in _expr_names(node))


def _mentions_ep(node: ast.AST) -> bool:
    return any(_EP.search(n) for n in _expr_names(node) if n)


@dataclass
class Finding:
    kind: str
    file: str
    line: int
    snippet: str

    def key(self) -> tuple:
        return (self.kind, self.file, self.line)


class Visitor(ast.NodeVisitor):
    def __init__(self, rel: str, src: str, key_positioned: set):
        self.rel, self.lines, self.found = rel, src.splitlines(), []
        self._fn_stack: list[ast.AST] = []
        self._key_positioned = key_positioned

    def _snip(self, node) -> str:
        i = getattr(node, "lineno", 1) - 1
        return self.lines[i].strip()[:140] if 0 <= i < len(self.lines) else ""

    def _add(self, kind, node):
        self.found.append(Finding(kind, self.rel, getattr(node, "lineno", 0),
                                  self._snip(node)))

    def visit_FunctionDef(self, node):
        self._fn_stack.append(node)
        self.generic_visit(node)
        self._fn_stack.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_AnnAssign(self, node):
        # `WAVELENGTH_CONVENTION: dict[str, X] = {...}` is a type annotation, not a key.
        if node.value is not None:
            self.visit(node.value)

    def _enclosing_has_ep(self) -> bool:
        return any(_mentions_ep(fn) for fn in self._fn_stack)

    def visit_Call(self, node):
        f = node.func
        fname = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")

        # round(wave, n) / np.round(wave, n) — ONLY when the result reaches a key position.
        if fname == "round" and node.args and _mentions_wave(node.args[0]):
            if node in self._key_positioned:
                self._add("ROUNDED_WAVE_KEY", node)

        # argmin / idxmin over a wavelength distance: "nearest wins"
        if fname in ("argmin", "idxmin") and _mentions_wave(node.func):
            self._add("ARGMIN_ON_WAVE", node)
        if fname in ("argmin", "idxmin") and isinstance(f, ast.Attribute) \
                and _mentions_wave(f.value):
            self._add("ARGMIN_ON_WAVE", node)

        # df.merge(...) / pd.merge(...) keyed on a wavelength
        if fname in ("merge", "join"):
            keys = [kw.value for kw in node.keywords
                    if kw.arg in ("on", "left_on", "right_on", "by")]
            if any(_mentions_wave(k) for k in keys) and not any(_mentions_ep(k) for k in keys):
                self._add("WAVE_ONLY_MERGE", node)

        self.generic_visit(node)

    def visit_Compare(self, node):
        # abs(a - b) <op> tol   with a wavelength on the left and no EP in the function
        left = node.left
        if (isinstance(left, ast.Call)
                and getattr(getattr(left, "func", None), "id", "") == "abs"
                and left.args and isinstance(left.args[0], ast.BinOp)
                and isinstance(left.args[0].op, ast.Sub)
                and _mentions_wave(left.args[0])
                # both operands must be data. `abs(wl - 6247.557)` names ONE known line for
                # a self-check; it is not a join between two tables.
                and not any(isinstance(o, ast.Constant)
                            for o in (left.args[0].left, left.args[0].right))
                and not self._enclosing_has_ep()):
            self._add("WAVE_ONLY_TOL", node)
        self.generic_visit(node)

    def visit_Subscript(self, node):
        # d[wave] used as a lookup/assignment key.
        sl = node.slice
        if not _mentions_wave(sl) or _mentions_ep(sl):
            return self.generic_visit(node)
        # df['wavelength_air_A'] is a COLUMN; dict[str, X] is an ANNOTATION;
        # out[np.isfinite(wave)] is a boolean MASK. None of those is a key.
        if any(isinstance(n, ast.Constant) for n in ast.walk(sl)):
            return self.generic_visit(node)
        if any(isinstance(n, ast.Call) for n in ast.walk(sl)):
            return self.generic_visit(node)
        # A READ (`r[wave_col]`) is ambiguous: the slice may hold a column NAME rather than
        # a wavelength value, and nothing static separates them. A WRITE (`d[wave] = x`) is
        # unambiguously a key being formed. Flag writes only — a check that cannot tell
        # signal from noise gets waived into uselessness.
        if not isinstance(node.ctx, ast.Store):
            return self.generic_visit(node)
        if any(isinstance(n, (ast.Compare, ast.BoolOp)) for n in ast.walk(sl)):
            return self.generic_visit(node)   # boolean mask, not a key
        if _expr_names(sl):
            self._add("WAVE_DICT_KEY", node)
        self.generic_visit(node)


def _key_positions(tree: ast.AST) -> set:
    """Call nodes whose VALUE lands somewhere a key is formed.

    This is the "and then used as a key" half of the rule. Without it the check flags every
    `round(wave, 1)` written for a report — 100 findings on the first pass, almost all of
    them display formatting.
    """
    out: set = set()

    def mark(node):
        for n in ast.walk(node):
            if isinstance(n, ast.Call):
                out.add(n)

    for n in ast.walk(tree):
        # key = ...  /  idx = ...  /  d[...] = ...
        if isinstance(n, ast.Assign):
            names = [t.id for t in n.targets if isinstance(t, ast.Name)]
            if any(re.search(r"(key|idx|index|id)$", x, re.I) for x in names):
                mark(n.value)
            for t in n.targets:
                if isinstance(t, ast.Subscript):
                    mark(t.slice)
        # dict literal keys, set members
        if isinstance(n, ast.Dict):
            for k in n.keys:
                if k is not None:
                    mark(k)
        if isinstance(n, ast.Set):
            for e in n.elts:
                mark(e)
        # x in <dict/set>, .map(), .merge(on=), .groupby()
        if isinstance(n, ast.Compare) and any(isinstance(o, (ast.In, ast.NotIn))
                                              for o in n.ops):
            mark(n.left)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute):
            if n.func.attr in ("map", "groupby", "set_index", "reindex", "isin"):
                for a in n.args:
                    mark(a)
            if n.func.attr in ("merge", "join"):
                for kw in n.keywords:
                    if kw.arg in ("on", "left_on", "right_on", "by"):
                        mark(kw.value)
    return out


def scan(root: Path) -> list[Finding]:
    out = []
    for d in SCAN_DIRS:
        for p in sorted((root / d).rglob("*.py")):
            rel = str(p.relative_to(root))
            if rel == CANONICAL:
                continue
            try:
                tree = ast.parse(p.read_text())
            except SyntaxError:
                continue
            v = Visitor(rel, p.read_text(), _key_positions(tree))
            v.visit(tree)
            out.extend(v.found)
    return out


def grep_count(root: Path) -> int:
    r = subprocess.run(
        ["grep", "-rEn", r"merge\(.*wavelength|round\(.*wave|np\.round\(.*wave",
         "pipeline", "scripts"], cwd=root, capture_output=True, text=True)
    return len([l for l in r.stdout.splitlines() if l.strip()])


def load_waivers(path: Path) -> dict[tuple[str, str], dict]:
    if not path.exists():
        return {}
    import yaml
    doc = yaml.safe_load(path.read_text()) or {}
    out = {}
    for w in doc.get("waivers", []):
        if not w.get("ticket") or not w.get("reason"):
            raise SystemExit(f"waiver without a ticket/reason is a suppression: {w}")
        out[(w["kind"], w["file"])] = w
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="CI mode: exit 1 on any finding")
    ap.add_argument("--root", default=str(ROOT))
    a = ap.parse_args()
    root = Path(a.root)

    found = scan(root)
    waivers = load_waivers(Path(a.root) / "config" / "line_key_waivers.yaml")
    seen_counts = collections.Counter((f.kind, f.file) for f in found)
    over = {k: (n, waivers[k].get("count", 0))
            for k, n in seen_counts.items()
            if k in waivers and n > int(waivers[k].get("count", 0))}
    waived = [f for f in found if (f.kind, f.file) in waivers and (f.kind, f.file) not in over]
    found = [f for f in found if (f.kind, f.file) not in waivers or (f.kind, f.file) in over]
    if over:
        print("\n  🔴 MORE findings than the waiver allows — a new occurrence in a waived "
              "file:")
        for (kind, fil), (now, allowed) in sorted(over.items()):
            print(f"    {kind}  {fil}: {now} now, waiver covers {allowed}")
    by_kind: dict[str, list[Finding]] = {}
    for f in found:
        by_kind.setdefault(f.kind, []).append(f)

    print(f"=== RYA-1037 line-key audit — {len(found)} finding(s) across "
          f"{len({f.file for f in found})} file(s) ===")
    for kind in sorted(by_kind):
        print(f"\n  {kind}  ({len(by_kind[kind])})")
        for f in by_kind[kind][:60]:
            print(f"    {f.file}:{f.line}  {f.snippet}")
        if len(by_kind[kind]) > 60:
            print(f"    ... {len(by_kind[kind]) - 60} more")

    if waived:
        print(f"\n  waived (known, owned, not fixed here): {len(waived)} across "
              f"{len({f.file for f in waived})} file(s)")
        seen = set()
        for f in waived:
            t = waivers[(f.kind, f.file)]["ticket"]
            if (f.kind, f.file) not in seen:
                seen.add((f.kind, f.file))
                print(f"    [{t}] {f.kind}  {f.file}")

    g = grep_count(root)
    print(f"\n  cross-check: the ticket's grep finds {g}; the AST finds {len(found)} "
          f"(RYA-1033's lesson — grep missed ten sites)")

    if not a.check:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "rya1037_line_key_inventory.json").write_text(json.dumps(
            {"n_findings": len(found), "grep_count": g,
             "by_kind": {k: len(v) for k, v in by_kind.items()},
             "n_waived": len(waived),
             "waived": [asdict(f) | {"ticket": waivers[(f.kind, f.file)]["ticket"]}
                        for f in waived],
             "findings": [asdict(f) for f in found]}, indent=2) + "\n")
        print(f"\n[out] {OUT}/rya1037_line_key_inventory.json")
    return 1 if (a.check and found) else 0


if __name__ == "__main__":
    raise SystemExit(main())
