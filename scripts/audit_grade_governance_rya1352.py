"""Read-only RYA-1352 classification/product regression against an explicit Git base.

No source guessing, line matching, measurement, regeneration or publication. The
small classifier is loaded from its AST so auditing cannot invoke engine loaders.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import platform
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROTECTED = ("pipeline/", "scripts/", "config/", "data/linelists/",
             "data/reference/", "data/products/", "data/catalog/")
SELF = "scripts/audit_grade_governance_rya1352.py"


def git(*args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=ROOT)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def classifier(source: str):
    """Extract only the current public predicate and its declared NIST classes."""
    tree = ast.parse(source)
    keep = [node for node in tree.body if
            (isinstance(node, ast.FunctionDef) and node.name == "is_gf_graded") or
            (isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and
             t.id == "GRADED_NIST_CLASSES" for t in node.targets))]
    if len(keep) != 2:
        raise ValueError("classifier structure changed; review the audit extractor")
    scope = {"pd": pd}
    exec(compile(ast.Module(body=keep, type_ignores=[]), "gf_grades audit", "exec"), scope)
    return scope["is_gf_graded"]


def audit(base: str) -> dict:
    base = git("rev-parse", "--verify", f"{base}^{{commit}}").decode().strip()
    paths = git("ls-tree", "-rz", "--name-only", base).decode().split("\0")
    protected = [p for p in paths if p.startswith(PROTECTED)]
    manifest = {}
    changed = []
    for p in protected:
        before = git("show", f"{base}:{p}")
        path = ROOT / p
        after = path.read_bytes() if path.is_file() else b""
        manifest[p] = digest(before)
        if not path.is_file() or before != after:
            changed.append(p)
    # Include untracked additions so a new runtime module cannot escape the proof.
    additions = git("ls-files", "--others", "--exclude-standard", "-z").decode().split("\0")
    additions += git("diff", "--name-only", "--diff-filter=A", base).decode().splitlines()
    changed += sorted({p for p in additions if p.startswith(PROTECTED) and p != SELF})
    canonical = "data/linelists/canonical_gf.csv"
    rows = pd.read_csv(io.BytesIO(git("show", f"{base}:{canonical}")), low_memory=False)
    before = classifier(git("show", f"{base}:pipeline/gf_grades.py").decode())(rows)
    current_rows = pd.read_csv(ROOT / canonical, low_memory=False)
    after = classifier((ROOT / "pipeline/gf_grades.py").read_text())(current_rows)
    # A different input population is a failed proof, never aligned by row position.
    same_ids = rows.line_id.tolist() == current_rows.line_id.tolist()
    mask_equal = same_ids and before.tolist() == after.tolist()
    accepted = rows[before]
    sigma = pd.to_numeric(rows.gf_sigma_dex, errors="coerce")
    invalid = before & sigma.notna() & ((sigma <= 0) | (sigma == float("inf")))
    products = {}
    for p in protected:
        if p.startswith("data/products/") and p.endswith(".json"):
            doc = json.loads(git("show", f"{base}:{p}"))
            if isinstance(doc, dict) and isinstance(doc.get("products"), list):
                products[p] = dict(sorted(Counter(str(r.get("grade")) for r in doc["products"]).items()))
    return {
        "schema": "rya1352.audit.v1", "base_sha": base,
        "host": platform.system(), "python": platform.python_version(), "pandas": pd.__version__,
        "scope": "Committed repository inputs/products and classifier; no external holdings or live website",
        "protected_file_sha256": manifest, "changed_protected_paths": sorted(set(changed)),
        "canonical_rows": len(rows), "accepted_current_predicate": int(before.sum()),
        "accepted_by_tier": dict(sorted(Counter(accepted.gf_tier.fillna("").astype(str)).items())),
        "accepted_not_exact_LAB": int((accepted.gf_tier != "LAB").sum()),
        "accepted_with_nonpositive_or_infinite_stored_sigma": int(invalid.sum()),
        "classification_unchanged": bool(mask_equal),
        "accepted_line_ids_sha256": digest("\n".join(rows.loc[before, "line_id"]).encode()),
        "product_grade_counts": products,
        "pass": not changed and mask_equal,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    report = audit(args.base)
    with args.output.open("x") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(json.dumps({k: v for k, v in report.items() if k not in
                      ("protected_file_sha256", "product_grade_counts")}, sort_keys=True))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
