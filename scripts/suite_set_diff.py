#!/usr/bin/env python3
"""RYA-1138 — the merge gate: diff a branch's failure set against a baseline's.

    python3 scripts/suite_set_diff.py --baseline-ref origin/main [--check]
    python3 scripts/suite_set_diff.py --baseline-root ~/codex/somewhere [--check]

"Run the suite and confirm the failure set is unchanged from main" IS the merge gate
in this repo. Nothing else asks whether a branch broke something: the suite has stood
at 8 red on clean main for weeks, so an absolute pass count is meaningless and only the
DIFF carries information.

🔴 THAT MAKES THE BASELINE A LOAD-BEARING MEASUREMENT, AND IT HAS BEEN COMPUTED BY HAND
EVERY TIME. On 2026-08-30 a baseline worktree was created under a scratch directory;
`ISPEC_DIR` resolves as `ROOT.parent / 'ispec'`, so from there `ispec` did not exist,
25 test modules failed to IMPORT, and the baseline reported 25 collection errors and no
failures. The branch had eleven real failures. Diffed against that baseline, ELEVEN
REGRESSIONS WOULD HAVE READ AS AN IMPROVEMENT -- and the output would have looked
entirely normal.

This harness removes the two hand steps that made that possible:

  1. IT CREATES THE BASELINE ITSELF, under the canonical worktree parent, so the
     baseline's location is never a caller's choice. Siting it was the error; a tool
     that sites it correctly cannot make that error (RYA-1138 item 1).
  2. IT REFUSES TO EMIT A VERDICT unless both checkouts have the same capability
     fingerprint (RYA-1138 item 2, RYA-869's loud-fail). A baseline that cannot import
     what the branch imports is not a baseline.

⚠️ IT COMPARES SKIPS, NOT ONLY FAILURES. A skipped test is in neither the pass set nor
the fail set, so a baseline that skips 60 tests the branch runs produces a clean diff
while measuring nothing -- the RYA-1035 lesson ("a skipped test is not a passing
test"), which a failure-set diff structurally cannot see.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import canonical_checkout as cc                        # noqa: E402

OUT = ROOT / "data" / "results" / "rya1138" / "suite_set_diff.json"


def run_suite(root: Path, target: str, python: str) -> dict:
    """Run the suite and return its failure set, skip set and counts."""
    r = subprocess.run([python, "-m", "pytest", target, "-q", "-p", "no:randomly",
                        "-rs"],
                       cwd=root, capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    failed = sorted({m.group(1) for m in
                     re.finditer(r"^(?:FAILED|ERROR)\s+(\S+)", out, re.M)})
    skipped = sorted({m.group(1) for m in
                      re.finditer(r"^SKIPPED\s+\[\d+\]\s+(\S+?):", out, re.M)})
    tail = re.search(r"^(\d+ failed.*|\d+ passed.*)$", out, re.M)
    return {"root": str(root), "failed": failed, "skipped": skipped,
            "n_failed": len(failed), "n_skipped_sites": len(skipped),
            "summary": tail.group(1) if tail else "(no pytest summary line)",
            "returncode": r.returncode, "raw_tail": out[-4000:]}


def ensure_baseline(ref: str) -> Path:
    """Create (or reuse) a baseline worktree UNDER THE CANONICAL PARENT.

    🔴 The location is not a parameter. Letting the caller choose it is precisely the
    degree of freedom that produced the false baseline, and a convention that has to be
    remembered at the moment of use is the thing this ticket exists to replace.
    """
    parent = cc.CANONICAL_WORKTREE_PARENT
    parent.mkdir(parents=True, exist_ok=True)
    slug = re.sub(r"[^a-zA-Z0-9]+", "-", ref).strip("-").lower()
    dest = parent / f"_setdiff-baseline-{slug}"
    if dest.exists():
        subprocess.run(["git", "worktree", "remove", "--force", str(dest)],
                       cwd=ROOT, capture_output=True, text=True)
    subprocess.run(["git", "fetch", "-q", "origin"], cwd=ROOT, check=False)
    r = subprocess.run(["git", "worktree", "add", "-q", "--detach", str(dest), ref],
                       cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"could not create baseline worktree at {dest}:\n{r.stderr}")
    return dest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--branch-root", type=Path, default=ROOT)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--baseline-ref", help="git ref; a worktree is created for it "
                                          "under the canonical worktree parent")
    g.add_argument("--baseline-root", type=Path,
                   help="an EXISTING checkout to use as the baseline (still probed "
                        "and still refused if its fingerprint differs)")
    ap.add_argument("--target", default="tests/")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--keep", action="store_true",
                    help="keep the created baseline worktree")
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if the branch adds a failure or a skip")
    a = ap.parse_args()

    created = None
    if a.baseline_ref:
        created = ensure_baseline(a.baseline_ref)
        baseline_root = created
    else:
        baseline_root = a.baseline_root.expanduser().resolve()

    branch_root = a.branch_root.expanduser().resolve()

    print("=" * 88)
    print("RYA-1138 — suite set-diff, with the baseline's competence checked FIRST")
    print("=" * 88)

    b_caps = cc.probe(branch_root, python=a.python)
    l_caps = cc.probe(baseline_root, python=a.python)
    b_caps.collected_tests, b_caps.collection_errors = cc.collect(
        branch_root, a.target, python=a.python)
    l_caps.collected_tests, l_caps.collection_errors = cc.collect(
        baseline_root, a.target, python=a.python)

    print("\nBRANCH");   print(cc.describe(b_caps))
    print("\nBASELINE"); print(cc.describe(l_caps))

    try:
        cc.assert_comparable(b_caps, l_caps)
    except cc.NonComparableCheckoutsError as e:
        print(f"\n🔴 {e}")
        return 2
    print("\n✅ fingerprints match — the two sides are comparable, so a difference in "
          "their failure sets is attributable to the branch.")

    branch = run_suite(branch_root, a.target, a.python)
    base = run_suite(baseline_root, a.target, a.python)

    new_fail = [t for t in branch["failed"] if t not in set(base["failed"])]
    fixed = [t for t in base["failed"] if t not in set(branch["failed"])]
    new_skip = [t for t in branch["skipped"] if t not in set(base["skipped"])]

    print("\n" + "-" * 88)
    print(f"  baseline : {base['summary']}")
    print(f"  branch   : {branch['summary']}")
    print("-" * 88)
    for label, items, mark in (("NEW FAILURE", new_fail, "🔴"),
                               ("newly passing", fixed, "✅"),
                               ("NEWLY SKIPPED", new_skip, "⚠️")):
        for t in items:
            print(f"  {mark} {label}: {t}")
    if not new_fail and not new_skip:
        print("  ✅ failure set IDENTICAL and nothing newly skipped — nothing moved.")

    doc = {"ticket": "RYA-1138", "branch": cc.as_dict(b_caps),
           "baseline": cc.as_dict(l_caps), "baseline_ref": a.baseline_ref,
           "new_failures": new_fail, "newly_passing": fixed,
           "newly_skipped": new_skip,
           "branch_summary": branch["summary"], "baseline_summary": base["summary"]}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")

    if created and not a.keep:
        subprocess.run(["git", "worktree", "remove", "--force", str(created)],
                       cwd=ROOT, capture_output=True, text=True)

    if a.check and (new_fail or new_skip):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
