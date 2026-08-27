#!/usr/bin/env python3
"""RYA-1087 — READ-ONLY worktree prune/archive plan. Deletes nothing.

🔴 THE DECISIVE CHECK IS THE ARTIFACT CHECK, NOT THE MERGE CHECK. A fully merged worktree
can still be the only place a gitignored product exists — RYA-1080 rescued four from
/private/tmp, RYA-1011 recorded that anchor CSVs "survived on luck". `git status` reports
a worktree with orphan artifacts as CLEAN, because they are ignored. So merge status
decides whether COMMITS are safe; only the artifact scan decides whether DATA is.

🔴 BRANCH STATUS IS DECIDED BY CONTENT, NOT BY GREPPING FOR AN RYA NUMBER. The ticket
names RYA-424 as a false negative: no trace of the number on main, yet the work is Done.
For every file the branch touched, this asks whether the branch's blob already matches
main's. All match -> the content is in main however it landed (merge, squash, rebase,
cherry-pick). That is the test; the ticket marks the grep version CRITICAL.

⚠️ SCOPE, STATED HONESTLY: this run happens AFTER a cleanup earlier the same day that
removed 20 worktrees and 317 merged remote branches at Ryan's request. The ticket's
handed-in map (440 branches) predates it. Those removals required merged + zero dirty
TRACKED files -- which is exactly the check this ticket says is insufficient, because it
does not see gitignored artifacts. Verified after the fact: all 93 feed rows across
data/products/solar/*.json still resolve to existing artifacts, so nothing the feeds
depend on was lost. That is evidence, not proof, and it is the reason this audit exists.
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "audit" / "worktree_audit"

#: Paths that hold load-bearing gitignored output. From CONVENTIONS.md's RYA-461 policy
#: (plots, data intermediates, diagnostics, grids, atlases) plus the product surfaces this
#: project has actually lost things from.
LOAD_BEARING = (
    "data/results", "debug/intake", "data/processed", "data/audit",
    "results/plots", "band_products",
)
#: The canonical artifact store — outside any worktree, per CONVENTIONS.md.
STORE = Path.home() / "Documents" / "Exoplanet Codex"


def sh(*args, cwd=None) -> str:
    r = subprocess.run(args, cwd=cwd, capture_output=True, text=True)
    return r.stdout.strip()


def worktrees() -> list[dict]:
    out, cur = [], {}
    for line in sh("git", "worktree", "list", "--porcelain", cwd=ROOT).splitlines():
        if line.startswith("worktree "):
            if cur:
                out.append(cur)
            cur = {"path": line[9:], "branch": None, "detached": False}
        elif line.startswith("branch "):
            cur["branch"] = line[7:].replace("refs/heads/", "")
        elif line == "detached":
            cur["detached"] = True
    if cur:
        out.append(cur)
    return out


def branch_status(branch: str | None) -> str:
    """merged / landed (content in main) / UNMERGED / detached — by CONTENT."""
    if not branch:
        return "detached"
    if sh("git", "merge-base", "--is-ancestor", branch, "origin/main", cwd=ROOT) == "" and \
       subprocess.run(["git", "merge-base", "--is-ancestor", branch, "origin/main"],
                      cwd=ROOT, capture_output=True).returncode == 0:
        return "merged"
    base = sh("git", "merge-base", branch, "origin/main", cwd=ROOT)
    if not base:
        return "UNMERGED(no-base)"
    files = [f for f in sh("git", "diff", "--name-only", base, branch, cwd=ROOT).splitlines() if f]
    if not files:
        return "merged(empty-diff)"
    differing = []
    for f in files:
        b = sh("git", "rev-parse", f"{branch}:{f}", cwd=ROOT)
        m = sh("git", "rev-parse", f"origin/main:{f}", cwd=ROOT)
        if b != m:
            differing.append(f)
    if not differing:
        return "landed(content-in-main)"
    return f"UNMERGED({len(differing)}/{len(files)} files differ)"


def md5(p: Path) -> str:
    h = hashlib.md5()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def store_hashes() -> set[str]:
    """md5s already preserved in the canonical artifact store."""
    man = STORE / "ARTIFACT_MANIFEST.csv"
    got = set()
    if man.exists():
        try:
            for row in csv.DictReader(man.open()):
                for k in ("md5", "checksum", "hash"):
                    if row.get(k):
                        got.add(row[k].strip())
        except Exception:
            pass
    return got


def ignored_artifacts(wt: Path) -> list[Path]:
    """Gitignored files under load-bearing paths — the ones a prune would destroy."""
    names = sh("git", "-C", str(wt), "ls-files", "--others", "--ignored",
               "--exclude-standard").splitlines()
    out = []
    for n in names:
        # ⚠️ MATCH PATH COMPONENTS, NOT SUBSTRINGS. A first version used
        # `any(seg in n ...)`, and "band_products" duly matched
        # `pipeline/__pycache__/band_products.cpython-312.pyc` — byte-compiled caches
        # counted as load-bearing orphans. That inflates the archive list and, worse,
        # keeps worktrees alive to protect a .pyc.
        if n.endswith(".pyc") or "__pycache__" in n.split("/"):
            continue
        parts = n.split("/")
        joined2 = {"/".join(parts[i:i + 2]) for i in range(len(parts) - 1)}
        if not (set(parts) & set(LOAD_BEARING) or joined2 & set(LOAD_BEARING)):
            continue
        p = wt / n
        if p.is_file() and p.stat().st_size > 0:
            out.append(p)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tracked_in_main = set(sh("git", "ls-tree", "-r", "--name-only", "origin/main",
                             cwd=ROOT).splitlines())
    preserved = store_hashes()
    rows = []
    wts = worktrees()
    print(f"RYA-1087 — auditing {len(wts)} worktrees (READ-ONLY)\n")

    for w in wts:
        p = Path(w["path"])
        if not p.exists():
            continue
        size_mb = int(sh("du", "-sm", str(p)).split()[0] or 0)
        dirty = len([l for l in sh("git", "-C", str(p), "status", "--porcelain").splitlines()
                     if not l.startswith("??")])
        status = branch_status(w["branch"])
        protected = str(p).startswith(str(STORE))

        arts = ignored_artifacts(p)
        orphans = []
        for a in arts:
            rel = str(a.relative_to(p))
            if rel in tracked_in_main:            # committed — safe
                continue
            try:
                if md5(a) in preserved:           # already in the artifact store
                    continue
            except Exception:
                pass
            orphans.append(rel)

        if w["branch"] and status.startswith("UNMERGED"):
            verdict, action = "ARCHIVE", "git bundle/tag the branch before removal"
        elif orphans:
            verdict, action = "ARCHIVE", f"save_artifact() {len(orphans)} orphan file(s) first"
        elif dirty:
            verdict, action = "KEEP", "uncommitted tracked changes present"
        elif protected:
            verdict, action = "KEEP", "inside the artifact store root"
        else:
            verdict, action = "PRUNE", "merged/landed, no orphan artifacts, clean"

        rows.append(dict(path=str(p), size_mb=size_mb, branch=w["branch"] or "(detached)",
                         branch_status=status, protected=protected, dirty_tracked=dirty,
                         n_ignored_artifacts=len(arts), n_orphan_artifacts=len(orphans),
                         orphan_examples="; ".join(orphans[:3]), verdict=verdict,
                         save_action=action))

    rows.sort(key=lambda r: (-r["size_mb"]))
    recl = sum(r["size_mb"] for r in rows if r["verdict"] == "PRUNE")
    arch = sum(r["size_mb"] for r in rows if r["verdict"] == "ARCHIVE")
    counts = {v: sum(1 for r in rows if r["verdict"] == v) for v in ("PRUNE", "ARCHIVE", "KEEP")}
    print(f"  reclaimable now (PRUNE): {recl/1024:.1f} GB across {counts['PRUNE']} worktrees")
    print(f"  reclaimable after archiving: {arch/1024:.1f} GB across {counts['ARCHIVE']}")
    print(f"  KEEP: {counts['KEEP']}")
    print(f"  total orphan artifacts found: {sum(r['n_orphan_artifacts'] for r in rows)}")

    with (OUT / "worktree_plan.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    md = [f"# Worktree prune/archive plan — RYA-1087 (READ-ONLY)\n",
          f"**Reclaimable now (PRUNE): {recl/1024:.1f} GB** across {counts['PRUNE']} worktrees. ",
          f"A further **{arch/1024:.1f} GB** in {counts['ARCHIVE']} worktrees needs archiving first. ",
          f"{counts['KEEP']} KEEP.\n",
          "\n> Nothing here was deleted, pruned, or unprotected. The plan is the deliverable.\n",
          "\n## Verdicts\n",
          "| path | GB | branch status | prot | dirty | orphans | verdict | action |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| `{Path(r['path']).name}` | {r['size_mb']/1024:.2f} | {r['branch_status']} | "
                  f"{'Y' if r['protected'] else ''} | {r['dirty_tracked']} | "
                  f"{r['n_orphan_artifacts']} | **{r['verdict']}** | {r['save_action']} |")
    (OUT / "worktree_plan.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"\n  wrote {OUT}/worktree_plan.md + .csv")


if __name__ == "__main__":
    main()
