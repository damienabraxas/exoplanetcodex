#!/usr/bin/env python3
"""
scripts/rescue_worktree_artifacts_rya461.py
===========================================
RYA-461 Part B — non-destructive rescue of gitignored artifacts from every git
worktree into the canonical store (pipeline.artifact_store). COPIES ONLY; never
deletes. De-duplicates by md5 and provenance-tags each file (source worktree +
branch + date), then prints a manifest: found / newly-rescued / already-present,
plus any expected grids found MISSING (re-download candidates).

Usage:  python scripts/rescue_worktree_artifacts_rya461.py [--apply]
        (default is a dry run; --apply performs the copies)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import artifact_store as A  # noqa: E402

# (glob relative to a worktree root, store kind). Only GITIGNORED artifacts that are
# valuable + expensive to regenerate: diagnostic plots, normalized-spectrum
# intermediates, and NLTE .grd grids. Tracked files (committed test fixtures,
# committed curation diagnostics, committed atlases) are intentionally excluded.
RESCUE_GLOBS = [
    ('results/plots/**/*.png', 'plots'),
    ('results/plots/**/*.pdf', 'plots'),
    ('results/plots/**/*.svg', 'plots'),
    ('data/processed/*.png', 'plots'),
    ('data/processed/*_normalized.csv', 'data'),
    ('data/nlte_grids/**/*.grd', 'grids'),
]


def _worktrees(repo: Path) -> list[tuple[Path, str]]:
    out = subprocess.run(['git', '-C', str(repo), 'worktree', 'list', '--porcelain'],
                         capture_output=True, text=True).stdout
    trees, cur = [], {}
    for line in out.splitlines():
        if line.startswith('worktree '):
            cur = {'path': line.split(' ', 1)[1]}
        elif line.startswith('branch '):
            cur['branch'] = line.split('/')[-1]
        elif line == '':
            if cur:
                trees.append((Path(cur['path']), cur.get('branch', 'detached')))
                cur = {}
    if cur:
        trees.append((Path(cur['path']), cur.get('branch', 'detached')))
    return trees


def _is_ignored(wt: Path, rel: str) -> bool:
    r = subprocess.run(['git', '-C', str(wt), 'check-ignore', '-q', rel],
                       capture_output=True)
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true', help='perform copies (default: dry run)')
    a = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    trees = _worktrees(repo)
    print(f"Store root: {A.store_root()}")
    print(f"Worktrees:  {len(trees)}\n")

    found = newly = already = 0
    grid_hits = 0
    per_kind: dict[str, int] = {}
    for wt, branch in trees:
        if not wt.exists():
            continue
        for glob, kind in RESCUE_GLOBS:
            for f in sorted(wt.glob(glob)):
                if not f.is_file():
                    continue
                rel = str(f.relative_to(wt))
                if not _is_ignored(wt, rel):
                    continue                       # only gitignored artifacts
                found += 1
                if kind == 'grids':
                    grid_hits += 1
                prov = {'worktree': wt.name, 'branch': branch}
                if not a.apply:
                    print(f"  [dry] {kind:11s} {wt.name}/{rel}  ({A.md5sum(f)[:8]})")
                    continue
                before = {r['md5'] for r in A._read_manifest()
                          if r.get('kind') == kind}
                dest = A.save_artifact(f, kind, provenance=prov)
                after = {r['md5'] for r in A._read_manifest() if r.get('kind') == kind}
                if after - before:
                    newly += 1
                    per_kind[kind] = per_kind.get(kind, 0) + 1
                    print(f"  RESCUED  {kind:11s} {wt.name}/{rel} -> {dest.name}")
                else:
                    already += 1

    print(f"\n── summary ───────────────────────────────")
    print(f"  found (gitignored, at-risk): {found}")
    if a.apply:
        print(f"  newly rescued:               {newly}  {per_kind}")
        print(f"  already in store (dedup):    {already}")
        assert newly + already == found, "rescued + already-present must equal found"
        print(f"  INVARIANT OK: rescued + already-present == found ({newly}+{already}=={found})")
    print(f"  NLTE .grd grids found:       {grid_hits}"
          + ("  <-- NONE on disk: re-download candidates (grid-audit)" if grid_hits == 0 else ""))


if __name__ == '__main__':
    main()
