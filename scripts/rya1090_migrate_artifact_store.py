#!/usr/bin/env python3
"""
scripts/rya1090_migrate_artifact_store.py
=========================================
RYA-1090 — reconcile the FORKED artifact store into the single canonical root.

WHY
---
`store_root()` used to be `Path(__file__).resolve().parents[2]` — the parent of
whichever worktree the code ran from. Once worktrees moved from
`~/Documents/Exoplanet Codex/` to `~/codex/`, that produced TWO stores. Artifacts
saved to one are invisible to the other, md5 de-duplication no longer holds across
them, and a save-before-clean rescue can report success while the check reads the
other manifest (this happened during RYA-1087).

WHAT THIS DOES
--------------
Moves every artifact in each legacy root into `CANONICAL_STORE_ROOT`, then retires
the legacy manifest so the divergence guard goes quiet.

Three properties this script is built around:

1. **Provenance is preserved, never restamped.** The migrated manifest row keeps the
   ORIGINAL `date`, `source_worktree`, `source_branch` and `source_relpath`. A
   migration that stamps today's date destroys the only record of where an artifact
   came from — which is the entire reason the store exists.
2. **Copy, verify by md5, and only then remove.** No `mv`. The source is deleted
   only after the destination is confirmed byte-identical.
3. **Files with no manifest row are migrated too.** A store can hold files its
   manifest never recorded; those are the easiest to lose and are exactly what a
   manifest-driven migration would silently drop. They are carried over with
   provenance marked `unmanifested`.

Dry-run by default. `--apply` to act.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline.artifact_store import (  # noqa: E402
    KINDS, _MANIFEST_HEADER, _MANIFEST_NAME, CANONICAL_STORE_ROOT,
    LEGACY_STORE_ROOTS, ensure_store, md5sum,
)

RETIRED_NAME = 'ARTIFACT_MANIFEST.migrated-rya1090.csv'


def read_manifest(root: Path) -> list[dict]:
    mp = root / _MANIFEST_NAME
    if not mp.is_file():
        return []
    with open(mp, newline='') as fh:
        return list(csv.DictReader(fh))


def write_manifest(root: Path, rows: list[dict]) -> None:
    with open(root / _MANIFEST_NAME, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=_MANIFEST_HEADER)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in _MANIFEST_HEADER})


def plan_root(legacy: Path, dest_root: Path, dest_rows: list[dict]) -> dict:
    """Classify every artifact under `legacy` against the destination store."""
    # md5 -> stored_name, per kind, on the destination side
    have: dict[tuple[str, str], str] = {}
    for r in dest_rows:
        if r.get('md5'):
            have[(r.get('kind', ''), r['md5'])] = r.get('stored_name', '')
    dest_names: dict[str, set[str]] = {
        k: {p.name for p in (dest_root / k).glob('*')} for k in KINDS
    }

    rows = read_manifest(legacy)
    by_file: dict[Path, dict] = {}
    for r in rows:
        k, name = r.get('kind', ''), r.get('stored_name', '')
        if k in KINDS and name:
            by_file[legacy / k / name] = r

    migrate, dedup, missing, corrupt, unmanifested = [], [], [], [], []

    for src, row in by_file.items():
        if not src.is_file():
            missing.append((src, row))
            continue
        actual = md5sum(src)
        if row.get('md5') and row['md5'] != actual:
            corrupt.append((src, row, actual))
            continue
        kind = row['kind']
        if (kind, actual) in have:
            dedup.append((src, row))
            continue
        migrate.append((src, row, actual, kind))

    # files present in the legacy store that the legacy manifest never recorded
    for k in KINDS:
        d = legacy / k
        if not d.is_dir():
            continue
        for p in sorted(d.glob('*')):
            if not p.is_file() or p in by_file:
                continue
            actual = md5sum(p)
            if (k, actual) in have:
                dedup.append((p, {'kind': k, 'stored_name': p.name, 'md5': actual}))
                continue
            unmanifested.append((p, actual, k))

    return dict(migrate=migrate, dedup=dedup, missing=missing,
                corrupt=corrupt, unmanifested=unmanifested,
                dest_names=dest_names, have=have)


def dest_path_for(dest_root: Path, kind: str, name: str, digest: str,
                  taken: set[str]) -> Path:
    """Non-clobbering destination, same disambiguation rule as save_artifact."""
    p = Path(name)
    if name in taken:
        cand = dest_root / kind / p.name
        if not (cand.is_file() and md5sum(cand) == digest):
            name = f"{p.stem}__{digest[:8]}{p.suffix}"
    return dest_root / kind / name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--apply', action='store_true',
                    help='actually migrate (default: dry-run)')
    ap.add_argument('--keep-sources', action='store_true',
                    help='do not delete legacy files after verified copy')
    a = ap.parse_args()

    dest_root = CANONICAL_STORE_ROOT
    print(f"canonical store : {dest_root}")
    print(f"mode            : {'APPLY' if a.apply else 'DRY-RUN'}\n")

    if a.apply:
        ensure_store(check_divergence=False)

    dest_rows = read_manifest(dest_root)
    print(f"destination holds {len(dest_rows)} manifest rows "
          f"({len({r['md5'] for r in dest_rows if r.get('md5')})} distinct md5)\n")

    total_new = total_dedup = total_fail = 0

    for legacy in LEGACY_STORE_ROOTS:
        legacy = legacy.expanduser()
        if legacy == dest_root or not legacy.is_dir():
            continue
        pl = plan_root(legacy, dest_root, dest_rows)
        n_mig, n_un = len(pl['migrate']), len(pl['unmanifested'])
        print(f"── legacy root: {legacy}")
        print(f"     migrate      : {n_mig}")
        print(f"     unmanifested : {n_un}   (files the legacy manifest never recorded)")
        print(f"     already held : {len(pl['dedup'])}   (same md5 both sides)")
        print(f"     manifest rows with no file : {len(pl['missing'])}")
        print(f"     md5 MISMATCH : {len(pl['corrupt'])}")
        for src, row, actual in pl['corrupt']:
            print(f"        🔴 {src.name}: manifest {row.get('md5','')[:8]} vs actual {actual[:8]}")

        if not a.apply:
            total_new += n_mig + n_un
            total_dedup += len(pl['dedup'])
            continue

        taken = {k: set(v) for k, v in pl['dest_names'].items()}
        new_rows, moved = [], []

        def carry(src: Path, kind: str, digest: str, prov: dict) -> bool:
            dp = dest_path_for(dest_root, kind, src.name, digest, taken[kind])
            dp.parent.mkdir(parents=True, exist_ok=True)
            if not (dp.is_file() and md5sum(dp) == digest):
                shutil.copy2(src, dp)
            if md5sum(dp) != digest:                       # verify BEFORE any delete
                print(f"        🔴 VERIFY FAILED, source kept: {src}")
                return False
            taken[kind].add(dp.name)
            new_rows.append({
                'date': prov.get('date', ''),
                'kind': kind,
                'stored_name': dp.name,
                'md5': digest,
                'bytes': prov.get('bytes', '') or str(src.stat().st_size),
                'source_worktree': prov.get('source_worktree', ''),
                'source_branch': prov.get('source_branch', ''),
                'source_relpath': prov.get('source_relpath', ''),
            })
            moved.append(src)
            return True

        for src, row, digest, kind in pl['migrate']:
            if not carry(src, kind, digest, row):
                total_fail += 1
        for src, digest, kind in pl['unmanifested']:
            prov = {'date': '', 'source_worktree': str(legacy.name),
                    'source_branch': 'unmanifested',
                    'source_relpath': str(src)}
            if not carry(src, kind, digest, prov):
                total_fail += 1

        write_manifest(dest_root, dest_rows + new_rows)
        dest_rows = read_manifest(dest_root)
        print(f"     ✅ carried over {len(new_rows)}; destination now {len(dest_rows)} rows")

        if not a.keep_sources:
            for src in moved:
                src.unlink()
            for src, _row in pl['dedup']:
                if src.is_file():
                    src.unlink()
        lm = legacy / _MANIFEST_NAME
        if lm.is_file():
            lm.rename(legacy / RETIRED_NAME)
            print(f"     ✅ retired legacy manifest → {RETIRED_NAME}")

        total_new += len(new_rows)
        total_dedup += len(pl['dedup'])

    print(f"\nnew in canonical store : {total_new}")
    print(f"deduped (already held) : {total_dedup}")
    print(f"failures               : {total_fail}")
    if not a.apply:
        print("\n(dry-run — re-run with --apply)")
    return 1 if total_fail else 0


if __name__ == '__main__':
    raise SystemExit(main())
