#!/usr/bin/env python3
"""
RYA-1228 -- verify the preserved VALD archive against the repository's LFS objects.

🔴 THIS IS PHASE 2's PRECONDITION. Nothing may be removed from history until every LFS
object has a hash-verified copy outside git. "I copied them" is not that check; this is.
Run it immediately before any rewrite, and again after.

⚠️ IT VERIFIES ALL REFS, NOT HEAD. 23 of the 53 objects are at HEAD; the other 30 exist
only in history, and a history rewrite destroys those too. Verifying only the working tree
would leave 413.7 MiB -- including every tau Boo extraction -- unprotected.

The LFS oid IS the sha256 of the content, so `archive file hash == oid` is an identity
check rather than a re-derivation. VALD3 has no API: none of these is re-fetchable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys

ARCHIVE = pathlib.Path.home() / "Documents" / "Exoplanet Codex" / "vald_raw_archive"


def sha256(path: pathlib.Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def lfs_objects(repo: pathlib.Path) -> list[dict]:
    """Every LFS object across ALL refs, with its oid and path."""
    out = subprocess.run(["git", "-C", str(repo), "lfs", "ls-files", "--all", "--long"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(f"git lfs ls-files failed: {out.stderr.strip()}")
    seen, rows = set(), []
    for line in out.stdout.splitlines():
        oid, _mark, path = line.split(None, 2)
        if oid in seen:
            continue
        seen.add(oid)
        rows.append({"oid": oid, "path": path})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archive", default=str(ARCHIVE))
    ap.add_argument("--repo", default=".")
    args = ap.parse_args()
    arch = pathlib.Path(args.archive).expanduser()
    if not arch.exists():
        print(f"ARCHIVE ABSENT: {arch}", file=sys.stderr)
        return 2

    manifest = json.loads((arch / "MANIFEST.json").read_text())
    by_oid = {m["sha256"]: m for m in manifest["objects"]}
    objects = lfs_objects(pathlib.Path(args.repo).resolve())

    missing, mismatched, ok = [], [], 0
    for obj in objects:
        entry = by_oid.get(obj["oid"])
        if entry is None:
            missing.append(f"{obj['path']} ({obj['oid'][:10]}): not in the archive manifest")
            continue
        copy = arch / entry["archive"]
        if not copy.exists():
            missing.append(f"{entry['archive']}: manifest entry has no file")
            continue
        got = sha256(copy)
        if got != obj["oid"]:
            mismatched.append(f"{entry['archive']}: sha256 {got[:16]} != oid {obj['oid'][:16]}")
            continue
        ok += 1

    print("RYA-1228 archive verification")
    print(f"  archive         : {arch}")
    print(f"  LFS objects     : {len(objects)} (all refs)")
    print(f"  hash-verified   : {ok}")
    print(f"  missing         : {len(missing)}")
    print(f"  hash MISMATCH   : {len(mismatched)}")
    for m in missing + mismatched:
        print(f"    {m}")
    if missing or mismatched:
        print("\nREFUSE: Phase 2 must not proceed. Every object needs a hash-verified copy.")
        return 1
    print("\nOK: every LFS object across all refs has a hash-verified preserved copy.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
