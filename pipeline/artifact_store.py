"""
pipeline/artifact_store.py
==========================
RYA-461 — the canonical local artifact store + the save-before-clean helper.

WHY
---
Gitignored artifacts produced inside a git worktree (diagnostic plots, downloaded
atlases, NLTE .grd grids, normalized-spectrum intermediates) live ONLY in that
worktree. When the worktree is cleaned/removed they are LOST — and for the big
NLTE grids that turns "we have the grid" into "re-download 26 GB." This module is
the durable home for anything gitignored that may be needed later, plus a one-call
helper so future runs copy their artifacts to the store as part of the run instead
of remembering by hand.

THE STANDING CONVENTION (see docs/CONVENTIONS.md)
-------------------------------------------------
Any gitignored artifact a pipeline/diagnostic run produces — plots, diagnostic
CSV/JSON, downloaded atlases, NLTE grids, normalized-spectrum intermediates — must
be copied to the matching store subfolder AS PART OF THE RUN, never left only in a
worktree. Call `save_artifact(path, kind)` right after you write the file.

STORE LAYOUT
------------
The store lives OUTSIDE any worktree, in the directory that CONTAINS the worktrees
(default: the parent of the repo root, i.e. ~/Documents/Exoplanet Codex/), or at
$CODEX_ARTIFACT_STORE if set. Subfolders, one per `kind`:
    plots/        diagnostic plots (png/pdf/svg)
    data/         normalized-spectrum + other reusable data intermediates
    diagnostics/  diagnostic CSV/JSON outputs (audit/proof/verdict tables)
    grids/        NLTE .grd / model grids (large, external, expensive to re-fetch)
    atlases/      downloaded reference atlases (Kitt Peak, CALSPEC, IRTF, ...)

Every saved file is de-duplicated by md5 and recorded in store/ARTIFACT_MANIFEST.csv
with its source worktree/branch + date (provenance).
"""
from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import os
import shutil
import subprocess
from pathlib import Path

KINDS = ('plots', 'data', 'diagnostics', 'grids', 'atlases')

_MANIFEST_NAME = 'ARTIFACT_MANIFEST.csv'
_MANIFEST_HEADER = ['date', 'kind', 'stored_name', 'md5', 'bytes',
                    'source_worktree', 'source_branch', 'source_relpath']


def store_root() -> Path:
    """The canonical store root: $CODEX_ARTIFACT_STORE, else the directory that
    contains this repo's worktree (the worktrees' shared parent)."""
    env = os.environ.get('CODEX_ARTIFACT_STORE')
    if env:
        return Path(env).expanduser()
    # this file is <store_root>/<worktree>/pipeline/artifact_store.py
    return Path(__file__).resolve().parents[2]


def ensure_store() -> Path:
    """Create the store root + all `KINDS` subfolders. Returns the root."""
    root = store_root()
    for k in KINDS:
        (root / k).mkdir(parents=True, exist_ok=True)
    return root


def md5sum(path) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _git(repo: Path, *args) -> str:
    try:
        out = subprocess.run(['git', '-C', str(repo), *args],
                             capture_output=True, text=True, timeout=15)
        return out.stdout.strip()
    except Exception:
        return ''


def git_provenance(path) -> dict:
    """Best-effort (worktree, branch) for the file's containing repo."""
    p = Path(path).resolve()
    repo = _git(p.parent, 'rev-parse', '--show-toplevel')
    worktree = Path(repo).name if repo else p.parent.name
    branch = _git(p.parent, 'rev-parse', '--abbrev-ref', 'HEAD') if repo else ''
    return {'worktree': worktree, 'branch': branch}


def _manifest_path() -> Path:
    return store_root() / _MANIFEST_NAME


def _read_manifest() -> list[dict]:
    mp = _manifest_path()
    if not mp.exists():
        return []
    with open(mp, newline='') as fh:
        return list(csv.DictReader(fh))


def _append_manifest(row: dict) -> None:
    mp = _manifest_path()
    new = not mp.exists()
    with open(mp, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=_MANIFEST_HEADER)
        if new:
            w.writeheader()
        w.writerow(row)


def save_artifact(path, kind: str, provenance: dict | None = None,
                  dedup: bool = True) -> Path:
    """Copy a gitignored artifact into the canonical store and record its provenance.

    Parameters
    ----------
    path        : the artifact to preserve (must exist).
    kind        : one of KINDS — selects the store subfolder.
    provenance  : optional {'worktree':..., 'branch':...}; auto-detected from git if None.
    dedup       : if True (default) and an identical-md5 file already lives in the
                  store for this kind, do NOT copy again — return the existing path.

    Returns the path in the store. Never deletes or mutates the source. Same-name
    different-content files are disambiguated with a short md5 suffix (no clobber).
    """
    src = Path(path)
    if not src.is_file():
        raise FileNotFoundError(f"artifact not found: {src}")
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}, got {kind!r}")
    ensure_store()
    digest = md5sum(src)
    nbytes = src.stat().st_size
    prov = provenance or git_provenance(src)
    dest_dir = store_root() / kind

    if dedup:
        for row in _read_manifest():
            if row.get('md5') == digest and row.get('kind') == kind:
                existing = dest_dir / row['stored_name']
                if existing.exists():
                    return existing            # already preserved — idempotent

    dest = dest_dir / src.name
    if dest.exists() and md5sum(dest) != digest:
        dest = dest_dir / f"{src.stem}__{digest[:8]}{src.suffix}"
    if not (dest.exists() and md5sum(dest) == digest):
        shutil.copy2(src, dest)

    _append_manifest({
        'date': _dt.date.today().isoformat(),
        'kind': kind,
        'stored_name': dest.name,
        'md5': digest,
        'bytes': nbytes,
        'source_worktree': prov.get('worktree', ''),
        'source_branch': prov.get('branch', ''),
        'source_relpath': str(src),
    })
    return dest
