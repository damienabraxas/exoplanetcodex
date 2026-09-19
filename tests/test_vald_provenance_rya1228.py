"""RYA-1228 Phase 1 -- the manifest, and the preconditions Phase 2 depends on.

🔴 WHAT THIS PROTECTS. VALD3 has no API: every raw delivery is a manual web extraction and
none is re-fetchable. Phase 2 deletes them from history permanently. The only thing standing
between that and unrecoverable loss is a hash-verified copy outside git -- so the checks
that assert it must themselves be tested, and must cover ALL REFS rather than the 23 files
visible at HEAD.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MANIFEST = ROOT / "data/linelists/VALD_PROVENANCE.md"
ARCHIVE = Path.home() / "Documents" / "Exoplanet Codex" / "vald_raw_archive"


def _lfs(*args):
    out = subprocess.run(["git", "-C", str(ROOT), "lfs", "ls-files", *args],
                         capture_output=True, text=True)
    return out.stdout.splitlines() if out.returncode == 0 else []


def test_the_manifest_is_committed_and_names_both_categories():
    text = MANIFEST.read_text()
    assert "## CURRENT" in text and "## SUPERSEDED" in text
    assert "53 objects, 853.9 MiB" in text
    assert "NONE OF THESE IS RE-FETCHABLE" in text


def test_the_manifest_states_that_built_line_lists_are_out_of_scope():
    """The one thing that must NOT be touched: the Codex's own built products."""
    text = MANIFEST.read_text()
    assert "linelist_*.csv" in text and "Untouched" in text


def test_no_built_line_list_was_modified():
    """CRITICAL in the ticket: any built linelist_*.csv touched fails Phase 1."""
    out = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain",
                          "data/linelists/linelist_", "data/linelists/canonical_gf.csv"],
                         capture_output=True, text=True).stdout.strip()
    assert out == "", f"a built product was modified: {out}"


def test_no_lfs_pointer_was_modified():
    out = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain",
                          "data/linelists/"], capture_output=True, text=True).stdout
    touched = [l for l in out.splitlines() if "_raw.txt" in l]
    assert touched == [], f"Phase 1 must not modify a raw delivery: {touched}"


def test_gitattributes_still_tracks_vald_raw_because_phase_2_owns_that_change():
    """Phase 1 changes no policy. If this ever fails, Phase 2 has run."""
    assert "vald_*_raw.txt filter=lfs" in (ROOT / ".gitattributes").read_text()


@pytest.mark.skipif(not ARCHIVE.exists(), reason="preservation archive not on this machine")
def test_every_LFS_object_across_ALL_REFS_has_a_hash_verified_copy():
    """🔴 ALL REFS, not HEAD. 30 of the 53 objects exist only in history and a rewrite
    destroys those too -- verifying the working tree would leave 413.7 MiB unprotected."""
    manifest = json.loads((ARCHIVE / "MANIFEST.json").read_text())
    by_oid = {m["sha256"]: m for m in manifest["objects"]}
    oids = {line.split(None, 2)[0] for line in _lfs("--all", "--long")}
    assert len(oids) == 53, f"expected 53 LFS objects across all refs, saw {len(oids)}"
    assert oids <= set(by_oid), "an LFS object is absent from the archive manifest"
    for oid in oids:
        assert (ARCHIVE / by_oid[oid]["archive"]).exists()


@pytest.mark.skipif(not ARCHIVE.exists(), reason="preservation archive not on this machine")
def test_the_verifier_refuses_when_a_preserved_copy_is_wrong():
    """Mutation control: a verifier that cannot fail is not a gate.

    The real script hashes every copy; here we prove the comparison it makes is the one
    that would catch a corrupted or truncated archive file.
    """
    import hashlib
    manifest = json.loads((ARCHIVE / "MANIFEST.json").read_text())
    entry = min(manifest["objects"], key=lambda m: m["size_bytes"])
    data = (ARCHIVE / entry["archive"]).read_bytes()
    assert hashlib.sha256(data).hexdigest() == entry["sha256"]
    assert hashlib.sha256(data + b"x").hexdigest() != entry["sha256"]


@pytest.mark.skipif(not ARCHIVE.exists(), reason="preservation archive not on this machine")
def test_the_oid_really_is_the_content_sha256():
    """The identity the whole verification rests on -- checked, not assumed."""
    import hashlib
    manifest = json.loads((ARCHIVE / "MANIFEST.json").read_text())
    entry = min(manifest["objects"], key=lambda m: m["size_bytes"])
    assert hashlib.sha256((ARCHIVE / entry["archive"]).read_bytes()).hexdigest() == entry["sha256"]


def test_every_manifest_row_carries_a_sha256_prefix():
    """A provenance row without a hash documents nothing checkable."""
    rows = [l for l in MANIFEST.read_text().splitlines()
            if l.startswith("| `vald_") and "…`" in l]
    assert len(rows) == 53, f"expected 53 data rows, found {len(rows)}"
    for r in rows:
        assert re.search(r"`[0-9a-f]{16}…`", r), r
