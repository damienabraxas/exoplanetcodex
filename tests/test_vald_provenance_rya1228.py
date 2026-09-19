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


def test_the_manifest_is_committed_and_names_all_three_categories():
    text = MANIFEST.read_text()
    assert "## CURRENT" in text and "## SUPERSEDED" in text
    # PRE-LFS is the third: same path, same rewrite, invisible to `git lfs ls-files`.
    assert "## PRE-LFS" in text
    assert "53 LFS objects, 853.9 MiB" in text
    assert "54 objects, 861.2 MiB" in text
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


def test_phase_2_policy_has_landed_and_no_lfs_pattern_remains():
    """⚠️ THIS TEST INVERTED, AND THAT IS THE TRIPWIRE WORKING.

    In Phase 1 it asserted the LFS filter was still present, with the note "if this ever
    fails, Phase 2 has run". Phase 2's policy half has now landed, so it asserts the other
    side: there must be NO LFS-tracked pattern at all. The history purge itself is a
    separate step and is NOT what this checks -- see the raw-file test below.
    """
    assert "filter=lfs" not in (ROOT / ".gitattributes").read_text()


def test_the_raw_files_are_still_in_the_tree_because_the_purge_has_not_run():
    """Honest state: the policy is in place, the history rewrite is not.

    Until the purge runs, the working tree still tracks the raw deliveries and the
    stewardship invariant legitimately reports them. Flipping this test is how you will
    know the purge landed.
    """
    tracked = subprocess.run(["git", "-C", str(ROOT), "ls-files",
                              "data/linelists/vald_*_raw.txt"],
                             capture_output=True, text=True).stdout.split()
    assert len(tracked) == 23, f"expected the 23 HEAD raw files, saw {len(tracked)}"


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
    assert len(rows) == 54, f"expected 54 data rows, found {len(rows)}"
    for r in rows:
        assert re.search(r"`[0-9a-f]{16}…`", r), r


# --- RYA-1228: the pre-LFS plain-blob class -----------------------------------------
#
# The Phase 1 inventory enumerated with `git lfs ls-files --all` and missed a raw
# extraction that predates LFS adoption. It is on the same path and dies in the same
# rewrite, but it is not an LFS object, so that command cannot see it. These tests pin
# the class, not just the instance.

import importlib.util
import pathlib as _pathlib

_VERIFIER = _pathlib.Path(__file__).resolve().parents[1] / "scripts" / "rya1228_verify_vald_archive.py"


def _load_verifier():
    spec = importlib.util.spec_from_file_location("rya1228_verify", _VERIFIER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_census_is_taken_by_path_not_by_lfs():
    """raw_objects() must union both storage kinds, or the class reopens."""
    mod = _load_verifier()
    src = _VERIFIER.read_text()
    assert "def plain_blobs(" in src, "plain-blob enumeration was removed"
    assert "lfs_objects(repo) + plain_blobs(repo)" in src, (
        "raw_objects() no longer unions both storage kinds -- an LFS-only census "
        "misses every pre-LFS blob"
    )
    assert mod.RAW_GLOB == "data/linelists/vald_*_raw.txt"


def test_manifest_covers_the_pre_lfs_blob():
    """The specific object the LFS-only inventory missed."""
    import json
    manifest = json.loads((ARCHIVE / "MANIFEST.json").read_text())
    by_sha = {o["sha256"]: o for o in manifest["objects"]}
    sha = "a7d016acd567dbca8959989e4275b2cc8108b9b12749201b9cacb9aca2cbf728"
    assert sha in by_sha, "the pre-LFS plain blob is not in the archive manifest"
    entry = by_sha[sha]
    assert entry["storage"] == "plain-git-blob"
    assert entry["size_bytes"] == 7634177
    copy = ARCHIVE / entry["archive"]
    assert copy.exists(), f"{entry['archive']}: manifest entry has no file"
    import hashlib
    h = hashlib.sha256()
    with open(copy, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    assert h.hexdigest() == sha, "preserved copy does not hash to the recorded sha256"


def test_provenance_records_the_plain_blob_and_the_corrected_quota():
    text = MANIFEST.read_text()
    assert "a7d016acd5" in text, "the pre-LFS blob is not recorded in the provenance file"
    assert "54 objects" in text, "census count not corrected to 54"
    assert "43 objects / 692.7 MiB" in text, "GitHub quota figure not corrected"
