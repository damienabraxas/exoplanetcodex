"""
tests/test_feed_repo_reconciliation_rya1080.py — RYA-1080
=========================================================
The feed publishes; `data/results/band_products/` is what people grep. When the two
disagree, a check of the tree answers a question about the feed with an answer about the
disk — which cost two wrong diagnoses in one session, the second ("Fe II has never been
re-run") nearly triggering a from-scratch re-run of finished work.

The positive control is the point of this file: a guard that has never been shown to fail
is not evidence of anything.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import feed_repo_reconciliation as frr  # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
BLOCKED = ROOT / "data" / "results" / "rya1080" / "rya1080_blocked.csv"

#: The ten products whose SOURCE file changed after the feed ingested it, so there is no
#: file left that matches the recorded sha256. Deliberately not committed: a file whose
#: checksum has moved is not the file the feed measured, whatever its numbers say.
#: Pinned as a SET so an eleventh is a hard failure rather than a bigger number.
KNOWN_BLOCKED = {
    ("MISSING_COPIED_TO", "Fe.json", "Fe", "I", "VIS", "1D-LTE"),
    ("MISSING_COPIED_TO", "Fe.json", "Fe", "I", "VIS", "ENGINE-A"),
    ("MISSING_COPIED_TO", "Fe.json", "Fe", "I", "near-UV", "1D-LTE"),
    ("MISSING_COPIED_TO", "Fe.json", "Fe", "I", "near-UV", "ENGINE-A"),
    ("MISSING_COPIED_TO", "Fe.json", "Fe", "I", "red-optical", "1D-LTE"),
    ("MISSING_COPIED_TO", "Fe.json", "Fe", "I", "red-optical", "ENGINE-A"),
}


@pytest.fixture(scope="module")
def findings():
    return frr.check()


# ── the guard's own contract ──────────────────────────────────────────────────

@pytest.mark.xfail(strict=True, reason=
    "KNOWN, REPORTED, NOT SILENT (RYA-1080). Ten live Fe.json products have no committed\n"
    "     artifact because their source file's sha256 no longer matches the one the feed\n"
    "     recorded — see data/results/rya1080/rya1080_blocked.csv. Nine differ in bytes\n"
    "     only; one also differs in sigma_stat (feed 0.0217, file 0.0218 — a round-half tie\n"
    "     on 0.02175). Committing any of them under the feed's hash would record a file the\n"
    "     feed never measured, which is the drift this ticket exists to close.\n"
    "     strict=True: when they are reconciled this flips to an unexpected PASS and fails,\n"
    "     so the marker cannot outlive the exception. Ryan's call — re-ingest or investigate.")
def test_every_live_product_has_a_committed_checksum_matching_artifact(findings):
    blocking = [f for f in findings if f.blocking]
    assert not blocking, frr.report(findings)


def test_the_blocking_set_is_exactly_the_known_ten(findings):
    """The xfail above says 'some drift is known'. THIS says WHICH — so an eleventh
    unreconciled product is a hard failure today, not a slightly larger number nobody
    reads."""
    got = {f.key() for f in findings if f.blocking}
    assert got == KNOWN_BLOCKED, (
        f"the blocking set has changed.\n  new: {sorted(got - KNOWN_BLOCKED)}\n"
        f"  gone: {sorted(KNOWN_BLOCKED - got)}")


def test_the_guard_exit_code_is_nonzero_while_drift_exists():
    """A guard that reports drift and exits 0 is decoration. CI reads the exit code."""
    r = subprocess.run([sys.executable, "pipeline/feed_repo_reconciliation.py"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 1


# ── POSITIVE CONTROLS — required by the ticket ────────────────────────────────

def _sandbox(tmp_path: Path) -> Path:
    """A throwaway repo-shaped tree with one reconciled product, so a control can break
    something without touching the real tree."""
    doc = json.loads(FEED.read_text())
    good = next(p for p in doc["products"]
                if p["provenance"].get("copied_to")
                and (ROOT / p["provenance"]["copied_to"]).exists())
    rel = good["provenance"]["copied_to"]
    (tmp_path / Path(rel).parent).mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / rel, tmp_path / rel)
    feed_dir = tmp_path / "data" / "products" / "solar"
    feed_dir.mkdir(parents=True, exist_ok=True)
    (feed_dir / "Fe.json").write_text(json.dumps({"products": [good]}, indent=2))
    # git-track it, since the guard requires tracking and band_products is gitignored
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "-f", rel], cwd=tmp_path, check=True)
    return tmp_path


def test_control_the_sandbox_passes_before_it_is_broken(tmp_path):
    """The control's own control. If the untouched sandbox already failed, breaking it
    would prove nothing."""
    s = _sandbox(tmp_path)
    assert not [f for f in frr.check(root=s) if f.blocking]


def test_control_a_renamed_artifact_makes_the_guard_fail(tmp_path):
    """POSITIVE CONTROL 1 — the exact drift this ticket exists to catch: the feed says the
    product is there, the file is not."""
    s = _sandbox(tmp_path)
    rel = json.loads((s / "data/products/solar/Fe.json").read_text())["products"][0][
        "provenance"]["copied_to"]
    (s / rel).rename(s / (rel + ".renamed"))
    kinds = {f.kind for f in frr.check(root=s) if f.blocking}
    assert "ARTIFACT_ABSENT" in kinds, kinds


def test_control_a_modified_artifact_makes_the_guard_fail(tmp_path):
    """POSITIVE CONTROL 2 — the file is present and is not the one the feed measured. This
    is the case that A-value equality would have waved through."""
    s = _sandbox(tmp_path)
    rel = json.loads((s / "data/products/solar/Fe.json").read_text())["products"][0][
        "provenance"]["copied_to"]
    (s / rel).write_bytes((s / rel).read_bytes() + b"\n")
    kinds = {f.kind for f in frr.check(root=s) if f.blocking}
    assert "SHA_MISMATCH" in kinds, kinds


def test_control_a_copied_to_outside_the_repo_is_not_accepted(tmp_path):
    """POSITIVE CONTROL 3 — 🔴 the trap that made the exposure look like 67 instead of 75.
    Eight rows carried a non-null copied_to pointing into /private/tmp/g3d/ and
    /private/tmp/sirius_orphans/, which READS as reconciled and is not committed at all."""
    s = _sandbox(tmp_path)
    f = s / "data/products/solar/Fe.json"
    doc = json.loads(f.read_text())
    doc["products"][0]["provenance"]["copied_to"] = "/private/tmp/g3d/somewhere.csv"
    f.write_text(json.dumps(doc, indent=2))
    kinds = {x.kind for x in frr.check(root=s) if x.blocking}
    assert "COPIED_TO_OUTSIDE_REPO" in kinds, kinds


def test_control_an_untracked_artifact_makes_the_guard_fail(tmp_path):
    """POSITIVE CONTROL 4 — band_products is gitignored (.gitignore:87), so a file being
    present on disk is NOT the same as it being committed. Present-but-untracked must
    fail, or a local run would satisfy the guard for everyone else."""
    s = _sandbox(tmp_path)
    rel = json.loads((s / "data/products/solar/Fe.json").read_text())["products"][0][
        "provenance"]["copied_to"]
    subprocess.run(["git", "rm", "--cached", "-q", rel], cwd=s, check=True)
    kinds = {f.kind for f in frr.check(root=s) if f.blocking}
    assert "ARTIFACT_UNTRACKED" in kinds, kinds


# ── the ticket's CRITICAL conditions ──────────────────────────────────────────

def test_no_published_value_was_edited_to_reconcile():
    """🔴 CRITICAL: reconcile the artifacts, not the numbers. Every published field of every
    product must be byte-for-byte what origin/main published; `copied_to` is the only key
    this ticket may write."""
    r = subprocess.run(["git", "show", "origin/main:data/products/solar/Fe.json"],
                       cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip("origin/main not available")
    before = {}
    for p in json.loads(r.stdout)["products"]:
        pr = dict(p["provenance"])
        pr.pop("copied_to", None)
        before[(p["element"], p["ion"], p["band"], p["instrument"], p.get("tier"),
                p["treatment"], pr["sha256"])] = {**{k: v for k, v in p.items()
                                                    if k != "provenance"}, "_prov": pr}
    after = {}
    for p in json.loads(FEED.read_text())["products"]:
        pr = dict(p["provenance"])
        pr.pop("copied_to", None)
        after[(p["element"], p["ion"], p["band"], p["instrument"], p.get("tier"),
               p["treatment"], pr["sha256"])] = {**{k: v for k, v in p.items()
                                                   if k != "provenance"}, "_prov": pr}
    assert before == after, "a published field changed — only copied_to may be written"


def test_regenerability_gaps_are_recorded_not_silent(findings):
    """Committing the artifact is the floor; regenerability is the net (RYA-1011). Where
    the holding is Mac-local the product cannot be rebuilt on the committed runner, and
    that must appear in the guard's own output rather than being discovered later."""
    gaps = [f for f in findings if f.kind == "REGENERABILITY_GAP"]
    assert len(gaps) >= 50
    assert all("RYA-1011" in f.detail for f in gaps)


def test_the_blocked_ledger_records_both_checksums():
    """A blocked product must say what the feed recorded AND what is on disk, or the next
    person cannot tell a cosmetic rewrite from a changed measurement."""
    if not BLOCKED.exists():
        pytest.skip("blocked ledger absent")
    import pandas as pd
    d = pd.read_csv(BLOCKED)
    assert len(d) == 10
    assert d.feed_sha256.notna().all() and d.actual_sha256.notna().all()
    assert (d.feed_sha256 != d.actual_sha256).all()
