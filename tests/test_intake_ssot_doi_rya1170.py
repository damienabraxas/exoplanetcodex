"""RYA-1170 — an intake must POINT at the reference SSOT, never copy a value out of it.

`Amarsi2019_Table1` carried DOI 10.1051/0004-6361/201936179 as the cited source of the
ENTIRE C I / O I atomic census. Crossref resolves it to Curran & Moss, "Quasi-stellar
object redshift estimates from optical, near-infrared, and ultraviolet colours", A&A 629
(2019) — a different paper, different authors, different volume. `bibliography.csv` held
the correct DOI the whole time, so the intake had made a second copy and the copy drifted:
RYA-355's single-source defect, not a typo.

The row's own CITATION string ("Amarsi, Nissen & Skuladottir 2019, A&A 630 A104") was
right and disagreed with its own DOI, which is what makes this detectable at all.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import check_intake_ssot_rya1170 as G  # noqa: E402

CNO_BIB = ROOT / "data/audit/rya1136_cno_intake/source_bibliography.csv"

#: Resolved LIVE at api.crossref.org on 2026-09-07 while implementing this ticket. Pinned
#: rather than re-fetched: CI has no network, and a test that silently skips when offline
#: is not a check. The wrong DOI is kept beside it so the pair stays legible.
CROSSREF_VERIFIED = {
    "10.1051/0004-6361/201936265": ("Carbon, oxygen, and iron abundances in disk and "
                                    "halo stars", "Amarsi", 630),
    "10.1051/0004-6361/201936179": ("Quasi-stellar object redshift estimates from "
                                    "optical, near-infrared, and ultraviolet colours",
                                    "Curran", 629),
}
WRONG_DOI = "10.1051/0004-6361/201936179"


def _rows(p: Path) -> dict[str, dict]:
    with p.open(newline="") as fh:
        return {r["source_id"]: r for r in csv.DictReader(fh)}


# ── the correction ───────────────────────────────────────────────────────────
def test_the_amarsi2019_doi_now_matches_the_ssot():
    ssot = G.ssot_rows()
    row = _rows(CNO_BIB)["Amarsi2019_Table1"]
    assert row["ssot_key"] == "amarsi2019"
    assert row["doi"] == ssot["amarsi2019"]["doi"] == "10.1051/0004-6361/201936265"
    assert row["doi"] != WRONG_DOI


def test_the_wrong_doi_is_gone_from_every_intake_artifact():
    """Not just from the row that was fixed — from anywhere it may have been copied."""
    hits = {f for doi, files in G.loose_dois().items() if doi == WRONG_DOI for f in files}
    assert not hits, f"the QSO-paper DOI survives in {sorted(hits)}"


def test_the_citation_string_and_the_doi_now_agree():
    """The tell that made this findable: the row named A&A 630 while its DOI was A&A 629."""
    row = _rows(CNO_BIB)["Amarsi2019_Table1"]
    _title, first_author, volume = CROSSREF_VERIFIED[row["doi"]]
    assert first_author in row["citation"]
    assert str(volume) in row["citation"]
    # and the old DOI would fail exactly that check
    assert CROSSREF_VERIFIED[WRONG_DOI][1] not in row["citation"]


# ── the sweep ────────────────────────────────────────────────────────────────
def test_no_linked_row_diverges_from_the_ssot():
    assert G.divergences() == []
    r = subprocess.run([sys.executable, "scripts/check_intake_ssot_rya1170.py", "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "0 divergence(s)" in r.stdout


def test_the_guard_fires_on_a_reintroduced_divergence(tmp_path, monkeypatch):
    """Mutation control. A guard only ever run against a clean tree has not been shown
    to fail — and this whole class survived a full QA pass unnoticed."""
    fake = tmp_path / "source_bibliography.csv"
    rows = list(csv.DictReader(CNO_BIB.open(newline="")))
    for r in rows:
        if r["source_id"] == "Amarsi2019_Table1":
            r["doi"] = WRONG_DOI
    with fake.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    monkeypatch.setattr(G, "INTAKE_BIBLIOGRAPHIES",
                        (str(fake.relative_to(ROOT)) if fake.is_relative_to(ROOT)
                         else str(fake),))
    monkeypatch.setattr(G, "ROOT", tmp_path if not fake.is_relative_to(ROOT) else G.ROOT)
    bad = G.divergences()
    assert any("Amarsi2019_Table1" in b and "!=" in b for b in bad), bad


def test_an_unlinked_row_is_reported_not_silently_accepted():
    """9 of the CNO intake's 11 sources are not in bibliography.csv at all. That is a real
    gap in the SSOT, and it must stay visible rather than pass as 'nothing to check'."""
    un = {sid for _f, sid, _d in G.unlinked()}
    assert "Li2015_CO" in un and "Brooke2013_C2" in un
    assert "Amarsi2019_Table1" not in un, "the corrected row must now be LINKED"
    assert "AGSS21" not in un


# ── the amendment: the held-document link ────────────────────────────────────
def test_agss21_points_at_the_held_document_with_a_checksum():
    """RYA-1170 amendment. The row carried asset='article' with an EMPTY sha256 while the
    SSOT records the paper as HELD (verified=extracted) — a broken link, not a missing
    source, and what left A1/A7 with no acquired referent."""
    ssot = G.ssot_rows()["asplund2021"]
    row = _rows(CNO_BIB)["AGSS21"]
    assert row["ssot_key"] == "asplund2021"
    assert ssot["local_file"] and row["asset"] == ssot["local_file"]
    assert len(row["sha256"]) == 64, "a key the SSOT records as HELD needs a checksum"
    assert row["status"] != "SOURCE_IDENTIFIED"


def test_an_empty_checksum_on_a_held_key_is_a_failure():
    """The amendment's rule, exercised rather than asserted."""
    rows = list(csv.DictReader(CNO_BIB.open(newline="")))
    for r in rows:
        if r["source_id"] == "AGSS21":
            r["sha256"] = ""
    import io
    bad = []
    ssot = G.ssot_rows()
    for row in rows:
        key = row["ssot_key"].strip()
        if key and ssot[key]["local_file"].strip() and not row["sha256"].strip():
            bad.append(row["source_id"])
    assert bad == ["AGSS21"]


def test_the_checksum_is_not_verified_when_the_library_is_absent(monkeypatch):
    """⚠️ The reference library sits above the repo and is a working directory, not a repo
    artifact — legitimately absent on CI. `generate_sources_page.audit_library` already
    treats absence as non-fatal; failing here would only teach people to disable it."""
    monkeypatch.setattr(G, "LIBRARY_ROOT", Path("/nonexistent-library-rya1170"))
    status = G.checksum_status("asplund2021", "0" * 64)
    assert status.startswith("UNVERIFIABLE-HERE")


def test_bibliography_csv_itself_was_not_touched():
    """The ticket's explicit constraint: the SSOT is correct and is not this ticket's to
    edit."""
    r = subprocess.run(["git", "diff", "--stat", "HEAD", "--", "data/refs/bibliography.csv"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.stdout.strip() == "", f"bibliography.csv was modified:\n{r.stdout}"
