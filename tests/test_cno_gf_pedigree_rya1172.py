"""RYA-1172 — "NIST grade" is not one authority across CNO.

The 2006 Wiese & Fuhr update is explicitly PARTIAL: C I, C II, N I, N II carry MCHF, while
O I / O II still rest on the 1996 WFD Monograph 7 / OPACITY Project values. Same label, ten
years and a different method apart, and nothing in the store said so.

Every string is quoted from `Reference documents/20060052476.pdf`. The assertions below
quote the paper back, so a future edit that paraphrases the pedigree fails here.
"""
from __future__ import annotations

import csv
import io
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from pipeline.cno_gf_pedigree import (  # noqa: E402
    MCHF_UPDATE, PEDIGREE_BY_SPECIES, UNESTABLISHED, WFD1996_OPACITY,
    is_cno, pedigree_for)

INTAKE = ROOT / "data/audit/rya1129_atomic_intake"
LEDGER = INTAKE / "cno_gf_pedigree_rya1172.csv"
MANIFESTS = {e: INTAKE / f"{e}_atomic_manifest.csv" for e in ("C", "N", "O")}


def _rows(p: Path) -> list[dict]:
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


# ── the two pedigrees ────────────────────────────────────────────────────────
def test_the_update_is_partial_and_oxygen_is_not_in_it():
    """🔴 The whole ticket. The paper's own title names Carbon, Nitrogen and IRON."""
    assert [s for s, p in PEDIGREE_BY_SPECIES.items() if p is MCHF_UPDATE] == \
        ["C I", "C II", "N I", "N II"]
    assert [s for s, p in PEDIGREE_BY_SPECIES.items() if p is WFD1996_OPACITY] == \
        ["O I", "O II"]
    assert "partial update" in MCHF_UPDATE.evidence
    assert "C I, C II, N I and N II" in MCHF_UPDATE.evidence
    assert "oxygen is not in it" in WFD1996_OPACITY.evidence


def test_the_vintages_and_methods_actually_differ():
    assert MCHF_UPDATE.vintage == 2006 and WFD1996_OPACITY.vintage == 1996
    assert "MCHF" in MCHF_UPDATE.method
    assert "OPACITY Project" in WFD1996_OPACITY.method
    assert MCHF_UPDATE.method != WFD1996_OPACITY.method


def test_the_strings_are_quoted_from_the_held_document():
    """No value from memory: each string must appear in the PDF we hold. Skipped where the
    reference library is absent (it lives above the repo and is not a repo artifact)."""
    pdf = Path.home() / "Documents/Exoplanet Codex/Reference documents/20060052476.pdf"
    if not pdf.is_file():
        pytest.skip("reference library not present on this machine")
    pypdf = pytest.importorskip("pypdf")
    text = "\n".join((p.extract_text() or "") for p in pypdf.PdfReader(str(pdf)).pages)
    # ⚠️ The extractor scatters spaces INSIDE words ("transition p robabilities"), so a
    # word-wise comparison fails on text that is verbatim. Strip all whitespace from both
    # sides and compare the letters.
    def squash(t: str) -> str:
        return "".join(t.split()).replace("ﬁ", "fi").replace("ﬀ", "ff").replace("ﬃ", "ffi")

    flat = squash(text)
    for phrase in ("Carbon, Nitrogen, and Iron",
                   "partial update for the transition probabilities of C I, C II",
                   "was primarily based on the very extensive calculational results of the OPACITY Project",
                   "Monograph No. 7"):
        assert squash(phrase) in flat, phrase


def test_neither_pedigree_is_a_laboratory_measurement():
    """⚠️ Both are calculations. Tiering these LAB would repeat RYA-1005's Al mistake."""
    assert not MCHF_UPDATE.is_laboratory and not WFD1996_OPACITY.is_laboratory
    assert not any(p.is_laboratory for p in PEDIGREE_BY_SPECIES.values())
    src = (ROOT / "pipeline/cno_gf_pedigree.py").read_text()
    assert "is_laboratory=True" not in src, "no code path may set a CNO pedigree LAB"


# ── the unestablished stages ─────────────────────────────────────────────────
def test_a_stage_no_source_covers_is_unknown_not_inherited():
    """🔴 The store holds C III-C V, N III-N V, O III-O VI and neither document states a
    pedigree for them. An unrecognised case must not be laundered into the recognised one
    (RYA-1072) — the 1996 volume covering 'oxygen' does not license O VI."""
    for sp in ("O III", "O VI", "C III", "C V", "N III", "N V"):
        p = pedigree_for(sp)
        assert p is UNESTABLISHED, sp
        assert p.compilation == "" and p.method == "" and p.vintage == 0
    assert pedigree_for("Fe I") is UNESTABLISHED, "non-CNO must not pick up a CNO pedigree"


def test_the_unestablished_group_is_not_described_as_theory():
    """`is_laboratory=False` there is the SAFE DEFAULT, not a claim about the method.
    Calling it 'theory' would be this ticket's own over-claim, one group down."""
    out = subprocess.run([sys.executable, "scripts/rya1172_cno_pedigree_report.py"],
                         cwd=ROOT, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    block = out.stdout.split("NOT_ESTABLISHED_BY_HELD_SOURCES")[1]
    assert "NOT established" in block
    assert "critically-evaluated THEORY" not in block.split("wrote")[0]


# ── carried on the rows ──────────────────────────────────────────────────────
def test_every_cno_manifest_row_carries_vintage_and_method():
    for element, path in MANIFESTS.items():
        rows = _rows(path)
        assert rows, element
        for r in rows:
            assert r["gf_authority"], f"{element}: a row with no authority"
            if r["gf_authority"] != UNESTABLISHED.key:
                assert r["gf_authority_vintage"] and r["gf_authority_method"]
            assert r["gf_authority_is_laboratory"] in ("False", "false")


def test_the_o_rows_read_1996_and_the_c_n_rows_read_mchf():
    o = {r["species"]: r for r in _rows(MANIFESTS["O"])}
    c = {r["species"]: r for r in _rows(MANIFESTS["C"])}
    assert o["O I"]["gf_authority_vintage"] == "1996"
    assert "OPACITY" in o["O I"]["gf_authority_method"]
    assert c["C I"]["gf_authority_vintage"] == "2006"
    assert "MCHF" in c["C I"]["gf_authority_method"]
    assert o["O I"]["gf_authority"] != c["C I"]["gf_authority"]


def test_the_ledger_groups_the_species_by_pedigree():
    rows = _rows(LEDGER)
    assert len(rows) == 16, "every held CNO species must appear"
    by = {}
    for r in rows:
        by.setdefault(r["gf_authority"], []).append(r["species"])
    assert sorted(by["NIST_MCHF_PARTIAL_UPDATE_2006"]) == ["C I", "C II", "N I", "N II"]
    assert sorted(by["WFD1996_MONOGRAPH7_OPACITY_PROJECT"]) == ["O I", "O II"]
    assert len(by) == 3
    assert all(r["source_document"].endswith("20060052476.pdf") for r in rows)


def test_no_gf_value_changed():
    """Provenance/schema only: the manifests gain columns and move none."""
    for path in MANIFESTS.values():
        before = subprocess.run(["git", "show", f"HEAD:{path.relative_to(ROOT)}"],
                                cwd=ROOT, capture_output=True, text=True)
        if before.returncode != 0 or not before.stdout:
            pytest.skip("pre-change revision unavailable")
        old = list(csv.DictReader(io.StringIO(before.stdout)))
        new = _rows(path)
        if set(old[0]) == set(new[0]):
            continue
        assert len(old) == len(new)
        assert set(old[0]) < set(new[0]), "a column was removed"
        for k in old[0]:
            assert all(a[k] == b[k] for a, b in zip(old, new)), f"{path.name}: {k} moved"


def test_canonical_gf_was_not_touched():
    r = subprocess.run(["git", "status", "--porcelain", "data/linelists/canonical_gf.csv"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.stdout.strip() == ""
