"""RYA-1171 — the O I 777 triplet claimed better precision than the source it cites.

RYA-1160's EP-aware control reproduced every stored CNO log gf to <= 0.001 dex against
NIST ASD — the values are right. The GRADES were not: the store carried **A+ (<=2%)** where
NIST publishes **A (<=3%)** on all three lines of the most-used oxygen diagnostic in the
programme, so every gf_sigma_dex derived from that grade was too small by 1.49x.

validate-don't-tune in its cleanest form: no log gf moves, only the claimed precision is
brought back to what the cited authority actually says.
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

from pipeline.gf_grades import NIST_ACC_PCT, nist_sigma_dex  # noqa: E402
import rya1171_grade_overclaim_sweep as S  # noqa: E402

CANONICAL = ROOT / "data/linelists/canonical_gf.csv"
NIST_REF = ROOT / "data/linelists/nist_reference.csv"
NIST_XC = ROOT / "data/linelists/nist_crosscheck.csv"
CONTROL = ROOT / "data/audit/rya1160_cno_nist_gf/control.csv"

#: (wavelength_air_A, EP_eV, log_gf) — the log gf values are PINNED so a future edit that
#: moves one fails here, which is the whole point of a grade-only correction.
TRIPLET = ((7771.944, 9.146, "0.369"),
           (7774.166, 9.146, "0.223"),
           (7775.388, 9.146, "0.002"))
SIGMA_A = "0.0128"


def _rows(p: Path) -> list[dict]:
    with p.open(newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def _triplet_rows() -> list[dict]:
    out = []
    for w, ep, _g in TRIPLET:
        hit = [r for r in _rows(CANONICAL)
               if r["species"].strip() == "O I"
               and abs(float(r["wavelength_air_A"]) - w) <= S.TOL_A
               and abs(float(r["excitation_potential_eV"]) - ep) <= S.TOL_EP]
        assert len(hit) == 1, f"{w} resolved {len(hit)} rows on the EP-aware key"
        out.append(hit[0])
    return out


# ── the correction ───────────────────────────────────────────────────────────
def test_the_triplet_is_grade_a_with_the_matching_sigma():
    for row, (_w, _ep, _lgf) in zip(_triplet_rows(), TRIPLET):
        assert row["nist_grade"].strip() == "A"
        assert row["gf_sigma_dex"].strip() == SIGMA_A
        assert row["loggf_reference"].strip().endswith("grade A")
        assert "A+" not in row["loggf_reference"]


def test_the_stored_sigma_is_the_one_the_grade_implies():
    """0.0128 is not a typed-in number: it is round(log10(1 + 3/100), 4), and it is what
    the other 137 grade-A rows in this store already carry."""
    assert f"{nist_sigma_dex('A'):.4f}" == SIGMA_A
    others = [r["gf_sigma_dex"].strip() for r in _rows(CANONICAL)
              if r["nist_grade"].strip() == "A" and r["gf_sigma_dex"].strip()]
    assert others.count(SIGMA_A) >= 130, "0.0128 is no longer the grade-A convention"


def test_no_log_gf_moved():
    """🔴 THE CONSTRAINT. A grade correction that moves a value is a tuning, not a fix."""
    for row, (_w, _ep, lgf) in zip(_triplet_rows(), TRIPLET):
        assert row["log_gf"].strip() == lgf


def test_no_log_gf_moved_anywhere_in_the_store():
    """Not just on the three rows — the edit was line surgery and must have touched
    nothing else (RYA-1084: a pandas round-trip once changed 48 rows when 2 were meant)."""
    before = subprocess.run(
        ["git", "show", f"HEAD:{CANONICAL.relative_to(ROOT)}"],
        cwd=ROOT, capture_output=True, text=True)
    if before.returncode != 0 or not before.stdout:
        pytest.skip("pre-correction revision not reachable in this checkout")
    old = list(csv.DictReader(io.StringIO(before.stdout)))
    new = _rows(CANONICAL)
    if len(old) != len(new):
        pytest.skip("row population changed — a different ticket is in between")
    moved = [(a, b) for a, b in zip(old, new) if a["log_gf"] != b["log_gf"]]
    assert not moved, f"{len(moved)} log_gf value(s) moved"
    cols = {c for a, b in zip(old, new) for c in a if a[c] != b[c]}
    if cols:                       # only when HEAD predates the correction
        assert cols <= {"nist_grade", "gf_sigma_dex", "loggf_reference"}, cols


def test_the_accuracy_classes_are_what_the_ticket_says():
    assert NIST_ACC_PCT["A+"] == 2.0 and NIST_ACC_PCT["A"] == 3.0
    assert nist_sigma_dex("A") > nist_sigma_dex("A+")
    assert round(nist_sigma_dex("A") / nist_sigma_dex("A+"), 2) == 1.49


# ── spec item 4 — the refusal stands ─────────────────────────────────────────
def test_oi_6300_is_still_refused_and_was_not_given_a_grade():
    """Two EP-aware candidates at NIST, so there is no published grade to adopt. Refusing
    is the correct answer to an ambiguity; assigning one would be a guess wearing a
    citation (RYA-1072)."""
    c = next(r for r in _rows(CONTROL) if r["wavelength_air_A"].startswith("6300"))
    assert c["match"] == "AMBIGUOUS(2)"
    assert c["nist_grade"].strip() == "" and c["nist_log_gf"].strip() == ""
    assert "refused" in c["verdict"]
    # and nothing in this ticket stamped a NIST grade onto it from that refusal
    live = [r for r in _rows(CANONICAL)
            if r["species"].strip() == "O I"
            and abs(float(r["wavelength_air_A"]) - 6300.304) <= S.TOL_A]
    assert len(live) == 1 and live[0]["nist_grade"].strip() == "A"


# ── spec item 5 — the sweep ──────────────────────────────────────────────────
def test_the_sweep_reports_no_confirmed_overclaim():
    over, _suspect, _notes = S.sweep()
    assert over == [], over
    r = subprocess.run([sys.executable, "scripts/rya1171_grade_overclaim_sweep.py", "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_the_sweep_still_names_the_two_suspected_siblings():
    """Li I 6707 and Fe I 5576 claim tighter than a LOCAL record. They are reported, not
    corrected: confirming them needs a NIST pull for those species, and RYA-1160's pull
    was CNO only. Silence here would be the same defect one element over."""
    _over, suspect, _notes = S.sweep()
    joined = " ".join(suspect)
    assert "Li I 6707.760" in joined and "Fe I 5576.090" in joined
    assert len(suspect) == 2, suspect


def test_the_sweep_fires_when_the_overclaim_is_restored(monkeypatch, tmp_path):
    """Mutation control — a sweep only ever run against a corrected store proves nothing."""
    rows = _rows(CANONICAL)
    for r in rows:
        if r["species"].strip() == "O I" and r["wavelength_air_A"].startswith("7771"):
            r["nist_grade"] = "A+"
    fake = tmp_path / "canonical_gf.csv"
    with fake.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    monkeypatch.setattr(S, "CANONICAL", fake)
    over, _s, _n = S.sweep()
    assert any("7771" in o and "TIGHTER" in o for o in over), over


# ── the other copies ─────────────────────────────────────────────────────────
def test_the_second_copy_agrees_with_the_store():
    """`nist_reference.csv` carried the same A+ over-claim. A corrected store beside an
    uncorrected copy is the RYA-355 divergence this project keeps re-finding."""
    ref = {r["wavelength_air_A"]: r for r in _rows(NIST_REF)
           if r["element"] == "O" and r["ion"] == "I"}
    for w, _ep, _g in TRIPLET:
        assert ref[f"{w:.3f}"]["nist_grade"].strip() == "A"


def test_the_crosscheck_file_was_deliberately_left_alone():
    """🔴 `nist_crosscheck.csv` exists to CROSS-CHECK canonical_gf. It grades [O I] 6300.304
    A+ where the store says A — a real disagreement. Editing it to agree would destroy the
    only thing it is for: a referee you harmonise stops being a referee (RYA-853)."""
    r = subprocess.run(["git", "diff", "--stat", "HEAD", "--",
                        str(NIST_XC.relative_to(ROOT))],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.stdout.strip() == "", "nist_crosscheck.csv was modified — it must stay independent"
    xc = {r["wavelength_air_A"]: r for r in _rows(NIST_XC) if r["element"] == "O"}
    assert xc["6300.304"]["nist_grade"].strip() == "A+", "the disagreement must remain visible"
