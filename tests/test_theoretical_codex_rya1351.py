"""RYA-1351 / RYA-728: the Theoretical Codex grade and its first set, CR_ELGUETA2026.

A Theoretical Codex set carries theoretical gf with no published uncertainty and is adopted
by name (Ryan, 2026-11-15). It must never be stamped Reference or Codex Grade, its gf must
never be called graded or laboratory, and every declared line is measured.
"""
import csv
from pathlib import Path

import pandas as pd
import pytest

from pipeline import line_sets
from pipeline.cno_grade import grade_for
from pipeline.coverage_report import _grade
from pipeline.run_descriptor import line_set, set_derivation_csv, star_graded_csv

ROOT = Path(__file__).resolve().parents[1]
SET = "CR_ELGUETA2026"


def _rows(path):
    with (ROOT / path).open(newline="") as fh:
        return list(csv.DictReader(fh))


def test_theoretical_set_is_its_own_grade_and_other_sets_are_unchanged():
    assert grade_for({"element": "Cr", "selector": f"SET-{SET}"}) == "Theoretical Codex"
    assert grade_for({"element": "Cr", "selector": "SET-CR_AGSS21"}) == "Reference Grade"
    assert _grade(f"set:{SET}") == "Theoretical Codex"
    assert _grade("set:CR_AGSS21") == "Reference Grade"


def test_every_declared_theoretical_set_carries_a_complete_theoretical_adoption():
    declared = [s["name"] for s in line_sets.load_sources()
                if s.get("grade") == line_sets.THEORETICAL_CODEX]
    assert SET in declared
    for name in declared:
        adoption = line_sets.theoretical_codex(name)
        assert all(adoption[k] for k in line_sets.ADOPTION_KEYS)
        assert adoption["gf_origin"] == "theoretical"


def test_incomplete_or_non_theoretical_adoption_is_refused(monkeypatch):
    base = {"name": "X_TEST", "grade": line_sets.THEORETICAL_CODEX}
    full = {k: "x" for k in line_sets.ADOPTION_KEYS} | {"gf_origin": "theoretical"}
    for adoption in ({}, {k: v for k, v in full.items() if k != "approved_by"},
                     full | {"gf_origin": "experimental"}):
        monkeypatch.setattr(line_sets, "load_sources", lambda a=adoption: [base | {"adoption": a}])
        with pytest.raises(ValueError):
            line_sets.theoretical_codex("X_TEST")


def test_theoretical_set_measures_every_declared_line_others_keep_their_graded_file():
    row = line_set(SET, "Cr", "I")
    declared = _rows(row["csv"])
    assert int(row["n_graded"]) == 0 and len(declared) == int(row["n_lines"])
    assert star_graded_csv(row, "solar") == row["csv"] == set_derivation_csv(row, "solar")
    agss = line_set("CR_AGSS21", "Cr", "I")
    assert star_graded_csv(agss, "solar") != agss["csv"]


def test_theoretical_lines_are_not_relabelled_graded():
    from pipeline.gf_grades import is_gf_graded
    declared = _rows(line_set(SET, "Cr", "I")["csv"])
    assert {r["gf_graded"] for r in declared} == {"False"}
    canon = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    cr = canon[canon.species == "Cr I"]
    graded = is_gf_graded(cr)
    for r in declared:
        hit = (cr.wavelength_air_A - float(r["wavelength_air_A"])).abs() < 0.01
        assert hit.sum() == 1 and not graded[hit].iloc[0], r["wavelength_air_A"]
