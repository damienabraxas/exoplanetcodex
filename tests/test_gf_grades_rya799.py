"""
tests/test_gf_grades_rya799.py — RYA-799
========================================
The gf grade must be honest in three specific ways, and each is a separate test because
each has its own way of going quietly wrong:

  * it never leaves a bar blank (the ratified rule: a larger bar is a result, a hidden
    bar is a defect);
  * it never wears NIST's letters (RYA-711 — our values are `GF-`/`systematic:` prefixed,
    and NIST's own letters stay unprefixed wherever they are quoted);
  * **it never certifies a number that was not used.** This is the one that matters. The
    repo holds a better gf than the pool was measured with for dozens of these lines, and
    attaching that gf's pedigree to an abundance derived from a different value would
    manufacture a grade. The SCALE-MISMATCH state exists for exactly that case and this
    file pins it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import gf_grades as gg  # noqa: E402
from config.constants import LINE_GRADE_THRESHOLDS  # noqa: E402

LAB = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"


def _lab_row():
    lab = gg.lab_lines()
    ir = lab[(lab.wavelength_air_A >= 6910) & (lab.wavelength_air_A <= 9199)]
    assert len(ir), "no primary-lab line in the RYA-799 window — fixture is empty"
    return ir.iloc[0]


# ── the vendored primary-lab table ────────────────────────────────────────────

def test_lab_table_uncertainties_are_positive_and_small():
    """A whitespace read of a fixed-width VizieR table lands the PUBLISHED log gf in the
    uncertainty column. That is finite, so a not-null check passes; only a positivity and
    magnitude check catches it. It did, on the first run of the fetch script."""
    lab = pd.read_csv(LAB)
    assert (lab.e_loggf_dex > 0).all()
    assert lab.e_loggf_dex.max() <= 0.5


def test_lab_table_levels_and_wavelength_describe_the_same_transition():
    """Eup - Elo must reproduce hc/lambda_vac. If a parser mis-associated a level column
    with the wrong row this is the only check that would notice."""
    lab = pd.read_csv(LAB)
    hc = 12398.419843320026
    s2 = (1e4 / lab.wavelength_air_A.values) ** 2
    n = 1.0 + (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2)) * 1e-8
    resid = np.abs((lab.eup_eV - lab.elo_eV).values - hc / (lab.wavelength_air_A.values * n))
    assert resid.max() < 5e-3


# ── the grade vocabulary ──────────────────────────────────────────────────────

def test_our_grade_values_never_wear_nists_letters():
    """RYA-711: the prefix is in the VALUES, not just the column, because a value
    outlives the header it was emitted under."""
    nist_letters = {"AAA", "AA", "A+", "A", "B+", "B", "C+", "C", "D+", "D", "E"}
    for value in (gg.GRADE_LAB, gg.GRADE_SYSTEMATIC, gg.GRADE_MISMATCH,
                  gg.GRADE_NIST):
        assert value not in nist_letters
        assert value.startswith(("GF-", "systematic:"))


def test_the_gf_grade_is_a_different_vocabulary_from_the_measurement_grade():
    """`MQ-A..D` is measurement quality; a gf grade is a claim about the atomic datum.
    Reusing MQ- here would be a category error dressed as single-sourcing."""
    mq = set(LINE_GRADE_THRESHOLDS)
    assert all(v.startswith("MQ-") for v in mq)
    for value in (gg.GRADE_LAB, gg.GRADE_SYSTEMATIC, gg.GRADE_MISMATCH,
                  gg.GRADE_NIST):
        assert value not in mq
        assert not value.startswith("MQ-")


# ── the grading decision ──────────────────────────────────────────────────────

def test_a_line_measured_on_its_lab_value_is_graded_with_the_papers_sigma():
    r = _lab_row()
    v = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV), float(r.loggf))
    assert v.gf_grade == gg.GRADE_LAB
    assert v.is_graded
    assert v.gf_sigma_dex == pytest.approx(float(r.e_loggf_dex))
    assert "DOI" in v.gf_grade_source


def test_the_same_line_measured_on_a_DIFFERENT_gf_is_NOT_graded():
    """The defect this module exists to prevent: a lab measurement covering the line does
    not certify an abundance derived from some other gf."""
    r = _lab_row()
    v = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV),
                      float(r.loggf) + 0.30)
    assert v.gf_grade == gg.GRADE_MISMATCH
    assert not v.is_graded
    assert v.gf_sigma_dex == gg.K07_SYSTEMATIC_DEX
    assert "does NOT describe this measurement" in v.note


def test_an_unknown_line_falls_to_the_systematic_and_never_to_a_blank():
    v = gg.grade_line(7123.4567, 4.321, -9.99)
    assert v.gf_grade.startswith("systematic:")
    assert v.gf_sigma_dex == gg.K07_SYSTEMATIC_DEX
    assert v.gf_grade_source.strip() != ""


def test_matching_needs_the_excitation_potential_not_wavelength_alone():
    """RYA-780 ratified matching on wavelength AND EP. Hand the same wavelength a wrong
    EP and the lab tie must be refused — otherwise a blend would inherit a neighbour's
    pedigree (RYA-785: a same-species neighbour does not cancel)."""
    r = _lab_row()
    v = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV) + 2.0, float(r.loggf))
    assert v.gf_grade != gg.GRADE_LAB


def test_grade_pool_reconciles_and_refuses_a_blank_bar():
    pool = pd.DataFrame({
        "wavelength_air_A": [7123.4567, 8123.4567, 6999.0001],
        "ep_eV": [4.321, 4.5, 2.2],
        "log_gf": [-9.99, -1.23, -3.0],
    })
    out = gg.grade_pool(pool)
    assert len(out) == len(pool)
    assert (out.gf_grade.astype(str).str.strip() != "").all()
    assert (out.gf_grade_source.astype(str).str.strip() != "").all()
    assert np.isfinite(out.gf_sigma_dex).all()


def test_every_bounded_line_carries_exactly_the_ratified_systematic():
    """RYA-161's 0.20 dex, never averaged down and never quietly widened either."""
    assert gg.K07_SYSTEMATIC_DEX == 0.20
    pool = pd.DataFrame({"wavelength_air_A": [7123.4567], "ep_eV": [4.321],
                         "log_gf": [-9.99]})
    out = gg.grade_pool(pool)
    assert float(out.gf_sigma_dex.iloc[0]) == 0.20


# ── RYA-822: the NIST accuracy class as a citable sigma ──────────────────────

def test_the_nist_ladder_is_monotonic_and_includes_the_plus_tiers():
    """RYA-592: omitting the '+' tiers once put B+ (<=7%) BELOW B (<=10%) — an inverted
    ladder that silently DEMOTED the better measurement. Assert the order, not the keys.
    """
    order = ["AAA", "AA", "A+", "A", "B+", "B", "C+", "C", "D+", "D", "E"]
    assert list(gg.NIST_ACC_PCT) == order
    pct = [gg.NIST_ACC_PCT[g] for g in order]
    assert pct == sorted(pct), "the ladder must widen monotonically, best to worst"
    sig = [gg.nist_sigma_dex(g) for g in order]
    assert sig == sorted(sig), "sigma must widen with the grade, or a B+ outranks an A"


def test_nist_sigma_is_the_percent_to_dex_bridge_and_never_guesses():
    import numpy as np
    # grade B = 10% on A_ki = log10(1.10) dex, exactly as error_budget documents
    assert gg.nist_sigma_dex("B") == pytest.approx(np.log10(1.10), abs=1e-12)
    assert gg.nist_sigma_dex("A") < gg.nist_sigma_dex("B") < gg.nist_sigma_dex("E")
    # an unmapped class must NOT become an accuracy
    assert np.isnan(gg.nist_sigma_dex("Z"))
    assert np.isnan(gg.nist_sigma_dex(""))


def test_gf_nist_is_not_gf_lab_because_a_compilation_is_not_a_measurement():
    """RYA-760: FMW *is* NIST and VALD copies it. A compiled value must not inherit a
    measured value's pedigree, so GF-NIST stays outside `is_graded`."""
    assert gg.GRADE_NIST != gg.GRADE_LAB
    v = gg.GradeVerdict(gg.GRADE_NIST, "src", 0.04, "tag", -1.0, 0.0, "note")
    assert not v.is_graded          # not a primary lab measurement
    assert v.has_cited_sigma        # but it does carry a citable per-line sigma
    lab = gg.GradeVerdict(gg.GRADE_LAB, "src", 0.02, "tag", -1.0, 0.0, "note")
    assert lab.is_graded and lab.has_cited_sigma
    sysm = gg.GradeVerdict(gg.GRADE_SYSTEMATIC, "src", 0.20, "tag", -1.0, 0.0, "note")
    assert not sysm.is_graded and not sysm.has_cited_sigma
