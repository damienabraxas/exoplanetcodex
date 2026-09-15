"""RYA-1214 — a species with no laboratory table can still be graded.

THE DEFECT. `gf_grades.grade_line` opened by reading `lab_lines(species)`, which raises
for an unregistered species. The raise is correct and stays (RYA-833: "we hold no table"
must never look like "no measurement exists"), but as the FIRST statement it also made
the NIST branch below unreachable for C/N/O — species that will never have a laboratory
table, because the accepted light-element standard is critically-evaluated theory
(RYA-1172). So 842 NIST accuracy classes adjudicated into `canonical_gf` were invisible
to every error budget, and the first C I product charged the 0.17 UNGRADED placeholder
over six lines every one of which was graded.

WHAT MUST NOT REGRESS, in both directions:

  * a CNO pool whose lines all carry a NIST class reaches rung 2 with a MEASURED sigma;
  * that sigma is charged through a term that says COMPILATION, never "cited lab";
  * rung 3 stays laboratory-only, so no CNO pool can ever wear a laboratory pedigree;
  * one ungraded line still returns the whole pool to rung 1 (the RYA-855 mixed rule);
  * Fe and Al — species that DO have tables — are bit-identical to before.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import error_budget, gf_grades, gf_rung  # noqa: E402

#: Six AGSS21 C I red-optical indicators, at the wavelength / EP / log gf the first
#: RYA-1214 product actually measured them on. Real lines, not a fixture: the point is
#: that THIS pool reaches rung 2.
CI_POOL = pd.DataFrame(
    [(8335.147, 7.685, -0.4365), (9061.430, 7.483, -0.3468), (9078.280, 7.483, -0.5834),
     (9088.508, 7.483, -0.4295), (9094.829, 7.488, 0.1461), (9111.799, 7.488, -0.2967)],
    columns=["wavelength_air_A", "ep_eV", "log_gf"]).assign(resolved=True)


def test_has_lab_table_answers_without_raising():
    assert gf_grades.has_lab_table("Fe I") is True
    assert gf_grades.has_lab_table("Fe II") is True
    assert gf_grades.has_lab_table("Al I") is True
    for sp in ("C I", "C II", "N I", "N II", "O I", "O II"):
        assert gf_grades.has_lab_table(sp) is False, (
            f"{sp} must NOT have a laboratory table — CNO gf is critically-evaluated "
            f"THEORY (Opacity Project / MCHF), and registering one would let a computed "
            f"value be graded as a measurement (RYA-1172/1005)")
    # The raising accessor is UNCHANGED: the two facts stay distinguishable.
    with pytest.raises(KeyError):
        gf_grades.lab_lines("C I")


def test_a_cno_line_grades_to_gf_nist_instead_of_raising():
    v = gf_grades.grade_line(9094.829, 7.488, 0.1461, species="C I")
    assert v.gf_grade == gf_grades.GRADE_NIST
    assert v.has_cited_sigma and np.isfinite(v.gf_sigma_dex) and v.gf_sigma_dex > 0
    assert "structurally unreachable" in v.note, (
        "the verdict must SAY that GF-LAB is unreachable for this species rather than "
        "merely unmet — that is the RYA-833 distinction the raise used to carry")


def test_a_cno_line_can_never_come_back_gf_lab():
    """The refusal the old gate existed for. A C I line must not acquire an Fe I
    laboratory pedigree by a wavelength-and-EP coincidence."""
    for w, ep, gf in CI_POOL[["wavelength_air_A", "ep_eV", "log_gf"]].itertuples(index=False):
        v = gf_grades.grade_line(float(w), float(ep), float(gf), species="C I")
        assert v.gf_grade != gf_grades.GRADE_LAB
        assert v.is_graded is False, "is_graded means PRIMARY LABORATORY and must stay so"


def test_the_cno_pool_reaches_rung_2_with_a_measured_sigma():
    r = gf_rung.decide("C", "I", CI_POOL)
    assert r.rung == 2, r.reason
    assert r.nist_class is True
    assert r.gf_graded is True
    assert r.cited_sigma_dex is not None and 0.0 < r.cited_sigma_dex < 0.1
    assert r.term_name == "gf scale (cited NIST class)"
    assert "COMPILATION, NOT LABORATORY" in r.reason
    kw = r.budget_kwargs()
    assert "nist_class_gf_sigma_dex" in kw
    assert "cited_gf_sigma_dex" not in kw, (
        "a NIST accuracy class routed through the laboratory channel would be described "
        "as a 'published per-line laboratory sigma' by cited_gf_term's own text")


def test_the_budget_charges_a_compilation_term_and_not_the_blanket():
    r = gf_rung.decide("C", "I", CI_POOL)
    b = error_budget.build("C", 8000.0, n_lines=len(CI_POOL), scatter_dex=0.0165,
                           harness_residual_dex=0.0, handler="SynthesisHandler",
                           **r.budget_kwargs())
    text = b.render() if hasattr(b, "render") else str(b)
    assert not error_budget.carries_ungraded_gf(text), (
        "the pool is fully NIST-graded; charging the 0.17 placeholder is what this "
        "ticket fixes, and the publication gate refuses it (RYA-1212)")
    assert "gf scale (cited NIST class)" in text
    assert "gf scale (cited lab)" not in text
    assert "COMPILATION" in text


def test_one_ungraded_line_still_sends_the_whole_pool_to_rung_1():
    """The RYA-855 mixed-pool rule is not relaxed by this change."""
    mixed = pd.concat([CI_POOL,
                       pd.DataFrame([(9999.999, 1.0, -5.0)],
                                    columns=["wavelength_air_A", "ep_eV", "log_gf"])
                       .assign(resolved=True)], ignore_index=True)
    r = gf_rung.decide("C", "I", mixed)
    assert r.rung == 1
    assert r.nist_class is False
    assert "every line in it is" in r.reason


def test_an_unresolved_line_still_forces_rung_1():
    pool = CI_POOL.copy()
    pool.loc[0, "resolved"] = False
    r = gf_rung.decide("C", "I", pool)
    assert r.rung == 1 and r.n_unresolved == 1


def test_the_three_gf_routes_stay_mutually_exclusive():
    """They describe ONE term. RYA-1214 added a third, and a pairwise check would have
    admitted (cited, NIST) silently."""
    common = dict(element="Fe", wavelength_A=5000.0, n_lines=10, scatter_dex=0.02,
                  harness_residual_dex=0.0, handler="h", gf_graded=True)
    for a, b in (("cited_gf_sigma_dex", "nist_class_gf_sigma_dex"),
                 ("cited_gf_sigma_dex", "empirical_gf_sigma_dex"),
                 ("nist_class_gf_sigma_dex", "empirical_gf_sigma_dex")):
        with pytest.raises(ValueError, match="same gf term"):
            error_budget.build(**common, **{a: 0.04, b: 0.05},
                               cited_gf_source="x", nist_class_gf_source="y",
                               empirical_gf_provenance="z")


def test_a_nist_class_term_is_never_unsourced():
    with pytest.raises(ValueError, match="nist_class_gf_source"):
        error_budget.build("C", 8000.0, n_lines=6, scatter_dex=0.02,
                           harness_residual_dex=0.0, handler="h", gf_graded=True,
                           nist_class_gf_sigma_dex=0.04)


def test_fe_and_al_take_the_unchanged_path():
    """Every species WITH a laboratory table must be untouched by this change: the new
    branch is entered only when `has_lab_table` is False."""
    for element, ion in (("Fe", "I"), ("Fe", "II"), ("Al", "I")):
        assert (element, ion) in gf_rung.LAB_GRADED_SPECIES
        # An empty pool still short-circuits before either branch, and a one-line pool of
        # a species WITH a table must never be flagged nist_class.
        r = gf_rung.decide(element, ion,
                           pd.DataFrame([(5000.0, 3.0, -1.0)],
                                        columns=["wavelength_air_A", "ep_eV", "log_gf"])
                           .assign(resolved=True))
        assert r.nist_class is False, (
            f"{element} {ion} has a laboratory table; the RYA-1214 branch must be "
            f"unreachable for it")


# ── the RYA-1191 validity bound, which was inert for everything but iron ─────
def test_the_validity_bound_now_reaches_every_element_with_a_published_centre():
    """RYA-1214. `fit_is_physical` returned True unconditionally for `element != "Fe"`,
    so the guard that exists to keep a non-convergent fit out of an aggregate did nothing
    for C, N, O or Al. C I 4890.653 fitted to A = 5.443 on solar_harps_molecfit_corrected
    — 3.0 dex below solar carbon — and entered the VIS C I pool with a blank
    excluded_reason."""
    from pipeline.fit_validity import (SOLAR_A_FE, VALIDITY_HALF_WIDTH_DEX,
                                       fit_is_physical, rejection_reason,
                                       solar_reference)

    # Fe is bit-identical: same constant, not routed through the table.
    assert solar_reference("Fe") == SOLAR_A_FE
    assert fit_is_physical(7.46) and fit_is_physical(7.0) and fit_is_physical(8.0)
    assert not fit_is_physical(4.539) and not fit_is_physical(10.988)

    # The case that motivated it.
    assert not fit_is_physical(5.443, "C"), (
        "A(C) = 5.443 is 3.0 dex below solar carbon — a factor of a thousand — and must "
        "not enter an aggregate")
    assert fit_is_physical(8.50, "C") and fit_is_physical(8.46, "C")
    for el, a_ok, a_bad in (("N", 7.83, 5.0), ("O", 8.69, 4.5), ("Al", 6.43, 9.5)):
        assert fit_is_physical(a_ok, el), f"{el}: solar must be inside the bound"
        assert not fit_is_physical(a_bad, el), f"{el}: {a_bad} must be outside it"

    # An element with no published centre keeps the ORIGINAL refusal rather than
    # inventing one (RYA-833: "we hold no reference" is not "any value is fine", but
    # guessing a centre is worse than declining to filter).
    assert solar_reference("Zz") is None
    assert fit_is_physical(5.443, "Zz") is True

    # The reason names the centre it measured against — never a borrowed one.
    r = rejection_reason(5.443, "C")
    assert "A(C)" in r and "8.46" in r and str(VALIDITY_HALF_WIDTH_DEX) in r
    assert "7.46" not in r, "a carbon rejection must not quote the iron reference"


def test_absence_and_nan_are_still_not_rejections():
    """A missing abundance is a different problem (RYA-833) and must not be recoded as a
    non-physical fit by the generalisation."""
    from pipeline.fit_validity import fit_is_physical
    for el in ("Fe", "C", "N", "O"):
        assert fit_is_physical(None, el)
        assert fit_is_physical(float("nan"), el)
