"""RYA-1173 — the RYA-946 freeze gate, and the Al census that discharges it.

RYA-1141 check D3: Al was stamped FROZEN on 168 manifest rows while the mandatory AGSS21
reference-line-set census did not exist for Al and no exception was recorded. These tests cover
the gate that makes that impossible, the reconstructed reference set, and the census.

🔴 THE POINT OF THE NEGATIVE TESTS. A gate that passes is worth nothing until you have watched it
fail. Three of these remove a different leg — the registered set, the file behind it, the rows
inside the file — and each must turn the gate red on its own.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from pipeline import reference_census_gate as gate
from pipeline import reference_lineset as rls
from pipeline.model_registry import LINE_SETS

ROOT = Path(__file__).resolve().parents[1]
CENSUS = ROOT / "data" / "audit" / "rya1173_al_agss21_census"
REF = ROOT / "data" / "reference" / "asplund2021_al"


# ── the gate ─────────────────────────────────────────────────────────────────────────────

def test_the_gate_actually_evaluates_something():
    """⚠️ NON-VACUITY, FIRST. `check_all()` returning [] is meaningless if it swept nothing.

    Pointing the gate at the element ledger instead would do exactly that: no row of
    `intake_status_ledger.csv` is FROZEN, so a ledger-only guard passes forever without ever
    testing an element. The gate reads the per-line manifests, where FROZEN is really written.
    """
    frozen = gate.frozen_elements_in_manifests()
    assert frozen, "the gate swept no frozen rows at all -- it is not testing anything"
    assert "Al" in frozen, f"expected Al among frozen elements, got {sorted(frozen)}"
    assert sum(w["n_frozen_rows"] for w in frozen["Al"]) > 0


def test_the_repo_satisfies_the_gate_today():
    assert gate.check_all() == []


def test_al_is_satisfied_by_a_reference_set_not_by_an_exception():
    st = gate.census_state("Al")
    assert st["satisfied"] and st["census_complete"]
    assert st["exception"] is None, "Al must pass on the census itself, not on a bypass"
    assert [s["line_set"] for s in st["reference_sets"]] == ["asplund-al"]


def test_the_gate_would_have_caught_rya1141_d3(monkeypatch):
    """🔴 THE REGRESSION. Remove the set RYA-1173 built and the gate must flag Al — the exact
    state the repo was in when 168 rows were stamped FROZEN and nothing objected."""
    without = {k: v for k, v in rls.SETS.items() if v.element != "Al"}
    monkeypatch.setattr(rls, "SETS", without)
    bad = gate.check_all()
    assert any(b.startswith("Al:") for b in bad), f"gate stayed green with no Al set: {bad}"
    assert "168" in " ".join(bad) or "manifest row" in " ".join(bad)


def test_an_element_with_no_reference_set_is_flagged(tmp_path):
    """A synthetic frozen manifest for an element nothing has a census for."""
    d = tmp_path / "data" / "audit" / "rya9999_zz_intake"
    d.mkdir(parents=True)
    pd.DataFrame({"species": ["Zz I", "Zz I"], "intake_status": ["FROZEN", "CROSSMATCH_REVIEW"],
                  "wavelength_air": [5000.0, 5001.0]}).to_csv(d / "zz_manifest.csv", index=False)
    found = gate.frozen_elements_in_manifests(tmp_path)
    assert found["Zz"][0]["n_frozen_rows"] == 1, "only the FROZEN row counts"
    bad = gate.check_all(tmp_path)
    assert any(b.startswith("Zz:") for b in bad)
    assert "no reference set and no exception" in " ".join(bad)


def test_a_registered_set_whose_file_is_gone_does_not_satisfy_the_gate(monkeypatch):
    """🔴 A REGISTRY ENTRY IS NOT A CENSUS. `census_state` LOADS the set, so a declared-but-absent
    file fails — which an existence check on the registry alone would not catch."""
    broken = dict(rls.SETS)
    broken["asplund-al"] = rls.ReferenceSet(
        **{**rls.SETS["asplund-al"].__dict__, "path": ROOT / "data" / "reference" / "nope.csv"})
    monkeypatch.setattr(rls, "SETS", broken)
    st = gate.census_state("Al")
    assert not st["satisfied"]
    assert st["reference_sets"][0]["loads"] is False


def test_an_empty_reference_set_does_not_satisfy_the_gate(monkeypatch, tmp_path):
    """🔴 AND NEITHER IS AN EMPTY FILE. RYA-1141's D3 asked only whether a DIRECTORY existed; an
    `mkdir` would have satisfied it. Zero measurable rows must not read as a discharged census."""
    src = pd.read_csv(REF / "asplund2021_al_lines.csv")
    empty = src.iloc[0:0]
    p = tmp_path / "empty_lines.csv"
    empty.to_csv(p, index=False)
    stub = dict(rls.SETS)
    stub["asplund-al"] = rls.ReferenceSet(**{**rls.SETS["asplund-al"].__dict__, "path": p})
    monkeypatch.setattr(rls, "SETS", stub)
    st = gate.census_state("Al")
    assert st["reference_sets"][0]["n_used"] == 0
    assert not st["satisfied"]


def test_every_freeze_spelling_is_a_freeze():
    """RYA-946 names FROZEN_READY_FOR_MEASUREMENT; RYA-1132 wrote FROZEN and
    FROZEN_SOURCE_CONTROL for the same act. Reading the gate narrowly is how it was missed."""
    for s in ("FROZEN", "FROZEN_READY_FOR_MEASUREMENT", "FROZEN_SOURCE_CONTROL", " frozen "):
        assert gate.is_frozen(s), s
    for s in ("CROSSMATCH_REVIEW", "BLOCKED_PIPELINE_COVERAGE", "", None, "NOT_STARTED"):
        assert not gate.is_frozen(s), s


def test_an_incomplete_exception_is_refused(monkeypatch, tmp_path):
    """⚠️ The bypass must name a ticket, an approver and a date, or it is not an approval."""
    p = tmp_path / "census_exceptions.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=gate.EXCEPTION_FIELDS)
        w.writeheader()
        w.writerow({"element": "Zz", "ticket": "", "approved_by": "", "approved_on": "",
                    "source_publication_problem": "", "scope": "", "notes": ""})
    monkeypatch.setattr(gate, "EXCEPTIONS", p)
    with pytest.raises(gate.CensusGateError, match="incomplete"):
        gate.read_exceptions()


def test_no_exception_is_recorded_today():
    """If this ever fails, someone took the bypass and the test should make them say so out loud."""
    assert gate.read_exceptions() == {}


def test_require_census_raises_for_an_element_with_none():
    with pytest.raises(gate.CensusGateError, match="RYA-946"):
        gate.require_census("Zz")


# ── the reference set ────────────────────────────────────────────────────────────────────

def test_asplund_al_is_a_separate_line_set_from_asplund():
    """🔴 They must never share a product-key value: `asplund` is AGSS21's OWN Fe table, while
    AGSS21 publishes no Al lines at all and `asplund-al` is reconstructed from its citations."""
    assert "asplund-al" in LINE_SETS and "asplund" in LINE_SETS
    assert rls.SETS["asplund"].element == "Fe"
    assert rls.SETS["asplund-al"].element == "Al"
    assert rls.SETS["asplund-al"].path != rls.SETS["asplund"].path


def test_the_set_is_six_used_plus_one_published_exclusion():
    ref = rls.load("asplund-al")
    assert len(ref) == 7
    assert len(rls.measurable(ref)) == 6
    excl = rls.excluded_by_source(ref)
    assert len(excl) == 1
    assert round(float(excl.iloc[0].wavelength_air_A), 3) == 10891.732
    assert "telluric" in excl.iloc[0].selection_reason.lower()


def test_measurable_refuses_the_source_excluded_line():
    """🔴 THE FOOTGUN THIS CLOSED. `gf_missing` cannot see a selection: 10891.732 A has a perfectly
    good published gf and was thrown out by the SOURCE. Before `excluded_by_source` existed,
    `measurable()` returned all 7 and a 'replication of AGSS21' would have measured a line
    AGSS21's own source rejected."""
    ref = rls.load("asplund-al")
    assert not ref["gf_missing"].any(), "no row here lacks a gf -- so gf_missing alone cannot help"
    assert 10891.732 not in set(rls.measurable(ref).wavelength_air_A.round(3))


def test_the_other_sets_are_unchanged_by_the_selection_column():
    """Sets whose file has no selection_status default to USED; their counts must not move."""
    assert len(rls.measurable(rls.load("asplund"))) == 53
    assert len(rls.measurable(rls.load("gbs"))) == 150
    for n in ("asplund", "gbs"):
        assert rls.excluded_by_source(rls.load(n)).empty


def test_the_published_count_conflict_is_recorded_not_resolved():
    """AGSS21 says five, the primary says six, Scott lists seven. Six is adopted on the primary's
    authority and the conflict must remain visible to anyone reading the artifact."""
    prov = json.loads((REF / "asplund2021_al_lines.prov.json").read_text())
    c = prov["🔴 published_line_count_conflict"]
    assert (c["scott2015b"], c["nordlander_lind2017"], c["agss21_prose"]) == (7, 6, 5)
    assert c["adopted_here"] == 6
    assert "reproduces from neither source it cites" in c["why_not_five"]


def test_every_gf_in_the_set_is_theory_and_the_artifact_says_so():
    """RYA-946: an AGSS21 abundance value is not a gf grade. This set is Opacity Project theory
    end to end, so a Codex line matching into it has gained nothing."""
    prov = json.loads((REF / "asplund2021_al_lines.prov.json").read_text())
    assert "theory, not laboratory" in prov["🔴 the_gf_are_theory"].lower()
    d = pd.read_csv(REF / "asplund2021_al_lines.csv")
    assert d.gf_source_per_line.str.contains("TOPbase").all()
    assert (d.gf_source_per_line.str.contains("Mendoza")).all()


def test_roles_are_not_guessed_for_the_wider_analysis_table():
    """RYA-946 forbids inferring a line's role from its wavelength. NL2017 Table A.1's note is a
    statement about the TABLE, so most rows must carry the unstated marker."""
    a = pd.read_csv(REF / "nordlander_lind_2017_analysis_lines.csv")
    assert len(a) == 55
    assert (a.role == "ANALYSIS_LINE_ROLE_NOT_STATED_PER_LINE").sum() == 43
    assert (a.role == "SOLAR_ABUNDANCE_USED").sum() == 6


# ── the census ───────────────────────────────────────────────────────────────────────────

def test_the_census_produces_every_rya946_deliverable():
    for name in ("per_line_join.csv", "band_coverage_matrix.csv", "four_way_comparison.csv",
                 "tolerance_plateau_and_null.csv", "lineage_note.md", "census_verdict.json"):
        assert (CENSUS / name).exists(), name


def test_the_per_line_join_carries_rya946s_eighteen_fields():
    j = pd.read_csv(CENSUS / "per_line_join.csv")
    for col in ("reference_line_set", "adopted_solar_value_lineage", "source_paper_table",
                "species_isotopologue", "wavelength_A", "wavelength_medium", "elo_eV",
                "transition_identity", "published_loggf", "published_gf_source",
                "published_ew_mA", "published_weight", "published_selection_status",
                "source_band", "codex_canonical_line_id", "codex_loggf", "codex_gf_tier",
                "codex_nist_grade", "codex_lab_source", "codex_gf_sigma_dex",
                "codex_gf_source_doi", "delta_loggf_codex_minus_published", "join_status",
                "ambiguity_status_note"):
        assert col in j.columns, col


def test_the_band_matrix_reports_every_band_including_the_empty_ones():
    """'The source used none in this band' is a RESULT, not a blank row to omit."""
    cov = pd.read_csv(CENSUS / "band_coverage_matrix.csv")
    assert set(cov.band) == set(
        ["FUV", "NUV", "near-UV", "VIS", "red-optical", "NIR", "J", "H", "K",
         "OUTSIDE_CURRENT_INSTRUMENT_REACH"])
    assert cov.published_used.sum() == 6
    assert cov.published_explicitly_rejected.sum() == 1
    # the six sit in VIS/red-optical/NIR only, and the empty bands are a published selection
    assert set(cov[cov.published_used > 0].band) == {"VIS", "red-optical", "NIR"}


def test_the_tolerance_plateau_has_a_measured_null():
    """🔴 A plateau must be qualified by its null. A flat count with a non-zero null is
    saturation, not agreement."""
    s = pd.read_csv(CENSUS / "tolerance_plateau_and_null.csv")
    v = json.loads((CENSUS / "census_verdict.json").read_text())
    lo, hi = v["match"]["plateau_A"]
    inside = s[(s.tol_A >= lo) & (s.tol_A <= hi)]
    assert len(inside) > 1, "a one-point 'plateau' is a coincidence of one window width"
    assert (inside.null_max == 0).all(), "the null is not zero across the adopted plateau"
    assert (inside.ambiguous == 0).all()
    assert inside.matched.nunique() == 1, "the count still moves inside the claimed plateau"
    # and the sweep must actually go wrong somewhere, or it proves nothing
    assert (s.null_max > 0).any(), "the null never fires -- the sweep does not reach chance"


def test_ep_contributes_nothing_to_this_sets_join():
    """⚠️ RYA-1151 in its sharpest form. The worst pair is not merely close in EP, it is IDENTICAL:
    6696.015 and 6698.672 A both leave the 4s 2S(1/2) lower level, so Delta EP is exactly 0.000 eV
    against EP_TOL_EV = 0.02. Only 2.657 A and the UPPER level separate them, and an EP tolerance
    is not a level identity."""
    v = json.loads((CENSUS / "census_verdict.json").read_text())
    ep = v["match"]["⚠️ ep_separation_margin"]
    assert ep["ep_would_separate_them"] is False
    assert ep["delta_EP_eV"] == 0.0, "the two most-weighted lines share a lower level exactly"
    assert sorted(ep["pair_A"]) == [6696.015, 6698.672]
    # and the wavelength gap that DOES separate them must sit far above the adopted window
    assert ep["delta_lambda_A"] > 100 * v["match"]["adopted_tol_A"]


def test_one_of_the_six_used_lines_is_absent_from_codex():
    """🔴 THE CENSUS'S OWN FINDING. Al I 10768.363 A carries part of AGSS21's adopted value and
    has no row in canonical_gf, so that value cannot be replicated on our line list today."""
    j = pd.read_csv(CENSUS / "per_line_join.csv")
    absent = j[j.join_status == "ABSENT_FROM_CODEX"]
    assert len(absent) == 1
    assert round(float(absent.iloc[0].wavelength_A), 3) == 10768.363
    assert absent.iloc[0].published_selection_status == "USED_BY_SOURCE_ANALYSIS"

    can = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    w = pd.to_numeric(can.wavelength_air_A, errors="coerce")
    # not an ingest hole: the region is well populated with other species
    assert ((w > 10700) & (w < 10850)).sum() > 20
    assert can[(w > 10760) & (w < 10775) & can.species.eq("Al I")].empty


def test_the_census_does_not_claim_to_unblock_measurement():
    """It discharges ONE gate. Saying more would be the RYA-1141 failure with the sign flipped."""
    v = json.loads((CENSUS / "census_verdict.json").read_text())
    assert v["census_complete"] is True
    assert "not" in v["does_not_unblock"].lower()
    assert "FROZEN_READY_FOR_MEASUREMENT" in v["does_not_unblock"]


def test_the_wall_does_not_apply_because_the_primaries_publish_the_rows():
    v = json.loads((CENSUS / "census_verdict.json").read_text())
    assert v["source_line_list_not_published_wall"] is False
    assert "AGSS21 publishes no Al line list" in v["wall_note"]
