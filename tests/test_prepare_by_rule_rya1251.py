"""RYA-1251: preparation is the run's own first step, applied by rule.

Ryan 2026-10-24: the orchestrator runs A to Z with little agent work; culls are applied by
rule (listed on the appendix); uncertainty is computed after the runs, so a missing sigma on
an otherwise-correct gf must never hold a measurement.
"""
import csv

import pytest

from pipeline import element_prepare as prep
from pipeline import run_matrix as rm


def _canon(tmp_path, rows):
    path = tmp_path / "canonical.csv"
    cols = ["line_id", "species", "wavelength_air_A", "log_gf", "gf_tier", "gf_sigma_dex",
            "loggf_reference", "gf_source_doi", "lab_source_tag", "adjudication_status"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    return path


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    from pipeline import line_sets
    monkeypatch.setattr(prep, "ROOT", tmp_path)
    monkeypatch.setattr(prep, "REGISTRY_PC", tmp_path / "problem_children.csv")
    (tmp_path / "problem_children.csv").write_bytes(b"header\r\n")
    monkeypatch.setattr(line_sets, "load_sources",
                        lambda: [{"name": "NA_AGSS21", "citation": "AGSS21", "doi": "x"}])
    monkeypatch.setattr(line_sets, "build", lambda *a, **k: None)
    return tmp_path


def test_protected_status_is_anything_another_ticket_decided():
    for open_ in ("", "nan", None, "single_source", "pending", "agreed", "adopted_rya1233_x"):
        assert not prep.is_protected(open_)
    for owned in ("nist", "isotope_rya1075", "adjudicated_kr1968_rya466"):
        assert prep.is_protected(owned)


def test_sigma_only_writes_the_sigma_and_nothing_else(sandbox, monkeypatch):
    path = _canon(sandbox, [dict(line_id="a", species="Na I", wavelength_air_A="6154.225",
                                 log_gf="-1.547", gf_tier="NIST-C+", loggf_reference="NIST",
                                 adjudication_status="nist")])
    monkeypatch.setattr(prep, "CANONICAL", path)
    done = prep.apply("solar", "Na", dict(
        adopt=[], missing_lines=[], cull_candidates=[],
        sigma_only=[dict(line_id="a", species="Na I", canonical_wavelength_A=6154.225,
                         authors_sigma_dex=0.0128, set="NA_AGSS21")]))
    row = next(csv.DictReader(path.open()))
    assert done["sigma_only"] == 1 and done["adopted"] == 0
    assert row["gf_sigma_dex"] == "0.0128"
    assert (row["log_gf"], row["gf_tier"], row["loggf_reference"], row["adjudication_status"]) == \
        ("-1.547", "NIST-C+", "NIST", "nist")


def test_an_adoption_never_rewrites_a_protected_row(sandbox, monkeypatch):
    path = _canon(sandbox, [dict(line_id="a", species="Cu I", wavelength_air_A="5218.201",
                                 log_gf="0.364", adjudication_status="adjudicated_kr1968_rya466")])
    monkeypatch.setattr(prep, "CANONICAL", path)
    prep.apply("solar", "Cu", dict(
        missing_lines=[], cull_candidates=[],
        adopt=[dict(line_id="a", species="Cu I", canonical_wavelength_A=5218.201,
                    authors_loggf=0.264, authors_sigma_dex=0.02, set="NA_AGSS21")]))
    row = next(csv.DictReader(path.open()))
    assert row["log_gf"] == "0.364" and row["adjudication_status"] == "adjudicated_kr1968_rya466"


def test_culls_are_recorded_for_the_appendix(sandbox, monkeypatch):
    monkeypatch.setattr(prep, "CANONICAL", _canon(sandbox, []))
    prep.apply("solar", "Na", dict(
        adopt=[], missing_lines=[], generated="2026-10-24",
        cull_candidates=[dict(species="Na I", wavelength_A=5889.951, problem_class="SATURATION_COG",
                              rew=-4.04, depth=0.95, holding="h")]))
    culls = list(csv.DictReader((sandbox / "data/audit/orchestrator_prepare/solar_Na_culls.csv").open()))
    assert [(c["species"], float(c["wavelength_A"])) for c in culls] == [("Na I", 5889.951)]
    assert b"applied by rule" in (sandbox / "problem_children.csv").read_bytes()


def test_resolve_applies_the_plan_then_reprepares(monkeypatch):
    calls = []
    reps = iter([{"pending_rule_plan": {"cull_candidates": 3}},
                 {"pending_rule_plan": {"cull_candidates": 0}, "ready": True}])
    monkeypatch.setattr(prep, "prepare", lambda s, e: (calls.append("prepare"), next(reps))[1])
    monkeypatch.setattr(prep, "apply", lambda s, e, r: (calls.append("apply"), {"culled": 3})[1])
    monkeypatch.setattr(prep, "has_published_feed", lambda s, e: False)
    rep, done = prep.resolve("solar", "Na")
    assert calls == ["prepare", "apply", "prepare"] and done == {"culled": 3} and rep["ready"]


def test_a_published_element_keeps_its_human_checkpoint(monkeypatch):
    monkeypatch.setattr(prep, "prepare", lambda s, e: {"pending_rule_plan": {"adopt": 2}})
    monkeypatch.setattr(prep, "apply", lambda *a: pytest.fail("finished work re-planned"))
    monkeypatch.setattr(prep, "has_published_feed", lambda s, e: True)
    assert prep.resolve("solar", "Fe")[1] == {}


def _report(**kw):
    base = dict(ready=False, scoped_readiness=True, nlte_missing=[], adopt=[], missing_lines=[],
                cull_candidates=[], holdings=[], literature={"litscan": True, "missing_bib_keys": []},
                blocking={})
    base.update(kw)
    return base


def _cell():
    return rm.RunDescriptor(element="Na", ion="I", instrument="harps", holding="h",
                            lo_A=5000.0, hi_A=6000.0, engine_deck="ts-lte",
                            method="synthesis", pool="graded", star="solar")


def test_an_unapplied_cull_still_holds_a_cell_that_bypasses_resolve():
    rep = _report(cull_candidates=[{"species": "Na I", "wavelength_A": 5889.951}],
                  blocking={"cull_candidates": 1})
    rm._RULE_PLAN_ASSUMED.clear()
    assert "cull_candidates" in rm.preparation_hold_for_cell(rep, _cell())


def test_a_dry_run_forecasts_the_plan_as_applied():
    rep = _report(cull_candidates=[{"species": "Na I", "wavelength_A": 5889.951}],
                  blocking={"cull_candidates": 1})
    rm._RULE_PLAN_ASSUMED.update(prep.RULE_KINDS)
    try:
        assert rm.preparation_hold_for_cell(rep, _cell()) == ""
    finally:
        rm._RULE_PLAN_ASSUMED.clear()


def test_a_sigma_alone_never_holds_a_measurement():
    rep = _report(sigma_only=[{"species": "Na I", "wavelength_A": 5500.0}])
    rm._RULE_PLAN_ASSUMED.clear()
    assert rm.preparation_hold_for_cell(rep, _cell()) == ""
