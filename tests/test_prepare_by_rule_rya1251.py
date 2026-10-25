"""RYA-1251: preparation is the run's own first step, applied by rule.

Ryan 2026-10-24: the orchestrator runs A to Z with little agent work; culls are applied by
Ryan 2026-10-24: the orchestrator runs A to Z with little agent work; uncertainty is computed
after the runs, so a missing sigma on an otherwise-correct gf must never hold a measurement;
"I didnt say cull anything" -- cull candidates are measured and reported, never culled.
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


def test_the_run_never_culls(sandbox, monkeypatch):
    monkeypatch.setattr(prep, "CANONICAL", _canon(sandbox, []))
    plan = dict(adopt=[], missing_lines=[], generated="2026-10-24",
                cull_candidates=[dict(species="Na I", wavelength_A=5889.951,
                                      problem_class="SATURATION_COG", rew=-4.04, depth=0.95,
                                      holding="h")])
    assert prep.apply("solar", "Na", plan)["culled"] == 0
    assert (sandbox / "problem_children.csv").read_bytes() == b"header\r\n"
    assert "cull_candidates" not in prep.RULE_KINDS
    # a human's --cull still records it, for the appendix
    assert prep.apply("solar", "Na", plan, cull=True)["culled"] == 1
    culls = list(csv.DictReader((sandbox / "data/audit/orchestrator_prepare/solar_Na_culls.csv").open()))
    assert [(c["species"], float(c["wavelength_A"])) for c in culls] == [("Na I", 5889.951)]


def test_resolve_applies_the_plan_then_reprepares(monkeypatch):
    calls = []
    reps = iter([{"pending_rule_plan": {"adopt": 3}},
                 {"pending_rule_plan": {"adopt": 0}, "ready": True}])
    monkeypatch.setattr(prep, "prepare", lambda s, e: (calls.append("prepare"), next(reps))[1])
    monkeypatch.setattr(prep, "apply", lambda s, e, r: (calls.append("apply"), {"adopted": 3})[1])
    monkeypatch.setattr(prep, "has_published_feed", lambda s, e: False)
    rep, done = prep.resolve("solar", "Na")
    assert calls == ["prepare", "apply", "prepare"] and done == {"adopted": 3} and rep["ready"]


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


def test_a_cull_candidate_never_holds_a_cell():
    rep = _report(cull_candidates=[{"species": "Na I", "wavelength_A": 5889.951}])
    rm._RULE_PLAN_ASSUMED.clear()
    assert rm.preparation_hold_for_cell(rep, _cell()) == ""


def test_an_unadopted_gf_still_holds_a_cell_that_bypasses_resolve():
    rep = _report(adopt=[{"species": "Na I", "wavelength_A": 5889.951}], blocking={"adopt": 1})
    rm._RULE_PLAN_ASSUMED.clear()
    assert "adopt" in rm.preparation_hold_for_cell(rep, _cell())


def test_a_dry_run_forecasts_the_plan_as_applied():
    rep = _report(adopt=[{"species": "Na I", "wavelength_A": 5889.951}], blocking={"adopt": 1})
    rm._RULE_PLAN_ASSUMED.update(prep.RULE_KINDS)
    try:
        assert rm.preparation_hold_for_cell(rep, _cell()) == ""
    finally:
        rm._RULE_PLAN_ASSUMED.clear()


def test_a_sigma_alone_never_holds_a_measurement():
    rep = _report(sigma_only=[{"species": "Na I", "wavelength_A": 5500.0}])
    rm._RULE_PLAN_ASSUMED.clear()
    assert rm.preparation_hold_for_cell(rep, _cell()) == ""


def test_an_nlte_cell_with_no_labelled_line_is_held_not_failed(tmp_path, monkeypatch):
    from dataclasses import replace
    from config import synth_bands as sb
    lst = tmp_path / "ispec_x" / "atomic_lines.tsv"
    lst.parent.mkdir()
    lst.write_text("element\twave_A\tnlte_label_low\tnlte_label_up\n"
                   "Na 1\t10746.4\tnone\tnone\nFe 1\t10750.0\t3d6\t4p\n")
    band = replace(next(iter(sb.SYNTH_BANDS.values())), lo_A=10000.0, hi_A=11000.0,
                   linelist_spec=str(lst))
    monkeypatch.setattr(sb, "SYNTH_BANDS", {"x": band})
    rm._LABEL_COUNTS.clear()
    d = rm.RunDescriptor(element="Na", ion="I", instrument="i", holding="h", lo_A=10000.0,
                         hi_A=11000.0, engine_deck="gerber-nlte", method="synthesis",
                         pool="graded", star="solar")
    assert "none of the 1 Na I lines" in rm._nlte_label_hold(d)
    assert rm._nlte_label_hold(replace(d, engine_deck="gerber-1d-lte")) == ""
    # a sidecar that labels the line clears the hold by itself
    (lst.parent / "nlte_labels.tsv").write_text("element\twave_A\tnlte\nNa 1\t10746.4\tT\n")
    rm._LABEL_COUNTS.clear()
    assert rm._nlte_label_hold(d) == ""


def test_gerber_departure_scratch_is_per_process(monkeypatch):
    """RYA-1250/1251: one shared /tmp/rya798_gerber let parallel jobs clobber each other's
    departure file; the default is now per process."""
    import os
    from pipeline import gerber_nlte as g
    assert g.default_workdir() == f"/tmp/rya798_gerber_{os.getpid()}"
    seen = {}
    # the deck's aux table lives on Sirius; the test is about the scratch dir, not the deck
    monkeypatch.setattr(g, "abundance_axis", lambda el: [7.5])
    monkeypatch.setattr(g, "has_abundance_axis", lambda el: False)
    monkeypatch.setattr(g, "_interpolate", lambda *a: seen.setdefault("wd", a[-1]) and (_ for _ in ()).throw(RuntimeError("stop")))
    g._CACHE.clear()
    with pytest.raises(RuntimeError):
        g.for_node("Fe", 5772.0, 4.44, 0.0)
    assert seen["wd"] == g.default_workdir()


def test_a_missing_generated_line_list_is_built_before_dispatch(tmp_path, monkeypatch):
    import subprocess
    calls = []
    monkeypatch.setattr(rm, "GENERATED_LINELISTS", {"ispec_x": ["build.py"]})
    monkeypatch.setattr(rm, "__file__", str(tmp_path / "pipeline" / "run_matrix.py"))

    def fake(argv, **kw):
        calls.append(argv[1:])
        p = tmp_path / "data/linelists/ispec_x/atomic_lines.tsv"
        p.parent.mkdir(parents=True)
        p.write_text("x")
        return subprocess.CompletedProcess(argv, 0, "", "")
    monkeypatch.setattr(subprocess, "run", fake)
    assert rm.ensure_generated_linelists(dry_run=True, echo=False)["ispec_x"].startswith("missing")
    assert calls == []
    assert rm.ensure_generated_linelists(dry_run=False, echo=False) == {"ispec_x": "built"}
    assert calls == [["build.py"]]
    assert rm.ensure_generated_linelists(dry_run=False, echo=False) == {"ispec_x": "present"}


def _hfs_list():
    import numpy as np
    rows = [("Na 1", 10746.439, -2.100), ("Na 1", 10746.440, -2.100), ("Na 1", 10746.440, -2.498),
            ("Na 1", 10746.455, -1.653), ("Na 1", 10746.455, -2.100), ("Na 1", 10746.455, -2.799),
            ("Fe 1", 10746.449, -1.0)]
    dt = [("element", "U8"), ("wave_A", "f8"), ("loggf", "f8"), ("lower_state_eV", "f8"),
          ("lower_state_cm1", "f8"), ("upper_state_cm1", "f8")]
    return np.array([(e, w, g, 3.1914 if e == "Na 1" else 1.0,
                      25740.377 if e == "Na 1" else 8000.0, 35043.175 if e == "Na 1" else 17300.0)
                     for e, w, g in rows], dtype=dt)


def test_a_set_line_is_found_as_the_centroid_of_its_hfs_components():
    """RYA-1251: Asplund's Na I 10746.449 is the gf-weighted centroid of VALD's components at
    .439/.440/.455; no single row is within the set's 5 mA, so the NIR set cells held empty."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import derive_band_products as dbp
    ll = _hfs_list()
    keep, missing = dbp._match_hfs_centroids([10746.449], [3.19135], ll, "Na 1", 0.005)
    assert missing == [] and ll[keep[0]]["wave_A"] == 10746.455 and ll[keep[0]]["loggf"] == -1.653
    # a different lower level is a different transition, never this line
    assert dbp._match_hfs_centroids([10746.449], [2.10], ll, "Na 1", 0.005) == ([], [10746.449])
    # a centroid outside the set's window stays missing
    assert dbp._match_hfs_centroids([10746.470], [3.19135], ll, "Na 1", 0.005)[0] == []


def test_a_held_set_cell_is_not_budgeted_as_the_graded_cell_beside_it(monkeypatch):
    """RYA-1251: the old-report alias SET -> GRADED budgeted 4 HELD Na NIR set cells (36 legs,
    all rc=1) because the GRADED cell with the same identity was built."""
    from dataclasses import replace
    from pipeline import element_publish as ep
    d = replace(_cell(), pool=rm.SET_POOL_PREFIX + "NA_AGSS21")
    band = type("B", (), {"name_declared": "VIS"})()
    monkeypatch.setattr(rm, "expand", lambda *a, **k: [d])
    monkeypatch.setattr(rm, "_band_of", lambda *a: band)
    monkeypatch.setattr(rm, "pool_tier", lambda d: "SET-NA_AGSS21")
    import pipeline.run_descriptor as rd
    step = {"name": "derive_products", "args": ["--element", "Na"]}
    monkeypatch.setattr(rd, "resolve", lambda *a, **k: type("R", (), {"steps": [step]})())
    graded = dict(band="VIS", holding="h", ion="I", engine="ts-lte", pool="GRADED", status="DONE")
    held = dict(graded, pool="SET-NA_AGSS21", status="HELD")
    st = next(iter(ep.BUILT))
    graded["status"] = st
    assert ep.units_for("solar", "Na", report={"cells": [graded, held]}, interpreter="") == []
    # an old report that never named the set cell still reaches it through the alias
    assert len(ep.units_for("solar", "Na", report={"cells": [graded]}, interpreter="")) == 1


def test_a_leg_with_no_line_in_common_holds_the_component_not_the_element():
    """RYA-1251: one Na unit's legs shared no in-aggregate line; UnpairableProducts escaped
    `_paired` and killed the whole element's budget assembly after 603 legs."""
    import sys
    from pathlib import Path
    import pandas as pd
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import rya1230_cno_budget as cb
    nom = pd.DataFrame({"wavelength_air_A": [6154.225], "abundance": [6.30], "in_aggregate": [True]})
    leg = pd.DataFrame({"wavelength_air_A": [6160.747], "abundance": [6.31], "in_aggregate": [True]})
    p = cb._paired(leg, nom, 1)
    assert p["moved"] and p["n_paired"] == 0


def test_an_abundance_dependent_nlte_fit_never_probes_off_its_deck(monkeypatch):
    """RYA-1251: the A(sun) -3/+5 bracket put Brent's second probe at ~8.2 for Na, above the
    deck's 7.29; one refused trial failed every Na I NLTE line (19/19 VIS products, n=0)."""
    import sys
    import numpy as np
    from config.constants import ISPEC_DIR
    sys.path.insert(0, str(ISPEC_DIR))
    pytest.importorskip("ispec")
    from pipeline import abundances_derive as ad
    from pipeline import gerber_nlte as g
    trials = []
    monkeypatch.setattr(g, "assert_linelist_supports_nlte", lambda *a, **k: None)
    monkeypatch.setattr(g, "has_abundance_axis", lambda el: True)
    monkeypatch.setattr(g, "abundance_axis", lambda el: [4.29, 7.29])

    def for_node(el, t, lg, fe, abundance=None):
        if not 4.29 <= abundance <= 7.29:
            raise g.GerberDeckError(f"A({el}) = {abundance:.3f} is outside this deck's axis")
        return "dep"
    monkeypatch.setattr(g, "for_node", for_node)
    monkeypatch.setattr(g, "as_ispec_tuple", lambda dep, a: dep)

    def synth(sw, trial_A, nlte_departures=None, **kw):
        trials.append(trial_A)
        return 1.0 - 0.5 * np.exp(-((sw - 615.42) / 0.01) ** 2) * (trial_A - 3.0) / 3.0
    monkeypatch.setattr(ad, "_synth_flux_at_abund", synth)
    w = np.linspace(615.3, 615.5, 201)
    obs = synth(w, 6.25)
    trials.clear()
    r = ad._fit_synth_flux(w, obs, None, 5772.0, 4.44, 0.0, 1.0, None, None, None, "Na", 11,
                           615.3, 615.5, 3.24, 11.24, 700000.0, 0.0, 0.0, nlte_deck="gerber")
    assert r["status"] == "ok" and abs(r["A_X"] - 6.25) < 0.01
    assert trials and min(trials) >= 4.29 and max(trials) <= 7.29


def test_a_cell_whose_every_line_failed_synthesis_is_not_current(tmp_path, monkeypatch):
    """RYA-1251: 17 Na NLTE cells failed every fit and were recorded as built -- current, so a
    fixed engine never re-ran them. A cell with one real measurement stays current."""
    monkeypatch.setattr(rm, "ARTIFACT_ROOT", tmp_path)
    hdr = "wavelength_air_A,abundance,in_aggregate,excluded_reason\n"
    dead = tmp_path / "dead_lines.csv"
    dead.write_text(hdr + "6154.225,,False,SYNTHESIS: failed synthesis error: off axis\n"
                          "6160.747,,False,SYNTHESIS: failed synthesis error: off axis\n")
    live = tmp_path / "live_lines.csv"
    live.write_text(hdr + "6154.225,6.30,True,\n6160.747,,False,SYNTHESIS: failed x\n")

    def rec(p):
        return {"inputs_hash": "h", "recorded_at": "t",
                "artifacts": [{"path": p.name,
                               "sha256": rm._file_fingerprint(p)}]}
    assert rm.is_current([], "h", rec(live))[0]
    ok, why = rm.is_current([], "h", rec(dead))
    assert not ok and "failed SYNTHESIS" in why
    assert rm._every_line_failed_synthesis([{"path": str(dead)}]) == str(dead)


def test_a_deep_hfs_split_line_is_modelled_when_its_components_are_synthesised(tmp_path, monkeypatch):
    """RYA-1251: HFS is part of the recipe. A split line synthesised component by component
    is modelled at ANY depth (Na I 5688 is far from thin); one the synthesis sees as a single
    line is not, and only that holds."""
    import sys
    from pathlib import Path
    import pandas as pd
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import rya1230_cno_budget as cb
    ll = tmp_path / "atomic_lines.tsv"
    pd.DataFrame({"element": ["Na 1"] * 3 + ["Na 1"], "wave_A": [5688.20, 5688.21, 5688.22, 6154.23],
                  "lower_state_eV": [2.104] * 3 + [2.102], "theoretical_depth": [0.55, 0.40, 0.20, 0.30]}
                 ).to_csv(ll, sep="\t", index=False)
    monkeypatch.setattr(cb, "_synthesis_list_for", lambda w: ll)
    row = lambda w, ep: {"wavelength_air_A": w, "species": "Na I", "excitation_potential_eV": ep,
                         "hfs_n_components": 3}
    ev = cb._hfs_modelled_evidence(pd.DataFrame([row(5688.21, 2.104)]))
    assert ev and ev["lines"]["5688.210"]["n_components_in_synthesis"] == 3
    assert ev["lines"]["5688.210"]["max_component_depth"] > 0.10
    assert cb._hfs_modelled_evidence(pd.DataFrame([row(5688.21, 2.104), row(6154.23, 2.102)])) is None


def test_every_corrected_holding_has_a_telluric_route():
    """RYA-1251: telluric correction is done for every holding; the budget must be able to
    price it. Reiners 2016 was missing (KeyError), K had no route to the J arm's measured
    ratio, and the H rya1094 lines past the Vesta IDP segments had no sky."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import rya1230_cno_budget as cb
    assert [i for _, i in cb._INDEPENDENT["solar_iag_reiners2016"]] == [
        "solar_kpno_kurucz2005_corrected", "solar_kpno_molecfit_corrected"]
    assert cb._RATIO_FROM["solar_crires_plus_k_rya1219"][1] == 0.1632
    assert "solar_crires_plus_h_rya1232" in cb._OWN_SKY


def test_na_cited_solar_tables_cover_the_asplund_lines_and_refuse_a_double_cover(tmp_path, monkeypatch):
    """RYA-1251: AGSS21 Table A.1 gives the five Asplund Na lines a cited 1D (MARCS NLTE -
    MARCS LTE) and 3D (3D NLTE - MARCS LTE) correction; a second table may cite other lines
    (or one leg only) but never the same line twice."""
    import yaml
    from pipeline import solar_nlte_tables as s
    one = s.deltas("Na", "I", [6154.225, 6160.747, 10746.449], "1D", star="solar")
    three = s.deltas("Na", "I", [6154.225, 6160.747, 10746.449], "3D", star="solar")
    assert one == {6154.225: -0.047, 6160.747: -0.058, 10746.449: -0.006}
    assert three == {6154.225: -0.012, 6160.747: -0.015, 10746.449: 0.027}
    reg = yaml.safe_load(s.REGISTRY.read_text())
    dup = dict(next(r for r in reg["tables"] if r["element"] == "Na"), legs={"3D": "A_3d_nlte"})
    reg["tables"].append(dup)
    p = tmp_path / "reg.yaml"
    p.write_text(yaml.safe_dump(reg))
    monkeypatch.setattr(s, "REGISTRY", p)
    assert s.deltas("Na", "I", [6154.225], "1D", star="solar") == {6154.225: -0.047}
    with pytest.raises(ValueError, match="Ambiguous"):
        s.deltas("Na", "I", [6154.225], "3D", star="solar")


def test_the_engine_a_extension_is_owed_where_pysme_is_absent(tmp_path, monkeypatch):
    """RYA-1251: the run extends the ENGINE-A table itself; on a host without PySME the
    missing lines are reported as owed on Sirius, never silently skipped."""
    import json
    rep = tmp_path / "data/results/orchestrator/prepare"
    rep.mkdir(parents=True)
    (rep / "solar_Na.json").write_text(json.dumps({"nlte_missing": [
        {"species": "Na I", "wavelength_A": 8183.256}, {"species": "Na I", "wavelength_A": 8194.824}]}))
    monkeypatch.setattr(rm, "ROOT", tmp_path)
    import config.constants as C
    monkeypatch.setattr(C, "codex_path", lambda k: tmp_path / "no_such_python")
    out = rm.ensure_nlte_table("solar", "Na", dry_run=False, echo=False)
    assert out == {"I": "2 line(s) OWED: the ENGINE-A extension needs PySME "
                        "(engines.venv_pysme_python) -- run on Sirius"}
