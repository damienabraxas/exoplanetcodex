"""RYA-1233: the front door -- `--all-elements` sweeps every canonical element, Fe first.

What is pinned: the element axis is derived (never typed), Fe runs first and gates
the rest, no element is ever absent from the sweep report, one element breaking
never stops the others, and the legacy whole-star path still runs -- announced.
`run_matrix.run` is stubbed in the sweep tests: its own behaviour is pinned by
tests/test_rya1222_run_matrix.py, and what is under test here is the sweep's
bookkeeping around it.
"""
from __future__ import annotations

import ast
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


@pytest.fixture(scope="module", autouse=True)
def _kp(tmp_path_factory):
    """Same RYA-1064 import-time guard as test_rya1222_run_matrix."""
    kp = tmp_path_factory.mktemp("kp")
    (kp / "lm0296").touch()
    os.environ.setdefault("CODEX_KP_ATLAS", str(kp))


@pytest.fixture(scope="module")
def rm():
    from pipeline import run_matrix
    return run_matrix


@pytest.fixture(scope="module")
def sw():
    from pipeline import run_sweep
    return run_sweep


def _master(rm):
    return json.loads(rm.ELEMENTS_MASTER.read_text(encoding="utf-8"))["elements"]


def _stub_run(rm, monkeypatch, verdict):
    """run_matrix.run replaced by `verdict(target, ions) -> {status: n}` (or raise).

    The graded pool is stubbed to hold I and II of everything, so these tests are
    about the sweep's bookkeeping and not about which elements are graded today.
    """
    monkeypatch.setattr(rm, "graded_ions", lambda sym: ["I", "II"])
    calls = []

    def fake(star, target, *, ions=None, dry_run=False, **kw):
        calls.append((target, tuple(ions or ())))
        got = verdict(target, ions)
        counts = {s: 0 for s in rm.STATUSES}
        counts.update(got)
        cells = [{"steps": ["derive_products:ok"]}] if not dry_run and sum(got.values()) else []
        return {"counts": counts, "cells_total": sum(counts.values()), "cells": cells,
                "_report_path": f"x/{star}_{target.replace(' ', '')}_latest.json"}
    monkeypatch.setattr(rm, "run", fake)
    return calls


# ── the element axis ─────────────────────────────────────────────────────────

def test_the_plan_is_every_canonical_target_exactly_once(rm, sw):
    plan = sw.sweep_plan()
    assert sorted(e.target for e in plan) == sorted(e["symbol"] for e in _master(rm))
    assert len(plan) == len(_master(rm))


def test_the_sweep_begins_fe_i_then_fe_ii(sw):
    labels = [e.label for e in sw.sweep_plan()]
    assert labels[:2] == ["Fe I", "Fe II"], labels[:4]


def test_after_fe_the_order_is_the_master_files_own(rm, sw):
    rest = [e.target for e in sw.sweep_plan() if e.symbol not in sw.FE_FIRST]
    assert rest == [e["symbol"] for e in _master(rm)
                    if e["symbol"].split()[0] not in sw.FE_FIRST]


def test_bare_fe_does_not_rebuild_the_fe_ii_cells(sw):
    """`run("Fe")` expands EVERY graded ion; `Fe II` is its own master target. Run
    unfiltered, every Fe II cell would be built twice under two report names."""
    fe = next(e for e in sw.sweep_plan() if e.target == "Fe")
    assert "II" not in fe.ions and fe.ions, fe


def test_a_bare_entry_with_no_unclaimed_ion_is_reported_not_run_unfiltered(
        rm, sw, monkeypatch, tmp_path):
    calls = _stub_run(rm, monkeypatch, lambda t, i: {rm.SKIP: 1})
    monkeypatch.setattr(rm, "canonical_elements", lambda: ["Fe", "Fe II"])
    monkeypatch.setattr(rm, "graded_ions", lambda sym: ["II"])
    doc = sw.run_sweep("solar", report_dir=tmp_path, echo=False)
    fe = next(r for r in doc["elements"] if r["element"] == "Fe")
    assert fe["status"] == sw.MATRIX_ERROR and "second time" in fe["reason"]
    assert calls == [("Fe II", ())]


def test_no_literal_element_list_in_the_sweep_or_the_driver(rm):
    """RYA-914: FE_FIRST is the ONE symbol literal this ticket may add."""
    symbols = {e["symbol"].split()[0] for e in _master(rm)}
    for path in (ROOT / "pipeline" / "run_sweep.py", ROOT / "run_pipeline.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
                hits = [e.value for e in node.elts
                        if isinstance(e, ast.Constant) and e.value in symbols]
                assert len(hits) <= 1, f"{path.name}:{node.lineno} lists elements {hits}"


# ── the Fe gate ──────────────────────────────────────────────────────────────

def test_fe_with_no_usable_cell_halts_the_sweep_and_names_every_other_element(
        rm, sw, monkeypatch, tmp_path):
    calls = _stub_run(rm, monkeypatch, lambda t, i: {rm.BLOCKED: 3, rm.UNPUBLISHED: 2})
    doc = sw.run_sweep("solar", report_dir=tmp_path, echo=False)
    assert doc["fe_gate"] == "HALT" and doc["fe_gate_enforced"]
    assert [c[0] for c in calls] == ["Fe", "Fe II"], "the sweep ran past a halted gate"
    rest = [r for r in doc["elements"] if r["symbol"] not in sw.FE_FIRST]
    assert rest and all(r["status"] == sw.NOT_RUN_FE_GATE for r in rest)
    assert all(r["reason"] == sw.FE_GATE_REASON for r in rest)
    assert len(doc["elements"]) == len(_master(rm))
    assert sw.failed(doc)


def test_one_usable_fe_cell_passes_the_gate(rm, sw, monkeypatch, tmp_path):
    _stub_run(rm, monkeypatch,
              lambda t, i: {rm.SKIP: 1, rm.BLOCKED: 4} if t == "Fe II" else {rm.BLOCKED: 5})
    doc = sw.run_sweep("solar", report_dir=tmp_path, echo=False)
    assert doc["fe_gate"] == "PASS"
    assert all(r["status"] == sw.RAN for r in doc["elements"])


def test_unpublished_is_not_usable(rm, sw, monkeypatch, tmp_path):
    """Stages exiting 0 with nothing published is NOT a cell other elements can lean on."""
    _stub_run(rm, monkeypatch, lambda t, i: {rm.UNPUBLISHED: 9})
    assert sw.run_sweep("solar", report_dir=tmp_path, echo=False)["fe_gate"] == "HALT"


def test_a_dry_run_reports_the_gate_but_still_shows_every_matrix(
        rm, sw, monkeypatch, tmp_path):
    calls = _stub_run(rm, monkeypatch, lambda t, i: {rm.WOULD_RUN: 2})
    doc = sw.run_sweep("solar", dry_run=True, report_dir=tmp_path, echo=False)
    assert doc["fe_gate"] == "HALT" and not doc["fe_gate_enforced"]
    assert len(calls) == len(_master(rm))
    rest = [r for r in doc["elements"] if r["symbol"] not in sw.FE_FIRST]
    assert all(r.get("fe_gate_would_block") for r in rest)
    assert doc["executed"] == 0


# ── loud-fail-continue ───────────────────────────────────────────────────────

def test_one_element_failing_or_crashing_does_not_stop_the_sweep(
        rm, sw, monkeypatch, tmp_path, capsys):
    def verdict(t, i):
        if t == "Si":
            raise RuntimeError("boom")
        if t == "C":
            raise rm.MatrixError("no graded line")
        return {rm.FAILED: 1} if t == "O" else {rm.SKIP: 1}
    _stub_run(rm, monkeypatch, verdict)
    doc = sw.run_sweep("solar", report_dir=tmp_path, echo=False)
    by = {r["element"]: r for r in doc["elements"]}
    assert by["Si"]["status"] == sw.CRASHED and "boom" in by["Si"]["reason"]
    assert by["C"]["status"] == sw.MATRIX_ERROR
    assert by["O"]["status"] == sw.RAN and by["O"]["counts"][rm.FAILED] == 1
    assert by[_master(rm)[-1]["symbol"]]["status"] == sw.RAN, "the sweep stopped early"
    assert doc["totals"][rm.FAILED] == 1 and sw.failed(doc)
    assert "CRASHED" in capsys.readouterr().err


# ── the report ───────────────────────────────────────────────────────────────

def test_the_sweep_report_is_written_with_the_ticket_schema(rm, sw, monkeypatch, tmp_path):
    _stub_run(rm, monkeypatch, lambda t, i: {rm.SKIP: 2, rm.BLOCKED: 1})
    doc = sw.run_sweep("solar", report_dir=tmp_path, echo=False)
    latest = json.loads((tmp_path / "solar_SWEEP_latest.json").read_text())
    assert len(list(tmp_path.glob("solar_SWEEP_2*Z.json"))) == 1
    for k in ("star", "generated_at", "code_commit", "fe_gate", "elements", "totals"):
        assert k in latest
    for r in latest["elements"]:
        assert {"element", "status", "report", "counts"} <= set(r)
        assert set(r["counts"]) == set(rm.STATUSES)
    assert latest["totals"][rm.SKIP] == 2 * len(_master(rm))
    assert sum(latest["totals"].values()) == latest["cells_total"]
    assert not sw.failed(doc)


def test_the_timestamped_sweep_report_is_gitignored():
    import fnmatch
    pats = [l for l in (ROOT / "data/results/orchestrator/.gitignore").read_text().splitlines()
            if l and not l.startswith("#")]
    assert any(fnmatch.fnmatch("solar_SWEEP_20261001T120000Z.json", p) for p in pats)
    assert not any(fnmatch.fnmatch("solar_SWEEP_latest.json", p) for p in pats)


def test_stdout_is_one_line_per_element_plus_a_verdict(rm, sw, monkeypatch, tmp_path, capsys):
    _stub_run(rm, monkeypatch, lambda t, i: {rm.SKIP: 1})
    sw.run_sweep("solar", report_dir=tmp_path)
    out = capsys.readouterr().out
    assert out.isascii()
    lines = [l for l in out.splitlines() if l.split()[:1] and "DONE" in l and "SKIP" in l
             and not l.startswith("VERDICT")]
    assert len(lines) == len(_master(rm))
    assert lines[0].startswith("Fe I ") and lines[1].startswith("Fe II")
    assert any(l.startswith("VERDICT") for l in out.splitlines())


# ── the driver ───────────────────────────────────────────────────────────────

def _main(monkeypatch, *argv):
    import run_pipeline
    monkeypatch.setattr(sys, "argv", ["run_pipeline.py", *argv])
    return run_pipeline.main


def test_element_and_all_elements_are_mutually_exclusive(monkeypatch):
    with pytest.raises(SystemExit) as e:
        _main(monkeypatch, "--star", "solar", "--element", "Si", "--all-elements")()
    assert e.value.code == 2


def test_ion_with_all_elements_is_refused(monkeypatch, capsys):
    with pytest.raises(SystemExit):
        _main(monkeypatch, "--star", "solar", "--all-elements", "--ion", "I")()
    assert "--ion is per-element" in capsys.readouterr().err


def test_engine_paths_without_a_matrix_mode_are_refused(monkeypatch, capsys):
    with pytest.raises(SystemExit):
        _main(monkeypatch, "--star", "solar", "--interpreter", "/x/python")()
    assert "matrix-mode" in capsys.readouterr().err


def test_no_engine_env_warns_once_and_continues(monkeypatch, capsys):
    import argparse
    import run_pipeline
    for n in (*run_pipeline._INTERPRETER_ENV, *run_pipeline._ISPEC_ENV):
        monkeypatch.delenv(n, raising=False)
    got = run_pipeline._resolve_engine_env(argparse.Namespace(interpreter=None, ispec_dir=None))
    assert got == ("", "")
    err = capsys.readouterr().err
    assert err.count("WARNING") == 1 and "BLOCKED" in err


def test_the_ticket_env_name_wins_and_the_existing_one_still_works(monkeypatch):
    import argparse
    import run_pipeline
    ns = argparse.Namespace(interpreter=None, ispec_dir=None)
    monkeypatch.setenv("ISPEC_DIR", "/i")
    monkeypatch.delenv("CODEX_SYNTH_PYTHON", raising=False)
    monkeypatch.setenv("CODEX_INTERPRETER", "/old")
    assert run_pipeline._resolve_engine_env(ns) == ("/old", "/i")
    monkeypatch.setenv("CODEX_SYNTH_PYTHON", "/new")
    assert run_pipeline._resolve_engine_env(ns) == ("/new", "/i")
    ns.interpreter = "/flag"
    assert run_pipeline._resolve_engine_env(ns) == ("/flag", "/i")


def test_the_legacy_path_runs_unchanged_behind_its_banner(monkeypatch, capsys):
    """The whole-star chain still runs its own stages -- after the banner, not instead."""
    import run_pipeline
    ran = []

    def stage(name):
        return lambda star: ran.append(name)
    monkeypatch.setattr(run_pipeline, "_whole_star_stages", lambda: tuple(
        stage(n) for n in ("normalize", "fit", "params", "derive", "unc", "interp")))
    _main(monkeypatch, "--star", "solar", "--validate-only")()
    out = capsys.readouterr().out
    assert run_pipeline.LEGACY_BANNER in out
    assert out.index(run_pipeline.LEGACY_BANNER) < out.index("spectra_normalize")
    assert ran == ["normalize", "fit"]


def test_a_dry_run_counts_would_run_cells_that_already_have_a_product(rm, sw):
    """A projected HALT must not read as "Fe has nothing" when 160 products exist."""
    doc = {"cells": [{"status": rm.WOULD_RUN, "A": 7.45}, {"status": rm.WOULD_RUN, "A": None},
                     {"status": rm.BLOCKED, "A": 7.5}]}
    assert sw.would_run_published(doc) == 1


def test_an_element_with_no_graded_line_is_an_honest_empty_not_a_failure(
        rm, sw, monkeypatch, tmp_path):
    _stub_run(rm, monkeypatch, lambda t, i: {rm.SKIP: 1})
    monkeypatch.setattr(rm, "graded_ions", lambda sym: [] if sym == "Ni" else ["I", "II"])
    doc = sw.run_sweep("solar", report_dir=tmp_path, echo=False)
    ni = next(r for r in doc["elements"] if r["element"] == "Ni")
    assert ni["status"] == sw.NO_GRADED_POOL and "RYA-945" in ni["reason"]
    assert not sw.failed(doc)


def test_a_pinned_ion_the_pool_lacks_is_an_honest_empty(rm, sw, monkeypatch, tmp_path):
    calls = _stub_run(rm, monkeypatch, lambda t, i: {rm.SKIP: 1})
    monkeypatch.setattr(rm, "graded_ions", lambda sym: ["I"])
    doc = sw.run_sweep("solar", report_dir=tmp_path, echo=False)
    fe2 = next(r for r in doc["elements"] if r["element"] == "Fe II")
    assert fe2["status"] == sw.NO_GRADED_POOL
    assert ("Fe II", ()) not in calls
