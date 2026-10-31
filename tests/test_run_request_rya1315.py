"""RYA-1315: one run request, two modes -- what the CLI and a later front end both submit."""
import json

import pytest

from pipeline import run_request as rr


def test_full_mode_runs_everything_and_takes_no_selection():
    assert rr.RunRequest("solar", "Fe").problems() == []
    why = rr.RunRequest("solar", "Fe", bands=("VIS",)).problems()
    assert any("needs mode 'select'" in w for w in why)
    assert any("never the whole element" in w
               for w in rr.RunRequest("solar", "Fe", force=True).problems())


def test_select_mode_needs_a_selection():
    assert rr.RunRequest("solar", "Fe", mode="select").problems()
    assert rr.RunRequest("solar", "Fe", mode="select", bands=("VIS",), force=True).problems() == []


def test_a_dry_run_cannot_publish():
    why = rr.RunRequest("solar", "Fe", dry_run=True, publish=True).problems()
    assert any("cannot publish" in w for w in why)


def test_a_selection_is_checked_against_the_options():
    opts = {"bands": ["VIS", "near-UV"], "instruments": ["harps"],
            "holdings": ["solar_harps_molecfit_corrected"], "engines": ["ts-lte"]}
    ok = rr.RunRequest("solar", "Fe", mode="select", bands=("VIS",),
                       instruments=("solar_harps_molecfit_corrected",))
    assert ok.problems(opts) == []
    bad = rr.RunRequest("solar", "Fe", mode="select", bands=("UV",), engines=("nope",))
    why = bad.problems(opts)
    assert any("bands" in w for w in why) and any("engines" in w for w in why)


def test_round_trips_as_json(tmp_path):
    req = rr.RunRequest("solar", "Na", mode="select", bands=("VIS",), engines=("ts-lte",),
                        force=True, publish=True)
    p = tmp_path / "req.json"
    p.write_text(json.dumps(req.to_dict()))
    assert rr.RunRequest.from_json(p) == req
    assert rr.RunRequest.from_dict({"star": "solar", "element": "Na", "mode": "select",
                                    "bands": "VIS"}).bands == ("VIS",)
    with pytest.raises(ValueError, match="unknown field"):
        rr.RunRequest.from_dict({"star": "solar", "element": "Na", "colour": "red"})


def test_execute_refuses_before_dispatching(monkeypatch):
    from pipeline import run_matrix as rm
    monkeypatch.setattr(rm, "run", lambda *a, **k: pytest.fail("dispatched a refused request"))
    with pytest.raises(ValueError, match="refused"):
        rr.execute(rr.RunRequest("solar", "Fe", mode="select"), check_options=False)


def test_execute_passes_the_selection_through(monkeypatch):
    from pipeline import run_matrix as rm
    seen = {}
    monkeypatch.setattr(rm, "run", lambda star, el, **k: seen.update(k) or {"cells": []})
    out = rr.execute(rr.RunRequest("solar", "Fe", mode="select", bands=("VIS",),
                                   routes=("synthesis",), force=True), check_options=False)
    assert seen["bands"] == ["VIS"] and seen["methods"] == ["synthesis"] and seen["force"]
    assert out["publish"] is None and out["report"]["request"]["mode"] == "select"
