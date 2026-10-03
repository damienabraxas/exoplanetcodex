"""RYA-1234: every element run ends with its literature check (step 12) and its problem
lines (step 13). Decides no abundance; every band comes from the litscan YAML."""
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
    kp = tmp_path_factory.mktemp("kp")
    (kp / "lm0296").touch()
    os.environ.setdefault("CODEX_KP_ATLAS", str(kp))


@pytest.fixture(scope="module")
def ev():
    from pipeline import element_verdict
    return element_verdict


class _Lit:
    central, min, max, scale = 7.46, 7.42, 7.50, "3D-NLTE"
    def contains(self, v): return self.min <= v <= self.max
    def offset(self, v): return v - self.central
    def is_deviate(self, v): return abs(v - self.central) > 0.08


def test_in_and_out_of_band_with_signed_delta(ev):
    p = {"treatment": "ENGINE-A-3DNLTE", "route": "SYNTH", "band": "VIS", "A": 7.47}
    assert ev.judge(p, _Lit(), "", True)["status"] == ev.IN_BAND
    hi = ev.judge({**p, "A": 7.56}, _Lit(), "", True)
    assert hi["status"] == ev.OUT_OF_BAND and hi["delta"] == pytest.approx(0.06)
    lo = ev.judge({**p, "A": 7.40}, _Lit(), "", True)
    assert lo["delta"] == pytest.approx(-0.02)


def test_a_scale_mismatch_is_flagged_never_hidden(ev):
    j = ev.judge({"treatment": "1D-LTE", "route": "SYNTH", "band": "VIS", "A": 7.47},
                 _Lit(), "", True)
    assert "scale_mismatch" in j and "1D-LTE" in j["scale_mismatch"]


def test_no_litscan_is_no_reference_and_a_non_benchmark_is_reported(ev):
    p = {"treatment": "1D-LTE", "route": "SYNTH", "band": "VIS", "A": 7.5}
    assert ev.judge(p, None, "no litscan", True)["status"] == ev.NO_REFERENCE
    assert ev.judge(p, _Lit(), "", False)["status"] == ev.REPORTED


def test_fe_ii_is_not_judged_against_the_fe_i_litscan(ev):
    lit, why = ev.literature_for("Fe", "II")
    assert lit is None and "Fe I" in why


def test_the_sun_is_the_benchmark_anchor(ev):
    assert ev.is_benchmark("solar")[0] is True


def test_reason_codes_map_the_products_own_exclusions(ev):
    assert ev.reason_code("CHI2-GATE: red_chi2 118 >= 10") == "FIT_INVALID"
    assert ev.reason_code("ENGINE-A-NOT-SERVED: no delta") == "NOT_SERVED"
    assert ev.reason_code("QUARANTINED-TELLURIC: inside O2") == "TELLURIC_SATURATED"
    assert ev.reason_code("REW -4.83 above the -4.9 saturation ceiling") == "EW_SATURATED"
    assert ev.reason_code("something new") == ev.OTHER


def test_no_literal_tolerance_in_the_module():
    """Every band comes from the litscan: no float literal may be a tolerance here."""
    tree = ast.parse((ROOT / "pipeline" / "element_verdict.py").read_text())
    floats = [n.value for n in ast.walk(tree)
              if isinstance(n, ast.Constant) and isinstance(n.value, float) and n.value != 0.0]
    assert floats == [], floats


def test_the_verdict_is_written_into_the_report_and_changes_no_value(ev, tmp_path, monkeypatch):
    from pipeline import run_matrix as rm
    cell = {"band": "VIS", "instrument": "harps", "holding": "solar_harps_molecfit_corrected",
            "ion": "I", "engine": "ts-lte", "route": "SYNTH", "pool": "REFERENCE",
            "status": rm.WOULD_RUN}
    doc = {"generated_at": "2026-10-03T00:00:00Z", "cells": [cell]}
    (tmp_path / "solar_Fe_latest.json").write_text(json.dumps(doc))
    prod = {"element": "Fe", "ion": "I", "band": "VIS", "instrument": "harps",
            "holding": "solar_harps_molecfit_corrected", "route": "SYNTH",
            "selector": "REFERENCE", "treatment": "1D-LTE", "A": 7.60}
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": [dict(prod)]})
    monkeypatch.setattr(rm, "load_ledger", lambda: {})
    out = ev.verdict("solar", "Fe", report_dir=tmp_path)
    assert out["verdict"]["element_verdict"] == ev.REVIEW
    assert prod["A"] == 7.60, "the verdict touched a value"
    written = json.loads((tmp_path / "solar_Fe_latest.json").read_text())
    assert written["verdict"]["cells"][0]["status"] == ev.OUT_OF_BAND
    assert written["verdict"]["cells"][0]["A"] == 7.60


def test_no_products_is_incomplete_not_pass(ev, tmp_path, monkeypatch):
    from pipeline import run_matrix as rm
    (tmp_path / "solar_Si_latest.json").write_text(json.dumps(
        {"generated_at": "2026-10-03T00:00:00Z", "cells": []}))
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": []})
    monkeypatch.setattr(rm, "load_ledger", lambda: {})
    assert ev.verdict("solar", "Si", report_dir=tmp_path)["verdict"]["element_verdict"] \
        == ev.INCOMPLETE


def test_all_no_reference_is_incomplete_not_pass(ev, tmp_path, monkeypatch):
    """Fe II read PASS with zero comparisons: step 12 never happened."""
    from pipeline import run_matrix as rm
    cell = {"band": "VIS", "instrument": "harps", "holding": "solar_harps_molecfit_corrected",
            "ion": "II", "engine": "ts-lte", "route": "SYNTH", "pool": "REFERENCE",
            "status": rm.HELD}
    (tmp_path / "solar_FeII_latest.json").write_text(json.dumps(
        {"generated_at": "2026-10-03T00:00:00Z", "cells": [cell]}))
    prod = {"element": "Fe", "ion": "II", "band": "VIS", "instrument": "harps",
            "holding": "solar_harps_molecfit_corrected", "route": "SYNTH",
            "selector": "REFERENCE", "treatment": "1D-LTE", "A": 7.47}
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": [prod]})
    monkeypatch.setattr(rm, "load_ledger", lambda: {})
    v = ev.verdict("solar", "Fe II", report_dir=tmp_path)["verdict"]
    assert v["counts"]["NO_REFERENCE"] == 1 and v["element_verdict"] == ev.INCOMPLETE


def test_step_5_is_per_ion(ev):
    from pipeline import run_matrix as rm
    rows = {r["ion"]: r for r in rm.process_steps("solar", "Fe", ["I", "II"]) if r["step"] == 5}
    assert rows["I"]["ok"] is True and rows["II"]["ok"] is False
    assert "covers Fe I, not Fe II" in rows["II"]["evidence"]
