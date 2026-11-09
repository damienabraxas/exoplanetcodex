"""RYA-1310: IAG is not used in the solar VIS band (Ryan 2026-10-28).

The orchestrator builds no IAG VIS cell, says so in its report, and keeps IAG's red-optical
and NIR cells. The live feeds hold no IAG VIS product; the withdrawn ones sit in `archive`.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="module", autouse=True)
def _kp(tmp_path_factory):
    """`measure_band_ew` resolves the Kitt Peak atlas at import (as test_rya1222_run_matrix)."""
    kp = tmp_path_factory.mktemp("kp")
    (kp / "lm0296").touch()
    os.environ.setdefault("CODEX_KP_ATLAS", str(kp))


def test_no_iag_vis_cell_and_the_other_iag_bands_stay():
    from pipeline import run_matrix as rm
    cells = rm.expand("solar", "Mg", engines=["ts-lte"])
    iag = {d.band for d in cells if d.instrument == "iag_fts_solar_atlas"}
    assert "VIS" not in iag
    assert {"red-optical", "NIR"} <= iag
    assert any(d.band == "VIS" and d.instrument == "harps" for d in cells)


def test_the_report_names_the_exclusion():
    from pipeline import run_matrix as rm
    ex = rm.excluded_holding_bands("solar")
    assert {(e["holding"], e["band"]) for e in ex} == {("solar_iag", "VIS"),
                                                       ("solar_iag_reiners2016", "VIS")}
    assert all("RYA-1310" in e["reason"] for e in ex)
    src = (ROOT / "pipeline/run_matrix.py").read_text()
    assert '"excluded_holding_bands": excluded_holding_bands(star)' in src


def test_no_live_iag_vis_product_in_any_solar_feed():
    for f in sorted((ROOT / "data/products/solar").glob("*.json")):
        doc = json.loads(f.read_text())
        live = [p for p in doc["products"]
                if p.get("instrument") == "iag_fts_solar_atlas" and p.get("band") == "VIS"]
        assert not live, f"{f.name}: {len(live)} IAG VIS products still live"


def test_the_withdrawn_products_are_archived_with_the_reason():
    n = 0
    for f in sorted((ROOT / "data/products/solar").glob("*.json")):
        n += sum(1 for p in json.loads(f.read_text()).get("archive") or []
                 if "RYA-1310" in (p.get("superseded_reason") or ""))
    assert n == 97        # Fe 33, Mg 26, Na 25, Si 6, Al 4, C 3


def test_the_withdrawal_is_idempotent():
    out = subprocess.run([sys.executable, str(ROOT / "scripts/rya1310_withdraw_iag_vis.py")],
                         capture_output=True, text=True, check=True).stdout
    assert out.strip().endswith("total 0 (report only; --apply to withdraw)")
