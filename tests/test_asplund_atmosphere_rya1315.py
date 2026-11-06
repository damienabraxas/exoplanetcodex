"""RYA-1315: Asplund+2021's atmosphere terms come from Scott+2015's per-line 3D/<3D>/HM74 table."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT)]

import rya1230_cno_budget as b  # noqa: E402
import rya1230_cno_budget_legs as legs  # noqa: E402


def test_na_pool_terms_are_half_the_scott_differences():
    got = b.scott_atmosphere("Na", "I", [4751.822, 5148.838, 6154.225, 6160.747])
    assert got["n_scott_lines"] == 4 and got["basis"] == "pool lines"
    # Scott+2015a Table 2: <3D> - HM = -0.042,-0.042,-0.041,-0.041; 3D - <3D> = -0.022,-0.022,-0.020,-0.016
    assert abs(got["temperature_structure_dex"] - 0.02075) < 1e-9
    assert abs(got["inhomogeneity_dex"] - 0.01) < 1e-9


def test_pool_without_scott_lines_takes_the_species_value_and_unknown_species_is_none():
    assert b.scott_atmosphere("Ti", "II", [9999.0])["basis"].startswith("species mean")
    assert b.scott_atmosphere("C", "I", [5380.3]) is None


def test_no_rya587_leg_is_run_by_default():
    assert legs.DEFAULT_LEGS == ("nominal",)
    assert legs.FE_EXTRA_LEGS == ("cscale",)


def test_restamp_replaces_only_model_atmosphere_and_passes_the_gate():
    import json
    import subprocess
    import rya1315_asplund_restamp as rs
    from pipeline.uncertainty_contract import publication_problems
    feed = json.loads(subprocess.run(["git", "show", "origin/main:data/products/solar/Al.json"],
                                     capture_output=True, text=True, cwd=ROOT).stdout or "{}")
    if not feed.get("products"):
        import pytest
        pytest.skip("origin/main Al feed not available")
    p = json.loads(json.dumps(feed["products"][0]))
    before = {c["name"]: c for c in p["uncertainty"]["components"]}
    assert rs.restamp(p, "Al") is None
    after = {c["name"]: c for c in p["uncertainty"]["components"]}
    assert "Scott+2015" in after["model_atmosphere"]["source"]
    assert {k: v for k, v in before.items() if k != "model_atmosphere"} == \
        {k: v for k, v in after.items() if k != "model_atmosphere"}
    assert p["A"] == feed["products"][0]["A"]
    assert publication_problems(p) == []
