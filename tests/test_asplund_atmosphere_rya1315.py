"""RYA-1315: Asplund+2021's atmosphere terms come from Scott+2015's per-line 3D/<3D>/HM74 table."""
import pytest
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
    assert b.scott_atmosphere("C", "I", [5380.3])["basis"].startswith("Asplund+2021 Table 3")
    assert b.scott_atmosphere("Xx", "I", [5380.3]) is None


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


def _with_transition_na(evidence):
    import numpy as np
    from pipeline import uncertainty_contract as uc
    from tests.test_uncertainty_contract_rya587 import complete_product
    p = complete_product()
    u = p["uncertainty"]
    comps = [c for c in u["components"] if c["name"] != "transition_data"]
    comps.append(uc.component("transition_data", None, state="N/A",
                              source="not in the Asplund+2021 Sect. 2.1 budget", evidence=evidence))
    numeric = [c for c in comps if c["state"] in {"MEASURED", "DEFINED"}]
    p["uncertainty"] = uc.assemble(u["scope"], comps,
                                   covariance=np.diag([c["sigma_dex"] ** 2 for c in numeric]).tolist(),
                                   covariance_source=u["covariance_source"], assumptions=u["assumptions"])
    p["sigma_reported"] = p["uncertainty"]["sigma_reported"]
    return p


def test_unpriced_gf_line_is_na_on_an_asplund_budget_not_a_hold():
    """RYA-1315: Asplund+2021 leaves atomic data out; a line with no gf sigma no longer holds."""
    from pipeline import uncertainty_contract as uc
    ids = ["physical-transition-1", "physical-transition-2"]
    pool = uc.pool_digest(ids)
    ok = {"asplund_budget": True, "pool_sha256": pool, "indicator_ids": ids,
          "unpriced_ids": ids[1:]}
    assert uc.publication_problems(_with_transition_na(ok)) == []
    assert uc.publication_problems(_with_transition_na({**ok, "pool_sha256": "other"}))
    assert uc.publication_problems(_with_transition_na({**ok, "asplund_budget": None}))


def test_atomic_cno_take_asplund2021_table3():
    """RYA-1315: Scott+2015 has no C/N/O; Asplund+2021 Table 3's own <3D>/HM/3D rows price them."""
    from pipeline.asplund_atmosphere import scott_atmosphere
    n = scott_atmosphere("N", "I", [8629.24])
    assert n["temperature_structure_dex"] == pytest.approx((7.864 - 7.816) / 2)
    assert n["inhomogeneity_dex"] == pytest.approx((7.816 - 7.767) / 2)
    o = scott_atmosphere("O", "I", [6300.304])              # [O I] takes its own row
    assert o["temperature_structure_dex"] == pytest.approx((8.705 - 8.688) / 2)
    assert "[O I]" in o["basis"]
    assert scott_atmosphere("O", "II", [7000.0]) is None
