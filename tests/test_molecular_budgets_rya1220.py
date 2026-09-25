"""RYA-1220 -- the molecular budgets and the ratified policy.

Pins the decisions that were expensive to establish, so a later pass cannot quietly
re-promote a rejected route or bury the start-dependence back inside sigma_fit.
"""
from __future__ import annotations

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
BUDGETS = ROOT / "data/output/rya1220/molecular_route_budgets.json"
POLICY = ROOT / "data/output/rya1220/cno_method_policy.json"

pytestmark = pytest.mark.skipif(not BUDGETS.exists(), reason="budgets not generated")


def _routes():
    """Keyed by (diagnostic, instrument) -- fine where one run per pair."""
    return {(r["diagnostic"], r["instrument"]): r
            for r in json.loads(BUDGETS.read_text())["routes"]}


def _by_route():
    """Keyed by route directory: IAG runs THREE CN_AX_IR configurations and they
    carry different legs, so a (diagnostic, instrument) key silently keeps only one."""
    return {r["route"]: r for r in json.loads(BUDGETS.read_text())["routes"]}


def test_iag_start_dependence_is_priced_not_hidden_in_sigma_fit():
    """Two runs of the same diagnostic on the same holding land 0.225 dex apart."""
    r = _by_route()["cn_iag_3d_anchor_v1"]
    m = r["components"]["measurement"]
    assert m["state"] == "MEASURED_BOUNDED"
    assert m["sigma_dex"] == pytest.approx(0.225, abs=1e-3)
    assert r["fit_start_dependence_dex"] == pytest.approx(0.225, abs=1e-3)


def test_one_sided_continuum_is_not_halved_into_a_fake_symmetry():
    r = _routes()[("CN_red", "HARPS")]
    c = r["components"]["continuum"]
    assert c["state"] == "MEASURED_ONE_SIDED"
    assert c["sigma_dex"] == pytest.approx(0.849, abs=1e-3)
    assert "ONE-SIDED" in c["note"]


def test_unmeasured_rho_is_bounded_not_held():
    r = _by_route()["cn_iag_3d_anchor_v1"]        # the run that carries the C/O legs
    c = r["components"]["molecular_coupling"]
    assert c["state"] == "MEASURED_BOUNDED"
    assert c["sigma_dex"] > 0
    assert "worst rho" in c["note"]


def test_within_holding_scatter_is_not_charged_to_the_instrument():
    """holding_instrument uses per-holding means, so run scatter is not double counted."""
    r = _routes()[("CN_AX_IR", "Kitt Peak")]
    h = r["components"]["holding_instrument"]
    assert "per-holding mean" in h["note"]
    assert h["sigma_dex"] == pytest.approx(0.407, abs=1e-3)


def test_single_holding_routes_hold_repeatability_with_a_reason():
    for key in [("CN_AX_J", "CRIRES+"), ("OH_H", "CRIRES+"), ("CO_K", "CRIRES+")]:
        h = _routes()[key]["components"]["holding_instrument"]
        assert h["sigma_dex"] is None
        assert "a second is needed" in h["note"]


def test_policy_rejects_the_two_offband_routes_with_cited_reasons():
    pol = json.loads(POLICY.read_text())
    rejected = {r["diagnostic"].split()[0]: r for r in pol["molecular_N"]["rejected"]}
    assert set(rejected) == {"CN_red", "NH"}
    for r in rejected.values():
        assert r["verdict"].startswith("REJECTED")
        assert len(r["why"]) > 80, "a rejection must carry its reason"
        assert r["retain_as"].startswith("diagnostic only")


def test_policy_keeps_kitt_peak_as_solar_only():
    pol = json.loads(POLICY.read_text())
    assert pol["molecular_N"]["secondary"]["status"] == "SOLAR CROSS-CHECK ONLY"


def test_the_literature_tension_is_recorded_as_irreducible():
    pol = json.loads(POLICY.read_text())
    assert pol["irreducible"]["atomic_vs_molecular_N_dex"] == 0.13


def test_the_route_with_no_co_legs_holds_coupling_with_a_reason():
    """cn_iag_stellar_v2 ran teff/xi legs, not C/O; it must say so, not borrow."""
    c = _by_route()["cn_iag_stellar_v2"]["components"]["molecular_coupling"]
    assert c["sigma_dex"] is None
    assert "did not run on this route" in c["note"]
