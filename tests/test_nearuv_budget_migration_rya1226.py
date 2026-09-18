"""RYA-1226 part A -- what the legacy near-UV Fe rows can and cannot honestly carry.

🔴 WHAT THIS PROTECTS. RYA-587 exists to stop a placeholder being published as a
measurement. Migrating a legacy row is therefore only done with evidence that already
exists in an artifact -- and where none exists the component must stay HOLD and the
product must stay unpublishable. These tests pin both halves: the evidence that IS real,
and the refusal that must not be argued away.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

pd = pytest.importorskip("pandas")

BANDP = ROOT / "data/results/band_products"
MIGRATION = ROOT / "data/results/rya1226/nearuv_budget_migration.json"
RYA1190 = ROOT / "data/results/rya1190/rya1190_frontier_opacity.json"

KUR = "solar_kpno_kurucz2005_corrected"
MOL = "solar_kpno_molecfit_corrected"
DH19_DOI = "10.3847/1538-4365/ab322e"


def _lines(ion, holding, tier, treat):
    return pd.read_csv(BANDP / f"Fe{ion}_3000_3780_kpno_solar_atlas_{holding}_SYNTH_"
                               f"{tier}_{treat}_lines.csv")


@pytest.fixture(scope="module")
def migration():
    return json.loads(MIGRATION.read_text())


def test_sigma_stat_is_raw_scatter_over_sqrt_N_and_reproduces_the_published_value():
    """The contract's measurement identity holds on the real artifact, not on a fixture."""
    d = _lines("II", KUR, "DEEPGRADED", "1D-LTE")
    a = d.loc[d["in_aggregate"] == True, "abundance"].astype(float)  # noqa: E712
    se = a.std(ddof=1) / math.sqrt(len(a))
    assert len(a) == 12
    assert se == pytest.approx(0.0523, abs=5e-5), se


def test_every_line_carries_a_published_gf_sigma_from_ONE_common_source():
    """🔴 One source means the gf error does NOT average down.

    This is the fact that makes the contract's number differ from the legacy budget's
    RMS convention, so it must not quietly become 'independent' later.
    """
    g = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    g = g[g["species"] == "Fe II"]
    d = _lines("II", KUR, "DEEPGRADED", "1D-LTE")
    j = d.merge(g, on="wavelength_air_A", how="left")
    assert j["gf_sigma_dex"].notna().all()
    assert set(j["gf_tier"].dropna()) == {"LAB"}
    assert set(j["gf_source_doi"].dropna()) == {DH19_DOI}, "more than one gf source"


def test_the_gf_join_is_EP_consistent_not_wavelength_only():
    """RYA-1037: a lambda-only join is how a blend inherits another line's data."""
    g = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    g = g[g["species"] == "Fe II"]
    d = _lines("II", KUR, "DEEPGRADED", "1D-LTE")
    j = d.merge(g, on="wavelength_air_A", how="left")
    assert ((j["ep_eV"] - j["excitation_potential_eV"]).abs() < 0.005).all()


def test_hfs_is_N_A_on_evidence_not_on_convenience():
    g = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    g = g[g["species"] == "Fe II"]
    d = _lines("II", KUR, "DEEPGRADED", "1D-LTE")
    j = d.merge(g, on="wavelength_air_A", how="left")
    assert set(j["hfs_n_components"].dropna()) == {1}


@pytest.mark.parametrize("holding", [KUR, MOL])
def test_the_contract_refuses_the_fe2_rows_and_names_the_same_four_components(migration, holding):
    """🔴 THE REFUSAL IS THE RESULT. Twelve of sixteen components are real; four are not,
    and the product must stay unpublishable until they are measured."""
    row = next(r for r in migration["rows"]
               if r["ion"] == "II" and r["holding"] == holding and r["treatment"] == "1D-LTE")
    assert row["components_supplied"] == 12
    assert set(row["HOLDS"]) == {"blends", "continuum", "model_atmosphere", "profile_ew"}
    assert row["validator"], "validate() must refuse a document carrying HOLDs"
    assert row["sigma_reported_if_complete"] is None


def test_a_hold_is_the_contracts_verdict_not_a_hardcoded_list():
    """Control: supply one of the held components and the contract drops it from HOLDS.

    Without this, the four names above could be a constant someone typed.
    """
    import rya1226_migrate_nearuv_budget as M
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(x for x in feed["products"] if x["band"] == "near-UV" and x["ion"] == "II"
             and x["tier"] == "DEEPGRADED" and x["treatment"] == "1D-LTE"
             and x["holding"] == KUR)
    base = M.build(p, xi_slope=-0.1100)
    assert "continuum" in base["budget"]["holds"]

    ev = M.evidence_for(p, xi_slope=-0.1100, delta_xi=0.2912)
    ev["components"].append(dict(name="continuum", sigma_dex=0.01, state="DEFINED",
                                 source="TEST ONLY -- not a measurement", evidence={}))
    import numpy as np
    numeric = [c for c in ev["components"] if c["state"] in {"MEASURED", "DEFINED"}]
    doc = M.assemble(base["scope"], ev["components"],
                     covariance=np.diag([c["sigma_dex"] ** 2 for c in numeric]).tolist(),
                     covariance_source="test", assumptions="test")
    assert "continuum" not in doc["holds"]
    assert set(doc["holds"]) == {"blends", "model_atmosphere", "profile_ew"}


def test_blends_cannot_be_priced_because_RYA1190_refused_to_price_it():
    """The near-UV deficit is measured, real, and explicitly NOT catalogueable.

    Its ~0.028 'payoff' is offered as an ORDER, not a number to plan on -- so using it as
    the blends term would be inventing exactly what that ticket declined to supply.
    """
    d = json.loads(RYA1190.read_text())
    verdict = d["verdict_per_band"]["near-UV"]
    assert "OPACITY-DOMINATED" in verdict
    assert "NOT catalogueable" in verdict
    assert "not a number to plan on" in d["part_A_payoff_estimate"]["reading"]


def test_published_A_is_untouched_by_the_migration_survey():
    """The survey is read-only: it prices uncertainty and never an abundance."""
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    live = {(p["ion"], p["holding"], p["treatment"]): p["A"] for p in feed["products"]
            if p.get("band") == "near-UV" and p.get("tier") == "DEEPGRADED"}
    for r in json.loads(MIGRATION.read_text())["rows"]:
        if "error" in r:
            continue
        assert live[(r["ion"], r["holding"], r["treatment"])] == r["A"]
