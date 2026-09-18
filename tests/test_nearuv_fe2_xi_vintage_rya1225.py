"""RYA-1225 — the near-UV Fe II dA/dxi split, and the rule that closed it.

🔴 WHAT THIS PROTECTS. A dA/dxi belongs to a SYNTHESIS, not just to a line pool. RYA-1207
turned molecular opacity on for the near-UV band and nothing re-derived the derivatives
measured before it, so near-UV Fe II published two slopes 2.0-2.8x apart on a byte-identical
12-line pool. These tests pin the diagnosis, the substitution rule, and — the part most
worth guarding — the rule's REFUSAL to act on cells nobody adjudicated.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.xi_vintage import (  # noqa: E402
    ADJUDICATED, apply_vintage_supersession, certified_same_pool, is_stale)

A1168 = ROOT / "data/results/rya1168/nearuv_xi_dadxi.json"
A1213 = ROOT / "data/results/rya1213/reference_xi_dadxi.json"
LEDGER = ROOT / "data/audit/rya1213_reference_matrix/same_pool_ledger.json"
FEED = ROOT / "data/products/solar/Fe.json"
RCA = ROOT / "data/results/rya1225/nearuv_fe2_xi_rca.json"

NEARUV_FE2 = [("solar_kpno_kurucz2005_corrected", "1D-LTE"),
              ("solar_kpno_kurucz2005_corrected", "ENGINE-A"),
              ("solar_kpno_molecfit_corrected", "1D-LTE"),
              ("solar_kpno_molecfit_corrected", "ENGINE-A")]


@pytest.fixture(scope="module")
def feed():
    return json.loads(FEED.read_text())


@pytest.fixture(scope="module")
def ledger():
    return json.loads(LEDGER.read_text())


def _band_idx():
    """The SAME set of band runs the emitter loads.

    ⚠️ Loading only the two files this ticket is about would make the scope-control test
    pass vacuously: the stale red-optical and NIR entries it must decline to touch live in
    the OTHER two runs, so a narrow fixture would have nothing to decline.
    """
    import importlib
    sys.path.insert(0, str(ROOT / "scripts"))
    emit = importlib.import_module("rya1178_emit_fe_schema")
    idx = {}
    for path in emit.XI_BAND_RUNS:
        doc = json.loads(path.read_text())
        species = (doc.get("species") or "").split()
        for e in doc["pools"]:
            ion = e.get("ion") or (species[-1] if species else None)
            k = (ion, e["holding"], e["tier"], e["treatment"], e["band"])
            idx[k] = {**e, "_source": str(path.relative_to(ROOT)), "_ticket": doc.get("ticket")}
    return idx


def test_the_two_runs_really_do_disagree_2_to_2_8x():
    """The premise. If this ever goes quiet the ticket is moot and these tests are noise."""
    a = {(p["holding"], p["treatment"]): p for p in json.loads(A1168.read_text())["pools"]
         if p["ion"] == "II"}
    b = {(p["holding"], p["treatment"]): p for p in json.loads(A1213.read_text())["pools"]
         if p["ion"] == "II" and p["band"] == "near-UV"}
    ratios = [abs(b[k]["dA_dxi"] / a[k]["dA_dxi"]) for k in NEARUV_FE2]
    assert all(2.0 <= r <= 2.8 for r in ratios), ratios


def test_the_pool_is_byte_identical_so_case_a_is_refuted():
    """Not (a) POOL MEMBERSHIP — measured from the committed nominal per-line files.

    Equal `n_lines` is NOT this check. RYA-1204's 40-line and 58-line pools shared two
    lines; identity is a lambda+EP set comparison plus equal per-line abundances.
    """
    pd = pytest.importorskip("pandas")
    B = ROOT / "data/results/band_products"
    for holding, treat in NEARUV_FE2:
        stem = f"FeII_3000_3780_kpno_solar_atlas_{holding}_SYNTH"
        dg = pd.read_csv(B / f"{stem}_DEEPGRADED_{treat}_lines.csv")
        rf = pd.read_csv(B / f"{stem}_REFERENCE_{treat}_lines.csv")
        key = ["wavelength_air_A", "ep_eV"]
        sa = {(round(a, 4), round(b, 4)) for a, b in zip(dg[key[0]], dg[key[1]])}
        sb = {(round(a, 4), round(b, 4)) for a, b in zip(rf[key[0]], rf[key[1]])}
        assert sa == sb, f"{holding} {treat}: pools differ"
        m = dg.merge(rf, on=key, suffixes=("_dg", "_rf"))
        d = (m["abundance_rf"] - m["abundance_dg"]).dropna()
        assert d.empty or d.abs().max() < 1e-9


def test_rya1168_recorded_the_abundance_it_ran_against_and_it_is_superseded(feed):
    """The staleness signal is MEASURED, not asserted: entry A vs the product's A today."""
    live = {(p["ion"], p["holding"], p["tier"], p["treatment"], p["band"]): p
            for p in feed["products"]}
    for p in json.loads(A1168.read_text())["pools"]:
        k = (p["ion"], p["holding"], p["tier"], p["treatment"], "near-UV")
        assert is_stale(p, live[k]), f"{k} no longer reads as stale"


def test_an_entry_with_no_recorded_A_is_not_guessed_either_way():
    """RYA-1213's entries carry no `A`. Absence must not read as 'fresh' OR as 'stale'."""
    assert is_stale({"dA_dxi": -0.11}, {"A": 7.617}) is False
    assert is_stale({"A": 7.617}, None) is False


def test_same_pool_certification_requires_BOTH_zero_delta_and_same_n(ledger):
    """Equal n alone must never certify identity — that is the RYA-1204 trap."""
    cert = certified_same_pool(ledger)
    for holding, treat in NEARUV_FE2:
        assert ("II", holding, treat, "near-UV") in cert
    bad = {"rows": [{"ion": "Fe I", "holding": "h", "treatment": "t", "band": "b",
                     "delta_dex": 0.006, "same_n": True}]}
    assert certified_same_pool(bad) == set()


def test_supersession_fires_on_the_four_adjudicated_cells_and_nowhere_else(feed, ledger):
    """🔴 THE SCOPE CONTROL. The staleness rule finds stale red-optical and NIR entries
    too, whose cause this ticket did NOT diagnose. Acting on them would be re-pricing
    published bars nobody adjudicated."""
    _, report = apply_vintage_supersession(_band_idx(), feed["products"], ledger)
    done = [r for r in report if r["action"] == "SUPERSEDED"]
    assert len(done) == 4
    assert {(r["key"]["holding"], r["key"]["treatment"]) for r in done} == set(NEARUV_FE2)
    assert all(r["key"]["ion"] == "II" and r["key"]["band"] == "near-UV" for r in done)
    others = [r for r in report if r["action"] != "SUPERSEDED"]
    assert others, "the rule must still REPORT the cells it declines to touch"
    assert all(r["action"].startswith("FLAGGED") for r in others)


def test_widening_the_adjudicated_set_is_what_would_have_changed_other_bands(feed, ledger):
    """Mutation control: the narrow scope is load-bearing, not incidental."""
    wide = frozenset(ADJUDICATED | {("I", "red-optical"), ("I", "NIR")})
    _, report = apply_vintage_supersession(_band_idx(), feed["products"], ledger,
                                           adjudicated=wide)
    assert len([r for r in report if r["action"] == "SUPERSEDED"]) > 4


def test_after_supersession_near_uv_fe2_carries_ONE_slope_across_both_tiers(feed, ledger):
    """The ticket's own smoke test: the 2.0-2.8x split is closed."""
    idx, _ = apply_vintage_supersession(_band_idx(), feed["products"], ledger)
    for holding, treat in NEARUV_FE2:
        dg = idx[("II", holding, "DEEPGRADED", treat, "near-UV")]["dA_dxi"]
        rf = idx[("II", holding, "REFERENCE", treat, "near-UV")]["dA_dxi"]
        assert dg == rf, f"{holding} {treat}: {dg} != {rf}"


def test_the_superseded_entry_records_what_it_replaced_and_why(feed, ledger):
    """A silent substitution would be the defect wearing a fix's clothes."""
    idx, _ = apply_vintage_supersession(_band_idx(), feed["products"], ledger)
    s = idx[("II", NEARUV_FE2[0][0], "DEEPGRADED", NEARUV_FE2[0][1], "near-UV")]["_superseded"]
    assert s["replaced_ticket"] == "RYA-1168"
    assert s["replaced_dA_dxi"] == -0.0475
    assert s["donor_tier"] == "REFERENCE"
    assert "synthesis changed" in s["cause"]


def test_the_rca_artifact_classifies_the_case_with_evidence():
    """STEP 1's deliverable exists and says (c) with per-line evidence, not an adjective."""
    d = json.loads(RCA.read_text())
    assert d["case"] == "c"
    assert "RYA-1207" in d["physical_cause"]
    assert len(d["rows"]) == 4
    for r in d["rows"]:
        assert r["pool_identity_deepgraded_vs_reference"]["identical"] is True
        assert r["rya1168_measured_against_a_superseded_A"] is True
        assert len(r["reproduced_from_rya1168_own_legs"]["per_line"]) == 12


def test_the_refuted_band_confound_is_recorded_so_it_is_not_re_run():
    """The VIS-window control refutes the obvious hypothesis; keep it written down."""
    d = json.loads(RCA.read_text())
    assert "REFUTED" in d["refuted_hypothesis_band_confound"]
    vis = [r["vis_window_control_dA_dxi"] for r in d["rows"]]
    assert all(v is None or abs(v + 0.11) > 0.02 for v in vis), vis
