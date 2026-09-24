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
def test_the_fe2_rows_now_carry_a_COMPLETE_validated_budget(migration, holding):
    """Ryan's REDIRECT: three of the four I held were already measured. They are wired.

    16 of 16 components, no HOLDs, and validate() returns clean -- which is what makes the
    row publishable at all.
    """
    row = next(r for r in migration["rows"]
               if r["ion"] == "II" and r["holding"] == holding and r["treatment"] == "1D-LTE")
    assert row["components_supplied"] == 16
    assert row["HOLDS"] == []
    assert row["validator"] is None, row["validator"]
    assert row["sigma_reported_if_complete"] > 0


def test_each_wired_term_cites_a_real_prior_measurement():
    """🔴 'Unpriceable' with no cited reason re-buries measured work. Each must name its own.

    Ryan's 2026-09-19 correction moved the normalization term onto RYA-846's MEASUREMENT and
    withdrew the observed-spread floor that stood in for it, so what this pins has changed:
    pseudo_continuum is now measured, and `continuum` is N/A because pricing both would count
    normalization twice.
    """
    import rya1226_migrate_nearuv_budget as M
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(x for x in feed["products"] if x["band"] == "near-UV" and x["ion"] == "II"
             and x["tier"] == "DEEPGRADED" and x["treatment"] == "1D-LTE"
             and x["holding"] == KUR)
    terms = {c["name"]: c for c in M.build(p, xi_slope=-0.1100)["budget"]["components"]}

    #: pseudo_continuum -- RYA-846's measurement, read from its artifact and not restated.
    pc = terms["pseudo_continuum"]
    assert pc["state"] == "MEASURED" and not pc["evidence"].get("bound")
    assert "RYA-846" in pc["source"]
    rya846 = json.loads((ROOT / "data/results/rya846/rya846_sigma_delta.json").read_text())
    assert pc["sigma_dex"] == pytest.approx(
        rya846["decomposition"]["term_net_dex"], abs=1e-6), "must be the artifact's value"
    #: 🔴 THE NET TERM, NOT THE RAW ONE -- the raw figure double-charges Wallace's own
    #: internal scatter, and the two differ by 1.8x so a swap would be visible here.
    assert pc["sigma_dex"] != pytest.approx(rya846["decomposition"]["term_raw_dex"], abs=1e-6)
    assert pc["evidence"]["supersedes_assumed_dex"] == rya846["assumed_term_dex"] == 0.1
    assert pc["evidence"]["n_lines_outside_wallace"] == 0

    #: continuum -- N/A, and explicitly BECAUSE pseudo_continuum already carries it.
    cont = terms["continuum"]
    assert cont["state"] == "N/A" and cont["sigma_dex"] is None
    assert cont["evidence"]["not_double_counted"] is True
    assert "pseudo_continuum" in cont["evidence"]["subsumed_by"]

    #: model_atmosphere -- Ryan asked for the grid-shift lever to be VERIFIED, not assumed.
    atm = terms["model_atmosphere"]
    assert atm["sigma_dex"] == 0.004 and "RYA-1032" in atm["source"]
    assert atm["evidence"]["same_pool"] is True
    assert atm["evidence"]["grid_shift_lever_exists"] is True
    #: and the STAGGER step is recorded as the DIFFERENT axis it is, not quietly omitted.
    stg = atm["evidence"]["stagger_step_is_a_different_axis"]
    assert stg["varied"] == "dim" and stg["delta_dex"] == 0.101

    pew = terms["profile_ew"]
    assert pew["sigma_dex"] == 0.2160 and "RYA-1220" in pew["source"]
    assert pew["evidence"]["cross_element"] is True


def test_every_bound_is_flagged_as_a_bound_and_not_as_a_measurement():
    """A ceiling that reads like a measurement is the defect this ticket is about."""
    import rya1226_migrate_nearuv_budget as M
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(x for x in feed["products"] if x["band"] == "near-UV" and x["ion"] == "II"
             and x["tier"] == "DEEPGRADED" and x["treatment"] == "1D-LTE"
             and x["holding"] == KUR)
    terms = {c["name"]: c for c in M.build(p, xi_slope=-0.1100)["budget"]["components"]}
    #: After Ryan's 2026-09-19 correction only TWO terms are bounds. `continuum` and
    #: `blends` are N/A on evidence, and pseudo_continuum is a measurement -- so a bound
    #: creeping back onto any of those three fails here.
    for name in ("model_atmosphere", "profile_ew"):
        assert terms[name]["evidence"].get("bound") is True, name
    for name in ("measurement", "transition_data", "stellar.xi", "pseudo_continuum"):
        assert terms[name]["state"] == "MEASURED", name
        assert not terms[name]["evidence"].get("bound"), name
    for name in ("continuum", "blends"):
        assert terms[name]["state"] == "N/A", name
        assert terms[name]["sigma_dex"] is None, name
        assert not terms[name]["evidence"].get("bound"), name


def test_the_bar_is_dominated_by_bounds_and_that_is_stated_not_hidden():
    """🔴 85% of the variance is bounds, 45% is profile_ew alone. If that ever stops being
    true the reporting must change with it, so it is pinned rather than narrated once."""
    import rya1226_migrate_nearuv_budget as M
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(x for x in feed["products"] if x["band"] == "near-UV" and x["ion"] == "II"
             and x["tier"] == "DEEPGRADED" and x["treatment"] == "1D-LTE"
             and x["holding"] == KUR)
    b = M.build(p, xi_slope=-0.1100)["budget"]
    total = b["sigma_reported"] ** 2
    bound = sum(c["sigma_dex"] ** 2 for c in b["components"]
                if c["sigma_dex"] and c["evidence"].get("bound"))
    assert bound / total > 0.80
    pew = next(c for c in b["components"] if c["name"] == "profile_ew")
    assert pew["sigma_dex"] ** 2 / total > 0.40


def test_stellar_xi_carries_real_perturbation_evidence_and_a_measured_asymmetry():
    """A derived response is not a perturbation -- the contract says so and it is right."""
    import rya1226_migrate_nearuv_budget as M
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(x for x in feed["products"] if x["band"] == "near-UV" and x["ion"] == "I"
             and x["tier"] == "DEEPGRADED" and x["treatment"] == "1D-LTE"
             and x["holding"] == KUR)
    xi = next(c for c in M.build(p, xi_slope=-0.17)["budget"]["components"]
              if c["name"] == "stellar.xi")
    ev = xi["evidence"]
    assert xi["state"] == "MEASURED"
    assert ev["signed_response_dex"] == pytest.approx((ev["delta_plus_dex"] - ev["delta_minus_dex"]) / 2)
    assert abs(ev["signed_response_dex"]) == pytest.approx(xi["sigma_dex"])
    assert ev["delta_parameter"] == 0.2912
    assert "ASYMMETRIC" in ev["response_assessment"], "the measured asymmetry must be recorded"


def test_stellar_teff_is_a_sourced_zero_with_its_omission_sized():
    """Declaring it zero is only honest if the size of the omission is on the record."""
    import rya1226_migrate_nearuv_budget as M
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(x for x in feed["products"] if x["band"] == "near-UV" and x["ion"] == "II"
             and x["tier"] == "DEEPGRADED" and x["treatment"] == "1D-LTE"
             and x["holding"] == KUR)
    t = next(c for c in M.build(p, xi_slope=-0.1100)["budget"]["components"]
             if c["name"] == "stellar.teff")
    assert t["state"] == "DEFINED" and t["sigma_dex"] == 0.0
    assert t["evidence"]["derived_but_not_perturbed_dex"] == 0.000665


def test_blends_is_NOT_a_component_and_carries_no_bound():
    """🔴 Ryan, 2026-09-19, superseding his own earlier "bound it as a ceiling" directive:
    blends is NOT an open uncertainty component and must not be bounded.

    The previous pass charged 0.1721 dex -- 28.7% of the published variance -- for something
    that does not exist. RYA-1190's own verdict is that this band's catalogued opacity is
    already COMPLETE to the VALD threshold (there is nothing to catalogue, so no un-modelled
    resolvable blends), and the deficit it did find was an OPACITY deficit that RYA-1204
    corrected and RYA-1207 wired into the very synthesis these products are measured on.
    Both candidate numbers were PRE-correction measurements.
    """
    import rya1226_migrate_nearuv_budget as M
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(x for x in feed["products"] if x["band"] == "near-UV" and x["ion"] == "II"
             and x["tier"] == "DEEPGRADED" and x["treatment"] == "1D-LTE"
             and x["holding"] == KUR)
    bl = next(c for c in M.build(p, xi_slope=-0.1100)["budget"]["components"]
              if c["name"] == "blends")
    assert bl["state"] == "N/A"
    assert bl["sigma_dex"] is None, "a bound must not come back"
    assert not bl["evidence"].get("bound")
    #: both withdrawn numbers stay on the record WITH their reason, so the withdrawal is
    #: auditable rather than a silent deletion.
    assert bl["evidence"]["withdrawn_bound_dex"] == 0.1721
    assert bl["evidence"]["withdrawn_alternative_dex"] == 0.028
    assert bl["evidence"]["opacity_deficit_corrected_by"] == ["RYA-1204", "RYA-1207"]
    assert bl["evidence"]["synthesis_is_post_correction"] is True
    #: and the citation is checked against RYA-1190's artifact, not paraphrased from memory.
    d = json.loads(RYA1190.read_text())
    verdict = d["verdict_per_band"]["near-UV"]
    assert "OPACITY-DOMINATED" in verdict
    assert "already complete to the VALD threshold" in verdict


def test_published_A_is_untouched_by_the_migration_survey():
    """The survey is read-only: it prices uncertainty and never an abundance."""
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    live = {(p["ion"], p["holding"], p["treatment"]): p["A"] for p in feed["products"]
            if p.get("band") == "near-UV" and p.get("tier") == "DEEPGRADED"}
    for r in json.loads(MIGRATION.read_text())["rows"]:
        if "error" in r:
            continue
        assert live[(r["ion"], r["holding"], r["treatment"])] == r["A"]


def test_the_feed_actually_carries_the_migrated_budgets():
    """Part D: the rows publish. This is the smoke test the ticket asked for."""
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    with_block = [p for p in feed["products"] if "uncertainty" in p]
    assert len(with_block) == 6, "six near-UV rows migrated; the rest stay legacy"
    for p in with_block:
        assert p["band"] == "near-UV" and p["tier"] == "DEEPGRADED"
        assert p["star"] == "solar"
        assert p["sigma_reported"] == p["uncertainty"]["sigma_reported"]


def test_the_two_incomplete_rows_did_NOT_publish_a_budget():
    """The migration's whole point: an incomplete budget stays unpublishable."""
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    held = [p for p in feed["products"]
            if p.get("band") == "near-UV" and p.get("ion") == "I"
            and p.get("tier") == "DEEPGRADED" and p["holding"] == MOL
            and p["treatment"] in ("1D-LTE", "ENGINE-A")]
    assert len(held) == 2
    for p in held:
        assert "uncertainty" not in p, "molecfit Fe I still owes a xi measurement"


def test_the_contract_accepts_the_migrated_feed():
    """RYA-587 is ASKED, not bypassed -- that acceptance is the result."""
    from pipeline.uncertainty_contract import assert_publication_feed
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    assert_publication_feed(feed, previous=feed)


def test_a_migrated_rows_total_is_the_contracts_not_the_legacy_quadrature():
    """🔴 Regression: update_xi_budget silently recomputed the two-term legacy quadrature
    over the contract total, discarding fourteen components and restoring the understated
    bar. It must defer once a row carries an uncertainty block."""
    import math, importlib, sys as _s
    _s.path.insert(0, str(ROOT / "scripts"))
    emit = importlib.import_module("rya1178_emit_fe_schema")
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(x for x in feed["products"] if "uncertainty" in x)
    before = p["sigma_reported"]
    legacy = math.sqrt(sum(t * t for t in (p.get("sigma_stat"), p.get("sigma_syst_complete")) if t))
    assert abs(before - legacy) > 0.05, "the two totals must actually differ for this to bite"
    #: RYA-1224's index contract refuses a hand-rolled partial dict, so build the real one.
    hold, inst, models, xi_doc = emit.load_sources()
    emit.update_xi_budget(p, emit.xi_index(xi_doc, feed))
    assert p["sigma_reported"] == before, "the contract total was clobbered"
