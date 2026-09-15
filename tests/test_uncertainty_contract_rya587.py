"""Synthetic science controls for the universal publication contract."""
import copy
import json
import math

import numpy as np
import pytest

from pipeline import uncertainty_contract as uc
from pipeline import product_eligibility as pe


def complete_product(element="Fe", star="solar"):
    p = dict(element=element, ion="I", band="VIS", instrument="harps",
             holding="solar_harps_molecfit_corrected", tier="GRADED", selector="GRADED",
             route="SYNTH", treatment="1D-LTE", star=star,
             uncertainty_indicator_ids=["physical-transition-1", "physical-transition-2"],
             A=7.5, n_lines=2, sigma_stat=0.1, sigma_syst=0.05,
             provenance={"artifact_mtime": "2026-09-15T00:00:00Z"})
    scope = uc.product_scope(p, star=star, indicator_ids=p["uncertainty_indicator_ids"])
    stat = uc.scatter_measurement([7.4, 7.6], indicator_ids=p["uncertainty_indicator_ids"],
                                  source="synthetic independent lines", independent=True)
    gf = uc.transition_data(p["uncertainty_indicator_ids"], [0.05, 0.05], [0.5, 0.5],
                            covariance=[[0.0025, 0.0025], [0.0025, 0.0025]],
                            sources=["synthetic common lab scale"] * 2,
                            covariance_source="synthetic perfectly correlated gf scale")
    terms = [stat, gf]
    for name in uc.COMPONENTS[2:]:
        # Fixture represents exact stellar parameters and explicitly absent other
        # effects. These are synthetic truths, never defaults for actual products.
        defined = name.startswith("stellar.")
        terms.append(uc.component(name, 0.0 if defined else None,
                                  state="DEFINED" if defined else "N/A",
                                  source="synthetic control: exact parameters, isolated lines"))
    cov = np.diag([c["sigma_dex"] ** 2 for c in terms if c["sigma_dex"] is not None]).tolist()
    p["uncertainty"] = uc.assemble(scope, terms, covariance=cov,
                                   covariance_source="synthetic component covariance",
                                   assumptions="measurement independent of gf; parameters exact")
    p["sigma_reported"] = p["uncertainty"]["sigma_reported"]
    return p


@pytest.mark.parametrize("element", ["Fe", "Al", "Si", "C", "N", "O", "Ba"])
def test_one_contract_for_current_and_future_elements(element):
    p = complete_product(element)
    assert uc.publication_problems(p) == []
    assert p["sigma_reported"] == pytest.approx(math.hypot(0.1, 0.05))
    p.pop("uncertainty")
    assert "UNCERTAINTY_INCOMPLETE" in {r.code for r in pe.evaluate(p)}


def test_common_gf_scale_does_not_average_down():
    p = complete_product()
    gf = p["uncertainty"]["components"][1]
    assert gf["sigma_dex"] == pytest.approx(0.05)
    assert gf["sigma_dex"] != pytest.approx(0.05 / math.sqrt(2))


@pytest.mark.parametrize("mutation", ["missing", "nan", "negative", "no_source", "false_complete", "wrong_pool", "wrong_total"])
def test_bad_evidence_cannot_publish(mutation):
    p = complete_product()
    c = p["uncertainty"]["components"]
    if mutation == "missing": c.pop()
    if mutation == "nan": c[0]["sigma_dex"] = float("nan")
    if mutation == "negative": c[0]["sigma_dex"] = -1
    if mutation == "no_source": c[0]["source"] = ""
    if mutation == "false_complete": c[-1]["state"] = "HOLD"
    if mutation == "wrong_pool": p["uncertainty_indicator_ids"] = ["other"]
    if mutation == "wrong_total": p["sigma_reported"] = 0.001
    assert uc.publication_problems(p)


@pytest.mark.parametrize("cov", [[[1, 2], [2, 1]], [[1, 0.5], [0.4, 1]], [[1, float('nan')], [0, 1]], [[1]]])
def test_nonphysical_covariance_refused(cov):
    with pytest.raises(uc.UncertaintyError):
        uc.covariance_variance([1, -1], cov)


def test_shared_scale_cancels_in_matched_differential():
    target, ref = complete_product(star="procyon"), complete_product()
    a, b = target["uncertainty"], ref["uncertainty"]
    n = len(a["covariance"])
    cross = np.zeros((n, n)); cross[1, 1] = 0.0025
    result = uc.differential_uncertainty(a, b, matched_ids=target["uncertainty_indicator_ids"],
        target_scope=a["scope"], reference_scope=b["scope"], cross_covariance=cross.tolist(),
        source="synthetic shared lab scale", assumptions="independent measurements; fully shared gf")
    assert result["sigma_differential"] == pytest.approx(math.sqrt(0.02))
    assert result["sigma_differential"] < math.hypot(target["sigma_reported"], ref["sigma_reported"])
    assert result["target"]["sigma_reported"] == target["sigma_reported"]
    target.update(X_H=0.1, sigma_differential=result["sigma_differential"],
                  differential_uncertainty=result)
    assert uc.publication_problems(target) == []
    target["differential_uncertainty"]["cross_covariance"][1][1] = 0
    assert uc.publication_problems(target)


def test_differential_cannot_borrow_another_route():
    a, b = complete_product(star="procyon"), complete_product()
    b["route"] = "PROFILEFIT"
    b["uncertainty"]["scope"] = uc.product_scope(b, star="solar", indicator_ids=b["uncertainty_indicator_ids"])
    with pytest.raises(uc.UncertaintyError, match="unmatched route"):
        uc.differential_uncertainty(a["uncertainty"], b["uncertainty"], matched_ids=a["uncertainty_indicator_ids"],
            target_scope=a["uncertainty"]["scope"], reference_scope=b["uncertainty"]["scope"],
            cross_covariance=np.zeros((6, 6)).tolist(), source="test", assumptions="test")


def test_single_line_and_correlated_scatter_remain_hold():
    for values, ids, independent in [([7.5], ["a"], True), ([7.4, 7.6], ["a", "b"], False)]:
        c = uc.scatter_measurement(values, indicator_ids=ids, source="synthetic", independent=independent)
        assert c["state"] == "HOLD" and c["sigma_dex"] is None


def test_profile_sigma_is_not_divided_by_pixel_count():
    c = uc.fit_measurement(0.03, source="synthetic fit", likelihood="profile chi2 at delta=1",
                           correlation_treatment="full resampling covariance", n_pixels=10000)
    assert c["sigma_dex"] == 0.03
    assert "n_lines" not in c["evidence"]


def test_perturbation_preserves_signed_response_and_asymmetry():
    r = uc.paired_response({"a": 7.5, "b": 7.6}, {"a": 7.4, "b": 7.5},
                           {"a": 7.7, "b": 7.8}, delta=100, source="synthetic runs",
                           parameter_source="synthetic +/-100 K record")
    assert r["signed_response_dex"] == pytest.approx(0.15)
    assert r["delta_minus_dex"] == pytest.approx(-0.1)
    assert r["delta_plus_dex"] == pytest.approx(0.2)
    with pytest.raises(uc.UncertaintyError, match="pools differ"):
        uc.paired_response({"a": 7.5, "b": 7.6}, {"a": 7.4}, {"a": 7.7, "b": 7.8},
                           delta=100, source="test", parameter_source="test")


def test_ir_telluric_and_molecular_coupling_are_required():
    for field, value in [("band", "H"), ("selector", "MOL-CN_red")]:
        p = complete_product("N"); p[field] = value
        p["uncertainty"]["scope"] = uc.product_scope(p, star=p["star"], indicator_ids=p["uncertainty_indicator_ids"])
        assert uc.publication_problems(p)


def test_ratio_propagates_covariance_and_log_to_linear_units():
    result = uc.ratio_uncertainty([8.4, 8.7], [[0.01, 0.008], [0.008, 0.01]],
        numerator=0, denominator=1, source="coupled CNO fit", assumptions="joint fitted C,O",
        linear_ratio=True)
    assert result["value"] == pytest.approx(10**-0.3)
    assert result["sigma"] == pytest.approx(math.log(10) * 10**-0.3 * math.sqrt(0.004))


def test_canonical_budget_round_trips_through_publisher():
    import pandas as pd
    from scripts.publish_product import normalise
    p = complete_product()
    frame = pd.DataFrame([{**p, "stat_dex": p["sigma_stat"], "syst_dex": p["sigma_syst"],
                           "n_excluded": 0, "uncertainty": json.dumps(p["uncertainty"]),
                           "uncertainty_indicator_ids": json.dumps(p["uncertainty_indicator_ids"])}])
    rows = normalise(frame, holding=p["holding"], tier=p["tier"], route=p["route"], selector=p["selector"])
    assert rows[0]["uncertainty"] == p["uncertainty"]
    assert rows[0]["uncertainty_indicator_ids"] == p["uncertainty_indicator_ids"]


def test_fe_legacy_arithmetic_is_unchanged():
    from pipeline.error_budget import ErrorBudget, Term
    b = ErrorBudget("Fe", "VIS", 4).add(Term("scatter", 0.2, True, "synthetic"))
    b.add(Term("gf", 0.05, False, "synthetic"))
    assert b.total() == pytest.approx((0.1, 0.05))
    assert b.component_records()[0]["sigma_dex"] == pytest.approx(0.1)


def test_write_boundary_cannot_republish_incomplete_legacy_record(tmp_path):
    from scripts.publish_product import write_feed
    p = complete_product()
    p.pop("uncertainty")
    path = tmp_path / "feed.json"
    path.write_text("original")
    with pytest.raises(uc.UncertaintyError):
        write_feed(path, {"products": [p]})
    assert path.read_text() == "original"
