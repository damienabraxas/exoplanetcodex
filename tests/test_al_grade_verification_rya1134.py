"""Physical-identity and uncertainty counterexamples for the Al handoff."""
import math

import numpy as np
import pandas as pd
import pytest

from pipeline.al_grade_verification import (ROOT, build, burheim_sources, component_audit,
    finalize, fractional_bound_dex, johnson_source, loggf_from_a, unique_match, vujnovic_sources)


@pytest.fixture(scope="module")
def products():
    ledger, sources, candidates = build()
    ledger, components = component_audit(ledger, ROOT)
    return finalize(ledger, sources, components), sources, candidates, components


def test_wrong_species_or_lower_energy_cannot_promote():
    source = pd.DataFrame([dict(species="Al I", wavelength_air=6696.015, lower_EP=3.1427212)])
    for species, ep in (("Al II", 3.1427212), ("Al I", 4.0215), ("Al I", np.nan)):
        row = dict(species=species, wavelength_air=6696.015, lower_EP=ep)
        assert unique_match(row, source, .008, .0001).empty


def test_coincident_fine_structure_is_not_nearest_selected():
    source = pd.DataFrame([dict(species="Al I", wavelength_air=w, lower_EP=4.0216) for w in (7836.133, 7836.135)])
    row = dict(species="Al I", wavelength_air=7836.134, lower_EP=4.0216)
    assert len(unique_match(row, source, .008, .0001)) == 2


def test_a_conversion_uses_statistical_weight_and_vacuum():
    # Independent SI-constant derivation catches Angstrom/metre and g-factor errors.
    from scipy.constants import m_e, epsilon_0, c, e
    gf = loggf_from_a(2669.95, 1., 3330.)
    expected = math.log10(m_e * epsilon_0 * c / (2*math.pi*e**2) * (2669.95e-10)**2 * 3 * 3330)
    assert gf == pytest.approx(expected, abs=.00001)
    assert loggf_from_a(2669.95, 0., 3330.) == pytest.approx(gf-math.log10(3))
    with pytest.raises(ValueError):
        loggf_from_a(2669.95, 1., 0.)


def test_burheim_theory_lifetimes_never_get_primary_lab():
    source = burheim_sources(ROOT)
    assert source.source_class.value_counts().to_dict() == {"PRIMARY_LABORATORY": 8, "MIXED_LAB_THEORY": 4}
    assert (source.recomputed_loggf-source.loggf).abs().max() < .006


def test_vujnovic_flags_and_missing_error_survive():
    source = vujnovic_sources(ROOT).set_index("source_record_id")
    assert len(source) == 106
    weak = source.loc["vuj2002_t2_010"]
    assert weak.source_class == "MIXED_LAB_THEORY"
    assert math.isnan(weak.gf_bound_dex)
    assert source.loc["vuj2002_t2_011", "gf_bound_dex"] == pytest.approx(math.log10(1.02))
    assert source[source.uncertainty_limit.ne("")].gf_bound_dex.isna().all()
    assert source[source.aki_limit.ne("")].loggf.isna().all()
    assert source[source.theory_flag.eq("*")].source_class.eq("MIXED_LAB_THEORY").all()
    assert source.other_reference.str.len().gt(0).any()


def test_johnson_confidence_is_not_silently_one_sigma():
    row = johnson_source().iloc[0]
    assert "90%" in row.uncertainty_basis
    assert row.gf_bound_dex == pytest.approx(math.log10(1+230/3330))
    assert math.isnan(fractional_bound_dex(float("nan")))


def test_full_denominator_and_no_placeholder_uncertainties(products):
    ledger, _, _, _ = products
    assert len(ledger) == 505 and ledger.canonical_line_id.is_unique
    held = ledger[ledger.verified_source_record_id.eq("")]
    assert held.gf_uncertainty_bound_dex.isna().all()
    assert held.verified_loggf.isna().all()
    assert not ledger.verified_source_class.isna().any()


def test_6696_neighbor_cannot_inherit_lab_gf(products):
    ledger = products[0].set_index("canonical_line_id")
    assert ledger.loc["alphys_I_6696.0150_0352", "verified_source_class"] == "PRIMARY_LABORATORY"
    assert ledger.loc["alphys_I_6696.1850_0353", "verified_source_class"] == "UNRESOLVED"


def test_11254_component_never_becomes_blend_total(products):
    ledger = products[0].set_index("canonical_line_id")
    row = ledger.loc["alphys_I_11254.9239_0407"]
    assert row.reference_membership == "HOLD"
    assert row.disposition_reason == "COMPONENT_NOT_BLEND_TOTAL"
    assert row.hfs_component_status == "MULTIPLE_FINE_STRUCTURE_TRANSITIONS"


def test_component_proof_recovers_physical_identity_and_closes(products):
    ledger, _, _, components = products
    for line_id in ("alphys_I_6696.0150_0352", "alphys_I_6698.6730_0355",
                    "alphys_I_3944.0060_0335", "alphys_I_3961.5200_0336"):
        part = components[components.canonical_line_id.eq(line_id)]
        assert len(part) > 1 and part.raw_path.str.len().gt(0).all()
        assert part.gf_fraction.sum() == pytest.approx(1.)
        assert math.log10(np.power(10., part.component_loggf).sum()) == pytest.approx(part.parent_sum_loggf.iloc[0])
        assert math.log10(np.power(10., part.rescaled_component_loggf).sum()) == pytest.approx(part.source_total_loggf.iloc[0])
        assert ledger.set_index("canonical_line_id").loc[line_id, "physical_identity_status"] == "VERIFIED_RAW_J_AND_SOURCE_IDENTITY"


def test_reference_overlaps_depth_pools_and_replication_is_distinct(products):
    ledger = products[0]
    assert (ledger.reference_membership.eq("MEMBER") & ledger.codex_membership.eq("MEMBER")).any()
    assert (ledger.reference_membership.eq("MEMBER") & ledger.deep_membership.eq("MEMBER")).any()
    assert ledger.replication_membership.eq("MEMBER").sum() == 6
    assert ledger.replication_membership.eq("EXCLUDED_BY_SOURCE").sum() == 1
    assert ledger.reference_membership.eq("MEMBER").sum() > 6


def test_raw_J_conflict_blocks_handoff(products):
    ledger, sources, _, components = products
    bad = components.copy()
    bad.loc[bad.canonical_line_id.eq("alphys_I_6696.0150_0352"), "upper_J"] = 7.5
    checked = finalize(ledger, sources, bad).set_index("canonical_line_id")
    assert checked.loc["alphys_I_6696.0150_0352", "reference_membership"] == "HOLD"
    assert checked.loc["alphys_I_6696.0150_0352", "atomic_handoff_status"] == "HOLD"


def test_wrong_upper_energy_blocks_even_when_J_matches(products):
    ledger, sources, _, components = products
    bad = components.copy()
    line_id = "alphys_I_6696.0150_0352"
    bad.loc[bad.canonical_line_id.eq(line_id), "raw_upper_EP"] += .1
    checked = finalize(ledger, sources, bad).set_index("canonical_line_id")
    assert checked.loc[line_id, "physical_identity_status"] == "SOURCE_RAW_UPPER_ENERGY_CONFLICT"
    assert checked.loc[line_id, "atomic_handoff_status"] == "HOLD"


def test_holding_and_engine_denominators_remain_explicit(products):
    from pipeline.al_eligibility import matrices
    # No staged readers: absence must become HOLD, never an inferred clean spectrum.
    lines, engines = matrices(products[0], ROOT, holding_specs={})
    registry = pd.read_csv(ROOT / "data/catalog/holdings_manifest_registry.csv")
    solar = registry[registry.system_id.eq("solar")]
    assert len(lines) == 505 * len(solar)
    assert not lines.duplicated(["canonical_line_id", "holding_id"]).any()
    assert set(lines.holding_id) == set(solar.holding_id)
    assert lines.observed_conditioning.eq("unknown").all()
    assert not lines.status.eq("RUN").any()
    assert engines.n_verified_eligible.eq(0).all()
    assert engines.reason.str.len().gt(0).all()
    # Corrected/uncorrected siblings preserve their own declaration and same gf.
    left = lines[lines.holding_id.eq("solar_harps")].set_index("canonical_line_id")
    right = lines[lines.holding_id.eq("solar_harps_molecfit_corrected")].set_index("canonical_line_id")
    assert left.telluric_applied.eq("not-applied").all()
    assert right.telluric_applied.eq("applied").all()
    assert left.verified_source_record_id.equals(right.verified_source_record_id)
