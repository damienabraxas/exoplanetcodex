"""Controls for the census's scientifically consequential refusals."""
import importlib.util
import csv
import gzip
import hashlib
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "census", Path(__file__).resolve().parents[1] / "scripts/rya1353_line_census.py")
census = importlib.util.module_from_spec(spec)
spec.loader.exec_module(census)


def candidate(wave, ep, rows):
    return census.tuple_candidates(wave, ep, rows, "wave", "ep")


def test_wrong_lower_level_never_matches_lone_wavelength_neighbor():
    assert candidate(6065.490, 4.956, [{"wave": "6065.490", "ep": "2.608"}])[1] is None


def test_missing_ep_is_identity_hold_not_absence():
    assert candidate(7771.94, None, [{"wave": "7771.94", "ep": "9.146"}])[0] == "MISSING_IDENTITY_HOLD"


def test_unresolved_components_never_nearest_selected():
    status, row = candidate(7771.94, 9.146, [
        {"wave": "7771.939", "ep": "9.146"}, {"wave": "7771.941", "ep": "9.146"}])
    assert status == "AMBIGUOUS_HOLD" and row is None


def test_identical_wavelengths_cannot_collapse_different_gf_or_components():
    rows = [{"wave": "7771.94", "ep": "9.146", "gf": x} for x in ("-1", "-2")]
    assert candidate(7771.94, 9.146, rows) == ("DUPLICATE_WAVELENGTH_HOLD", None)


def test_tuple_is_only_candidate_not_identity_proof():
    rows = [{"wave": "7771.94", "ep": "9.146"}]
    status, row = candidate(7771.941, 9.146, rows)
    assert status == "TUPLE_CANDIDATE_COMPONENT_REVIEW" and row == rows[0]


def test_invalid_number_fails_and_blank_does_not_mean_zero():
    assert census.number("") is None
    assert census.number("0") == 0
    with pytest.raises(ValueError):
        census.number("inf")


def test_ordinal_review_keeps_experiment_as_benchmark_and_depth_as_prediction():
    row = dict(species="Fe I", gf_availability="FINITE_CANONICAL_GF",
               experimental_gf_state="RECORDED_LAB", predicted_depth_state="PREDICTED_STRONG")
    assert census.queue(row) == "FE_STRESS_TEST:WITHHELD_LAB_BENCHMARK"


def test_committed_census_conserves_every_canonical_row_without_science_mutation():
    out = census.ROOT / "data/audit/rya1353_line_census/v4"
    with gzip.open(out / "canonical_census.csv.gz", "rt") as stream:
        actual = {r['line_id']: r for r in csv.DictReader(stream)}
    expected = [r for r in census.read(census.ROOT / "data/linelists/canonical_gf.csv")
                if r['species'] in census.SPECIES]
    assert set(actual) == {r['line_id'] for r in expected}
    for row in expected:
        assert all(actual[row['line_id']][k] == v for k, v in row.items())
        assert actual[row['line_id']]['scientific_use'].startswith('HOLD')


def test_committed_input_manifest_authenticates_all_inputs():
    manifest = census.read(census.ROOT / "data/audit/rya1353_line_census/v4/input_manifest.csv")
    for row in manifest:
        data = (census.ROOT / row['path']).read_bytes()
        assert len(data) == int(row['bytes'])
        assert hashlib.sha256(data).hexdigest() == row['sha256'], row['path']


def test_iag_products_keep_distinct_runtime_spans_and_do_not_inherit_manifest_conditioning():
    rows = census.read(census.ROOT / 'data/audit/rya1353_line_census/v4/holding_spans.csv')
    exact = {r['holding_id']: r for r in rows if r['association'] == 'EXACT_RUNTIME_HOLDING_ID'}
    baker, reiners = exact['solar_iag'], exact['solar_iag_reiners2016']
    assert float(baker['wave_min_A']) > float(reiners['wave_min_A'])
    assert float(baker['wave_max_A']) > float(reiners['wave_max_A'])
    assert all(r['telluric_applied'] == 'UNKNOWN_EXACT_PRODUCT_MAPPING_PENDING'
               for r in rows if r['association'].startswith('MANIFEST_INSTRUMENT'))
