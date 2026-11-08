"""Fictional fixtures test arithmetic, identity, leakage and failure controls."""
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.rya1354_benchmark import (angularly_allowed, correlated_score, decay_summary,
                                      e1_gf, file_hash, score_frozen)
from scripts.rya1354_summarize import parse_grasp


def test_multipoles_have_distinct_parity_and_triangles():
    assert angularly_allowed(4, 2, -1, 1, "E1")
    assert not angularly_allowed(0, 0, -1, 1, "E1")
    assert not angularly_allowed(4, 2, 1, 1, "E1")
    assert angularly_allowed(4, 2, 1, 1, "M1")
    assert angularly_allowed(0, 4, 1, 1, "E2")
    assert not angularly_allowed(0, 2, 1, 1, "E2")
    assert not angularly_allowed(1, 2, -1, 1, "E1")


@pytest.mark.parametrize("multipole", ["M1", "E2"])
def test_e1_rate_helper_refuses_forbidden_multipoles(multipole):
    with pytest.raises(ValueError, match="E1 conversion refuses"):
        e1_gf(1., 1., 2, multipole=multipole)


@pytest.mark.parametrize("bad", [0., -1., float("inf"), float("nan")])
def test_nonphysical_numbers_never_become_valid_rates(bad):
    with pytest.raises(ValueError):
        e1_gf(bad, 1000., 2, multipole="E1")


def test_e1_wavelength_and_degeneracy_scaling():
    base = e1_gf(1e7, 1000., 0, multipole="E1")
    assert e1_gf(1e7, 2000., 2, multipole="E1") == pytest.approx(12*base)


def test_lifetime_requires_all_decay_channels():
    with pytest.raises(ValueError, match="complete decay"):
        decay_summary([2., 3.], complete_inventory=False)
    result = decay_summary([2., 3.], complete_inventory=True)
    assert result == {"lifetime_s": .2, "branching_fractions": [.4, .6]}


def test_shared_multiplet_covariance_does_not_average_like_independent_lines():
    result = correlated_score([10**.1, 10**.1], [1., 1.], [[.01, .009], [.009, .01]])
    assert result["chi_square"] == pytest.approx(2/1.9)
    independent = correlated_score([10**.1, 10**.1], [1., 1.], [[.01, 0.], [0., .01]])
    assert independent["chi_square"] == pytest.approx(2.)


@pytest.mark.parametrize("cov", [[[1., 2.], [2., 1.]], [[1., 1.], [1., 1.]],
                                  [[1., 0.], [.1, 1.]], [[1., 0.], [0., float("nan")]]])
def test_invalid_covariance_refuses_without_diagonal_fallback(cov):
    with pytest.raises(ValueError):
        correlated_score([1., 1.], [1., 1.], cov)


def fixture_bundle(tmp_path):
    p, b = tmp_path / "pred.json", tmp_path / "bench.json"
    model_hash = "a" * 64
    common = {"quantity": "A_ki", "unit": "s^-1", "rows": [{"transition_id": "fictional:t1", "value": 1.}]}
    p.write_text(json.dumps({**common, "origin": "theoretical", "model_manifest_sha256": model_hash,
                             "solar_calibrated": False}))
    b.write_text(json.dumps({**common, "origin": "experimental", "primary_source": "fictional:test-only",
                             "uncertainty_convention": "one_sigma_normal_log10", "reference_independent_of_model": True}))
    lock = {"prediction_sha256": file_hash(p), "benchmark_sha256": file_hash(b),
            "model_manifest_sha256": model_hash, "curator_receipt": "fictional:test-only",
            "exposure_status": "CURATOR_ATTESTED_UNEXPOSED", "ordered_transition_ids": ["fictional:t1"],
            "quantity": "A_ki", "unit": "s^-1", "combined_covariance_log10": [[.01]]}
    return lock, p, b


def test_frozen_score_still_requires_review(tmp_path):
    lock, p, b = fixture_bundle(tmp_path)
    score = score_frozen(lock, p, b)
    assert score["chi_square"] == 0
    assert score["scientific_admission"] == "HOLD_REQUIRES_REVIEW"
    assert "NOT_SOFTWARE_CERTIFIED" in score["blindness"]


def test_changed_predictions_and_exposed_benchmarks_are_refused(tmp_path):
    lock, p, b = fixture_bundle(tmp_path)
    exposed = {**lock, "exposure_status": "ALREADY_EXPOSED"}
    with pytest.raises(ValueError, match="previously exposed"):
        score_frozen(exposed, p, b)
    p.write_text(p.read_text() + " ")
    with pytest.raises(ValueError, match="changed after freeze"):
        score_frozen(lock, p, b)


@pytest.mark.parametrize("change,expected", [
    ({"unit": "ns"}, "unit mismatch"),
    ({"origin": "evaluated"}, "primary experimental"),
    ({"uncertainty_convention": "NIST-B"}, "uncertainty convention"),
    ({"reference_independent_of_model": False}, "training overlap"),
    ({"rows": [{"transition_id": "wrong-state", "value": 1.}]}, "identity/order"),
])
def test_reveal_never_matches_by_wavelength_or_inherits_a_false_pedigree(tmp_path, change, expected):
    lock, p, b = fixture_bundle(tmp_path)
    doc = json.loads(b.read_text()); doc.update(change); b.write_text(json.dumps(doc))
    lock["benchmark_sha256"] = file_hash(b)
    with pytest.raises(ValueError, match=expected):
        score_frozen(lock, p, b)


def test_parser_refuses_partial_transition_tables():
    with pytest.raises(ValueError, match="three explicitly"):
        parse_grasp("Execution complete. No usable table.")


def test_committed_capsule_is_complete_and_candidates_reproduce(tmp_path):
    import tarfile
    from scripts.rya1354_summarize import summarize
    archive = ROOT / "data/audit/rya1354/solver_smoke_v10.tar.gz"
    provenance = json.loads((ROOT / "docs/audit/rya1354/toolchain_provenance.json").read_text())
    assert file_hash(archive) == provenance["capsule"]["sha256"]
    with tarfile.open(archive) as tar:
        tar.extractall(tmp_path, filter="data")
    run = tmp_path / "solver_smoke_v10"
    receipt = json.loads((run / "receipt.json").read_text())
    assert receipt["status"] == "COMPLETED_WITH_SCIENTIFIC_FAILURES"
    for relative, expected in receipt["files"].items():
        assert file_hash(run / relative) == expected["sha256"], relative
    generated = summarize(run)
    stored = json.loads((ROOT / "data/audit/rya1354/exploratory_candidates_v3.json").read_text())
    assert generated == stored
    assert len(stored["candidates"]) == 3
    assert not stored["independent_method"]["gf_comparison_valid"]
    assert stored["modern_biorthogonal_attempt"]["gf_comparison_valid"]
    assert all(r["uncertainty"] is None and r["lifetime_s"] is None for r in stored["candidates"])
    for row in stored["candidates"]:
        for gauge in ("B", "C"):
            converted = e1_gf(row[f"A_{gauge}_s1"], row["wavelength_vac_A_native"], row["two_J_upper"], multipole="E1")
            assert converted == pytest.approx(row[f"gf_{gauge}"], rel=2e-5)  # native table/constant rounding
    from scripts.rya1354_solver_smoke import modern_coupling_layout
    for stem, directory in (("lower", "atsp_3s_5S"), ("upper", "atsp_3p_5P")):
        book = (run / directory / "cfg.out").read_text()
        exported = modern_coupling_layout(book)
        assert exported == (run / "atsp2k_biorthogonal_E1" / f"{stem}.c").read_text()
        assert exported[:30] == book.splitlines()[0].ljust(30)
    control = json.loads((ROOT / "docs/audit/rya1354/atsp2k_native_control.json").read_text())
    assert not control["semantic_failure"] and control["native_equals_expected_bytes"]
    assert control["native_ls_sha256"] == control["expected_ls_sha256"]


def test_legacy_zero_exit_does_not_hide_domain_refusal(tmp_path):
    from scripts.rya1354_solver_smoke import run
    fake = tmp_path / "fictional_solver"
    fake.write_text('#!/bin/sh\nprintf "CONFIGURATION DATA WRONG, PROGRAM HALTED\\n"\n')
    fake.chmod(0o755)
    result = run(fake, "", tmp_path)
    assert result["returncode"] == 0
    assert result["semantic_errors"] == ["CONFIGURATION DATA WRONG", "PROGRAM HALTED"]
