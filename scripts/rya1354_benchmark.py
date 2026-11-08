"""Atomic benchmark controls and diagnostics; no selection or publication imports.

The score CLI requires exact frozen prediction/benchmark byte hashes. It reports
curator attestations, never certifies that software can prove human blindness.
Synthetic tests use fictional transition IDs, not invented oxygen measurements.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

# NIST ASD lines help: lambda in vacuum Angstrom, A in s^-1.
# Rounded published coefficient; this helper is deliberately scoped to E1.
GF_A_COEFFICIENT = 1.49919e-16


def positive(value: float, subject: str) -> float:
    if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{subject} must be finite and positive")
    return value


def e1_gf(a_s1: float, wavelength_vac_A: float, two_j_upper: int, *, multipole: str) -> float:
    if multipole != "E1":
        raise ValueError("E1 conversion refuses M1/E2; retain native rates and multipole-specific physics")
    if type(two_j_upper) is not int or two_j_upper < 0:
        raise ValueError("two_j_upper must be a nonnegative integer")
    return GF_A_COEFFICIENT * (two_j_upper + 1) * positive(a_s1, "A") * positive(wavelength_vac_A, "wavelength")**2


def angularly_allowed(two_j_lower: int, two_j_upper: int, parity_lower: int,
                      parity_upper: int, multipole: str) -> bool:
    """Triangle/parity necessary conditions, not a guarantee of nonzero amplitude.

    No hard spin-forbidden rejection: relativistic mixing can produce amplitude.
    Parity is +1/-1; J is represented exactly as integer 2J.
    """
    if multipole not in ("E1", "M1", "E2"):
        raise ValueError("unsupported multipole")
    if any(type(j) is not int or j < 0 for j in (two_j_lower, two_j_upper)):
        raise ValueError("J must be supplied as nonnegative integer 2J")
    if any(type(p) is not int or p not in (-1, 1) for p in (parity_lower, parity_upper)):
        raise ValueError("parity must be +/-1")
    if (two_j_lower - two_j_upper) % 2:
        return False
    rank = 2 if multipole == "E2" else 1
    triangle = abs(two_j_lower - two_j_upper) <= 2*rank <= two_j_lower + two_j_upper
    parity_change = -1 if multipole == "E1" else 1
    return triangle and parity_lower * parity_upper == parity_change


def decay_summary(rates_s1: list[float], *, complete_inventory: bool) -> dict:
    if complete_inventory is not True:
        raise ValueError("A lifetime requires the complete decay inventory, including other multipoles")
    if not rates_s1:
        raise ValueError("empty decay inventory")
    rates = [positive(x, "decay rate") for x in rates_s1]
    total = math.fsum(rates)
    return {"lifetime_s": 1 / total, "branching_fractions": [x/total for x in rates]}


def correlated_score(predicted: list[float], measured: list[float], covariance_log10: list) -> dict:
    """Residual log10(predicted/measured), full-covariance chi-square; no fit.

    Covariance must already include the declared prediction/experiment cross
    covariance. Units are log10(quantity)^2 in the exact supplied row order.
    Unknown/singular/invalid covariance refuses calibration; no diagonal fallback.
    """
    if not predicted or len(predicted) != len(measured):
        raise ValueError("nonempty matched prediction/measurement arrays required")
    residual = np.log10([positive(x, "prediction") for x in predicted]) - np.log10([positive(x, "measurement") for x in measured])
    cov = np.asarray(covariance_log10, dtype=float)
    n = len(residual)
    if cov.shape != (n, n) or not np.isfinite(cov).all():
        raise ValueError("finite square covariance matching row order required")
    if not np.allclose(cov, cov.T, rtol=1e-12, atol=0):
        raise ValueError("covariance must be symmetric")
    try:
        factor = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as exc:
        raise ValueError("positive definite covariance required; singular shared-scale cases need an explicit reduced model") from exc
    white = np.linalg.solve(factor, residual)
    return {"residual_log10": residual.tolist(), "marginal_standardized_residuals": (residual/np.sqrt(np.diag(cov))).tolist(),
            "chi_square": float(white @ white), "n_comparisons": n,
            "fitted_parameters": 0, "quantity": "log10_ratio", "scientific_admission": "HOLD_REQUIRES_REVIEW"}


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def score_frozen(lock: dict, predictions: Path, benchmarks: Path) -> dict:
    for key, path in (("prediction_sha256", predictions), ("benchmark_sha256", benchmarks)):
        if file_hash(path) != lock.get(key):
            raise ValueError(f"{key} changed after freeze")
    if not lock.get("model_manifest_sha256") or not lock.get("curator_receipt"):
        raise ValueError("frozen model hash and independent-curator receipt required")
    if lock.get("exposure_status") != "CURATOR_ATTESTED_UNEXPOSED":
        raise ValueError("previously exposed data cannot be reported as a blind benchmark")
    pred = json.loads(predictions.read_text())
    bench = json.loads(benchmarks.read_text())
    if pred.get("origin") != "theoretical" or pred.get("model_manifest_sha256") != lock["model_manifest_sha256"]:
        raise ValueError("theoretical prediction must name the exact frozen model manifest")
    if pred.get("solar_calibrated") is not False:
        raise ValueError("explicit absence of Solar calibration required")
    ids = lock["ordered_transition_ids"]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("unique ordered transition IDs required")
    for name, doc in (("prediction", pred), ("benchmark", bench)):
        if [r["transition_id"] for r in doc["rows"]] != ids:
            raise ValueError(f"{name} identity/order mismatch; no wavelength-only alignment")
        if doc.get("quantity") != lock["quantity"] or doc.get("unit") != lock["unit"]:
            raise ValueError(f"{name} quantity/unit mismatch")
    if bench.get("origin") != "experimental" or not bench.get("primary_source"):
        raise ValueError("independent primary experimental benchmark required, not a NIST/theory label")
    if bench.get("uncertainty_convention") != "one_sigma_normal_log10":
        raise ValueError("explicit compatible uncertainty convention required; NIST letters are not Gaussian sigmas")
    if bench.get("reference_independent_of_model") is not True:
        raise ValueError("benchmark dependence or training overlap must be resolved")
    result = correlated_score([r["value"] for r in pred["rows"]], [r["value"] for r in bench["rows"]], lock["combined_covariance_log10"])
    result.update({"blindness": "CURATOR_ATTESTATION_RECORDED_NOT_SOFTWARE_CERTIFIED",
                   "prediction_sha256": file_hash(predictions), "benchmark_sha256": file_hash(benchmarks),
                   "curator_receipt": lock["curator_receipt"], "ordered_transition_ids": ids})
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--benchmarks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = score_frozen(json.loads(args.lock.read_text()), args.predictions, args.benchmarks)
    with args.out.open("x") as fh:
        json.dump(result, fh, indent=2, allow_nan=False)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
