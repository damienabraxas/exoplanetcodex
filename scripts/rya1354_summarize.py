"""Summarize retained exploratory solver outputs; never fabricates missing validation.

Requires the diagnostic receipt and native text. Raw wavefunctions and
full inputs remain in the receipt's Sirius run directory. No experimental gf is read.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re


def h(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(token: str) -> float:
    x = float(token.replace("D", "E"))
    if not math.isfinite(x) or x <= 0:
        raise ValueError(f"invalid native transition value: {token}")
    return x


def parse_grasp(text: str) -> list[dict]:
    pattern = (r"f2\s+1\s+([123])\s+\+\s+f1\s+1\s+2\s+-\s+([\d.]+)\s+C\s+([\d.D+\-]+)\s+([\d.D+\-]+)\s+([\d.D+\-]+)"
               r"\s+B\s+([\d.D+\-]+)\s+([\d.D+\-]+)\s+([\d.D+\-]+)")
    matches = re.findall(pattern, text)
    if len(matches) != 3 or {m[0] for m in matches} != {"1", "2", "3"}:
        raise ValueError("expected three explicitly identified J=1,2,3 -> J=2 E1 rows; no first-hit parse")
    rows = []
    for j, wave, a_c, gf_c, s_c, a_b, gf_b, s_b in matches:
        energy, ac, gc, sc, ab, gb, sb = map(number, (wave, a_c, gf_c, s_c, a_b, gf_b, s_b))
        rows.append({"candidate_id": f"OI_2p3-4S-3s_5So_J2__2p3-4S-3p_5P_J{j}",
                     "atomic_number": 8, "ion_stage": 1, "multipole": "E1",
                     "two_J_lower": 4, "two_J_upper": 2*int(j), "parity_lower": -1, "parity_upper": 1,
                     "transition_energy_cm1_native": energy, "wavelength_vac_A_native": 1e8/energy,
                     "A_B_s1": ab, "A_C_s1": ac, "gf_B": gb, "gf_C": gc,
                     "S_B_au": sb, "S_C_au": sc,
                     "gauge_fraction_A": abs(ab-ac)/max(ab, ac),
                     "origin": "theoretical", "method": "GRASP_Dirac_Fock_single_configuration",
                     "isotope": "UNSPECIFIED_POINT_STATIC_NUCLEUS",
                     "level_identity": "CONFIGURED_STATES_WITH_NATIVE_LS_COMPOSITION_NOT_CANONICAL_ADOPTION",
                     "observed_level_substitution": False, "uncertainty": None, "covariance": None,
                     "lifetime_s": None, "branching_fraction": None,
                     "status": "CANDIDATE_HOLD_NO_ACTIVE_SPACE_CONVERGENCE_OR_BLIND_VALIDATION"})
    return rows


def summarize(root: Path) -> dict:
    receipt_path = root / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    if not (receipt["status"].startswith("SOLVER_STEPS_COMPLETED") or
            receipt["status"] == "COMPLETED_WITH_SCIENTIFIC_FAILURES"):
        raise ValueError("solver run did not complete")
    for relative, record in receipt["files"].items():
        path = root / relative
        if path.exists() and h(path) != record["sha256"]:
            raise ValueError(f"changed native artifact: {path}")
    # These specific records are indispensable; omitted wavefunction files may
    # remain on Sirius, but required native outputs cannot silently disappear.
    for relative in ("grasp_E1/lower.upper.t", "grasp_3s_J2/state.lsj.lbl",
                     "grasp_3p_J123/state.lsj.lbl", "atsp_E1/Mltpol.stdout", "atsp_E1/Lstr.stdout"):
        if relative not in receipt["files"] or not (root / relative).is_file():
            raise ValueError(f"missing required native artifact {relative}")
    lower = (root / "grasp_3s_J2/state.lsj.lbl").read_text()
    upper = (root / "grasp_3p_J123/state.lsj.lbl").read_text()
    if "2p(3)4S.3s_5S" not in lower or upper.count("2p(3)4S.3p_5P") != 3:
        raise ValueError("LS composition does not establish the configured target states")
    rows = parse_grasp((root / "grasp_E1/lower.upper.t").read_text())
    for row in rows:
        row["calculation_receipt_sha256"] = h(receipt_path)
        row["result_source"] = {"kind": "THIS_UNPUBLISHED_CALCULATION", "doi": None,
                                "doi_status": "NO_DOI_ASSIGNED_TO_CANDIDATE", "method_reference_doi": "10.1016/j.cpc.2018.10.032"}
        row["native_table"] = {"path": "grasp_E1/lower.upper.t", "sha256": h(root / "grasp_E1/lower.upper.t")}
        row["version"] = "rya1354.exploratory.v3"
    angular = (root / "atsp_E1/Mltpol.stdout").read_text()
    refused = "TOO MANY ELECTRONS FOR THIS CODE" in angular
    modern = root / "atsp2k_biorthogonal_E1/biotr.stdout"
    if modern.exists() and str(modern.relative_to(root)) not in receipt["files"]:
        raise ValueError("unrecorded independent-method output")
    modern_refused = modern.is_file() and "CONFIGURATION DATA WRONG" in modern.read_text()
    modern_result = {"status": "FAILED_COUPLING_INPUT_REFUSAL" if modern_refused else "NOT_VALIDATED",
                     "gf_comparison_valid": False,
                     "reason": "Book and native modern coupling/eigenvector formats require explicit value-preserving export."}
    ls_path = root / "atsp2k_biorthogonal_E1/lower.upper.ls"
    if modern.is_file() and not modern_refused and ls_path.is_file():
        if str(ls_path.relative_to(root)) not in receipt["files"]:
            raise ValueError("unrecorded independent transition table")
        native = ls_path.read_text()
        values = re.findall(r"S\s*=\s*([\d.D+\-]+)\s+GF\s*=\s*([\d.D+\-]+)\s+AKI\s*=\s*([\d.D+\-]+)", native)
        if len(values) != 2 or "E1  length" not in native or "2p(3)4S3.3s_5S" not in native or "2p(3)4S3.3p_5P" not in native:
            raise ValueError("independent native result does not identify the target E1 multiplet and both gauges")
        sl, gl, al = map(number, values[0]); sv, gv, av = map(number, values[1])
        grasp_sum = math.fsum(r["gf_B"] for r in rows)
        modern_result = {"status": "EXPLORATORY_NATIVE_TRANSITION_COMPLETED", "gf_comparison_valid": True,
                         "scope": "MCHF_LS_multiplet_vs_sum_of_GRASP_fine_structure_gf_NOT_EXPERIMENTAL_VALIDATION",
                         "gf_length_multiplet": gl, "gf_velocity_multiplet": gv,
                         "A_length_term_average_s1": al, "A_velocity_term_average_s1": av,
                         "S_length_au": sl, "S_velocity_au": sv,
                         "grasp_gf_B_component_sum": grasp_sum,
                         "relative_gf_difference_to_mchf": (grasp_sum-gl)/gl,
                         "caveats": "Different small configuration models, relativistic Hamiltonians and native Rydberg conventions; no uncertainty, convergence or blind validation."}
    return {"schema": "rya1354.exploratory-candidates.v1", "receipt_sha256": h(receipt_path),
            "scientific_admission": "HOLD", "blindness": "NOT_BLIND_EXPLORATORY",
            "candidates": rows,
            "independent_method": {"code": "ATSP_book_MCHF", "status": "FAILED_ANGULAR_DOMAIN_REFUSAL" if refused else "REVIEW_REQUIRED",
                                   "gf_comparison_valid": False,
                                   "reason": "Relaxed open-shell nonorthogonal transition is refused by book angular code; exit 0 and printed near-zero gf are not a measurement." if refused else "Independent method has not been scientifically validated"},
            "modern_biorthogonal_attempt": modern_result,
            "energy_scope": "Native computed levels only; point static nucleus; no observed substitution or fitted shift",
            "missing": ["active-space_and_correlation_convergence", "Breit_QED_sensitivity", "radial_grid_sensitivity",
                        "independent_converged_transition_model", "withheld_primary_lab_benchmark", "uncertainty_calibration",
                        "multiplet_covariance", "complete_decay_inventory_for_lifetime_and_BF", "UV_IR_and_forbidden_line_calculations"],
            "publication_allowed": False}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    result = summarize(args.run)
    with args.out.open("x") as fh:
        json.dump(result, fh, indent=2, allow_nan=False)
        fh.write("\n")
    print(f"{len(result['candidates'])} exploratory candidate rows; admission HOLD; modern limited comparison {result['modern_biorthogonal_attempt']['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
