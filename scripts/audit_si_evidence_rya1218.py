#!/usr/bin/env python3
"""Reconcile recovered Si measurements, physical gf identities and RYA-587 holds."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import si_evidence as si
from pipeline.gf_grades import grade_line


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def audit(root=ROOT):
    evidence = root / "data/audit/rya1218_continuation"
    canonical = pd.read_csv(evidence / "original_si_atomic_snapshot.csv")
    current = pd.read_csv(root / "data/linelists/canonical_gf.csv", low_memory=False).set_index("line_id")
    reference = pd.read_csv(root / "data/reference/si_agss21/si_agss21_lines.csv").set_index("canonical_line_id")
    manifest = json.loads((evidence / "recovery_manifest.json").read_text())
    for rec in manifest["artifacts"]:
        if sha(root / rec["path"]) != rec["sha256"]:
            raise ValueError(f"recovered artifact changed: {rec['path']}")
    locations = {
        "solar_harps_molecfit_corrected": root / "data/results/rya1218/recovered_20260916/harps",
        "solar_iag": root / "data/results/rya1218/recovered_20260916/iag",
        "solar_kpno_molecfit_corrected": root / "data/results/rya1218/reference_synthesis",
    }
    products, line_records, pools, inputs = [], [], {}, {}
    for holding, folder in locations.items():
        for path in sorted(folder.glob("*products.csv")):
            row = pd.read_csv(path).iloc[0].to_dict()
            if row["ion"] != "I":
                continue
            prefix = path.name.removesuffix("_products.csv")
            suffix = "_lines.csv" if row["treatment"].startswith("ENGINE-") else "_1D-LTE_lines.csv"
            line_path = path.with_name(prefix + suffix)
            lines = si.attach_identity(pd.read_csv(line_path), canonical)
            used = si.accepted(lines)
            stats = si.diagnostic_statistics(lines, row)
            product = dict(element="Si", ion="I", band=row["band"], instrument=row["instrument"],
                           holding=holding, tier="ALL", selector="FROMEW", route="SYNTH",
                           treatment=row["treatment"], line_set="si-agss21", star="solar",
                           A=float(row["A"]), n_lines=int(row["n_lines"]),
                           uncertainty_indicator_ids=sorted(used.line_id), sigma_reported=None)
            rel = str(path.relative_to(root))
            product["uncertainty"] = si.hold_budget(product, product["uncertainty_indicator_ids"],
                                                   source=rel, statistics=stats)
            product.update(artifact=rel, line_artifact=str(line_path.relative_to(root)),
                           statistics=stats, legacy_stat_dex=float(row["stat_dex"]),
                           legacy_syst_dex=float(row["syst_dex"]), disposition="DIAGNOSTIC_HOLD",
                           line_set_meaning="External-list membership on canonical gf, not source-scale replication")
            products.append(product)
            pools[(holding, row["treatment"])] = lines
            inputs[rel] = sha(path)
            inputs[str(line_path.relative_to(root))] = sha(line_path)
            for rec in lines.to_dict("records"):
                cid = rec["line_id"]
                ref = reference.loc[cid]
                now = current.loc[cid]
                verdict = grade_line(rec["wavelength_air_A"], rec["ep_eV"], rec["snapshot_loggf"], "Si I")
                line_records.append(dict(
                    holding=holding, treatment=row["treatment"], line_id=cid,
                    wavelength_air_A=rec["wavelength_air_A"], ep_eV=rec["ep_eV"],
                    abundance=rec["abundance"], in_aggregate=rec["in_aggregate"],
                    excluded_reason=rec.get("excluded_reason", ""),
                    snapshot_loggf=rec["snapshot_loggf"], current_loggf=float(now.log_gf),
                    reference_loggf=float(ref.published_loggf),
                    delta_reference_minus_snapshot=float(ref.published_loggf) - rec["snapshot_loggf"],
                    current_gf_unchanged=bool(float(now.log_gf) == rec["snapshot_loggf"]),
                    reference_membership=ref.selection_status,
                    current_grade_for_snapshot_value=verdict.gf_grade,
                    grade_sigma_is_publishable=verdict.gf_grade == "GF-LAB",
                    gf_source_note="Snapshot reconstructed from source checkout; run-time gf not serialized per line",
                    sigma_A=rec.get("sigma_A"), red_chi2=rec.get("red_chi2")))
    comparisons = []
    for left, right in itertools.combinations(locations, 2):
        comparisons.append(dict(kind="holding", left=left, right=right, treatment="1D-LTE",
                                **si.paired_difference(pools[(left, "1D-LTE")], pools[(right, "1D-LTE")])))
    for holding in locations:
        comparisons.append(dict(kind="NLTE_on_matched_lines", left="ENGINE-A", right="1D-LTE", holding=holding,
                                **si.paired_difference(pools[(holding, "ENGINE-A")], pools[(holding, "1D-LTE")])))
        comparisons.append(dict(kind="engine_label_check", left="ENGINE-B", right="1D-LTE", holding=holding,
                                **si.paired_difference(pools[(holding, "ENGINE-B")], pools[(holding, "1D-LTE")])))
    ion_profiles = []
    for path in sorted((root / "data/results/rya1218/recovered_20260916/si_ii_profile").glob("*ew.csv")):
        for row in si.attach_identity(pd.read_csv(path), canonical).to_dict("records"):
            ion_profiles.append({"artifact": str(path.relative_to(root)), "line_id": row["line_id"],
                                 "instrument": row["instrument"], "ew_mA": row["ew_mA"],
                                 "in_aggregate": bool(row["in_aggregate"]),
                                 "A": None, "status": "PROFILE_MEASURED_ABUNDANCE_NOT_INVERTED",
                                 "snapshot_loggf": row["snapshot_loggf"],
                                 "reference_loggf": float(reference.loc[row["line_id"], "published_loggf"])})
    grades = pd.read_csv(root / "data/audit/rya1218_si_protocol/si_full_grading_matrix.csv")
    # Read the current complete matrix, not historical ticket counts.
    summary = {"ticket": "RYA-1218", "related": ["RYA-725", "RYA-754", "RYA-587"],
               "products": products, "comparisons": comparisons, "si_ii_profiles": ion_profiles,
               "grade_counts": grades.recipe_grade.value_counts().to_dict(),
               "input_sha256": inputs, "publication_admitted": 0,
               "limits": ["No full uncertainty supplied by recovered artifacts.",
                          "Median products do not acquire SE(mean) as their validated uncertainty.",
                          "No gf values adjusted; source-scale offsets are diagnostic evidence only.",
                          "ENGINE-B LTE rows duplicate the 1D-LTE solution; not an independent engine validation."]}
    return summary, pd.DataFrame(line_records)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        ap.error("output already exists; preserve previous evidence")
    summary, lines = audit()
    a.out.mkdir(parents=True)
    (a.out / "product_evidence.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    lines.to_csv(a.out / "line_grade_abundance.csv", index=False)
    print(json.dumps({"products": len(summary["products"]), "line_records": len(lines),
                      "admitted": summary["publication_admitted"], "grade_counts": summary["grade_counts"]}, indent=2))


if __name__ == "__main__":
    main()
