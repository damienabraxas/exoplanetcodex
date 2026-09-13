"""RYA-1218 evidence inventory. No abundance computation or grade promotion.

Coverage bounds are declarations, never a claim that pixels or corrected flux
were served. All unresolved work stays HOLD; absence from an intake is not
evidence that laboratory data do not exist.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Refusing empty inventory: {path}")
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def physical_matches(rows, species, wave, ep):
    return [r for r in rows if r["species"] == species
            and abs(float(r["wavelength_air_A"]) - wave) <= .03
            and abs(float(r["ep_eV"]) - ep) <= .01]


def holding_specs(path):
    """Read literal declarations without triggering the driver's atlas I/O.

    Fail if that declaration changes shape instead of silently losing holdings.
    This is an inventory adapter, not a second runtime coverage implementation.
    """
    tree = ast.parse(path.read_text())
    constants = {n.targets[0].id: n.value for n in tree.body
                 if isinstance(n, ast.Assign) and len(n.targets) == 1
                 and isinstance(n.targets[0], ast.Name)}

    def literal(node):
        if isinstance(node, ast.Name):
            return ast.literal_eval(constants[node.id])
        if isinstance(node, ast.Tuple):
            return tuple(literal(x) for x in node.elts)
        return ast.literal_eval(node)
    node = next(n.value for n in tree.body if isinstance(n, ast.AnnAssign)
                and isinstance(n.target, ast.Name) and n.target.id == "_INSTRUMENT_HOLDINGS")
    result = {}
    for key, value in zip(node.keys, node.values):
        for call in value.elts:
            if not isinstance(call, ast.Call) or call.func.id != "HoldingSpec":
                raise ValueError("Unrecognised holding declaration")
            fields = {k.arg: k.value for k in call.keywords}
            result[ast.literal_eval(call.args[0])] = {
                "instrument_id": ast.literal_eval(key),
                "reader": ast.literal_eval(fields["reader"]),
                "span": literal(fields["span_A"]) if "span_A" in fields else None,
            }
    return result


def build(root=ROOT):
    import yaml
    out = root / "data/audit/rya1218_si_protocol"
    out.mkdir(parents=True, exist_ok=True)
    paths = ["data/linelists/canonical_gf.csv", "data/linelists/linelist_solar.csv",
             "data/catalog/holdings_manifest_registry.csv", "data/catalog/model_registry.csv",
             "data/catalog/engine_coverage.csv", "config/synth_bands.yaml",
             "scripts/measure_band_ew.py",
             "data/audit/rya1169_si_intake/si_agss21_reference_lines.csv",
             "data/audit/rya1169_si_intake/si_primary_lab_gf_census.csv"]
    canonical = [r for r in read_csv(root / paths[0]) if r["species"] in ("Si I", "Si II")]
    refs, labs = (read_csv(root / p) for p in paths[-2:])
    bands = yaml.safe_load((root / "config/synth_bands.yaml").read_text())["bands"]

    def band(wave):
        hits = [name for name, b in bands.items() if b["lo_A"] <= wave < b["hi_A"]]
        return ";".join(hits) or "outside_configured_synthesis_bands"

    census = []
    for row in canonical:
        wave, ep = float(row["wavelength_air_A"]), float(row["excitation_potential_eV"])
        rr = physical_matches(refs, row["species"], wave, ep)
        ll = physical_matches(labs, row["species"], wave, ep)
        census.append({**row, "band": band(wave), "wavelength_frame": "air",
            "external_membership": ";".join(r["reference_status"] for r in rr) or "not_in_RYA1169_source_set",
            "reference_gf_source": ";".join(r["published_gf_source"] for r in rr),
            "reference_gf_sigma_dex": ";".join(r["published_gf_sigma_dex"] for r in rr),
            "DH23_matches": len(ll),
            "physical_identity_status": "HOLD_levels_and_J_not_in_canonical_schema",
            "reference_grade": "HOLD_best_available_per_line_adjudication",
            "codex_grade": "HOLD_primary_lab_and_observed_depth_audit",
            "deep_grade": "HOLD_primary_lab_and_observed_depth_audit",
            "measurement_eligibility": "HOLD", "abundance": ""})
    write_csv(out / "canonical_si_census.csv", census)

    gaps = []
    adapted = [{**r, "ep_eV": r["excitation_potential_eV"]} for r in canonical]
    for row in labs:
        hits = physical_matches(adapted, row["species"], float(row["wavelength_air_A"]), float(row["ep_eV"]))
        gaps.append({**row, "current_match_count": len(hits),
                     "current_candidate_ids": ";".join(r["line_id"] for r in hits),
                     "current_disposition": "candidate_identity_review" if len(hits) == 1 else
                         "canonical_gap" if not hits else "ambiguous_identity",
                     "gap_is_holding_absence": "NOT_INFERRED_FROM_CANONICAL_ABSENCE"})
    write_csv(out / "dh23_gap_recheck.csv", gaps)

    holdings = [r for r in read_csv(root / paths[2]) if r["system_id"] == "solar"]
    specs = holding_specs(root / "scripts/measure_band_ew.py")
    for row in holdings:
        spec = specs.get(row["holding_id"])
        row.update(reader=spec["reader"] if spec else "",
                   declared_min_A=spec["span"][0] if spec and spec["span"] else "",
                   declared_max_A=spec["span"][1] if spec and spec["span"] else "",
                   coverage_status="PIXELS_NOT_TESTED" if spec else "HOLD_no_band_harness_reader",
                   observed_conditioning="HOLD_exact_served_product_not_verified")
    write_csv(out / "solar_holdings.csv", holdings)
    # Include source rows missing from canonical rather than losing their cells.
    lines = [{"id": r["line_id"], "species": r["species"],
              "wave": float(r["wavelength_air_A"])} for r in canonical]
    lines += [{"id": f"DH23:{r['species']}:{r['wavelength_air_A']}", "species": r["species"],
               "wave": float(r["wavelength_air_A"])} for r in gaps if r["current_match_count"] == 0]
    matrix = []
    ew_path = "data/results/rya1218/source_set_ew_diagnostic/si_asplund_ew_per_line.csv"
    ew_rows = read_csv(root / ew_path)
    paths.append(ew_path)
    holding_aliases = {"solar_kpno_molecfit_corrected": "kpno_1984_molecfit_corrected",
                       "solar_kpno_kurucz2005_corrected": "kpno_kurucz2005_residual"}
    for line in lines:
        for h in holdings:
            lo, hi = h["declared_min_A"], h["declared_max_A"]
            outside = lo != "" and not lo <= line["wave"] <= hi
            cell = {"line_id": line["id"], "species": line["species"],
                "wavelength_air_A": line["wave"], "band": band(line["wave"]),
                "holding_id": h["holding_id"], "instrument_id": h["instrument_id"],
                "coverage": "outside_declared_span" if outside else "HOLD_pixel_coverage_not_tested",
                "telluric_applied_registry": h["telluric_applied"],
                "normalization_state_registry": h["normalization_state"],
                "observed_conditioning": h["observed_conditioning"],
                "depth": "", "blend_disposition": "HOLD_not_assessed",
                "route": "neither" if outside else "HOLD_EW_and_synthesis_not_adjudicated",
                "eligibility": "N/A" if outside else "HOLD",
                "reason": "outside_declared_span" if outside else
                    "line_identity_gf_grade_and_exact_product_evidence_pending"}
            # Join the runner evidence through the already physical-matched
            # reference set, including EP, never wavelength alone.
            source = next((r for r in canonical if r["line_id"] == line["id"]), None)
            matches = physical_matches(ew_rows, line["species"], line["wave"],
                                       float(source["excitation_potential_eV"])) if source else []
            matches = [r for r in matches if r["holding"] == holding_aliases.get(h["holding_id"])]
            cell.update(ew_mA="", evidence="")
            if len(matches) == 1:
                measured = matches[0]
                cell.update(evidence=ew_path, depth=measured["observed_core_depth"],
                            ew_mA=measured["ew_mA"])
                if measured["status"] == "served":
                    cell.update(coverage="in_range_pixels_served",
                                observed_conditioning="corrected_reader_pre_normalised",
                                route="EW_measured_abundance_HOLD",
                                reason="EW_diagnostic_only; gf_grade_and_model_adjudication_pending")
                else:
                    cell.update(coverage="HOLD_corrected_window_unserved",
                                reason=measured["reason_or_concern"])
            matrix.append(cell)
    write_csv(out / "line_holding_eligibility.csv", matrix)
    write_csv(out / "historical_grade_delta.csv", [{
        "line_id": r["canonical_line_id"], "species": r["species"],
        "historical_external_membership": r["reference_status"],
        "historical_asplund_label": "asplund" if r["reference_status"] == "used" else "excluded",
        "current_external_membership": r["reference_status"],
        "current_reference_grade": "HOLD_independent_best_available_adjudication",
        "reason": "external_membership_is_not_Reference_Grade"} for r in refs])
    engines = [r for r in read_csv(root / "data/catalog/engine_coverage.csv") if r.get("element") == "Si"]
    write_csv(out / "historical_engine_reach.csv", engines)
    models = read_csv(root / "data/catalog/model_registry.csv")
    model_cells = []
    for model in models:
        for species in ("Si I", "Si II"):
            for h in holdings:
                for b in bands:
                    model_cells.append({"model_id": model["model_id"],
                        "model_name": model["display_name"], "stored_token": model["stored_token"],
                        "species": species, "holding_id": h["holding_id"], "band": b,
                        "grade_scope": "Reference;Codex;Deep;external_replication (all unresolved)",
                        "status": "HOLD", "executed": False,
                        "reason": "grade_and_line_identity_preflight_incomplete; per_line_model_applicability_not_validated",
                        "abundance": "", "sigma_total": ""})
    write_csv(out / "model_inventory_cells.csv", model_cells)
    summary = {"ticket": "RYA-1218", "status": "PREFLIGHT_INCOMPLETE_MEASUREMENT_HOLD",
        "canonical_counts": dict(Counter(r["species"] for r in canonical)),
        "source_used_counts": dict(Counter(r["species"] for r in refs if r["reference_status"] == "used")),
        "DH23_canonical_gaps": sum(r["current_match_count"] == 0 for r in gaps),
        "solar_holdings": len(holdings), "line_holding_cells": len(matrix),
        "census_scope": "current canonical plus RYA1169 sources; not an exhaustive literature census",
        "reportable_products": 0,
        "source_sha256": {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in paths}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "source_sha256"}, indent=2))


if __name__ == "__main__":
    build()
