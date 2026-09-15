#!/usr/bin/env python3
"""Read-only uncertainty inventory and Reference xi recovery; never writes live feeds.

The Reference importer belongs to RYA-1213. Invoke that exact checkout to reuse
its run-stamp verification and pairing, and record its commit and file hashes.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.uncertainty_contract import COMPONENTS, publication_problems


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict): return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list): return [clean(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value): return None
    return value


def write(path, value):
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def write_legacy_table(path, inventory_doc):
    """Historical values only; absent derivatives never become measured zero."""
    with path.open("w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["element", "ion", "raw_sigma_legacy", "N_legacy", "sigma_stat_legacy",
                         "sigma_stellar_legacy_if_supported", "dominant_parameter_legacy", "canonical_state"])
        for element in inventory_doc["per_element"]:
            old = element["legacy_diagnostic_rows"]
            if not old:
                writer.writerow([element["element"], "", "", "", "", "", "", "HOLD: no measured product budget"])
            for row in old:
                supported = all(isinstance(row.get(k), (int, float)) and math.isfinite(row[k])
                                for k in ("dA_dTeff_per100K", "dA_dvmic_per_kms"))
                params = {k: row[k] for k in ("sigma_B_Teff", "sigma_B_logg", "sigma_B_vmic", "sigma_B_FeH")}
                sigma = math.sqrt(sum(v*v for v in params.values())) if supported else ""
                writer.writerow([row["element"], row.get("ion"), row.get("raw_sigma"), row.get("n_lines"),
                                 row.get("sigma_SE"), sigma, max(params, key=params.get) if supported else "",
                                 "HOLD: full component evidence absent"])


def inventory():
    legacy_path = ROOT / "data/audit/uncertainty/solar_uncertainty_rya158.json"
    legacy = json.loads(legacy_path.read_text())
    by_element = collections.defaultdict(list)
    for row in legacy["per_element"]:
        by_element[row["element"]].append(row)
    with (ROOT / "data/audit/element_status_tracker.csv").open() as f:
        tracker = list(csv.DictReader(line for line in f if not line.startswith("#")))
    elements = sorted({r["element"] for r in tracker})
    rows = []
    for element in elements:
        feed = ROOT / f"data/products/solar/{element}.json"
        products = json.loads(feed.read_text()).get("products", []) if feed.exists() else []
        problems = [{"product": {k: p.get(k) for k in
                     ("ion", "band", "holding", "tier", "route", "treatment")},
                     "problems": publication_problems(p)} for p in products]
        rows.append({"element": element, "legacy_diagnostic_rows": by_element[element],
                     "live_product_count": len(products), "product_audit": problems,
                     "canonical_complete_products": sum(not p["problems"] for p in problems),
                     "state": "HOLD", "reason": "No complete canonical component evidence supplied",
                     "required_components": list(COMPONENTS)})
    return {"ticket": "RYA-587", "base_sha": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "legacy_artifact_sha256": sha(legacy_path),
            "legacy_fe_rows": by_element["Fe"], "per_element": rows,
            "complete_budget_before": [], "complete_budget_after": [],
            "complete_budget_basis": "New full RYA-587 publication contract; historical Fe-only recipe is narrower",
            "reconciliation_module_exists": (ROOT / "pipeline/reconciliation_band.py").exists(),
            "legacy_warning": "Zero Type-B entries with absent derivatives are not proof of measured zero; preserved as historical diagnostics only"}


def reference_audit(repo, runs):
    # Subprocess isolates the unmerged importer's imports from this checkout.
    code = """import json,sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1])/'scripts'))
sys.path.insert(0, sys.argv[1])
from rya1213_xi_campaign import build
print(json.dumps(build(Path(sys.argv[2])), allow_nan=False))
"""
    fresh = json.loads(subprocess.check_output(
        [sys.executable, "-c", code, str(repo), str(runs)], cwd=repo, text=True))
    committed_path = repo / "data/results/rya1213/reference_xi_dadxi.json"
    old = json.loads(committed_path.read_text())
    feed_path = repo / "data/products/solar/Fe.json"
    feed = json.loads(feed_path.read_text())
    keys = ("ion", "holding", "tier", "treatment", "band")
    index = {tuple(p[k] for k in keys): p for p in fresh["pools"]}
    products, stamps = [], []
    for original in feed["products"]:
        if original.get("tier") != "REFERENCE": continue
        p = dict(original)
        q = index.get(tuple(p[k] for k in keys))
        verified = q is not None and q["xi_state"] == "MEASURED"
        stamps.append({**{k: p[k] for k in keys}, "previous_xi_state": p.get("xi_state"),
                       "previous_sigma_xi": p.get("sigma_xi"),
                       "recovered_sigma_xi": q["sigma_xi"] if verified else None,
                       "previous_stamp_matches": p.get("sigma_xi") == q["sigma_xi"] if verified else None})
        if verified:
            p["sigma_xi"] = q["sigma_xi"]
            p["xi_state"] = q["xi_state"]
            p["dA_dxi_dex_per_kms"] = q["dA_dxi"]
            p["xi_source"] = "reference_xi_recomputed.json"
            p["xi_note"] = q["xi_note"]
            p["sigma_syst_components"] = {"published_syst": p["sigma_syst"], "sigma_xi": p["sigma_xi"]}
            p["sigma_syst_complete"] = round(math.hypot(p["sigma_syst"], p["sigma_xi"]), 6)
            # Retain the legacy subtotal explicitly; not a full-contract total.
            p["legacy_sigma_reported"] = round(math.hypot(p["sigma_stat"], p["sigma_syst_complete"]), 6)
        p["previous_sigma_reported"] = original.get("sigma_reported")
        p["sigma_reported"] = None
        p["uncertainty_state"] = "HOLD"
        p["uncertainty_hold_reason"] = "Recovered xi is a component; full RYA-587 evidence still required"
        products.append(p)
    return fresh, {"ticket": "RYA-587", "reference_sha": subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
                  "reference_artifact_sha256": sha(committed_path), "reference_feed_sha256": sha(feed_path),
                  "importer_sha256": sha(repo / "scripts/rya1213_xi_campaign.py"),
                  "committed_pool_count": old["n_pools"], "committed_measured_count": old["n_measured"],
                  "recomputed_pool_count": fresh["n_pools"], "recomputed_measured_count": fresh["n_measured"],
                  "completed_run_directories": len(list(runs.glob("*/DONE"))), "product_stamps": stamps,
                  "files": [{"path": str(p.relative_to(runs)), "sha256": sha(p)}
                            for p in sorted(runs.rglob("*")) if p.is_file()],
                  "limitations": ["Raw runs use +/-0.1 km/s, not adopted +/-0.2912 km/s; nonlinear-response assessment remains owed",
                                   "Legacy nearest-wavelength pairing is preserved; exact nominal physical pool proof remains owed",
                                   "Mean-3D Gerber excluded from completion requirements per RYA-587 directive",
                                   "Reference branch is not integrated into main"]}, products


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference-repo", type=Path)
    parser.add_argument("--xi-runs", type=Path)
    args = parser.parse_args()
    if args.output.exists(): parser.error("output directory exists; preserve prior audit, choose a new directory")
    if bool(args.reference_repo) != bool(args.xi_runs): parser.error("supply both Reference checkout and xi runs")
    doc = inventory()
    reference = reference_audit(args.reference_repo, args.xi_runs) if args.xi_runs else None
    args.output.mkdir(parents=True)
    write(args.output / "element_inventory.json", doc)
    write_legacy_table(args.output / "legacy_budget_table.csv", doc)
    if reference:
        fresh, audit, products = reference
        write(args.output / "reference_xi_recomputed.json", fresh)
        write(args.output / "reference_xi_audit.json", audit)
        write(args.output / "reference_products_review.json", {"ticket": "RYA-587", "publication_state": "HOLD", "products": [], "review_products": products})
        print(f"Reference: {audit['committed_measured_count']} -> {fresh['n_measured']} measured xi pools")
    print(f"Audited {len(doc['per_element'])} elements; full-contract complete set: {doc['complete_budget_after']}")
    print(f"Wrote audit only: {args.output}; no live feeds or spectra modified")


if __name__ == "__main__":
    main()
