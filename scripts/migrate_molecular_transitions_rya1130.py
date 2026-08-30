#!/usr/bin/env python3
"""Audit or stream Turbospectrum assets into codex.molecular_transition/1 JSONL."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.molecular_canonical import PROVENANCE_FIELDS, canonical_rows  # noqa: E402
from pipeline.molecular_lists import VENDORED_DIR, load_manifest  # noqa: E402


def metadata_for(entry: dict) -> dict:
    # Existing v1 prose is retained as a source label, never mistaken for typed provenance.
    return {
        "transition_probability_source": entry.get("source"),
        "database": entry.get("database"), "release": entry.get("release"),
        "table": entry.get("table"), "doi": entry.get("doi"),
        "ads_bibcode": entry.get("ads_bibcode"), "catalog_id": entry.get("catalog_id"),
        "retrieved_at": entry.get("retrieved_at"), "license": entry.get("license"),
        "wavelength_medium": entry.get("wavelength_medium", "unknown"),
        "wavelength_frame": entry.get("wavelength_frame", "unknown"),
        "conversion_provenance": entry.get("conversion_provenance"),
        "probability_kind": entry.get("probability_kind", "unknown"),
        "partition_function_source": entry.get("partition_function_source"),
        "dissociation_energy_source": entry.get("dissociation_energy_source"),
        "component_convention": entry.get("component_convention"),
        "isotopologue_convention": entry.get("isotopologue_convention"),
        "precedence_decision": entry.get("precedence_decision"),
        "conflict_note": entry.get("conflict_note"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, help="optional transition JSONL output")
    parser.add_argument("--report", type=Path, help="write the machine-readable audit report")
    parser.add_argument("--allow-blocked", action="store_true",
                        help="write BLOCKED_PROVENANCE rows for review; never frozen-ready")
    args = parser.parse_args()

    manifest = load_manifest()
    counts, missing_by_molecule = Counter(), {}
    output = args.output.open("w") if args.output else None
    try:
        for molecule, entry in manifest["molecules"].items():
            metadata = metadata_for(entry)
            missing = [field for field in PROVENANCE_FIELDS if not metadata.get(field)]
            missing_by_molecule[molecule] = missing
            for filename in entry["files"]:
                path = VENDORED_DIR / entry["vendored_subdir"] / filename
                for row in canonical_rows(path, molecule, metadata):
                    counts[(molecule, row["intake_status"])] += 1
                    if output:
                        if row["intake_status"] == "BLOCKED_PROVENANCE" and not args.allow_blocked:
                            raise RuntimeError("refusing to emit blocked rows without --allow-blocked")
                        output.write(json.dumps(row, sort_keys=True) + "\n")
    finally:
        if output:
            output.close()

    report = {"schema": "codex.molecular_migration_audit/1",
              "counts": {f"{molecule}:{status}": n for (molecule, status), n in sorted(counts.items())},
              "missing_provenance": missing_by_molecule,
              "frozen_ready": not any(missing_by_molecule.values())}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered)
    print(rendered, end="")
    return 0 if report["frozen_ready"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
