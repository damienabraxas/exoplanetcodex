#!/usr/bin/env python3
"""Build the RYA-1129 intake ledger and atomic C/N/O manifests.

This is an inventory operation only.  It never runs an abundance engine and it
never changes ``canonical_gf.csv``.  Molecular rows are deliberately excluded:
their incomplete representation is tracked by RYA-1130.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ELEMENTS = ROOT / "data/config/elements_master.json"
GF = ROOT / "data/linelists/canonical_gf.csv"
OUT = ROOT / "data/audit/rya1129_atomic_intake"
ATOMIC_NUMBERS = {"C": "6", "N": "7", "O": "8"}


def _text(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def _source_class(row: pd.Series) -> str:
    tier = _text(row.gf_tier).upper()
    grade = _text(row.nist_grade)
    if tier == "LAB":
        return "PRIMARY_LAB"
    if grade or tier.startswith("NIST"):
        return "EVALUATED_NIST"
    if tier in {"KURUCZ", "VALD3"}:
        return "FALLBACK"
    return "OTHER_COMPILED_OR_THEORETICAL"


def _status(row: pd.Series) -> str:
    adjudication = _text(row.adjudication_status).lower()
    if "ambig" in adjudication or "review" in adjudication:
        return "CROSSMATCH_REVIEW"
    if not _text(row.loggf_reference):
        return "SOURCE_REVIEW"
    return "INVENTORIED_NOT_FROZEN"


def _manifest(frame: pd.DataFrame) -> pd.DataFrame:
    result = pd.DataFrame()
    result["canonical_line_id"] = frame.line_id
    result["physical_transition_id"] = frame.physical_id
    result["species"] = frame.species
    result["ion"] = frame.ion
    result["wavelength_air_A"] = frame.wavelength_air_A
    result["medium"] = "air"
    result["excitation_potential_eV"] = frame.excitation_potential_eV
    result["loggf"] = frame.log_gf
    result["gf_sigma_dex"] = frame.gf_sigma_dex
    result["nist_grade"] = frame.nist_grade
    result["gf_tier"] = frame.gf_tier
    result["source_class"] = frame.apply(_source_class, axis=1)
    result["adopted_source"] = frame.loggf_reference
    result["source_doi"] = frame.gf_source_doi
    result["source_set_membership"] = frame.apply(
        lambda r: "|".join(
            name
            for name, present in (
                ("synth", r.in_synth),
                ("regions", r.in_regions),
                ("linelist", r.in_linelist),
            )
            if pd.notna(present) and bool(present)
        ),
        axis=1,
    )
    result["hfs_isotope_status"] = frame.hfs_n_components.apply(
        lambda n: "HFS_OR_ISOTOPE_COMPONENTS" if pd.notna(n) and float(n) > 1 else "NONE_RECORDED"
    )
    result["intake_status"] = frame.apply(_status, axis=1)
    result["notes"] = "Atomic inventory only; not an abundance-quality selection."
    return result.sort_values(["species", "wavelength_air_A", "excitation_potential_eV"])


def build(out_dir: Path = OUT) -> None:
    registry = json.loads(ELEMENTS.read_text())
    targets = registry["elements"]
    symbols = list(dict.fromkeys(e["symbol"].split()[0] for e in targets))
    gf = pd.read_csv(GF, low_memory=False)

    out_dir.mkdir(parents=True, exist_ok=True)
    for element, number in ATOMIC_NUMBERS.items():
        rows = gf[gf.key_z.astype(str) == number].copy()
        _manifest(rows).to_csv(out_dir / f"{element}_atomic_manifest.csv", index=False)

    ledger_rows = []
    for symbol in symbols:
        rows = gf[gf.key_z.astype(str) == str(next(e["atomic_number"] for e in targets if e["symbol"].split()[0] == symbol))]
        classes = rows.apply(_source_class, axis=1) if len(rows) else pd.Series(dtype=str)
        status = "CROSSMATCH_REVIEW" if symbol in ATOMIC_NUMBERS else "NOT_STARTED"
        ledger_rows.append(
            {
                "element": symbol,
                "species": "|".join(sorted(rows.species.dropna().astype(str).unique())),
                "intake_status": status,
                "primary_source_status": f"PRIMARY_LAB={int((classes == 'PRIMARY_LAB').sum())};EVALUATED_NIST={int((classes == 'EVALUATED_NIST').sum())};FALLBACK={int((classes == 'FALLBACK').sum())}",
                "canonical_gf_status": f"{len(rows)} atomic rows inventoried",
                "HFS_status": f"{int((pd.to_numeric(rows.hfs_n_components, errors='coerce').fillna(0) > 1).sum())} component rows" if len(rows) else "NO_ROWS",
                "band_coverage_status": "REPORT_PENDING" if symbol not in ATOMIC_NUMBERS else "ATOMIC_WAVELENGTH_INVENTORY_COMPLETE",
                "source_ticket": "RYA-1129" if symbol in ATOMIC_NUMBERS else "RECONCILE_EXISTING_TICKET",
                "last_verified": "2026-08-29",
                "notes": "CNO blocked by molecular provenance schema RYA-1130" if symbol in ATOMIC_NUMBERS else "Phase A wave not yet executed",
            }
        )
    pd.DataFrame(ledger_rows).to_csv(out_dir / "intake_status_ledger.csv", index=False)

    metadata = {
        "ticket": "RYA-1129",
        "generated_from": [str(ELEMENTS.relative_to(ROOT)), str(GF.relative_to(ROOT))],
        "registry_version": registry["version"],
        "registry_n_targets_declared": registry["n_targets"],
        "registry_n_rows_observed": len(targets),
        "registry_n_unique_atomic_symbols": len(symbols),
        "registry_symbols": symbols,
        "known_delta": "Zn is described elsewhere as the 28th canonical element but is absent from this authoritative registry; no Zn ledger row was invented.",
        "molecular_scope": "STOPPED: RYA-1130",
        "abundances_generated": False,
    }
    (out_dir / "manifest_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=OUT)
    args = parser.parse_args()
    build(args.out_dir)


if __name__ == "__main__":
    main()
