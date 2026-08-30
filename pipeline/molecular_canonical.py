"""RYA-1130 transition-level molecular canonicalisation.

Molecular transition lists are not atomic laboratory-gf measurements.  This module
therefore never reads or writes ``canonical_gf.csv`` and never assigns Codex/Deep
grades.  It converts a Turbospectrum row into the dedicated, provenance-first schema;
missing source semantics remain explicit and prevent a frozen-ready verdict.
"""
from __future__ import annotations

import hashlib
import re
import shlex
from pathlib import Path

SCHEMA = "codex.molecular_transition/1"
PROVENANCE_FIELDS = (
    "database", "release", "retrieved_at", "license",
    "partition_function_source", "dissociation_energy_source",
    "component_convention", "isotopologue_convention",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def isotopologue_from_filename(path: Path) -> str:
    name = path.stem
    if name == "CO_IR_Li2015":
        return "12C16O"
    explicit = re.match(r"([0-9]+[A-Za-z]+)-([0-9]+[A-Za-z]+)__", name)
    if explicit:
        return "".join(explicit.groups())
    match = re.match(r"([0-9]+[A-Za-z]+(?:[0-9]+[A-Za-z]+)?)_", name)
    if not match:
        raise ValueError(f"cannot derive isotopologue from {path.name}")
    return match.group(1)


def bands_for_wavelength(wavelength_A: float) -> list[str]:
    # Reporting bins only; the physical wavelength remains authoritative.
    if wavelength_A < 2000:
        return ["FUV"]
    if wavelength_A < 4000:
        return ["NUV"]
    if wavelength_A < 7000:
        return ["VIS"]
    if wavelength_A < 10000:
        return ["RED_OPTICAL"]
    if wavelength_A < 25000:
        return ["NIR"]
    return ["IR"]


def _transition_quantum_fields(raw_label: str) -> dict:
    # Preserve the source label losslessly. Parse only the unambiguous ExoMol CO form;
    # heterogeneous PGopher labels remain raw until their release-specific decoder exists.
    match = re.search(r"v([^_]+)_J([0-9.]+)-([0-9.]+)_", raw_label)
    if match:
        vib, ju, jl = match.groups()
        vu, vl = (vib.split("-", 1) + [None])[:2]
        return {"electronic_system": None, "v_upper": vu, "v_lower": vl,
                "j_upper": ju, "j_lower": jl, "raw_label": raw_label}
    return {"electronic_system": None, "v_upper": None, "v_lower": None,
            "j_upper": None, "j_lower": None, "raw_label": raw_label}


def parse_bsyn_rows(path: Path):
    """Yield the numeric core and raw source identity from a .bsyn/.dat file."""
    with path.open(errors="replace") as fh:
        for line_no, line in enumerate(fh, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("'"):
                continue
            parts = shlex.split(stripped)
            if len(parts) < 15:
                raise ValueError(f"{path}:{line_no}: malformed molecular row")
            try:
                wavelength_A = float(parts[0])
                lower_energy_eV = float(parts[1])
                log_gf = float(parts[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: invalid numeric core") from exc
            # Some legacy C2 rows end in an empty quoted label.  The source omission is
            # data, not permission to synthesize an identity or discard the transition.
            raw_label = parts[-1].strip() or "SOURCE_LABEL_NOT_SUPPLIED"
            yield line_no, wavelength_A, lower_energy_eV, log_gf, raw_label


def canonical_rows(path: Path, molecule: str, metadata: dict):
    digest = sha256(path)
    isotopologue = isotopologue_from_filename(path)
    missing = [field for field in PROVENANCE_FIELDS if not metadata.get(field)]
    status = "BLOCKED_PROVENANCE" if missing else "SOURCE_VERIFIED"
    for line_no, wavelength_A, lower_energy_eV, log_gf, raw_label in parse_bsyn_rows(path):
        stable = f"{molecule}|{isotopologue}|{digest}|{line_no}|{raw_label}"
        transition_id = "mt_" + hashlib.sha256(stable.encode()).hexdigest()[:16]
        row_status = ("CROSSMATCH_REVIEW" if raw_label == "SOURCE_LABEL_NOT_SUPPLIED"
                      and status == "SOURCE_VERIFIED" else status)
        yield {
            "schema": SCHEMA,
            "transition_id": transition_id,
            "molecule": molecule,
            "isotopologue": isotopologue,
            "wavelength": {"value_A": wavelength_A, "medium": metadata.get("wavelength_medium", "unknown"),
                           "frame": metadata.get("wavelength_frame", "unknown"),
                           "conversion_provenance": metadata.get("conversion_provenance")},
            "lower_energy_eV": lower_energy_eV,
            "transition_probability": {"log_gf": log_gf, "kind": metadata.get("probability_kind", "unknown"),
                                       "source_id": metadata.get("transition_probability_source", "UNRESOLVED")},
            "transition_identity": _transition_quantum_fields(raw_label),
            "source": {key: metadata.get(key) for key in
                       ("database", "release", "table", "doi", "ads_bibcode", "catalog_id",
                        "retrieved_at", "license")}
                      | {"source_file": str(path), "source_sha256": digest},
            "thermochemistry": {"partition_function_source": metadata.get("partition_function_source"),
                                "dissociation_energy_source": metadata.get("dissociation_energy_source")},
            "normalization": {"component_convention": metadata.get("component_convention"),
                              "isotopologue_convention": metadata.get("isotopologue_convention")},
            "precedence": {"decision": metadata.get("precedence_decision"),
                           "conflict_note": metadata.get("conflict_note")},
            "band_reach": bands_for_wavelength(wavelength_A),
            "intake_status": row_status,
        }
