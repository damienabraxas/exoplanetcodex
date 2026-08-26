#!/usr/bin/env python3
"""Acquire and normalize Nandakumar et al. (2024) IGRINS CDS tables.

The CDS catalogue contains line-by-line abundances, not a transition list: it
does not publish excitation energies, level identities, log(gf), HFS, or the
gf references.  Accordingly, wavelength matches emitted here are discovery
candidates only and can never promote a canonical gf grade.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import statistics
import urllib.request
from dataclasses import dataclass
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "literature" / "igrins_nandakumar_2024"
RAW = OUT / "raw"
BASE_URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/684/A15"
DOI = "10.1051/0004-6361/202348462"
BIBCODE = "2024A&A...684A..15N"
ELEMENTS = ("f", "na", "mg", "al", "si", "s", "k", "ca", "sc", "ti", "v", "cr", "mn", "co", "ni", "cu", "zn", "y", "ce", "nd", "yb")
FILES = ("ReadMe", "stars.dat", *(f"table_{e}.dat" for e in ELEMENTS))


@dataclass(frozen=True)
class Column:
    start: int
    end: int
    label: str
    explanation: str


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def acquire() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        destination = RAW / name
        with urllib.request.urlopen(f"{BASE_URL}/{name}") as response:
            payload = response.read()
        destination.write_bytes(payload)


def read_schemas(readme: str) -> dict[str, list[Column]]:
    schemas: dict[str, list[Column]] = {}
    current: str | None = None
    row_re = re.compile(r"^\s*(\d+)(?:-\s*(\d+))?\s+\S+\s+\S+\s+(\S+)\s+(.*)$")
    for line in readme.splitlines():
        match = re.match(r"Byte-by-byte Description of file: (\S+)", line)
        if match:
            current = match.group(1)
            schemas[current] = []
            continue
        if current:
            match = row_re.match(line)
            if match:
                start, end = int(match.group(1)), int(match.group(2) or match.group(1))
                schemas[current].append(Column(start, end, match.group(3), match.group(4)))
    return schemas


def parse_fixed(path: Path, columns: list[Column]) -> list[dict[str, str]]:
    rows = []
    for line in path.read_text(encoding="ascii").splitlines():
        rows.append({col.label: line[col.start - 1 : col.end].strip() for col in columns})
    return rows


def line_columns(columns: list[Column]) -> list[tuple[Column, str, float]]:
    found = []
    pattern = re.compile(r"^A([A-Z][a-z]?)(\d+(?:\.\d+)?)$")
    for col in columns:
        match = pattern.match(col.label)
        if match:
            found.append((col, match.group(1), float(match.group(2))))
    return found


def band(wavelength: float) -> str:
    if 14900 <= wavelength < 18000:
        return "H"
    if 19500 <= wavelength <= 24000:
        return "K"
    return "H/K-gap-or-edge"


def load_canonical() -> list[dict[str, str]]:
    path = REPO / "data" / "linelists" / "canonical_gf.csv"
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def candidate_matches(element: str, wavelength: float, canonical: list[dict[str, str]]) -> list[dict[str, str]]:
    matches = []
    for row in canonical:
        if row.get("species", "").split()[0].lower() != element.lower():
            continue
        try:
            gap = abs(float(row["wavelength_air_A"]) - wavelength)
        except (KeyError, TypeError, ValueError):
            continue
        if gap <= 0.20:
            matches.append((gap, row))
    return [row for _, row in sorted(matches, key=lambda item: item[0])]


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def normalize() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    schemas = read_schemas((RAW / "ReadMe").read_text(encoding="ascii"))
    canonical = load_canonical()
    observations: list[dict[str, object]] = []
    lines: list[dict[str, object]] = []
    for element in ELEMENTS:
        filename = f"table_{element}.dat"
        columns = schemas[filename]
        records = parse_fixed(RAW / filename, columns)
        for col, symbol, wavelength in line_columns(columns):
            values = []
            for record in records:
                raw_value = record[col.label]
                value = None if not raw_value or raw_value == "-9.99" else float(raw_value)
                observations.append({
                    "species": symbol, "wavelength_air_A": f"{wavelength:.2f}",
                    "star": record["Star"], "abundance_dex": "" if value is None else f"{value:.2f}",
                    "used": value is not None, "source_table": filename,
                })
                if value is not None:
                    values.append(value)
            candidates = candidate_matches(symbol, wavelength, canonical)
            lines.append({
                "species": symbol, "element": symbol, "ion": "",
                "wavelength_air_A": f"{wavelength:.2f}", "wavelength_system": "air (CDS label; paper convention)",
                "band": band(wavelength), "n_stars_used": len(values), "n_stars_total": len(records),
                "median_abundance_dex": "" if not values else f"{statistics.median(values):.3f}",
                "source_table": filename, "doi": DOI, "bibcode": BIBCODE,
                "excitation_potential_eV": "", "lower_level": "", "upper_level": "",
                "log_gf_igrins": "", "gf_source_igrins": "", "hfs_isotopic_detail": "",
                "canonical_wavelength_candidates": ";".join(r["line_id"] for r in candidates),
                "crossmatch_status": "HOLD_TRANSITION_IDENTITY_INCOMPLETE" if candidates else "NO_CANONICAL_WAVELENGTH_CANDIDATE",
                "promotion_allowed": False,
            })
    write_csv(OUT / "igrins_nandakumar_2024_lines.csv", list(lines[0]), lines)
    write_csv(OUT / "igrins_nandakumar_2024_observations.csv", list(observations[0]), observations)
    write_delta(lines, canonical)
    write_al_audit(lines)
    write_provenance()


def write_delta(lines: list[dict[str, object]], canonical: list[dict[str, str]]) -> None:
    rows = []
    for element in sorted({str(row["element"]) for row in lines}):
        selected = [row for row in lines if row["element"] == element]
        with_candidates = sum(bool(row["canonical_wavelength_candidates"]) for row in selected)
        bands = sorted({str(row["band"]) for row in selected})
        rows.append({
            "element": element, "igrins_lines": len(selected), "canonical_wavelength_candidates": with_candidates,
            "transition_confirmed_matches": 0, "newly_discovered": len(selected) - with_candidates,
            "primary_lab_gradeable": 0, "weaker_gf_candidate": len(selected),
            "crires_bands_reachable": ";".join(bands), "telluric_risk": "REQUIRES_PER_SPECTRUM_VERIFICATION",
            "verdict": "DISCOVERY_ONLY_TRANSITION_IDENTITY_AND_GF_PROVENANCE_REQUIRED",
        })
    write_csv(OUT / "igrins_element_delta.csv", list(rows[0]), rows)


def write_al_audit(lines: list[dict[str, object]]) -> None:
    path = REPO / "data" / "results" / "rya1001" / "rya1001_al_line_census.csv"
    with path.open(newline="") as stream:
        census = list(csv.DictReader(stream))
    rows = []
    for line in (row for row in lines if row["element"] == "Al"):
        wavelength = float(line["wavelength_air_A"])
        near = []
        for existing in census:
            try:
                gap = abs(float(existing["wave_air_A"]) - wavelength)
            except (KeyError, TypeError, ValueError):
                continue
            if gap <= 0.20:
                near.append((gap, existing))
        near.sort(key=lambda item: item[0])
        best = near[0][1] if near else {}
        rows.append({
            "wavelength_air_A": f"{wavelength:.2f}", "igrins_n_stars_used": line["n_stars_used"],
            "rya1001_wavelength_candidate_A": best.get("wave_air_A", ""),
            "rya1001_tier": best.get("tier", ""), "rya1001_hfs_carried": best.get("hfs_carried", ""),
            "burheim_transition": best.get("burheim_transition", ""),
            "classification": "WAVELENGTH_CANDIDATE_ALREADY_IN_CENSUS_HOLD_IDENTITY" if best else "NEW_DISCOVERY_HOLD_IDENTITY_GF_HFS",
            "promotion_allowed": False,
        })
    write_csv(OUT / "igrins_al_completeness_audit.csv", list(rows[0]), rows)


def write_provenance() -> None:
    files = {name: {"sha256": sha256(RAW / name), "bytes": (RAW / name).stat().st_size} for name in FILES}
    payload = {
        "source": "CDS/VizieR J/A+A/684/A15", "article_doi": DOI, "bibcode": BIBCODE,
        "source_url": BASE_URL, "retrieved_utc_date": date.today().isoformat(),
        "raw_files_immutable": True, "files": files,
        "scientific_boundary": "CDS supplies empirical line usage/abundances but no EP, levels, log(gf), gf provenance, or HFS. Wavelength-only candidates are never grade promotions.",
    }
    (OUT / "provenance.json").write_text(json.dumps(payload, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquire", action="store_true", help="download and overwrite the raw CDS package")
    args = parser.parse_args()
    if args.acquire:
        acquire()
    missing = [name for name in FILES if not (RAW / name).exists()]
    if missing:
        raise SystemExit(f"missing raw files ({', '.join(missing)}); run with --acquire")
    normalize()


if __name__ == "__main__":
    main()
