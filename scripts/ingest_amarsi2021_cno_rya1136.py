#!/usr/bin/env python3
"""Normalize the 408-line Amarsi et al. (2021) Solar molecular CNO table."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HOLDING = ROOT / "data/reference/amarsi2021_cno"
RAW = HOLDING / "raw/table2.dat"
OUT = HOLDING / "derived/amarsi2021_cno_molecular_lines.csv"
MANIFEST = HOLDING / "manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    rows = []
    for line_number, line in enumerate(RAW.read_text().splitlines(), 1):
        if not line.strip():
            continue
        rows.append(
            {
                "reference_line_set": "AmarsiEtAl2021_AA656_A113_Table2",
                "adopted_solar_value_lineage": "AGSS21_CNO",
                "use_status": "USED",
                "source_row": line_number,
                "species": line[0:6].strip(),
                "system": line[7:11].strip(),
                "delta_nu": line[12:13].strip(),
                "band": line[14:19].strip(),
                "wavelength_vac_nm": line[20:29].strip(),
                "lower_energy_eV": line[30:35].strip(),
                "published_loggf": line[36:42].strip(),
                "equivalent_width_pm": line[43:49].strip(),
                "element_parameter": line[50:57].strip(),
                "abundance_3d": line[58:63].strip(),
                "abundance_mean_3d": line[64:69].strip(),
                "abundance_atmo": line[70:75].strip(),
                "abundance_marcs": line[76:81].strip(),
                "abundance_hm74": line[82:87].strip(),
                "transition_identity_status": "SOURCE_TABLE_LACKS_ROTATIONAL_QUANTUM_IDENTITY",
                "canonical_join_status": "CROSSMATCH_REVIEW",
            }
        )

    if len(rows) != 408:
        raise SystemExit(f"expected 408 rows, found {len(rows)}")
    counts = Counter(row["species"] for row in rows)
    expected = {"12C16O": 80, "C2": 39, "CH": 54, "CN": 59, "NH": 31, "OH": 145}
    if dict(counts) != expected:
        raise SystemExit(f"species counts differ: {dict(counts)}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)

    readme = HOLDING / "raw/ReadMe"
    manifest = {
        "schema": "codex.reference_holding/1",
        "ticket": "RYA-1136",
        "reference": "AmarsiEtAl2021_AA656_A113_Table2",
        "citation": "Amarsi, Grevesse, Asplund & Collet 2021, A&A, 656, A113",
        "doi": "10.1051/0004-6361/202141384",
        "ads_bibcode": "2021A&A...656A.113A",
        "catalog_id": "J/A+A/656/A113",
        "article_url": "https://doi.org/10.1051/0004-6361/202141384",
        "catalog_url": "https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A+A/656/A113",
        "retrieved_at": str(date.today()),
        "row_count": len(rows),
        "species_counts": dict(counts),
        "source_files": {
            "raw/ReadMe": sha256(readme),
            "raw/table2.dat": sha256(RAW),
        },
        "limitations": [
            "Table 2 identifies system and vibrational band but omits rotational quantum labels.",
            "A physical-identity join to primary molecular releases remains CROSSMATCH_REVIEW.",
            "The table contains used lines only; rejected synthesis features require upstream-source review.",
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
