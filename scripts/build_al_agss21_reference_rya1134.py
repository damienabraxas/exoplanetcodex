#!/usr/bin/env python3
"""Build the Scott -> Nordlander -> AGSS21 Solar Al reference lineage."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/audit/rya1132_al_intake/al_line_manifest.csv"
OUT = ROOT / "data/linelists/reference_sets/agss21_solar_al_rya1134.csv"
PROVENANCE = ROOT / "data/reference/asplund2021_al/provenance_rya1134.json"

# Scott et al. (2015), Table 1. Nordlander & Lind (2017) adopt this selection
# and these weights, except for the telluric-contaminated 10891 A feature.
SCOTT_LINES = (
    ("alphys_I_6696.0150_0352", 6696.023, 3.143, -1.569, 35.80, 2, 6.465, True),
    ("alphys_I_6698.6730_0355", 6698.673, 3.143, -1.870, 21.00, 3, 6.451, True),
    ("alphys_I_7835.3090_0386", 7835.309, 4.022, -0.689, 41.00, 1, 6.372, True),
    ("alphys_I_8912.9000_0397", 8912.900, 4.085, -1.963, 3.00, 1, 6.385, True),
    ("alphys_I_10768.3650_0401", 10768.365, 4.085, -2.020, 4.00, 1, 6.452, True),
    ("alphys_I_10872.9730_0404", 10872.973, 4.085, -1.326, 15.30, 1, 6.379, True),
    ("alphys_I_10891.7360_0405", 10891.736, 4.087, -1.027, 31.60, 1, 6.427, False),
)

FIELDS = (
    "lineage_set", "species", "scott2015_wavelength_air_A", "scott2015_elow_eV",
    "scott2015_loggf", "scott2015_ew_mA", "scott2015_weight",
    "scott2015_3d_nlte_abundance", "nordlander2017_solar_selected",
    "agss21_adoption_status", "canonical_line_id", "manifest_wavelength_air_A",
    "match_distance_mA", "manifest_elow_eV", "manifest_loggf_adopted",
    "manifest_gf_source", "manifest_gf_grade", "manifest_intake_status",
    "grade_verification_status", "lineage_note",
)


def main() -> None:
    with MANIFEST.open(newline="") as stream:
        manifest = {row["canonical_line_id"]: row for row in csv.DictReader(stream)}

    rows = []
    for line_id, wave, ep, loggf, ew, weight, abundance, selected in SCOTT_LINES:
        source = manifest[line_id]
        assert source["species"] == "Al I"
        distance = abs(float(source["wavelength_air"]) - wave) * 1000.0
        assert distance <= 10.0, (line_id, distance)
        assert abs(float(source["lower_EP"]) - ep) <= 0.001
        status = "PENDING_RYA1134_GF_GRADE_VERIFICATION" if selected else "EXCLUDED_TELLURIC_10891"
        note = (
            "Physical match uses wavelength plus lower excitation/transition identity; "
            "Scott wavelength differs by 8 mA from the newer 6696.015 value."
            if line_id == "alphys_I_6696.0150_0352" else
            "Excluded by Nordlander & Lind (2017) and AGSS21 because of telluric contamination."
            if not selected else
            "Adopted by lineage: Scott 2015 selection/weight -> Nordlander 2017 3D-NLTE -> AGSS21 abundance."
        )
        rows.append({
            "lineage_set": "Scott2015_Nordlander2017_AGSS21",
            "species": "Al I",
            "scott2015_wavelength_air_A": f"{wave:.3f}",
            "scott2015_elow_eV": f"{ep:.3f}",
            "scott2015_loggf": f"{loggf:.3f}",
            "scott2015_ew_mA": f"{ew:.2f}",
            "scott2015_weight": weight,
            "scott2015_3d_nlte_abundance": f"{abundance:.3f}",
            "nordlander2017_solar_selected": "Y" if selected else "N",
            "agss21_adoption_status": "ADOPTED_VIA_NORDLANDER2017" if selected else "EXCLUDED_TELLURIC",
            "canonical_line_id": line_id,
            "manifest_wavelength_air_A": source["wavelength_air"],
            "match_distance_mA": f"{distance:.1f}",
            "manifest_elow_eV": source["lower_EP"],
            "manifest_loggf_adopted": source["loggf_adopted"],
            "manifest_gf_source": source["gf_source"],
            "manifest_gf_grade": source["gf_grade"],
            "manifest_intake_status": source["intake_status"],
            "grade_verification_status": status,
            "lineage_note": note,
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    provenance = {
        "ticket": "RYA-1134",
        "purpose": "Solar Al grade-verification reference lineage; not a new abundance measurement",
        "agss21_abundance": {"A_Al": 6.43, "sigma_dex": 0.03, "scale": "3D-NLTE"},
        "counts": {"scott2015_retained": 7, "nordlander2017_after_telluric_exclusion": 6},
        "published_count_discrepancy": {
            "status": "OPEN_DOCUMENTED_SOURCE_INCONSISTENCY",
            "agss21_prose_count": 5,
            "traceable_selection_count": 6,
            "explanation": "AGSS21 says Nordlander used Scott's lines except 1089.1 nm and then calls the result five lines. Scott Table 1 contains seven, so the stated rule leaves six. No unreported second exclusion was found.",
        },
        "sources": [
            {"id": "Scott2015", "doi": "10.1051/0004-6361/201424109", "arxiv": "1405.0279", "role": "seven lines, weights, EW, line data and 3D+NLTE results"},
            {"id": "NordlanderLind2017", "doi": "10.1051/0004-6361/201730427", "arxiv": "1708.01949", "role": "3D-NLTE analysis; Scott selection/weights excluding 10891 A"},
            {"id": "AGSS21", "doi": "10.1051/0004-6361/202140445", "arxiv": "2105.01661", "role": "adopts Nordlander & Lind solar abundance"},
        ],
        "input_manifest": str(MANIFEST.relative_to(ROOT)),
        "output_lineset": str(OUT.relative_to(ROOT)),
    }
    PROVENANCE.parent.mkdir(parents=True, exist_ok=True)
    PROVENANCE.write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
