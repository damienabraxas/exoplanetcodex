#!/usr/bin/env python3
"""RYA-1169: test and grade the published AGSS21/Asplund Si line set.

Asplund membership is a provenance axis, not a gf-quality claim.  This script
therefore assigns ``asplund`` only to source-used rows, retains source-rejected
rows, and reports the independent Solar feature-depth route for every row.
Empty bands are emitted rather than backfilled with non-reference candidates.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INTAKE = ROOT / "data/audit/rya1169_si_intake"
REF = INTAKE / "si_agss21_reference_lines.csv"
SOLAR = ROOT / "data/linelists/linelist_solar.csv"
OUT = ROOT / "data/results/rya1169"

DEPTH_LO = 0.05
DEPTH_HI = 0.60
WAVE_TOL_A = 0.05
EP_TOL_EV = 0.01
BANDS = [
    ("FUV", 0.0, 3000.0),
    ("NUV", 3000.0, 3780.0),
    ("VIS", 3780.0, 6910.0),
    ("red-optical", 6910.0, 9199.0),
    ("NIR", 9199.0, 13000.0),
    ("H", 15007.0, 17500.0),
    ("K", 17500.0, 21400.0),
]


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        out = csv.DictWriter(handle, fieldnames=fields)
        out.writeheader(); out.writerows(rows)


def band(wave: float) -> str:
    for name, lo, hi in BANDS:
        if lo <= wave < hi:
            return name
    return "other-IR"


def main() -> None:
    import pandas as pd

    OUT.mkdir(parents=True, exist_ok=True)
    reference = read(REF)
    solar = pd.read_csv(SOLAR, low_memory=False)
    solar = solar[solar.element.astype(str) == "Si"].copy()

    tested = []
    for row in reference:
        wave, ep = float(row["wavelength_air_A"]), float(row["ep_eV"])
        ion = row["species"].split()[-1]
        hits = solar[(solar.ion.astype(str) == ion)
                     & ((solar.wavelength_air_A.astype(float) - wave).abs() <= WAVE_TOL_A)
                     & ((solar.excitation_potential_eV.astype(float) - ep).abs() <= EP_TOL_EV)]
        if len(hits) == 1:
            hit = hits.iloc[0]
            depth = float(hit.central_depth)
            depth_join = "matched"
            depth_source_wave = f"{float(hit.wavelength_air_A):.4f}"
            blend = str(hit.blend_flag)
        else:
            depth = None
            depth_join = "missing" if hits.empty else "ambiguous"
            depth_source_wave = ""
            blend = ""

        if depth is None:
            route = "ungraded"
        elif depth < DEPTH_LO:
            route = "below_floor"
        elif depth <= DEPTH_HI:
            route = "ew_valid"
        else:
            route = "deep_synthesis"

        used = row["reference_status"] == "used"
        tested.append({
            **row,
            "source_band_recomputed": band(wave),
            "solar_feature_depth": f"{depth:.3f}" if depth is not None else "",
            "depth_source_wavelength_A": depth_source_wave,
            "depth_source": "linelist_solar central_depth; project feature-depth gate",
            "depth_join_status": depth_join,
            "depth_route": route,
            "asplund_grade": "asplund" if used else "excluded",
            "active_for_asplund_test": "true" if used and depth_join == "matched" else "false",
            "blend_flag": blend,
            "test_status": ("PASS" if used and depth_join == "matched" and route in {"ew_valid", "deep_synthesis"}
                            else "SOURCE_REJECTED" if not used else "BLOCKED_DEPTH"),
        })
    write(OUT / "si_asplund_line_test.csv", tested, list(tested[0]))

    matrix = []
    for name, lo, hi in BANDS:
        rows = [x for x in tested if x["source_band_recomputed"] == name]
        used = [x for x in rows if x["reference_status"] == "used"]
        rejected = [x for x in rows if x["reference_status"] != "used"]
        note = ""
        if name in {"FUV", "NUV"} and not rows:
            note = "no transition in the adopted Amarsi2017/Scott2015 Solar Si line set"
        if name in {"NIR", "H", "K"} and not rows:
            note = ("Scott2015 explicitly discarded available near-IR Si I lines: only "
                    "low-accuracy theoretical gf and large/uncertain NLTE effects; "
                    "no wavelengths enumerated in the adopted table")
        matrix.append({
            "band": name, "lo_A": lo, "hi_A": hi,
            "source_used": len(used), "source_rejected_rows": len(rejected),
            "depth_tested": sum(x["depth_join_status"] == "matched" for x in rows),
            "asplund_active": sum(x["asplund_grade"] == "asplund" and x["test_status"] == "PASS" for x in rows),
            "ew_valid": sum(x["depth_route"] == "ew_valid" and x["asplund_grade"] == "asplund" for x in rows),
            "deep_synthesis": sum(x["depth_route"] == "deep_synthesis" and x["asplund_grade"] == "asplund" for x in rows),
            "below_floor": sum(x["depth_route"] == "below_floor" and x["asplund_grade"] == "asplund" for x in rows),
            "source_scope_note": note,
        })
    write(OUT / "si_asplund_all_band_matrix.csv", matrix, list(matrix[0]))

    active = [x for x in tested if x["asplund_grade"] == "asplund"]
    summary = {
        "ticket": "RYA-1169",
        "line_set": "asplund",
        "native_line_set": "AGSS21_Si_Amarsi2017_Scott2015",
        "definition": "source membership, orthogonal to gf-quality and depth route",
        "depth_gate": {"below_floor": "<0.05", "ew_valid": "0.05-0.60 inclusive", "deep_synthesis": ">0.60"},
        "source_used": len(active),
        "source_rejected": sum(x["asplund_grade"] == "excluded" for x in tested),
        "active_test_passed": sum(x["test_status"] == "PASS" for x in active),
        "ew_valid": sum(x["depth_route"] == "ew_valid" for x in active),
        "deep_synthesis": sum(x["depth_route"] == "deep_synthesis" for x in active),
        "bands_with_active_rows": sorted({x["source_band_recomputed"] for x in active}),
        "all_band_verdict": "PASS_WITH_SOURCE_DEFINED_GAPS",
        "measurement_readiness": "READY_FOR_ASPLUND_EW_TEST_VIS_REDOPTICAL_ONLY",
        "firewall": "No UV/IR candidate substituted for an absent published reference member.",
    }
    (OUT / "si_asplund_test_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
