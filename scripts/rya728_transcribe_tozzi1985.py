#!/usr/bin/env python3
"""RYA-728 -- Cr I lab gf from Tozzi, Brunner & Huber 1985 (MNRAS 217, 423), Table 1, transcribed by
hand from the scanned paper (pp. 427-429; local "Reference documents/Tozzi_1985_MNRAS_217_423_CrI.pdf")
into data/reference/lab_gf_tables/transcribed/tozzi1985_cr1.csv.

Branching fractions measured on the Kitt Peak FTS with laser lifetimes of Kwiatkowski et al.; the
table prints air wavelength (nm), the upper level and its energy (eV), and log gf with its
uncertainty in dex. The three upper limits ("<") are not transcribed. The table prints no lower-level
energy, so the lower and upper levels come from the single VALD solar line at that wavelength; the
upper-level energy Tozzi prints is the check on each transcribed row (|dE| <= 0.01 eV), and a row that
fails it, or has no single VALD line, is dropped and counted. A is derived from log gf for the
schema only. Only 3000 A and redder is kept.
Usage: rya728_transcribe_tozzi1985.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rya728_parse_gurell2010 import COLS, LO_A, VALD  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reference" / "lab_gf_tables" / "transcribed" / "tozzi1985_cr1.csv"
DOI = "10.1093/mnras/217.2.423"
VALD_TOL_A = 0.03
UPPER_TOL_EV = 0.01

# (upper level, upper energy eV as printed): [(air wavelength nm, log gf, sigma dex), ...]
TABLE1 = {
    ("z5P1", 3.323): [(520.451, -0.22, 0.02), (524.757, -1.65, 0.06), (526.416, -1.31, 0.05),
                      (529.669, -1.41, 0.04)],
    ("z5P2", 3.322): [(373.080, -2.78, 0.05), (520.604, -0.01, 0.02), (526.572, -1.74, 0.03),
                      (529.828, -1.16, 0.05), (534.580, -0.98, 0.03)],
    ("z5P3", 3.321): [(373.202, -2.63, 0.05), (520.842, 0.15, 0.02), (530.074, -2.12, 0.04),
                      (534.831, -1.29, 0.03), (540.977, -0.74, 0.04)],
    ("y5P1", 3.648): [(458.004, -1.65, 0.04), (461.336, -1.70, 0.04), (462.617, -1.31, 0.03),
                      (465.128, -1.44, 0.03)],
    ("y5P2", 3.668): [(337.910, -3.07, 0.06), (454.594, -1.36, 0.03), (459.139, -1.74, 0.03),
                      (461.612, -1.18, 0.02), (465.215, -1.03, 0.02)],
    ("y5P3", 3.698): [(335.196, -2.86, 0.03), (449.684, -1.13, 0.03), (456.550, -2.04, 0.08),
                      (460.074, -1.24, 0.03), (464.615, -0.70, 0.02)],
    ("z5G2", 5.271): [(454.460, -0.58, 0.02), (454.533, -1.28, 0.04), (895.576, -1.92, 0.15)],
    ("z5G3", 5.274): [(453.976, -1.16, 0.04), (454.049, -0.52, 0.03), (454.106, -1.16, 0.04),
                      (893.921, -2.08, 0.15)],
    ("z5G4", 5.277): [(453.513, -0.99, 0.04), (453.509, -0.38, 0.03), (453.575, -1.16, 0.05)],
    ("z5G5", 5.280): [(452.984, -1.41, 0.05), (453.068, -0.91, 0.05), (453.073, -0.28, 0.03)],
    ("z5G6", 5.282): [(452.844, -0.15, 0.02), (452.733, -1.07, 0.03), (891.624, -1.41, 0.12)],
    ("x3I5", 6.904): [(314.844, -0.04, 0.03), (316.061, -1.35, 0.18), (405.876, 0.26, 0.03)],
    ("x3I6", 6.911): [(315.515, 0.04, 0.03), (316.957, -1.37, 0.20), (404.879, 0.31, 0.03),
                      (404.977, -0.75, 0.11), (454.865, -0.95, 0.20), (503.468, -1.08, 0.25)],
    ("x3I7", 6.918): [(316.375, 0.13, 0.03), (403.908, 0.39, 0.03), (403.930, -0.66, 0.11),
                      (503.256, -0.94, 0.18)],
}


def build(vald: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rows, counts = [], {"transcribed": 0, "bluer_than_3000": 0, "no_single_vald_line": 0,
                        "upper_level_mismatch": 0}
    for (_, e_up), lines in TABLE1.items():
        for nm, loggf, sig in lines:
            counts["transcribed"] += 1
            w = nm * 10.0
            if w < LO_A:
                counts["bluer_than_3000"] += 1
                continue
            lv = vald[(vald.wavelength_air_A - w).abs() <= VALD_TOL_A]
            lv = lv[(lv.ep_upper_eV - e_up).abs() <= UPPER_TOL_EV]
            if len(lv) != 1:
                counts["no_single_vald_line" if lv.empty or len(lv) > 1 else "upper_level_mismatch"] += 1
                continue
            v = lv.iloc[0]
            g_up = 2.0 * float(v.j_upper) + 1.0
            wa = float(v.wavelength_air_A)
            A = 10 ** loggf / (1.49919e-16 * wa ** 2 * g_up) / 1e6
            rows.append({"paper": "tozzi1985_cr1", "species": "Cr I", "wavelength_air_A": round(wa, 3),
                         "ep_eV": float(v.ep_eV), "j_lower": float(v.j_lower),
                         "ep_upper_eV": float(v.ep_upper_eV), "A_1e6_s": round(A, 4),
                         "dA_1e6_s": round(A * (10 ** sig - 1), 4), "loggf": loggf,
                         "loggf_sigma_dex": sig, "doi": DOI, "vizier": ""})
    counts["written"] = len(rows)
    return pd.DataFrame(rows, columns=COLS), counts


def main() -> None:
    vald = pd.read_csv(VALD)
    out, counts = build(vald[vald.species == "Cr I"])
    out.to_csv(OUT, index=False)
    print(counts)


if __name__ == "__main__":
    main()
