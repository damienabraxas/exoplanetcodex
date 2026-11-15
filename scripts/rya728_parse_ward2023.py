#!/usr/bin/env python3
"""RYA-728 -- Cr II lab gf from Ward et al. 2023 (ApJ 959, 8), NIST FTS branching fractions with
lifetimes, parsed from the journal's machine-readable Table 5 (vendored at
data/reference/lab_gf_tables/journal/ApJ_959_8/apjacfafet5_mrt.txt) into
data/reference/lab_gf_tables/transcribed/ward2023_cr2.csv.

The table gives Ritz air wavelength (nm), both level energies (cm-1) and designations with J, the
experimental A-value and its one-sigma uncertainty, so no line-list join is needed for the levels.
log gf = log10(1.49919e-16 * lambda_A^2 * g_up * A), g_up = 2 J_up + 1. sigma(log gf) =
log10(1 + dA/A), the convention of the Gurell and Nilsson parsers. Only 3000 A and redder is kept.
The authors' semiempirical A (calc column) is theory and is NOT used.
Usage: rya728_parse_ward2023.py
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rya728_parse_gurell2010 import COLS, LO_A  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MRT = ROOT / "data" / "reference" / "lab_gf_tables" / "journal" / "ApJ_959_8" / "apjacfafet5_mrt.txt"
OUT = ROOT / "data" / "reference" / "lab_gf_tables" / "transcribed" / "ward2023_cr2.csv"
DOI = "10.3847/1538-4357/acfafe"
CM1_TO_EV = 1.0 / 8065.543937


def _j(designation: str) -> float:
    return float(Fraction(designation.split()[-1]))


def parse(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        if len(line) < 116 or not line[:10].strip().replace(".", "").isdigit():
            continue
        lam_nm, low, up = float(line[0:10]), line[22:44].strip(), line[45:69].strip()
        lo_cm, up_cm = float(line[70:80]), float(line[81:91])
        A, dA = float(line[105:111]), float(line[112:116])
        lam = lam_nm * 10.0
        g_up = 2.0 * _j(up) + 1.0
        rows.append({"paper": "ward2023_cr2", "species": "Cr II", "wavelength_air_A": round(lam, 4),
                     "ep_eV": round(lo_cm * CM1_TO_EV, 4), "j_lower": _j(low),
                     "ep_upper_eV": round(up_cm * CM1_TO_EV, 4), "A_1e6_s": A, "dA_1e6_s": dA,
                     "loggf": round(math.log10(1.49919e-16 * lam ** 2 * g_up * A * 1e6), 3),
                     "loggf_sigma_dex": round(math.log10(1.0 + dA / A), 4), "doi": DOI, "vizier": ""})
    return pd.DataFrame(rows, columns=COLS)


def main() -> None:
    t = parse(MRT.read_text())
    n = len(t)
    t = t[t.wavelength_air_A >= LO_A]
    t.to_csv(OUT, index=False)
    print({"table5": n, "at_or_redder_than_3000": len(t), "written": len(t)})


if __name__ == "__main__":
    main()
