#!/usr/bin/env python3
"""RYA-1268 -- per-line lab gf AND uncertainty from the lab papers' machine-readable tables.

The Asplund-cited solar tables print no gf error; the lab papers they cite do, as a transition
probability A +/- dA. sigma_dex = log10(1 + dA/A), the bridge pipeline.gf_grades.nist_sigma_dex
uses for a NIST percent class. The tables are the papers' own machine-readable versions, mirrored
free by CDS VizieR (vendored under data/reference/lab_gf_tables/vizier/ with their ReadMe; the
arXiv copies carry only a 5-row stub). Not on VizieR: Lawler+2019 Sc, Lawler+2013 Ti I,
Lawler+2001 Eu II -- their gf reach us through VALD (codes LHSNWC / LGWSC / LWHS).
Output: data/reference/lab_gf_tables/lab_gf_lines.csv
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "reference" / "lab_gf_tables"
CM_PER_EV = 8065.544
# catalog dir, table, bib key, doi, species (None: Mn I/II by J), wavelength col, nm?, Elow, Jlow, A, eA, gf
# (the upper level is the same table's E1/UpLev/Lev1/EU column, see UPPER)
CATS = [
    ("J_ApJS_208_27", "table4.dat", "wood2013_ti2", "10.1088/0067-0049/208/2/27", "Ti II", "lamAir", False, "E0", "J0", "A", "e_A", "log(gf)"),
    ("J_ApJS_215_20", "table3.dat", "lawler2014_v1", "10.1088/0067-0049/215/2/20", "V I", "lamAir", False, "E0", "J0", "A", "e_A", "log(gf)"),
    ("J_ApJS_214_18", "table3.dat", "wood2014_v2", "10.1088/0067-0049/214/2/18", "V II", "lambda", False, "E0", "J0", "A", "e_A", "log.gf"),
    ("J_ApJS_228_10", "table4.dat", "lawler2017_cr2", "10.3847/1538-4365/228/1/10", "Cr II", "lamAir", True, "LowLev", "LowJ", "A", "e_A", "log(gf)"),
    ("J_ApJS_220_13", "table3.dat", "lawler2015", "10.1088/0067-0049/220/1/13", "Co I", "lamAir", False, "Lev0", "J0", "A", "e_A", "log(gf)"),
    ("J_ApJS_211_20", "table3.dat", "wood2014", "10.1088/0067-0049/211/2/20", "Ni I", "lambda", False, "E0", "J0", "A", "e_A", "log(gf)"),
    # Den Hartog+2011 Table 3: integral J = Mn II, half-integral J = Mn I (the table caption)
    ("J_ApJS_194_35", "table3.dat", "denhartog2011", "10.1088/0067-0049/194/2/35", None, "lamAir", False, "E0", "J0", "TranP", "e_TranP", "log(gf)"),
    ("J_ApJ_667_1267", "table3.dat", "sobeck2007", "10.1086/519987", "Cr I", "lambda", False, "EL", "JL", "Aexp", "e_Aexp", "log(gf)"),
]


UPPER = {"J_ApJS_228_10": "UpLev", "J_ApJS_220_13": "Lev1", "J_ApJ_667_1267": "EU"}


def main() -> None:
    from astropy.table import Table
    rows = []
    for d, f, bib, doi, sp, wc, nm, ec, jc, ac, eac, gc in CATS:
        t = Table.read(BASE / "vizier" / d / f, format="ascii.cds", readme=str(BASE / "vizier" / d / "ReadMe"))
        for r in t:
            a, da = float(r[ac]), float(r[eac])
            jl = float(r[jc])
            rows.append({"paper": bib, "species": sp or ("Mn I" if jl % 1 else "Mn II"),
                         "wavelength_air_A": round(float(r[wc]) * (10 if nm else 1), 4),
                         "ep_eV": round(float(r[ec]) / CM_PER_EV, 4), "j_lower": jl,
                         "ep_upper_eV": round(float(r[UPPER.get(d, "E1")]) / CM_PER_EV, 4),
                         "A_1e6_s": a, "dA_1e6_s": da, "loggf": float(r[gc]),
                         "loggf_sigma_dex": round(math.log10(1 + da / a), 4),
                         "doi": doi, "vizier": d.replace("_", "/")})
    with (BASE / "lab_gf_lines.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(Counter((r["paper"], r["species"]) for r in rows))


if __name__ == "__main__":
    main()
