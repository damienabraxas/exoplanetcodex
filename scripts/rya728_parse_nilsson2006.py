#!/usr/bin/env python3
"""RYA-728 -- Cr II lab gf from Nilsson et al. 2006 (A&A 445, 1165), parsed from the paper's online
tables into data/reference/lab_gf_tables/transcribed/nilsson2006_cr2.csv.

Table 3 prints vacuum wavelength, wavenumber and this work's log gf first. Table 1 prints the air
wavelength and each line's A-value uncertainty (%). The tables print term labels rather than level
energies, so the lower and upper levels come from the VALD row at that wavelength. A line with no
stated uncertainty, or with no single VALD transition there, is dropped and counted. Only 3000 A
and redder is kept: nothing bluer is observed.

Scott+2015 Sect 6.4.1: the few good solar lines here "return abundances that are highly inconsistent
with each other". AGSS21 moved to Lawler+2017, which outranks this table wherever both have a line.
Usage: rya728_parse_nilsson2006.py <aa4152-05.pdf>
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import pandas as pd
import pypdf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rya728_parse_gurell2010 import COLS, JOIN_TOL_A, LO_A, VALD, vac_to_air  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reference" / "lab_gf_tables" / "transcribed" / "nilsson2006_cr2.csv"
DOI = "10.1051/0004-6361:20054152"
VALD_TOL_A = 0.03


def read_tables(pdf: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    t1, t3 = [], []
    for page in pypdf.PdfReader(pdf).pages:
        text = (page.extract_text() or "").replace("−", "-")
        if "Online Material" not in text[:200]:
            continue
        t1 += re.findall(r"\d{2} \d{3}\.\d{2} (\d{4}\.\d{3}) \d\.\d{2} \d\.\d{2} (\d+\.\d+) (\d+)\s*$", text, re.M)
        t3 += re.findall(r"^(\d{4}\.\d{3}) (\d{2} \d{3}\.\d{2}) ([-+]\d\.\d{3})", text, re.M)
    table1 = pd.DataFrame(t1, columns=["air", "A_1e7", "unc_pct"]).astype(float)
    table3 = pd.DataFrame(t3, columns=["vac", "sigma", "loggf"])
    table3["sigma"] = table3.sigma.str.replace(" ", "")
    return table1, table3.astype(float)


def build(table1: pd.DataFrame, table3: pd.DataFrame, vald: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    counts = {"table3": len(table3)}
    t3 = table3.assign(air=lambda d: vac_to_air(d.vac)).sort_values("air")
    t3 = t3[t3.air >= LO_A]
    counts["at_or_redder_than_3000"] = len(t3)
    m = pd.merge_asof(t3, table1.sort_values("air").rename(columns={"air": "air2"}), left_on="air",
                      right_on="air2", direction="nearest", tolerance=JOIN_TOL_A)
    counts["no_stated_uncertainty"] = int(m.unc_pct.isna().sum())
    m = m[m.unc_pct.notna()]
    rows, no_transition = [], 0
    for r in m.itertuples():
        lv = vald[(vald.wavelength_air_A - r.air).abs() <= VALD_TOL_A]
        if len(lv) != 1:
            no_transition += 1
            continue
        A, dA = r.A_1e7 * 10.0, r.A_1e7 * 10.0 * r.unc_pct / 100.0
        rows.append({"paper": "nilsson2006_cr2", "species": "Cr II", "wavelength_air_A": round(r.air, 3),
                     "ep_eV": float(lv.ep_eV.iloc[0]), "j_lower": float(lv.j_lower.iloc[0]),
                     "ep_upper_eV": float(lv.ep_upper_eV.iloc[0]), "A_1e6_s": round(A, 3),
                     "dA_1e6_s": round(dA, 3), "loggf": r.loggf,
                     "loggf_sigma_dex": round(math.log10(1 + r.unc_pct / 100.0), 4),
                     "doi": DOI, "vizier": ""})
    counts["no_single_vald_transition"] = no_transition
    counts["written"] = len(rows)
    return pd.DataFrame(rows, columns=COLS), counts


def main() -> None:
    table1, table3 = read_tables(Path(sys.argv[1]))
    vald = pd.read_csv(VALD)
    out, counts = build(table1, table3, vald[vald.species == "Cr II"])
    out.to_csv(OUT, index=False)
    print(counts, "table1", len(table1))


if __name__ == "__main__":
    main()
