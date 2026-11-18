#!/usr/bin/env python3
"""RYA-728 -- Cr II lab gf from Gurell et al. 2010 (A&A 511, A68; arXiv 1001.2139), parsed from
the paper's own online tables into data/reference/lab_gf_tables/transcribed/gurell2010_cr2.csv.

Table 3 prints vacuum wavelength, transition wavenumber, lower-level energy and log gf. Table 2
prints the air wavelength and each line's A-value uncertainty (%); the two join on air wavelength.
A line without a stated uncertainty cannot be graded and is dropped. A transition wavenumber
Table 3 prints against more than one lower level (27328.36 cm-1) names no single transition and
is dropped. Both are counted. Only 3000 A and redder is kept: nothing bluer is observed.

Table 3 prints lower levels 32603.40 and 32355.68 cm-1 against four lines (3617.27, 3633.99, 3665.97,
3684.22) whose upper levels would then be no known Cr II level. Table 2 calls their lower levels
3d5 2F, and VALD has each line at the same wavelength from 4.93-4.94 eV, reaching a known level. Those
rows take VALD's levels and are counted.

Scott+2015 Sect 6.4.1: this paper "measured no useful solar lines". It enters as lab-referenced
lines under the Reference Grade definition, never as a reason to expect a better Cr II value.
Usage: rya728_parse_gurell2010.py <gurell2010.pdf>
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import pandas as pd
import pypdf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reference" / "lab_gf_tables" / "transcribed" / "gurell2010_cr2.csv"
VALD = ROOT / "data" / "audit" / "rya1268_line_coverage" / "vald_solar_gf_refs.csv"
DOI = "10.1051/0004-6361/200913672"
CM_PER_EV = 8065.544
COLS = ["paper", "species", "wavelength_air_A", "ep_eV", "j_lower", "ep_upper_eV",
        "A_1e6_s", "dA_1e6_s", "loggf", "loggf_sigma_dex", "doi", "vizier"]
JOIN_TOL_A = 0.02
SAME_LINE_TOL_A = 0.01
LEVEL_TOL_EV = 0.002
VALD_TOL_A = 0.06
LO_A = 3000.0


def vac_to_air(vac: pd.Series) -> pd.Series:
    """Morton (2000), the standard IAU conversion."""
    s2 = (1e4 / vac) ** 2
    return vac / (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))


def read_tables(pdf: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    reader = pypdf.PdfReader(pdf)
    t2, t3 = [], []
    for page in reader.pages:
        text = (page.extract_text() or "").replace("−", "-")
        if "Online Material" not in text[:200]:
            continue
        if "Table 3" in text[:300] or "Line list" in text[:300]:
            t3 += re.findall(r"^(\d{4}\.\d{3}) (\d{5}\.\d{2}) (\d{4,5}\.\d{2}) (-?\d\.\d{3})", text, re.M)
        elif "Table 2" in text[:300]:
            t2 += re.findall(r"(\d{4}\.\d{3}) \d\.\d{3} \d\.\d{3} (\d+\.\d+) (\d+)\s*$", text, re.M)
    table2 = pd.DataFrame(t2, columns=["air", "A_1e7", "unc_pct"]).astype(float)
    table3 = pd.DataFrame(t3, columns=["vac", "sigma", "elow", "loggf"]).astype(float)
    return table2, table3


def build(table2: pd.DataFrame, table3: pd.DataFrame, vald: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    counts = {"table3": len(table3)}
    ambiguous = table3.duplicated("sigma", keep=False)
    counts["ambiguous_lower_level"] = int(ambiguous.sum())
    t3 = table3[~ambiguous].assign(air=lambda d: vac_to_air(d.vac)).sort_values("air")
    t3 = t3[t3.air >= LO_A]
    counts["at_or_redder_than_3000"] = len(t3)
    m = pd.merge_asof(t3, table2.sort_values("air").rename(columns={"air": "air2"}), left_on="air",
                      right_on="air2", direction="nearest", tolerance=JOIN_TOL_A)
    counts["no_stated_uncertainty"] = int(m.unc_pct.isna().sum())
    m = m[m.unc_pct.notna()]
    rows = []
    counts["lower_level_from_vald"] = 0
    for r in m.itertuples():
        ep = r.elow / CM_PER_EV
        ep_up = (r.elow + r.sigma) / CM_PER_EV
        same = vald[(vald.wavelength_air_A - r.air).abs() <= SAME_LINE_TOL_A]
        if not (vald.ep_upper_eV - ep_up).abs().le(LEVEL_TOL_EV).any() and len(same) == 1:
            ep, ep_up = float(same.ep_eV.iloc[0]), float(same.ep_upper_eV.iloc[0])
            counts["lower_level_from_vald"] += 1
        lv = vald[((vald.wavelength_air_A - r.air).abs() <= VALD_TOL_A) & ((vald.ep_eV - ep).abs() <= 0.02)]
        A, dA = r.A_1e7 * 10.0, r.A_1e7 * 10.0 * r.unc_pct / 100.0
        rows.append({"paper": "gurell2010_cr2", "species": "Cr II", "wavelength_air_A": round(r.air, 3),
                     "ep_eV": round(ep, 4), "j_lower": float(lv.j_lower.iloc[0]) if len(lv) else "",
                     "ep_upper_eV": round(ep_up, 4),
                     "A_1e6_s": round(A, 3), "dA_1e6_s": round(dA, 3), "loggf": r.loggf,
                     "loggf_sigma_dex": round(math.log10(1 + r.unc_pct / 100.0), 4),
                     "doi": DOI, "vizier": ""})
    counts["written"] = len(rows)
    return pd.DataFrame(rows, columns=COLS), counts


def main() -> None:
    table2, table3 = read_tables(Path(sys.argv[1]))
    vald = pd.read_csv(VALD)
    vald = vald[vald.species == "Cr II"]
    out, counts = build(table2, table3, vald)
    out.to_csv(OUT, index=False)
    print(counts)


if __name__ == "__main__":
    main()
