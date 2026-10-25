#!/usr/bin/env python3
"""RYA-1268 -- lab gf + uncertainty from lab papers that have NO VizieR table, parsed from the
papers' own free PDFs (arXiv / author copies) into data/reference/lab_gf_tables/transcribed/.

  bw2006_ti1      Blackwell-Whitehead et al. 2006, MNRAS 373, 1603 (arXiv astro-ph/0612656), Table 3
                  -- Ti I 0.5-3.9 um, log gf + Unc. (dex) as printed; VALD code BLNP
  bw2011_mn1      Blackwell-Whitehead et al. 2011, A&A 525, A44, Table 1 -- Mn I 1.29-1.40 um,
                  log gf + Unc. (dex) as printed; VALD code B-WPNP
  denhartog2021_ca1  Den Hartog et al. 2021, ApJS 255, 27 (arXiv 2105.10534), Table 8 MRT --
                  ONLY rows with an experimental reference (ExptRef); the rest are the M17 theory
                  values and are not lab. sigma from ExptUnc (%)
  lawler2018_co2  Lawler et al. 2018, ApJS 238, 7 (arXiv 1806.00581), Table 3 -- Co II A +/- dA

The two Blackwell-Whitehead tables print no level energies: EP, upper EP and J_lower come from
the VALD row at the same wavelength (+/-0.06 A), the row carrying the paper's own code first.
Unmatched rows (outside our VALD extracts), and rows whose upper level no other row of the
table shares (a different transition at the same wavelength), are dropped and counted.
Usage: rya1268_parse_lab_pdfs.py <dir holding the four PDFs>
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import pandas as pd
import pypdf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reference" / "lab_gf_tables" / "transcribed"
VALD = ROOT / "data" / "audit" / "rya1268_line_coverage" / "vald_solar_gf_refs.csv"
CM_PER_EV = 8065.544
COLS = ["paper", "species", "wavelength_air_A", "ep_eV", "j_lower", "ep_upper_eV",
        "A_1e6_s", "dA_1e6_s", "loggf", "loggf_sigma_dex", "doi", "vizier"]


def _text(pdf: Path, pages: range) -> str:
    r = pypdf.PdfReader(pdf)
    return "\n".join(r.pages[i].extract_text() or "" for i in pages).replace("−", "-")


def _air(vac_A: pd.Series) -> pd.Series:
    s2 = (1e4 / vac_A) ** 2  # Morton 2000 / IAU vacuum -> air
    return vac_A / (1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2))


def _levels_from_vald(d: pd.DataFrame, species: str, code: str) -> pd.DataFrame:
    v = pd.read_csv(VALD)
    v = v[v.species == species]
    out = []
    for r in d.itertuples():
        c = v[(v.wavelength_air_A - r.wavelength_air_A).abs() <= 0.06]
        # the paper's own code first; the levels are physics, so any code's row serves after it
        c = c[c.gf_ref_code == code] if (c.gf_ref_code == code).any() else c
        if c.empty:
            out.append((math.nan, math.nan, math.nan))
            continue
        c = c.loc[(c.loggf - r.loggf).abs().idxmin()]  # the component closest in gf
        out.append((c.ep_eV, c.j_lower, c.ep_upper_eV))
    d[["ep_eV", "j_lower", "ep_upper_eV"]] = out
    # Each table block is one upper level with several branches; a row whose upper level no other
    # row shares matched a different transition at the same wavelength (Ti I 24713: 5.70 eV, not 2.66).
    up = d.ep_upper_eV.round(3)
    d.loc[up.map(up.value_counts()) < 2, "ep_eV"] = math.nan
    n = d.ep_eV.isna().sum()
    print(f"  {species} {code}: {len(d) - n}/{len(d)} matched to VALD levels")
    return d.dropna(subset=["ep_eV"])


def bw2006(pdf: Path) -> pd.DataFrame:
    pat = re.compile(r"([a-z] ?\d[A-Z]\d)\s+(\d{3,4}\.\d{4})\s+(\d+\.\d{3})\s+(\d\.\d{3})\s+(\d+)\s+"
                     r"([\d.]+)\s+(-?\d\.\d\d)\s+(0\.\d\d)")
    d = pd.DataFrame([m.groups() for m in pat.finditer(_text(pdf, range(4, 7)))],
                     columns=["lower", "vac_nm", "sigma", "BF", "BFunc", "A_1e7", "loggf", "unc"])
    d = d.assign(**{c: pd.to_numeric(d[c]) for c in d.columns[1:]})
    d = d.assign(wavelength_air_A=_air(d.vac_nm * 10).round(3), A_1e6_s=(d.A_1e7 * 10).round(4),
                 dA_1e6_s=math.nan, loggf_sigma_dex=d.unc)
    return _levels_from_vald(d, "Ti I", "BLNP").assign(
        paper="bw2006_ti1", species="Ti I", doi="10.1111/j.1365-2966.2006.11161.x", vizier="")


def bw2011(pdf: Path) -> pd.DataFrame:
    pat = re.compile(r"([a-z]\s?\d[A-Z]\d+/2)\s+(\d{1,2}\s?\d{3}\.\d\d)\s+(\d{1,2}\s?\d{3}\.\d\d)\s+"
                     r"(\d\.\d{3})\s+(\d+)\s+([\d.]+)\s+\(([\d.]+)\)\s+(-?\d\.\d\d)\s+(0\.\d\d)")
    d = pd.DataFrame([m.groups() for m in pat.finditer(_text(pdf, range(2, 5)))],
                     columns=["lower", "sigma", "wl", "BF", "BFunc", "A_1e7", "eA_1e7", "loggf", "unc"])
    d = d.assign(**{c: d[c].str.replace(" ", "").astype(float) for c in d.columns[1:]})
    d = d.assign(wavelength_air_A=d.wl, A_1e6_s=(d.A_1e7 * 10).round(4), dA_1e6_s=(d.eA_1e7 * 10).round(4),
                 loggf_sigma_dex=d.unc)
    return _levels_from_vald(d, "Mn I", "B-WPNP").assign(
        paper="bw2011_mn1", species="Mn I", doi="10.1051/0004-6361/201015006", vizier="")


def denhartog2021(pdf: Path) -> pd.DataFrame:
    pat = re.compile(r"^\s*(\d{4,5}\.\d{4})\s+(\d+\.\d{3})\s+(\d+\.\d{3})\s+(\d\.\d{3})\s*(-?\d+\.\d{3})"
                     r"\s+(-?\d+)?\s*(-?\d+\.\d+)\s*(\S+)?\s*(\d+(?:\.\d+)?)?\s*$", re.M)
    d = pd.DataFrame([m.groups() for m in pat.finditer(_text(pdf, range(53, 58)))],
                     columns=["wl", "Eup", "Elow", "ep", "gf_M17", "LV", "loggf", "ref", "unc_pct"])
    assert len(d) == 202, len(d)  # the caption: "202 lines of Ca I"
    d = d[d.ref.notna()].copy()  # no ExptRef = the M17 theory value
    d = d.assign(**{c: pd.to_numeric(d[c]) for c in ["wl", "Eup", "Elow", "loggf", "unc_pct"]})
    return d.assign(paper="denhartog2021_ca1:" + d.ref, species="Ca I", wavelength_air_A=d.wl,
                    ep_eV=(d.Elow / CM_PER_EV).round(4), j_lower=math.nan,
                    ep_upper_eV=(d.Eup / CM_PER_EV).round(4), A_1e6_s=math.nan, dA_1e6_s=math.nan,
                    loggf_sigma_dex=(1 + d.unc_pct / 100).map(math.log10).round(4),
                    doi="10.3847/1538-4365/ac04b1", vizier="")


def lawler2018(pdf: Path) -> pd.DataFrame:
    pat = re.compile(r"^\s*(\d{4,5}\.\d{3,4})\s+(\d+\.\d{3})\s+(\d)\s+(\d+\.\d{3})\s+(\d)\s+([\d.]+)\s*"
                     r"±\s*([\d.]+)\s+(-?\d+\.\d+)", re.M)
    r = pypdf.PdfReader(pdf)
    d = pd.DataFrame([m.groups() for m in pat.finditer(_text(pdf, range(45, len(r.pages))))],
                     columns=["wl", "Eup", "Jup", "Elow", "Jlow", "A", "dA", "loggf"]).apply(pd.to_numeric)
    assert len(d) == 84, len(d)  # the caption: "84 lines of Co II"
    return d.assign(paper="lawler2018_co2", species="Co II", wavelength_air_A=d.wl,
                    ep_eV=(d.Elow / CM_PER_EV).round(4), j_lower=d.Jlow,
                    ep_upper_eV=(d.Eup / CM_PER_EV).round(4), A_1e6_s=d.A, dA_1e6_s=d.dA,
                    loggf_sigma_dex=(1 + d.dA / d.A).map(math.log10).round(4),
                    doi="10.3847/1538-4365/aac773", vizier="")


PAPERS = {"bw2006_ti1": ("BlackwellWhitehead_2006_MNRAS_373_1603_TiI_IR.pdf", bw2006),
          "bw2011_mn1": ("BlackwellWhitehead_2011_AA_525_A44_MnI_IR.pdf", bw2011),
          "denhartog2021_ca1": ("DenHartog_2021_ApJS_255_27_CaI.pdf", denhartog2021),
          "lawler2018_co2": ("Lawler_2018_CoII_1806.00581.pdf", lawler2018)}


def main() -> None:
    src = Path(sys.argv[1])
    OUT.mkdir(parents=True, exist_ok=True)
    for key, (fn, parse) in PAPERS.items():
        d = parse(src / fn)[COLS]
        d.to_csv(OUT / f"{key}.csv", index=False)
        print(f"{key}: {len(d)} rows, {d.wavelength_air_A.min():.0f}-{d.wavelength_air_A.max():.0f} A")


if __name__ == "__main__":
    main()
