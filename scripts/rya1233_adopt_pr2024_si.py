#!/usr/bin/env python3
"""RYA-1233 -- grade the Si I lines beyond Garz 1973's reach (>= 8000 A) with Pehlivan Rhodin
et al. 2024 (A&A 682, A184; CDS J/A+A/682/A184), the gf source Lodders, Bergemann & Palme
2025 call "preferred in recent studies (Magg et al. 2022, Deshmukh et al. 2022)" and the one
Deshmukh et al. 2022 use for every IR line of their solar Si set.

    python3 scripts/rya1233_adopt_pr2024_si.py [--check]

Order of preference per canonical Si I row at >= 8000 A (Garz measured 2500-8000 A only):
  1. Table B3 -- EXPERIMENTAL log gf (FTS branching fractions x lifetimes), its own dex error.
     gf_tier LAB, lab_source_tag PR2024_exp.
  2. Table B4 -- calculated gf RESCALED to the experimental lifetimes (`gfresc`). B4's own
     `funcert` u is the calculation's INTERNAL spread (0.2% on 10627.648, which sits 0.64 dex
     from VALD) -- not an accuracy. The accuracy is MEASURED from the paper itself: Table B3
     prints experimental and calculated log gf for the same 17 lines, rms(exp - calc) = 0.046
     dex (derived at run time, `calc_floor_dex`). sigma_dex = sqrt(log10(1+u)^2 + floor^2).
     The floor comes from strong low-lying lines, so for weak lines it is a LOWER bound.
     gf_tier PR2024-CALC: priced, not LAB-graded.
Rows already LAB (a lab value we adopted on purpose) are left alone. Optical rows (< 8000 A)
are left on Asplund's Garz values -- Asplund 2021 keeps Garz there and this does not overrule
that. WAVELENGTHS, MEASURED not read from the ReadMe (which says vacuum for 2000-20000 A
for both): Table B3 is AIR -- 14 of 17 rows sit within 0.003 A of a canonical air line as
printed (1210.354 nm = Bergemann's 12103.54 A) and 0 after a vac->air shift; Table B4 is
VACUUM -- 3 matches as printed, 417 after conversion.
A row is matched on air wavelength, unique within 0.03 A with no second candidate within
0.10 A; anything else is reported, never guessed. Every replaced value is written to
data/audit/rya1233_si_pr2024/adoption.csv.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/linelists/primary_gf/pehlivan_rhodin2024"
GF = ROOT / "data/linelists/canonical_gf.csv"
AUDIT = ROOT / "data/audit/rya1233_si_pr2024/adoption.csv"
DOI = "10.1051/0004-6361/202245686"
REF = "Pehlivan Rhodin et al. 2024, A&A 682, A184"
GARZ_LIMIT_A = 8000.0
TOL_A, CLEAR_A = 0.03, 0.10
STATUS = "adjudicated_rya1233_pr2024"


def vac_to_air(w):
    s2 = (1e4 / np.asarray(w, float)) ** 2
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
    return np.asarray(w, float) / n


def _vac_A(nm):
    w = np.asarray(nm, float) * 10.0
    return np.where((w > 2000) & (w < 20000), vac_to_air(w), w)


def read_b3() -> pd.DataFrame:
    rows = []
    for l in (SRC / "tableb3.dat").read_text().splitlines():
        if not l.strip():
            continue
        exp, eexp = l[76:82].strip(), l[83:87].strip()
        if not exp:
            continue
        rows.append({"lam_nm": float(l[39:47]), "loggf": float(exp), "sigma": float(eexp),
                     "upper": l[0:17].strip(), "lower": l[18:38].strip()})
    d = pd.DataFrame(rows)
    d["air_A"] = d.lam_nm * 10.0          # B3 is printed in AIR (measured; see module doc)
    return d


def calc_floor_dex() -> float:
    """rms(experimental - calculated) log gf over Table B3 -- the calculation's measured accuracy."""
    d = [float(l[76:82]) - float(l[88:94]) for l in (SRC / "tableb3.dat").read_text().splitlines()
         if l.strip() and l[76:82].strip() and l[88:94].strip()]
    return float(np.sqrt(np.mean(np.square(d))))


def read_b4() -> pd.DataFrame:
    floor = calc_floor_dex()
    rows = []
    for l in (SRC / "tableb4.dat").read_text().splitlines():
        if not l.strip():
            continue
        gfr, u = l[148:157].strip(), l[158:163].strip()
        if not gfr or not u or float(gfr) <= 0:
            continue
        rows.append({"lam_nm": float(l[104:113]), "loggf": math.log10(float(gfr)),
                     "sigma": math.hypot(math.log10(1 + float(u)), floor), "u": float(u),
                     "upper": l[0:46].strip(), "lower": l[47:93].strip()})
    d = pd.DataFrame(rows)
    d["air_A"] = _vac_A(d.lam_nm)
    return d


def match(cg_si: pd.DataFrame, w: float):
    dw = (cg_si.wavelength_air_A - w).abs()
    near = dw[dw <= CLEAR_A].sort_values()
    if len(near) == 0 or near.iloc[0] > TOL_A:
        return None, "no canonical row within 0.03 A"
    if len(near) > 1:
        return None, f"{len(near)} canonical rows within 0.10 A -- ambiguous"
    return near.index[0], ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="report only, write nothing")
    a = ap.parse_args(argv)
    cg = pd.read_csv(GF, low_memory=False)
    si = cg[(cg.species == "Si I") & (cg.wavelength_air_A >= GARZ_LIMIT_A)]
    audit, done = [], set()
    for kind, tab in (("exp", read_b3()), ("calc", read_b4())):
        for _, r in tab[tab.air_A >= GARZ_LIMIT_A].iterrows():
            i, why = match(si, float(r.air_A))
            rec = {"source_table": "B3" if kind == "exp" else "B4", "air_A": round(r.air_A, 3),
                   "upper": r.upper, "lower": r.lower, "pr_loggf": round(r.loggf, 3),
                   "pr_sigma_dex": round(r.sigma, 4)}
            if i is None:
                audit.append({**rec, "action": "unmatched", "why": why})
                continue
            if i in done:
                audit.append({**rec, "action": "skipped", "why": "already set from Table B3"})
                continue
            old = cg.loc[i]
            if str(old.gf_tier) == "LAB" and str(old.lab_source_tag) != "PR2024_exp":
                audit.append({**rec, "action": "kept", "why": f"LAB ({old.lab_source_tag})"})
                continue
            audit.append({**rec, "action": "adopted", "why": "",
                          "canonical_A": old.wavelength_air_A, "ep_eV": old.excitation_potential_eV,
                          "old_loggf": old.log_gf, "old_tier": old.gf_tier,
                          "old_reference": old.loggf_reference,
                          "d_loggf": round(r.loggf - float(old.log_gf), 3)})
            done.add(i)
            if a.check:
                continue
            cg.at[i, "log_gf"] = round(float(r.loggf), 3)
            cg.at[i, "gf_sigma_dex"] = round(float(r.sigma), 4)
            cg.at[i, "loggf_reference"] = REF + (" Table B3 (experimental)" if kind == "exp"
                                                 else " Table B4 (calculated, rescaled)")
            cg.at[i, "gf_source_doi"] = DOI
            cg.at[i, "gf_tier"] = "LAB" if kind == "exp" else "PR2024-CALC"
            cg.at[i, "lab_source_tag"] = "PR2024_exp" if kind == "exp" else "PR2024_calc"
            # A fixed token (RYA-834's convention); the value it replaced is in the audit.
            cg.at[i, "adjudication_status"] = STATUS
    au = pd.DataFrame(audit)
    print(au.groupby(["source_table", "action"]).size().to_string())
    if not a.check:
        cg.to_csv(GF, index=False)
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        au.to_csv(AUDIT, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
