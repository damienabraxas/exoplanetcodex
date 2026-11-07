#!/usr/bin/env python3
"""RYA-1311 -- grade the H- and K-band Ca I lines with Civis et al. 2013 (A&A 554, A24; Table 5),
quantum-defect-theory f-values for Ca I from 800 to 9000 cm-1.

    python3 scripts/rya1311_adopt_civis2013_ca.py [--check]

Why this source: no laboratory gf exists for any Ca line in H or K, and every non-Kurucz gf
there (APOGEE, Shetrone 2015 / Smith 2021) is fitted to the solar spectrum, so circular for a
solar measurement. Civis 2013 is a calculation, not a solar fit, and it is TESTED: on the 28
H/K/Y/J Mg I features where we hold a lab, NIST or published-set gf, its Mg values (same paper,
Table 6) differ from ours by a median -0.025 dex, MAD 0.04. Ryan 2026-11-07: include it.

Per-line sigma (MEASURED, not assumed): the RMS of Civis minus Elgueta et al. 2024 on the two
Ca I lines both cover (11955.955, 12105.841), after Elgueta's measured 0.141 dex scale offset
(gf_error_model E2024_astro) = 0.105 -> SIGMA_DEX 0.11. The paper itself states no f-value
uncertainty; its authors note Ca agrees with NIST less well than Mg.

A Civis feature is the set of its components sharing a lower level within 1.5 A. A canonical
Ca I row in H/K (14000-25000 A) with a KURUCZ gf adopts log10(sum gf) of the components within
0.6 A of it, but only when no other canonical row of the same lower level shares them (a
multiplet we split differently is reported, never guessed). Graded rows are never touched.
Every replaced value is written to data/audit/rya1311_civis2013/adoption.csv.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from rya1311_adopt_elgueta2024_ca import _splice  # noqa: E402

SRC = ROOT / "data/audit/rya1311_civis2013/civis2013_table5_CaI.tsv"
GF = ROOT / "data/linelists/canonical_gf.csv"
AUDIT = ROOT / "data/audit/rya1311_civis2013/adoption.csv"
DOI = "10.1051/0004-6361/201321052"
REF = "Civis et al. 2013, A&A 554, A24, Table 5 (quantum defect theory)"
TIER, TAG = "C2013-QDT", "C2013_qdt"
SIGMA_DEX = 0.11
LO_A, HI_A, TOL_A, ELO_TOL_EV = 14000.0, 25000.0, 0.6, 0.02
STATUS = "adjudicated_rya1311_civis2013"


def read_table5() -> pd.DataFrame:
    d = pd.read_csv(SRC, sep="\t", comment="#", dtype=str)
    d = d[d.recno.str.strip().str.isdigit()].copy()
    for c in ("Elow", "lambda", "fik"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["Jlo"] = d.Trans.map(lambda s: int(re.findall(r"_(\d+)_", s.split("--")[0])[-1]))
    d["gf"] = (2 * d.Jlo + 1) * d.fik
    d["air_A"] = d["lambda"] * 10.0
    d["elo_eV"] = d.Elow / 8065.544
    return d[(d.air_A >= LO_A) & (d.air_A <= HI_A)]


def _multiplet(civ: pd.DataFrame, comp: pd.DataFrame) -> pd.DataFrame:
    """Every Civis component chained to `comp` within 1.5 A on the same lower level."""
    same = civ[(civ.elo_eV - comp.elo_eV.iloc[0]).abs() <= ELO_TOL_EV].sort_values("air_A")
    lo, hi = comp.air_A.min(), comp.air_A.max()
    grew = True
    while grew:
        sel = same[(same.air_A >= lo - 1.5) & (same.air_A <= hi + 1.5)]
        grew = sel.air_A.min() < lo or sel.air_A.max() > hi
        lo, hi = sel.air_A.min(), sel.air_A.max()
    return sel


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="report only, write nothing")
    a = ap.parse_args(argv)
    cg = pd.read_csv(GF, low_memory=False)
    ca = cg[(cg.species == "Ca I") & (cg.wavelength_air_A >= LO_A) & (cg.wavelength_air_A <= HI_A)]
    civ = read_table5()
    audit, edits = [], {}
    for i, row in ca.iterrows():
        w, e = float(row.wavelength_air_A), float(row.excitation_potential_eV)
        comp = civ[((civ.air_A - w).abs() <= TOL_A) & ((civ.elo_eV - e).abs() <= ELO_TOL_EV)]
        if comp.empty:
            continue
        rec = {"canonical_A": w, "elo_eV": round(e, 4), "n_components": len(comp),
               "civis_trans": comp.Trans.iloc[0], "old_loggf": row.log_gf, "old_tier": row.gf_tier}
        rivals = ca[(ca.index != i) & ((ca.excitation_potential_eV - e).abs() <= ELO_TOL_EV)
                    & ca.wavelength_air_A.map(lambda x: ((comp.air_A - x).abs() <= TOL_A).any())]
        if len(rivals):
            audit.append({**rec, "action": "ambiguous",
                          "why": f"{len(rivals)} other canonical row(s) share these components"})
            continue
        # The whole Civis multiplet (same lower level, components chained within 1.5 A) must be
        # carried by our rows: 20494-20497 A has 6 components, our 3 rows catch 3 of them, and
        # adopting a lone weak component read 0.9-1.1 dex off.
        left = _multiplet(civ, comp)
        mates = ca[((ca.excitation_potential_eV - e).abs() <= ELO_TOL_EV)]
        caught = left[left.air_A.map(lambda x: ((mates.wavelength_air_A - x).abs() <= TOL_A).any())]
        if len(caught) < len(left):
            audit.append({**rec, "action": "ambiguous",
                          "why": f"our rows carry {len(caught)} of the multiplet's {len(left)} components"})
            continue
        if str(row.gf_tier) != "KURUCZ":
            audit.append({**rec, "action": "kept", "why": f"already {row.gf_tier}"})
            continue
        new = float(np.log10(comp.gf.sum()))
        audit.append({**rec, "action": "adopted", "why": "", "civis_loggf": round(new, 3),
                      "d_loggf": round(new - float(row.log_gf), 3)})
        edits[row.line_id] = {"log_gf": f"{new:.3f}", "gf_sigma_dex": f"{SIGMA_DEX:.3f}",
                              "loggf_reference": f"{REF} | prior: {row.loggf_reference}",
                              "gf_source_doi": DOI, "gf_tier": TIER, "lab_source_tag": TAG,
                              "adjudication_status": STATUS}
    au = pd.DataFrame(audit)
    print(au.to_string(index=False))
    print(f"adopted {len(edits)}; ambiguous {(au.action == 'ambiguous').sum()}")
    if not a.check:
        _splice(edits)
        au.to_csv(AUDIT, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
