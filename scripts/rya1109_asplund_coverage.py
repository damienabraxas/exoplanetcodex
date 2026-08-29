#!/usr/bin/env python3
"""RYA-1109 - coverage of AGSS21's Fe line set against our VIS holdings and our graded pool.

    python3 scripts/rya1109_asplund_coverage.py

Answers the two questions RYA-1109 step 4 asks:

  1. How many AGSS21 lines fall inside each VIS holding's wavelength span?
  2. What is the ACTUAL overlap with our own graded (lab-gf) pool? Ryan recalled "~1".

🔴 THE "~1" IS ABOUT A DIFFERENT LIST. The 1-of-50 figure comes from RYA-1104/1106 and is
AP2002's set (`amarsi2022_solar_control_lines.csv`) against our graded pool. AGSS21's own
Table A.2 is a different list, and it is built on Den Hartog (2014) and Belmonte (2017) -
laboratory sources that OUR graded tier also draws on - so there is no reason to expect the
same answer. Measured here rather than assumed.

MATCHING IS THE CANONICAL MATCHER, IN STRICT MODE. `pipeline.line_match.match(...,
require_ep=True)` on the lambda+EP dual key (RYA-1037). Not a rounded wavelength join: a
2-dp wavelength is not a line identity (RYA-1033), and one in-window candidate is not the
same as the right transition (RYA-853). A wavelength that cannot be resolved WITH an EP is
reported unresolved rather than counted.

🔴 THE TOLERANCE IS DERIVED FROM THE SOURCE'S PRINTED PRECISION, AND THE FIRST RUN OF THIS
SCRIPT GOT IT WRONG. `line_match.MATCH_TOL_A` is 0.005 A, which suits a table printed to
0.01 A. AGSS21 Table A.2 prints lambda in NANOMETRES to two decimals -- 0.1 A resolution,
half-width 0.05 A -- so the default window is 20x too tight FOR THAT SOURCE and throws away
real matches. Measured: the AGSS21 Fe I overlap with our graded pool reads 2/40 at 0.005 A
and 19/40 at 0.05 A. The first number was an artifact of the tolerance, not a fact about
the line lists.

So each reference list declares the precision it was PRINTED at, and the tolerance follows
from it. The counts are not tuned by it either -- `--plateau` sweeps the tolerance and shows
the count STOPS MOVING well before the window is wide enough to sweep in a neighbouring
line (AGSS21 Fe I: 19 at 0.05 A and still 19 at 0.50 A, zero ambiguous). A count that keeps
climbing with the window is a count of coincidences; one that plateaus is the answer.

AP2002 prints lambda in ANGSTROMS to two decimals, so the module default is already right
for it, and its 2/50 is stable from 0.005 A to 0.1 A. The two lists genuinely differ.
"""
from __future__ import annotations

import csv
from collections import Counter
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import line_match  # noqa: E402

ASPLUND = ROOT / "data" / "reference" / "asplund2021_fe" / "asplund2021_fe_lines.csv"
AP2002 = (ROOT / "data" / "reference" / "amarsi2022_training"
          / "amarsi2022_solar_control_lines.csv")
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"
OUT = ROOT / "data" / "results" / "rya1109" / "asplund_coverage.json"

#: The VIS solar holdings the Asplund replication targets (RYA-1106's scope).
VIS_HOLDINGS = ["solar_kpno_kurucz2005_corrected", "solar_kpno_molecfit_corrected",
                "solar_harps_molecfit_corrected", "solar_iag"]


#: Half the printed wavelength resolution of each reference table. DERIVED, not chosen:
#: AGSS21 Table A.2 prints nm to 2 dp (0.1 A), AP2002 Table 2 prints A to 2 dp (0.01 A).
TOL_A = {"asplund_agss21": 0.05, "ap2002": 0.005}


def fe_pool(ion: str, *, lab_only: bool) -> pd.DataFrame:
    d = pd.read_csv(CANON, low_memory=False)
    d = d[d["species"].astype(str) == f"Fe {ion}"]
    if lab_only:
        d = d[d["gf_tier"].astype(str) == "LAB"]
    return d.dropna(subset=["wavelength_air_A", "excitation_potential_eV"]).reset_index(
        drop=True)


def _match(ref: pd.DataFrame, pool: pd.DataFrame, tol: float):
    return line_match.match(ref["wavelength_air_A"].to_numpy(float),
                            pool["wavelength_air_A"].to_numpy(float),
                            want_ep=ref["elo_eV"].to_numpy(float),
                            src_ep=pool["excitation_potential_eV"].to_numpy(float),
                            require_ep=True, tol_A=tol)


def overlap(ref: pd.DataFrame, pool: pd.DataFrame, tol: float) -> dict:
    """Strict lambda+EP match of a reference list against a pool, at a DERIVED tolerance."""
    if ref.empty or pool.empty:
        return {"n_ref": int(len(ref)), "n_pool": int(len(pool)), "matched": 0,
                "matched_lambda": [], "unresolved": int(len(ref)), "tol_A": tol}
    r = _match(ref, pool, tol)
    idx = np.asarray(r.index)
    hit = idx >= 0
    return {
        "n_ref": int(len(ref)),
        "n_pool": int(len(pool)),
        "matched": int(hit.sum()),
        "tol_A": tol,
        "matched_lambda": [round(float(x), 3)
                           for x in ref["wavelength_air_A"].to_numpy(float)[hit]],
        "unresolved": int((~hit).sum()),
    }


def plateau(ref: pd.DataFrame, pool: pd.DataFrame,
            tols=(0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.25, 0.50)) -> dict:
    """The count as a function of the window. A real overlap plateaus; coincidences climb."""
    return {str(t): int((np.asarray(_match(ref, pool, t).index) >= 0).sum()) for t in tols}


def main() -> int:
    asp = pd.read_csv(ASPLUND)
    doc: dict = {"ticket": "RYA-1109", "source": asp["source"].iloc[0], "holdings": {},
                 "graded_overlap": {}, "ap2002_comparison": {}}

    # ── 1. per-holding wavelength coverage ────────────────────────────────────
    # ⚠️ THE SPAN IS THE INSTRUMENT'S, NOT THE HOLDING'S. holdings_manifest_registry maps
    # holding -> instrument and carries NO wavelength range, so the range comes from
    # instrument_catalog. A holding can cover less than its instrument's full range (a
    # single setting, a trimmed product), so these counts are an UPPER BOUND on what each
    # holding can actually serve. The measured per-holding coverage is a run, not a
    # catalogue lookup, and it belongs to RYA-1106.
    hold = pd.read_csv(HOLDINGS)
    inst = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    print("AGSS21 Fe line set vs our VIS holdings")
    print("  (span = the INSTRUMENT's declared range; an upper bound on the holding)")
    for h in VIS_HOLDINGS:
        row = hold[hold["holding_id"].astype(str) == h]
        if row.empty:
            doc["holdings"][h] = {"status": "holding not in the registry"}
            print(f"  {h:34} NOT IN holdings_manifest_registry")
            continue
        iid = str(row.iloc[0]["instrument_id"])
        irow = inst[inst["instrument_id"].astype(str) == iid]
        if irow.empty or pd.isna(irow.iloc[0].get("wavelength_min_nm")):
            doc["holdings"][h] = {"instrument_id": iid, "status": "no span declared"}
            print(f"  {h:34} instrument {iid}: no span declared")
            continue
        lo = float(irow.iloc[0]["wavelength_min_nm"]) * 10.0
        hi = float(irow.iloc[0]["wavelength_max_nm"]) * 10.0
        n = {}
        for ion in ("I", "II"):
            s = asp[asp.ion == ion]
            n[ion] = int(((s.wavelength_air_A >= lo) & (s.wavelength_air_A <= hi)).sum())
        doc["holdings"][h] = {"instrument_id": iid, "span_A": [lo, hi],
                              "span_basis": "instrument_catalog (upper bound)",
                              "Fe I": n["I"], "Fe II": n["II"],
                              "total": n["I"] + n["II"]}
        print(f"  {h:34} {lo:7.1f}-{hi:7.1f} A   Fe I {n['I']:3}/40   "
              f"Fe II {n['II']:3}/13")

    # ── 2. overlap with OUR pools, at the DERIVED tolerance ───────────────────
    tol = TOL_A["asplund_agss21"]
    print(f"\nOverlap with OUR line list, strict lambda+EP (RYA-1037), tol {tol} A")
    print("  (derived from AGSS21's printed 0.1 A resolution, NOT the 0.005 A default)")
    doc["match_tolerance_A"] = {"asplund_agss21": tol, "basis":
                                "half the printed wavelength resolution (nm to 2 dp)"}
    for ion in ("I", "II"):
        s = asp[asp.ion == ion]
        allp = overlap(s, fe_pool(ion, lab_only=False), tol)
        lab = overlap(s, fe_pool(ion, lab_only=True), tol)
        doc["graded_overlap"][f"Fe {ion}"] = {"in_canonical_gf": allp, "in_LAB_tier": lab}
        print(f"  Fe {ion:3} n={allp['n_ref']:3}  in canonical_gf {allp['matched']:3}"
              f"/{allp['n_ref']:<3}   in GRADED (LAB) tier {lab['matched']:3}"
              f"/{lab['n_ref']:<3}")
        doc["graded_overlap"][f"Fe {ion}"]["plateau_LAB"] = plateau(
            s, fe_pool(ion, lab_only=True))

    print("\n  PLATEAU (Fe I vs LAB tier) - a real overlap stops moving, a coincidence "
          "keeps climbing:")
    pl = doc["graded_overlap"]["Fe I"]["plateau_LAB"]
    print("    " + "  ".join(f"{t}A:{n}" for t, n in pl.items()))

    # ── 3. WHY the non-LAB lines are missing: tier, not absence ───────────────
    fe1 = asp[asp.ion == "I"].reset_index(drop=True)
    allp = fe_pool("I", lab_only=False)
    idx = np.asarray(_match(fe1, allp, tol).index)
    low = (fe1.elo_eV < 2.85).to_numpy()
    held = [str(allp.iloc[j]["gf_tier"]) for j in idx[low] if j >= 0]
    doc["low_elo"] = {
        "n_below_2p85_eV": int(low.sum()),
        "present_in_canonical_gf": int((idx[low] >= 0).sum()),
        "in_LAB_tier": 0 if not held else sum(t == "LAB" for t in held),
        "tiers_we_hold_them_at": dict(Counter(held)),
    }
    d = doc["low_elo"]
    print(f"\n  LOW-Elo (<2.85 eV, below our graded floor): n={d['n_below_2p85_eV']}")
    print(f"    present in canonical_gf: {d['present_in_canonical_gf']}   "
          f"in LAB tier: {d['in_LAB_tier']}")
    print(f"    we hold them at: {d['tiers_we_hold_them_at']}")

    # ── 4. the same measurement on AP2002, at ITS OWN precision ───────────────
    if AP2002.exists():
        ap = pd.read_csv(AP2002)
        t2 = TOL_A["ap2002"]
        print(f"\nSame measurement on AP2002 (a DIFFERENT paper's list), tol {t2} A "
              f"(it prints A to 2 dp):")
        for ion in ("I", "II"):
            ov = overlap(ap[ap.ion == ion], fe_pool(ion, lab_only=True), t2)
            ov["plateau_LAB"] = plateau(ap[ap.ion == ion], fe_pool(ion, lab_only=True))
            doc["ap2002_comparison"][f"Fe {ion}"] = ov
            print(f"  Fe {ion:3} AP2002 {ov['n_ref']:3} vs graded pool {ov['n_pool']:4}  "
                  f"-> MATCHED {ov['matched']:3}  {ov['matched_lambda']}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
