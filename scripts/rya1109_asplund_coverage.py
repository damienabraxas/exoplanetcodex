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
"""
from __future__ import annotations

import csv
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


def graded_fe(ion: str) -> pd.DataFrame:
    d = pd.read_csv(CANON, low_memory=False)
    d = d[(d["species"].astype(str) == f"Fe {ion}") & (d["gf_tier"].astype(str) == "LAB")]
    return d.dropna(subset=["wavelength_air_A", "excitation_potential_eV"])


def overlap(ref: pd.DataFrame, pool: pd.DataFrame) -> dict:
    """Strict lambda+EP match of a reference list against a pool."""
    if ref.empty or pool.empty:
        return {"n_ref": int(len(ref)), "n_pool": int(len(pool)), "matched": 0,
                "matched_lambda": [], "unresolved": int(len(ref))}
    r = line_match.match(ref["wavelength_air_A"].to_numpy(float),
                         pool["wavelength_air_A"].to_numpy(float),
                         want_ep=ref["elo_eV"].to_numpy(float),
                         src_ep=pool["excitation_potential_eV"].to_numpy(float),
                         require_ep=True)
    idx = np.asarray(r.index)
    hit = idx >= 0
    return {
        "n_ref": int(len(ref)),
        "n_pool": int(len(pool)),
        "matched": int(hit.sum()),
        "matched_lambda": [round(float(x), 3)
                           for x in ref["wavelength_air_A"].to_numpy(float)[hit]],
        "unresolved": int((~hit).sum()),
    }


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

    # ── 2. overlap with OUR graded pool ───────────────────────────────────────
    print("\nOverlap with our GRADED (lab-gf) pool, strict lambda+EP (RYA-1037):")
    for ion in ("I", "II"):
        pool = graded_fe(ion)
        ov = overlap(asp[asp.ion == ion], pool)
        doc["graded_overlap"][f"Fe {ion}"] = ov
        print(f"  Fe {ion:3} AGSS21 {ov['n_ref']:3} vs graded pool {ov['n_pool']:4}  "
              f"-> MATCHED {ov['matched']:3}")
        if ov["matched"]:
            print(f"        {ov['matched_lambda']}")

    # ── 3. the same measurement on AP2002, so the two are comparable ──────────
    if AP2002.exists():
        ap = pd.read_csv(AP2002).rename(columns={"elo_eV": "elo_eV"})
        print("\nSame measurement on AP2002 (the list we already held), for comparison:")
        for ion in ("I", "II"):
            ov = overlap(ap[ap.ion == ion], graded_fe(ion))
            doc["ap2002_comparison"][f"Fe {ion}"] = ov
            print(f"  Fe {ion:3} AP2002 {ov['n_ref']:3} vs graded pool {ov['n_pool']:4}  "
                  f"-> MATCHED {ov['matched']:3}  {ov['matched_lambda']}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
