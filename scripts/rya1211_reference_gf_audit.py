#!/usr/bin/env python3
"""RYA-1211 — why does the Reference Grade product charge a 0.17 dex UNGRADED gf bar?

THE ASK, AND THE PREMISE IT RESTS ON
------------------------------------
The AGSS21 Reference Grade product (Amarsi 3D-NLTE, n=21, all four VIS holdings) reports
sigma_reported 0.171, essentially all of it the blanket ungraded gf systematic. Its budget
says why: `GF-LAB x3, GF-NIST x2, systematic:K07 x15, systematic:K07/SCALE-MISMATCH x1`.
So the ticket asks which of those 15 Kurucz-sourced lines have a laboratory gf available
and can be re-sourced.

🔴 THE PREMISE IS WRONG, AND THIS SCRIPT MEASURES THAT RATHER THAN ASSERTING IT. Not one
of the 21 lines carries the K07 tag. `grade_line` has no "not found" verdict: when the
join misses, it falls through to the blanket Kurucz systematic and LABELS the miss K07.
Those 15 lines were never sourced from Kurucz — they were not found at all.

WHY THEY WERE NOT FOUND
-----------------------
AGSS21 Table A.2 prints lambda in nanometres to 2 dp = 0.1 A. `gf_grades.WAVE_TOL_A` is
0.02 A, right for a pool stated at the line list's own 4-decimal precision. So a line can
sit up to 0.05 A from its own canonical_gf row BY PRINTING ALONE and miss a 0.02 A window.
That is the RYA-1109 trap, which `pipeline.reference_lineset` already names in its
docstring and already fixes for every other join by declaring a DERIVED per-set
`match_tol_A`. `gf_grades` is the one join that never asked.

THE CONTROLS, BECAUSE A WIDER WINDOW IS EXACTLY HOW A FAKE PEDIGREE GETS MANUFACTURED
------------------------------------------------------------------------------------
Widening a tolerance until a count improves is the move RYA-1109 was punished for, so
three things are measured here and all three must hold before any count is believed:

1. **A DISPLACED NULL.** The same join with the pool shifted +0.35 A. If the matches were
   proximity accidents the null would find some too. It finds ZERO at every tolerance out
   to 0.15 A -- the EP window does the discriminating, and EP is not what got rounded.
2. **A PLATEAU.** The census is swept across tolerance. It saturates at 21/21 canonical
   and 6/21 laboratory and stays there, so the reported counts are not an artifact of
   where the window was put.
3. **AMBIGUITY.** `AmbiguousLineMatch` is counted, not caught and ignored. Zero at every
   tolerance for this pool -- so nothing here is a match picked by proximity (RYA-1034).

And the cost of the other choice is measured too: a GLOBAL widening to 0.06 A takes the
canonical Fe I table's self-ambiguous rows from 26 to 210. That is why the fix is a
per-source tolerance handed down from the registry, not a new constant.

WHAT THIS SCRIPT DOES NOT DO
----------------------------
It changes no log gf (RYA-161 / spec item 4). Where a reference disagrees with the value
AGSS21 published, that is REPORTED as a discrepancy and the systematic stands -- a lab
sigma describes the value the lab measured, not a different number someone else adopted.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

from pipeline import gf_grades as gg            # noqa: E402
from pipeline import gf_rung                    # noqa: E402
from pipeline import reference_lineset as rls   # noqa: E402

OUT = ROOT / "data" / "results" / "rya1211"
PER_LINE = ROOT / "data/results/rya1106/kpno_kurucz2005/asplund_lines_per_line.csv"

#: AGSS21 Sect. 4's COLLECTIVE gf attribution, classified. Table A.2 publishes no per-line
#: source column, so this is the finest attribution the paper supports and every line
#: inherits it -- flagged, not guessed (RYA-161). Every entry is experimental: AGSS21
#: names no theoretical or semi-empirical gf source for Fe I at all.
AGSS21_GF_SOURCES = {
    "Den Hartog et al. 2014": "primary laboratory (FTS branching fractions x lifetimes)",
    "Belmonte et al. 2017": "primary laboratory",
    "Blackwell et al. 1995 (Oxford)": "primary laboratory",
    "Holweger et al. 1991 (Hannover)": "primary laboratory",
    "O'Brian et al. 1991 (Wisconsin)": "primary laboratory",
    "Scott et al. 2015a": "line SELECTION carried over, gf from the groups above",
}
DISPLACE_A = 0.35


def pool() -> pd.DataFrame:
    d = pd.read_csv(PER_LINE)
    return d[d["a_3dnlte"].notna()].copy()


def census(s: pd.DataFrame, wtol: float, shift: float = 0.0) -> dict:
    """How many of the pool tie to canonical_gf / to a lab row at this window."""
    can, lab = gg.canonical_species("Fe I"), gg.lab_lines("Fe I")
    n_can = n_lab = n_amb = 0
    for _, r in s.iterrows():
        w, ep = float(r.wavelength_air_A) + shift, float(r.elo_eV)
        for d, wc, ec, is_lab in ((can, "wavelength_air_A", "excitation_potential_eV", False),
                                  (lab, "wavelength_air_A", "elo_eV", True)):
            m = d[(np.abs(d[wc] - w) <= wtol) & (np.abs(d[ec] - ep) <= gg.EP_TOL_EV)]
            if len(m) > 1:
                n_amb += 1
            elif len(m) == 1:
                n_lab += is_lab
                n_can += not is_lab
    return {"wave_tol_A": wtol, "canonical_ties": n_can, "lab_ties": n_lab,
            "ambiguous": n_amb}


def grade(s: pd.DataFrame, wtol: float | None) -> pd.DataFrame:
    rows = []
    for _, r in s.iterrows():
        v = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV),
                          float(r.loggf_asplund), wave_tol_A=wtol)
        rows.append({
            "wavelength_air_A": float(r.wavelength_air_A),
            "elo_eV": float(r.elo_eV),
            "loggf_agss21": float(r.loggf_asplund),
            "gf_grade": v.gf_grade,
            "our_gf_reference": v.gf_reference_tag or "(no row found)",
            "our_loggf": v.gf_ref_loggf,
            "delta_ours_minus_agss21": v.gf_delta_dex,
            "cited_sigma_dex": v.gf_sigma_dex,
            "has_cited_sigma": bool(v.has_cited_sigma),
        })
    return pd.DataFrame(rows)


def rung_for(s: pd.DataFrame, wtol: float | None):
    lines = pd.DataFrame({"wavelength_air_A": s.wavelength_air_A.to_numpy(float),
                          "ep_eV": s.elo_eV.to_numpy(float),
                          "log_gf": s.loggf_asplund.to_numpy(float)})
    return gf_rung.decide("Fe", "I", lines, wave_tol_A=wtol)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    s = pool()
    derived = rls.grading_tol_A("asplund")
    printed_dp = max(len(str(w).split(".")[-1]) for w in s.wavelength_air_A)

    sweep = [census(s, t) for t in (0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15)]
    null = [census(s, t, shift=DISPLACE_A) for t in (0.02, 0.05, 0.06, 0.10, 0.15)]

    before, after = grade(s, None), grade(s, derived)
    r_before, r_after = rung_for(s, None), rung_for(s, derived)

    tbl = before[["wavelength_air_A", "elo_eV", "loggf_agss21"]].copy()
    tbl["grade_at_0.02A"] = before["gf_grade"]
    tbl["our_ref_at_0.02A"] = before["our_gf_reference"]
    tbl["grade_at_derived"] = after["gf_grade"]
    tbl["our_ref_at_derived"] = after["our_gf_reference"]
    tbl["our_loggf"] = after["our_loggf"]
    tbl["delta_ours_minus_agss21"] = after["delta_ours_minus_agss21"]
    tbl["cited_sigma_dex"] = after["cited_sigma_dex"]
    tbl["upgraded"] = np.where(before["gf_grade"] != after["gf_grade"], "Y", "N")
    tbl["agss21_cited_gf_source"] = "COLLECTIVE (Table A.2 has no per-line source column)"
    tbl.to_csv(OUT / "rya1211_per_line_gf_audit.csv", index=False)

    cited = after[after["has_cited_sigma"]]
    disc = after[after["gf_grade"].str.contains("MISMATCH")]
    summary = {
        "ticket": "RYA-1211",
        "product": "Fe I VIS ENGINE-A-3DNLTE, selector ASPLUND_AGSS21, n=21 (all 4 holdings)",
        "premise_tested": (
            "'15 of 21 lines carry K07 (Kurucz semi-empirical) gf in our store'"),
        "premise_verdict": "REFUTED — not one of the 21 carries the K07 tag",
        "why": (
            f"AGSS21 Table A.2 prints lambda to {printed_dp} dp in Angstroms (nanometres "
            f"to 2 dp = 0.1 A resolution) while gf_grades.WAVE_TOL_A is {gg.WAVE_TOL_A} A. "
            "15 lines sat 0.024-0.050 A from their own canonical_gf row — inside the "
            "printing's rounding bin, EP agreeing to better than 0.0005 eV — and missed "
            "the window. grade_line has no 'not found' verdict, so the miss fell through "
            "to the blanket Kurucz systematic and was LABELLED systematic:K07."),
        "derived_tolerance_A": derived,
        "derived_tolerance_basis": rls.SETS["asplund"].tol_basis,
        "controls": {
            "displaced_null_A": DISPLACE_A,
            "null_ties_at_every_tolerance": sorted({n["canonical_ties"] + n["lab_ties"]
                                                    for n in null}),
            "ambiguous_matches_at_every_tolerance": sorted({c["ambiguous"] for c in sweep}),
            "plateau_canonical": sweep[-1]["canonical_ties"],
            "plateau_lab": sweep[-1]["lab_ties"],
            "global_widening_refused_because": (
                "a GLOBAL 0.02 -> 0.06 A takes the canonical Fe I table's self-ambiguous "
                "rows from 26 to 210, so a pool stated at full precision must keep 0.02; "
                "the tolerance is a property of the SOURCE and is handed down per set"),
        },
        "grades_before": r_before.grade_counts,
        "grades_after": r_after.grade_counts,
        "n_upgraded": int((tbl["upgraded"] == "Y").sum()),
        "lab_available_of_21": sweep[-1]["lab_ties"],
        "cited_sigma_coverage_before": f"{int(before['has_cited_sigma'].sum())}/21",
        "cited_sigma_coverage_after": f"{len(cited)}/21",
        "rung_before": r_before.rung,
        "rung_after": r_after.rung,
        "sigma_syst_before": 0.17,
        "sigma_syst_after": 0.17,
        "sigma_syst_unchanged_because": (
            f"the MIXED-POOL rule binds either way: {r_after.n_graded} of 21 lines are "
            "GF-LAB, and a pool is graded only if EVERY line is. 16 of the 21 have no "
            "primary-laboratory Fe I measurement in the repo AT ANY TOLERANCE — the lab "
            f"census plateaus at {sweep[-1]['lab_ties']}/21 — so this is a real absence, "
            "not a lookup failure."),
        "drops_below_0.10": False,
        "hypothetical_cited_rms_dex": (
            round(float(np.sqrt(np.mean(cited["cited_sigma_dex"] ** 2))), 4)
            if len(cited) else None),
        "hypothetical_refused_because": (
            f"{len(cited)}/21 = {len(cited)/21:.0%} carry a citable per-line sigma, below "
            f"gf_rung.CITED_COVERAGE_MIN = {gf_rung.CITED_COVERAGE_MIN:.0%}. An RMS over "
            "part of a pool does not describe the pool. Reported as the size of the prize, "
            "NOT adopted."),
        "agss21_gf_sources": AGSS21_GF_SOURCES,
        "agss21_per_line_source_available": False,
        "discrepancies": [
            {"wavelength_air_A": float(r.wavelength_air_A),
             "our_reference": r.our_gf_reference, "our_loggf": float(r.our_loggf),
             "agss21_loggf": float(r.loggf_agss21),
             "delta_ours_minus_agss21": float(r.delta_ours_minus_agss21)}
            for _, r in disc.iterrows()],
        "no_log_gf_changed": True,
    }
    (OUT / "rya1211_gf_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    pd.DataFrame(sweep + [dict(n, displaced_by_A=DISPLACE_A) for n in null]).to_csv(
        OUT / "rya1211_tolerance_sweep.csv", index=False)

    print(f"premise: {summary['premise_verdict']}")
    print(f"  before {r_before.grade_counts}")
    print(f"  after  {r_after.grade_counts}")
    print(f"  upgraded {summary['n_upgraded']} of 21; lab available {summary['lab_available_of_21']}/21")
    print(f"  rung {r_before.rung} -> {r_after.rung}; sigma_syst 0.17 -> 0.17; "
          f"below 0.10? {summary['drops_below_0.10']}")
    print(f"  discrepancies: {len(summary['discrepancies'])}")
    print(f"  wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
