#!/usr/bin/env python3
"""RYA-1233 -- adjudicate the NIST ASD Si I / Si II pull into canonical_gf (governing
process steps 6-7: secure the correct lines, grade them).

For every canonical Si I / Si II row, find the NIST line of the same species on
RYA-1037's dual key: wavelength within TOL_A AND lower-level EP within TOL_EP.
  * exactly one candidate graded C+ or better  -> the row takes NIST's log gf, its OWN
    grade (RYA-1171: never a tighter class than the source), gf_sigma_dex from
    gf_grades.nist_sigma_dex, the source and DOI, and gf_tier NIST-C+;
  * one candidate graded below C+               -> recorded in the audit, row untouched
    (a grade must describe the gf the row USES, so it is not stamped on a Kurucz value);
  * two or more candidates                      -> UNRESOLVED, row untouched (RYA-1072);
  * none                                        -> no NIST record, row untouched.
LAB rows (Den Hartog 2023) are never overwritten. Every decision is written to
data/audit/rya1233_si_nist/adjudication.csv. Re-running is idempotent.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.gf_grades import NIST_ACC_PCT, nist_sigma_dex  # noqa: E402

CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
PULLS = {"Si I": ROOT / "data/linelists/primary_gf/nist_asd_SiI_3000_25000.tsv",
         "Si II": ROOT / "data/linelists/primary_gf/nist_asd_SiII_3000_25000.tsv"}
AUDIT = ROOT / "data" / "audit" / "rya1233_si_nist" / "adjudication.csv"
TOL_A, TOL_EP = 0.005, 0.002            # RYA-1171 / RYA-1160 EP-aware match
C_PLUS_PCT = NIST_ACC_PCT["C+"]
SOURCE = "NIST ASD (Kelleher & Podobedova 2008, JPCRD 37, 1285)"
DOI = "10.1063/1.2734566"   # Crossref-verified 2026-10-03 (bibliography: kelleher2008_si)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    cg = pd.read_csv(CANONICAL, low_memory=False)
    audit = []
    n_promoted = 0
    for sp, path in PULLS.items():
        nist = pd.read_csv(path, sep="\t")
        nist = nist[nist.log_gf.notna() & nist.wavelength_A.notna() & nist.ei_eV.notna()]
        rows = cg.index[cg.species.astype(str) == sp]
        for i in rows:
            r = cg.loc[i]
            w, ep = float(r.wavelength_air_A), float(r.excitation_potential_eV)
            cand = nist[(abs(nist.wavelength_A - w) <= TOL_A) & (abs(nist.ei_eV - ep) <= TOL_EP)]
            rec = {"species": sp, "wavelength_air_A": w, "ep_eV": ep,
                   "old_log_gf": r.log_gf, "old_tier": r.gf_tier, "n_candidates": len(cand)}
            if len(cand) == 0:
                rec["decision"] = "NO_NIST_RECORD"
            elif len(cand) > 1:
                rec["decision"] = "UNRESOLVED_AMBIGUOUS"
            else:
                c = cand.iloc[0]
                g = str(c.nist_grade).strip() if pd.notna(c.nist_grade) else ""
                rec.update(nist_log_gf=round(float(c.log_gf), 3), nist_grade=g,
                           delta_log_gf=round(float(c.log_gf) - float(r.log_gf), 3)
                           if pd.notna(r.log_gf) else np.nan)
                if str(r.gf_tier) == "LAB":
                    rec["decision"] = "KEPT_LAB"
                elif g in NIST_ACC_PCT and NIST_ACC_PCT[g] <= C_PLUS_PCT:
                    rec["decision"] = "PROMOTED_NIST_C+"
                    n_promoted += 1
                    if not a.dry_run:
                        cg.at[i, "log_gf"] = round(float(c.log_gf), 3)
                        cg.at[i, "nist_grade"] = g
                        cg.at[i, "gf_sigma_dex"] = round(nist_sigma_dex(g), 4)
                        cg.at[i, "gf_tier"] = "NIST-C+"
                        cg.at[i, "loggf_reference"] = SOURCE
                        cg.at[i, "gf_source_doi"] = DOI
                elif str(r.gf_tier) == "NIST-C+":
                    # RYA-1171 class: the store claims a NIST grade C+ or better that NIST's
                    # own record does not support (Si I 5793.073: stored B, NIST ASD E).
                    rec["decision"] = f"OVERCLAIM_DEMOTED(stored {r.nist_grade}, NIST {g or 'ungraded'})"
                    if not a.dry_run:
                        cg.at[i, "gf_tier"] = "OTHER"
                        cg.at[i, "nist_grade"] = np.nan
                        cg.at[i, "gf_sigma_dex"] = np.nan
                        cg.at[i, "adjudication_status"] = (
                            f"RYA-1233: NIST ASD grades this line {g or 'ungraded'} "
                            f"(log gf {float(c.log_gf):.3f}); the stored NIST grade "
                            f"{r.nist_grade} was an over-claim and is withdrawn")
                else:
                    rec["decision"] = f"NIST_GRADE_BELOW_C+({g or 'ungraded'})"
            audit.append(rec)
    out = pd.DataFrame(audit)
    print(out.groupby(["species", "decision"]).size().to_string())
    print(f"promoted to NIST-C+: {n_promoted}")
    if a.dry_run:
        return
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(AUDIT, index=False)
    cg.to_csv(CANONICAL, index=False)
    print(f"wrote {CANONICAL.relative_to(ROOT)} and {AUDIT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
