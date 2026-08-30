#!/usr/bin/env python3
"""RYA-1111 - where the gf for the 21 GBS lines Jofre never published could come from.

    python3 scripts/rya1111_gbs_missing_gf_options.py

RYA-1110 built the GBS reference set and found 21 of its 159 lines carry no log gf: their
`gf_provenance_gbs` reads "NOT PUBLISHED IN TABLES 4/5". RYA-1111 refuses to measure them
rather than substituting a value. This script answers the next question - WHAT COULD
SUPPLY THEM, and can that source be VALIDATED rather than merely asserted.

🔴 HEITER+2021 IS NOT A SUBSTITUTE SOURCE. IT IS THE PUBLICATION OF THE LIST GBS USED.
Jofre+2014's line list is GES-v3, cited there as "Heiter et al., in prep". Heiter et al.
2021 (A&A 645, A106; `heiter2021` in bibliography.csv) is that work published. So taking a
gf from it is not reaching for a different table - it is reading the same list from its
citable version. That is the whole reason this is answerable at all.

⚠️ BUT "same list, later version" IS A CLAIM, AND VERSION DRIFT IS REAL. GES-v3 is 2014;
Heiter+2021 is seven years of revisions later. Whether a given line's gf survived that is
not something to assume.

🔴 THE VALIDATION IS FREE, AND IT IS THE POINT OF THIS SCRIPT. GBS publishes a gf for the
OTHER 138 lines. Those are a held-out set: for each we can ask whether Heiter+2021
reproduces the value Jofre actually printed. That converts "Heiter is probably the right
source" into a measured per-reference-code reproduction rate, and it is measured on the
same table, the same species and the same list version boundary as the 21 in question.

The result is counterintuitive and worth stating plainly. Splitting the 138 by GES's OWN
quality flag:

    flag U   57 lines   57/57 reproduce Jofre EXACTLY   max deviation 0.000 dex
    flag Y   81 lines   36/81 reproduce exactly          max deviation 0.180 dex

The lines GES itself marks UNDECIDED are the ones that transfer perfectly, and the ones it
marks GOOD are the ones that moved. That is not a paradox: a `U` line is one nobody
revisited, so it carries the same old reference (mostly MRW) unchanged across versions,
while a `Y` line is one that GOT updated with newer laboratory data between GES-v3 and
Heiter+2021 - which is exactly what makes it disagree with what Jofre printed in 2014.

So the flag is not a quality signal here. THE REFERENCE CODE IS, and this script grades
each of the 21 by the measured reproduction rate of its own code.

⚠️ EVIDENCE STRENGTH IS GRADED, because "2 of 2 exact" is not validation. A code seen on
fewer than `MIN_N` published lines is reported THIN, not adoptable - the same refusal to
read a small-n agreement as a result that RYA-282 applies to abundances.

NOTHING HERE ADOPTS ANYTHING. It writes an options report. Which gf the GBS replication
runs on is RYA-1110's open decision and Ryan's to make (RYA-161).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import line_match  # noqa: E402

GBS = ROOT / "data" / "linelists" / "reference_sets" / "gbs_solar_fe_rya1110.csv"
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
OUT = ROOT / "data" / "results" / "rya1111" / "gbs_missing_gf_options.json"

#: Below this many held-out lines, a code's agreement is THIN, not validation.
MIN_N = 5
#: RYA-1110's derived tolerance for this set (re-checked in RYA-1117).
TOL_A = 0.015


def main() -> int:
    d = pd.read_csv(GBS)
    pub = d[d.log_gf_gbs.notna()].copy()
    miss = d[d.log_gf_gbs.isna()].copy()
    pub["dev"] = (pub.log_gf_gbs - pub.heiter2021_log_gf).abs()

    by_flag = {
        str(f): {"n": int(len(g)), "exact": int((g.dev < 1e-9).sum()),
                 "median_dev_dex": round(float(g.dev.median()), 4),
                 "max_dev_dex": round(float(g.dev.max()), 4)}
        for f, g in pub.groupby("heiter2021_gfflag")}

    rate = {}
    for code, g in pub.groupby("heiter2021_r_loggf"):
        rate[str(code)] = {"n": int(len(g)), "exact": int((g.dev < 1e-9).sum()),
                           "median_dev_dex": round(float(g.dev.median()), 4),
                           "max_dev_dex": round(float(g.dev.max()), 4)}

    print("HELD-OUT VALIDATION: does Heiter+2021 reproduce the gf Jofre PRINTED?")
    print("  (measured on the 138 GBS lines that DO carry a published gf)\n")
    print(f"  by GES's own flag:")
    for f, s in sorted(by_flag.items()):
        print(f"    flag {f}  n={s['n']:3}  exact {s['exact']:3}/{s['n']:<3} "
              f"median {s['median_dev_dex']:.4f}  max {s['max_dev_dex']:.4f} dex")
    print("\n  -> the UNDECIDED lines transfer perfectly and the GOOD ones moved: a `U`"
          "\n     line is one nobody revisited, so its old reference carried across"
          "\n     versions unchanged. The flag is not the quality signal here.\n")

    rows = []
    print("THE 21 MISSING LINES, graded by their OWN reference code's measured rate:\n")
    print(f"  {'species':7} {'lambda':>9} {'flag':>4} {'code':<18} {'heiter gf':>9}  "
          f"{'code n':>6} {'exact':>6} {'maxdev':>7}  verdict")
    for r in miss.itertuples():
        code = str(r.heiter2021_r_loggf)
        s = rate.get(code)
        if s is None:
            verdict, why = "UNVALIDATABLE", "code appears on none of the 138 held-out lines"
        elif s["exact"] == s["n"] and s["n"] >= MIN_N:
            verdict, why = "ADOPTABLE", f"reproduces Jofre exactly on {s['n']}/{s['n']}"
        elif s["exact"] == s["n"]:
            verdict, why = "THIN", f"exact but only n={s['n']} held-out lines (< {MIN_N})"
        else:
            verdict, why = "RISKY", (f"{s['exact']}/{s['n']} exact, deviations to "
                                     f"{s['max_dev_dex']} dex")
        rows.append({"species": r.species,
                     "wavelength_air_A": float(r.wavelength_air_A),
                     "elo_eV": float(r.excitation_potential_eV),
                     "ges_flag": str(r.heiter2021_gfflag),
                     "reference_code": code,
                     "heiter2021_log_gf": float(r.heiter2021_log_gf),
                     "our_held_log_gf": float(r.log_gf_ours),
                     "code_stats": s, "verdict": verdict, "why": why})
        n = s["n"] if s else 0
        ex = s["exact"] if s else 0
        mx = s["max_dev_dex"] if s else float("nan")
        print(f"  {r.species:7} {r.wavelength_air_A:9.2f} {r.heiter2021_gfflag:>4} "
              f"{code:<18} {r.heiter2021_log_gf:9.3f}  {n:6} {ex:6} {mx:7.3f}  "
              f"{verdict}")

    tally = {v: sum(1 for x in rows if x["verdict"] == v) for v in
             ("ADOPTABLE", "THIN", "RISKY", "UNVALIDATABLE")}
    print(f"\n  {tally}")

    # ── what else we can reach for these lines ───────────────────────────────
    other = {}
    can = pd.read_csv(CANON, low_memory=False)
    for sp in sorted(miss.species.unique()):
        s = miss[miss.species == sp]
        pool = can[can.species == sp].dropna(
            subset=["wavelength_air_A", "excitation_potential_eV"]).reset_index(drop=True)
        idx = np.asarray(line_match.match(
            s.wavelength_air_A.to_numpy(float), pool.wavelength_air_A.to_numpy(float),
            want_ep=s.excitation_potential_eV.to_numpy(float),
            src_ep=pool.excitation_potential_eV.to_numpy(float),
            require_ep=True, tol_A=TOL_A).index)
        got = pool.iloc[[j for j in idx if j >= 0]]
        other[sp] = {
            "matched_into_canonical_gf": int((idx >= 0).sum()), "n": int(len(s)),
            "gf_tier": {k: int(v) for k, v in got.gf_tier.astype(str).value_counts().items()},
            "with_a_nist_grade": int(got.nist_grade.notna().sum()),
            "with_a_vald_value": int(got.gf_linelist_vald.notna().sum()),
        }
    print("\nOTHER SOURCES WE ALREADY REACH FOR THESE LINES:")
    for sp, o in other.items():
        print(f"  {sp}: {o['matched_into_canonical_gf']}/{o['n']} in canonical_gf; "
              f"tiers {o['gf_tier']}; NIST-graded {o['with_a_nist_grade']}; "
              f"VALD {o['with_a_vald_value']}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "ticket": "RYA-1111",
        "question": ("where could the gf come from for the 21 GBS lines Jofre+2014 "
                     "Tables 4/5 do not publish?"),
        "candidate_source": ("Heiter et al. 2021, A&A 645, A106 -- the PUBLISHED version "
                             "of the GES line list that Jofre+2014 used as GES-v3 "
                             "('Heiter et al., in prep'). Not a substitute table: the "
                             "same list, later version."),
        "validation_method": ("held-out: GBS publishes a gf for the other 138 lines, so "
                              "Heiter+2021 can be scored against what Jofre actually "
                              "printed, per reference code, on the same table and the "
                              "same version boundary."),
        "min_n_for_adoptable": MIN_N,
        "held_out_by_ges_flag": by_flag,
        "held_out_by_reference_code": rate,
        "per_line": rows, "tally": tally,
        "other_sources_reachable": other,
        "decision": ("NOT TAKEN HERE. Which gf the GBS replication runs on is RYA-1110's "
                     "open decision and Ryan's (RYA-161). This is an options report."),
    }, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
