#!/usr/bin/env python3
"""RYA-853 scope 1 — how much of the graded pool the audit has actually refereed.

    python3 scripts/rya853_graded_pool_coverage.py

The ticket says *"For every line carrying a stored grade in canonical_gf, verify the stored
grade MATCHES the actual source grade ... Sweep the whole pool."* The executed pass
refereed the two HAND-MAINTAINED extracts (60 rows) and found 70% of their grades wrong.
This reports what that leaves, because "we audited the graded pool" and "we audited 3% of
the graded pool" are different claims and only one of them is true.

🔴 THE POOL GREW. RYA-853 was written against 843 graded lines. RYA-945 and RYA-822 have
since bulk-ingested Fe I/Fe II from NIST ASD, so the count is now far larger — which means
a coverage figure quoted from the ticket text is stale by construction and has to be
re-measured rather than carried.

⚠️ AND THE REMAINDER IS NOT THE SAME KIND OF RISK, WHICH IS THE POINT OF SPLITTING IT.
The 70%-wrong rows were TRANSCRIBED BY HAND into `nist_reference.csv` /
`nist_crosscheck.csv`, where a grade can simply be typed wrong — that is what happened
(`A` stored where ASD says `D`). The bulk rows carry the grade the ASD QUERY ITSELF
RETURNED: `scripts/rya822_pull_nist_nearuv.py` reads it straight from the response's
`Acc.` column. Fabrication in the hand-transcription sense cannot occur there.

That is NOT a clean bill of health, and this script does not give one. A machine pull has
its own failure mode — matching the WRONG LINE and stamping a correct grade from it — which
is precisely the defect the RYA-853/RYA-1037 EP guard found live in `crosscheck_nist()`
(Fe I 6065.490, NIST-row EP 2.608 against our line's 4.956). So the two classes are counted
separately and named, and neither is reported as the other.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CG = ROOT / "data/linelists/canonical_gf.csv"
EXTRACTS = {
    "nist_reference": ROOT / "data/linelists/nist_reference.csv",
    "nist_crosscheck": ROOT / "data/linelists/nist_crosscheck.csv",
}
LABPOOL = ROOT / "data/results/rya853/rya853_fe1_labpool_referee.json"


def _read_extract(p: Path) -> pd.DataFrame:
    import io
    # ⚠️ `newline=''` on OPEN, not on Path.read_text -- nist_crosscheck.csv has MIXED line
    # endings (21 CRLF, 5 LF) and read_text() normalises them (RYA-853/PR #448). This only
    # counts rows, but the habit is the point: never round-trip these two files.
    with open(p, newline="") as fh:
        body = [l for l in fh.read().splitlines() if not l.startswith("#")]
    return pd.read_csv(io.StringIO("\n".join(body)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/results/rya853/rya853_graded_pool_coverage.json")
    a = ap.parse_args(argv)

    cg = pd.read_csv(CG, low_memory=False)
    grade = cg.nist_grade.astype(str).str.strip()
    graded = cg[cg.nist_grade.notna() & grade.ne("") & grade.ne("nan")]
    ref = graded.loggf_reference.astype(str)

    # The provenance split. `NIST-C+` and `NIST ASD (astroquery)` are the machine pulls;
    # anything else carrying a grade is hand-maintained or inherited.
    machine = ref.str.startswith(("NIST-C+", "NIST ASD (astroquery)"))
    n_machine, n_other = int(machine.sum()), int((~machine).sum())

    refereed_rows = sum(len(_read_extract(p)) for p in EXTRACTS.values())
    labpool = json.loads(LABPOOL.read_text()) if LABPOOL.exists() else {}

    out = {
        "ticket": "RYA-853 scope 1 — graded-pool coverage",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_rows_with_a_stored_grade": int(len(graded)),
        "ticket_text_said": 843,
        "pool_grew_note": (
            "RYA-853 was written against 843 graded lines; RYA-822/RYA-945 have since bulk-"
            "ingested Fe I/Fe II from NIST ASD. Any coverage fraction quoted from the "
            "ticket text is stale — re-measure it, do not carry it."),
        "by_provenance_class": {
            "machine_pull_from_asd": {
                "rows": n_machine,
                "sources": ref[machine].str.split(" T").str[0].value_counts().to_dict(),
                "grade_origin": ("the ASD response's own `Acc.` column, read directly by "
                                 "scripts/rya822_pull_nist_nearuv.py:115"),
                "refereed_by_rya853": False,
                "risk": ("NOT hand-transcription. The failure mode here is matching the "
                         "WRONG LINE and stamping a correct grade off it — the defect the "
                         "EP guard caught live on Fe I 6065.490 (NIST EP 2.608 vs our "
                         "4.956). Un-refereed, and not claimed clean."),
            },
            "hand_maintained_extracts": {
                "rows_in_the_two_extract_files": refereed_rows,
                "refereed_by_rya853": True,
                "result": "70% of uniquely-matched grades disagreed with live NIST ASD",
                "corrections_tabulated": 45,
            },
        },
        "lab_tier_pool": {
            "rows_at_gf_tier_LAB": int(cg.gf_tier.astype(str).eq("LAB").sum()),
            "refereed_against_source_papers": labpool.get("n_refereed"),
            "loggf_mismatches": labpool.get("n_loggf_mismatch"),
            "sigma_mismatches": labpool.get("n_sigma_mismatch"),
            "note": ("scope 2. Refereed against Ruffoni 2014 / Den Hartog 2014 / Belmonte "
                     "2017's own tables — the pool RYA-850 promotes. ⚠️ `rows_at_gf_tier_"
                     "LAB` and `refereed_against_source_papers` are DIFFERENT SETS and "
                     "must not be read as a fraction of each other: the referee runs on "
                     "the curated `data/reference/fe_gf_lab/fe1_lab_loggf.csv` (465 Fe I "
                     "lines), while gf_tier==LAB counts canonical_gf rows across Fe I, "
                     "Fe II and Al I."),
        },
        "verdict": (
            f"The grade audit has refereed the {refereed_rows} rows of the two hand-"
            f"maintained extracts. {n_machine} rows carrying a machine-pulled ASD grade "
            f"and {n_other - refereed_rows if n_other > refereed_rows else 0} other graded "
            f"rows are NOT refereed. 'The graded pool was audited' is therefore FALSE as "
            f"stated; what is true is that the pool where a grade can be typed wrong was "
            f"audited, and it was 70% wrong."),
    }
    p = ROOT / a.out
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2) + "\n")

    print(f"rows with a stored grade : {out['n_rows_with_a_stored_grade']}"
          f"   (ticket text said 843)")
    print(f"  machine-pulled from ASD : {n_machine:>5}   refereed: NO")
    print(f"  hand-maintained extracts: {refereed_rows:>5}   refereed: YES (70% wrong)")
    lt = out["lab_tier_pool"]
    print(f"\n  [separate pool] gf_tier == LAB : {lt['rows_at_gf_tier_LAB']:>5} rows in "
          f"canonical_gf")
    print(f"                  scope-2 referee : {lt['refereed_against_source_papers']:>5} "
          f"lines from data/reference/fe_gf_lab/fe1_lab_loggf.csv, "
          f"{lt['loggf_mismatches']} log gf + {lt['sigma_mismatches']} sigma mismatches")
    print("                  ⚠️ these two counts are DIFFERENT SETS -- the referee runs on "
          "the\n                     curated Fe I lab list, not on the gf_tier column.")
    print(f"\n{out['verdict']}")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
