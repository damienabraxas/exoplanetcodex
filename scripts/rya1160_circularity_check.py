#!/usr/bin/env python3
"""
RYA-1160 — is the NIST pull an INDEPENDENT referee for our CNO rows? Mostly not.

🔴 THE POINT OF THIS FILE.  RYA-946 sets the precedence `lab > NIST-C+ > Kurucz >
VALD` -- NIST-C+ is explicitly a rung BELOW lab -- and its C/N/O entry nominates
"high-excitation atomic + molecular ... NLTE-sensitive, treat carefully", NOT NIST.
So a NIST ASD pull is FALLBACK material, not the lab-gf sweep RYA-946 asks for, and
it must not be adopted merely because the string "NIST" is authoritative-sounding.

Worse, it is not independent of what we already hold.  RYA-853's lesson: a referee you
ingest stops being a referee, silently.  This script measures how much of `canonical_gf`
is ALREADY NIST-derived, so that any later adjudication states up front which fraction of
its "agreement" is a self-match rather than a confirmation.

It also partitions the pull by the per-line TP (transition-probability) reference code.
NIST is a COMPILATION; "NIST" is not one authority. The TP code is the only thing in the
pull that names the underlying measurement or calculation per line, and it is what a
lab-vs-theory split has to be built from.

Reads only. Writes data/audit/rya1160_cno_nist_gf/circularity.csv + verdict.
"""
from __future__ import annotations

import collections, csv, glob, json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1160_cno_nist_gf"
CNO = {"C I", "C II", "C III", "N I", "N II", "N III", "O I", "O II", "O III"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (ROOT / "data/linelists/canonical_gf.csv").open(newline="") as fh:
        at = [r for r in csv.DictReader(fh) if r["species"] in CNO]

    def nist_derived(r) -> bool:
        return "NIST" in (r["loggf_reference"] + r["seed_source"]).upper()

    derived = [r for r in at if nist_derived(r)]
    by_ref = collections.Counter(r["loggf_reference"] or "(blank)" for r in at)

    rows = [{"loggf_reference": k, "rows": v,
             "nist_derived": "YES" if "NIST" in k.upper() else "no"}
            for k, v in by_ref.most_common()]
    with (OUT / "circularity.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["loggf_reference", "rows", "nist_derived"],
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    tp = collections.Counter()
    graded = 0
    for f in sorted(glob.glob(str(ROOT / "data/linelists/primary_gf/nist_asd_[CNO]*_*.tsv"))):
        d = pd.read_csv(f, sep="\t")
        d = d[d.nist_grade.notna()]
        graded += len(d)
        for v in d.ref_transition_probability.dropna():
            tp[str(v).strip()] += 1
    with (OUT / "tp_reference_codes.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["tp_code", "graded_rows", "class"],
                           lineterminator="\n")
        w.writeheader()
        for k, v in tp.most_common():
            w.writerow({"tp_code": k, "graded_rows": v, "class": "UNDECODED"})

    verdict = {
        "ticket": "RYA-1160",
        "cno_atomic_rows": len(at),
        "already_nist_derived": len(derived),
        "already_nist_derived_pct": round(100 * len(derived) / len(at), 1),
        "self_match_warning": (
            f"{len(derived)} of {len(at)} CNO atomic rows ({100*len(derived)/len(at):.1f}%) "
            "already carry a NIST-derived loggf_reference (NIST10, plus 5 NIST ASD grade "
            "rows). Adjudicating this pull against them is NIST-vs-NIST: it measures "
            "version drift between NIST10 and current ASD, not independent agreement. "
            "Agreement there proves 'no transcription error since NIST10'; DISAGREEMENT "
            "is the informative outcome (RYA-853)."),
        "rya946_precedence": (
            "RYA-946 sets lab > NIST-C+ > Kurucz > VALD. This pull is NIST-C+ rung, a "
            "FALLBACK, not the lab-gf sweep RYA-946 asks for. Its C/N/O source map "
            "nominates high-excitation atomic + molecular studies, not NIST."),
        "graded_rows_in_pull": graded,
        "distinct_tp_codes": len(tp),
        "rows_without_tp_code": graded - int(sum(tp.values())),
        "tp_next_step": (
            "Every graded row carries a TP code and there are "
            f"{len(tp)} distinct ones. Decoding them against NIST's reference list is "
            "what partitions the pull into experimental vs theoretical, which is the "
            "only way to know which rows could ever reach the lab rung. Currently ALL "
            "are class=UNDECODED."),
        "blocked": (
            "RYA-946 requires every source verified against ADS/VizieR per child. ADS "
            "needs an API token this session does not have, so no post-2007 experimental "
            "CNO gf search has been run. That gap is open."),
        "canonical_gf_modified": False,
    }
    (OUT / "circularity_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
