#!/usr/bin/env python3
"""
RYA-1160 — adjudicate the acquired CNO gf evidence into canonical_gf.

Populates the grade evidence that was blank on all 2,606 CNO atomic rows:
`nist_grade`, `gf_sigma_dex`, `gf_source_doi`, and `gf_tier` where a rung is earned.

🔴 VALIDATE-DON'T-TUNE (RYA-161). This attaches EVIDENCE; it never rewrites a
`log_gf`. Where our value and the source disagree the delta is RECORDED, not
applied -- adopting a new value is a separate, deliberate decision.

🔴 MATCHING IS EP-AWARE (RYA-1034/1037): wavelength AND excitation potential,
single-match-or-refuse, no argmin. A wavelength-only key is what put an unearned
lab sigma on a blend in RYA-1009.

RUNG, per RYA-946 (lab > NIST-C+ > Kurucz > VALD):
  * NIST ASD graded rows earn `NIST-C+`. They are critically evaluated but largely
    THEORETICAL for C/N/O (Opacity Project, then MCHF) -- never `LAB`.
  * `gf_sigma_dex` is derived from the NIST accuracy letter via NIST_ACC_PCT:
    sigma_dex = |log10(1 +/- pct/100)|, taking the larger side.
"""
from __future__ import annotations

import csv, glob, json, math
from bisect import bisect_left, bisect_right
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data/linelists/canonical_gf.csv"
OUT = ROOT / "data/audit/rya1160_cno_nist_gf"
CNO = {"C I", "C II", "C III", "N I", "N II", "N III", "O I", "O II", "O III"}

NIST_ACC_PCT = {"AAA": 0.3, "AA": 1.0, "A+": 2.0, "A": 3.0, "B+": 7.0, "B": 10.0,
                "C+": 18.0, "C": 25.0, "D+": 40.0, "D": 50.0, "E": 100.0}
DOI_NIST = "10.18434/T4W30F"          # NIST ASD, the citable dataset DOI
QUARANTINE: list[dict] = []
W_TOL, EP_TOL = 0.05, 0.02            # Angstrom, eV

#: 🔴 PHYSICAL-VALIDITY GATE. The oscillator-strength sum rule bounds f_ik at order
#: unity; a few strong resonance lines reach ~2, nothing reaches 10. Some NIST ASD rows
#: come back through astroquery with f_ik of 1e5-1e7 -- e.g. C II 5207.977 returns
#: fik=8.0e5 alongside Aki=2.0e5, which is not an f-value at all. Whatever that column
#: holds for those rows, it is not what tidy() assumed, so log gf = log10(gi*fik) is
#: meaningless there. They are QUARANTINED, not silently used and not silently dropped:
#: attaching a +6.68 dex "grade" to a real line would be far worse than attaching none.
FIK_MAX = 10.0


def sigma_dex(grade: str) -> float | None:
    pct = NIST_ACC_PCT.get(str(grade).strip())
    if pct is None:
        return None
    return abs(math.log10(1 - pct / 100)) if pct < 100 else 1.0


def load_nist() -> dict[str, list[tuple]]:
    idx: dict[str, list[tuple]] = {}
    for f in sorted(glob.glob(str(ROOT / "data/linelists/primary_gf/nist_asd_[CNO]*_*.tsv"))):
        slug = Path(f).name.split("nist_asd_")[1].rsplit("_", 2)[0]
        sp = slug[0] + " " + slug[1:]
        d = pd.read_csv(f, sep="\t")
        d = d[d.nist_grade.notna() & d.wavelength_A.notna() & d.ei_eV.notna()]
        bad = d[d.fik.notna() & (d.fik > FIK_MAX)]
        for _, r in bad.iterrows():
            QUARANTINE.append({"species": sp, "wavelength_A": f"{r.wavelength_A:.4f}",
                               "fik": f"{r.fik:.4E}", "aki_s-1": f"{r['aki_s-1']:.4E}",
                               "gi": r.gi, "bogus_log_gf": f"{r.log_gf:.4f}",
                               "nist_grade": str(r.nist_grade),
                               "reason": "f_ik exceeds the sum-rule bound; not an "
                                         "oscillator strength -- column meaning unknown"})
        d = d[~(d.fik.notna() & (d.fik > FIK_MAX))]
        for _, r in d.iterrows():
            idx.setdefault(sp, []).append(
                (float(r.wavelength_A), float(r.ei_eV),
                 None if pd.isna(r.log_gf) else float(r.log_gf),
                 str(r.nist_grade).strip(), str(r.ref_transition_probability)))
    for sp in idx:
        idx[sp].sort(key=lambda t: t[0])
    return idx


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    nist = load_nist()
    with CANON.open(newline="") as fh:
        reader = csv.DictReader(fh)
        fields = list(reader.fieldnames or [])
        rows = list(reader)

    ledger, stats = [], {"matched": 0, "ambiguous": 0, "absent": 0, "no_gf": 0}
    for r in rows:
        if r["species"] not in CNO:
            continue
        src = nist.get(r["species"], [])
        if not src:
            stats["absent"] += 1
            continue
        w = float(r["wavelength_air_A"]); ep = float(r["excitation_potential_eV"])
        waves = [t[0] for t in src]
        lo, hi = bisect_left(waves, w - W_TOL), bisect_right(waves, w + W_TOL)
        cand = [t for t in src[lo:hi] if abs(t[1] - ep) <= EP_TOL]
        if len(cand) != 1:
            stats["ambiguous" if cand else "absent"] += 1
            ledger.append({"line_id": r["line_id"], "species": r["species"],
                           "wavelength_air_A": r["wavelength_air_A"],
                           "EP_eV": r["excitation_potential_eV"],
                           "match": f"AMBIGUOUS({len(cand)})" if cand else "ABSENT",
                           "nist_grade": "", "gf_sigma_dex": "", "our_log_gf": r["log_gf"],
                           "nist_log_gf": "", "delta_dex": "",
                           "action": "no evidence attached"})
            continue
        lam, ei, nlg, grade, tp = cand[0]
        sd = sigma_dex(grade)
        delta = "" if nlg is None or not r["log_gf"] else f"{float(r['log_gf']) - nlg:+.4f}"
        if nlg is None:
            stats["no_gf"] += 1
        stats["matched"] += 1
        # EVIDENCE ONLY -- log_gf is never touched.
        r["nist_grade"] = grade
        r["gf_sigma_dex"] = f"{sd:.4f}" if sd is not None else ""
        r["gf_source_doi"] = DOI_NIST
        if r["gf_tier"] in ("KURUCZ", "VALD3", "OTHER") and nlg is not None:
            r["gf_tier"] = "NIST-C+"
        ledger.append({"line_id": r["line_id"], "species": r["species"],
                       "wavelength_air_A": r["wavelength_air_A"],
                       "EP_eV": r["excitation_potential_eV"], "match": "UNIQUE",
                       "nist_grade": grade, "gf_sigma_dex": r["gf_sigma_dex"],
                       "our_log_gf": r["log_gf"],
                       "nist_log_gf": "" if nlg is None else f"{nlg:.4f}",
                       "delta_dex": delta,
                       "action": "evidence attached; log_gf UNCHANGED"})

    with CANON.open("w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w_.writeheader(); w_.writerows(rows)

    if QUARANTINE:
        with (OUT / "nist_fik_quarantine.csv").open("w", newline="") as fh:
            w_ = csv.DictWriter(fh, fieldnames=list(QUARANTINE[0]), lineterminator="\n")
            w_.writeheader(); w_.writerows(QUARANTINE)

    with (OUT / "adjudication_ledger.csv").open("w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=list(ledger[0]), lineterminator="\n")
        w_.writeheader(); w_.writerows(ledger)

    deltas = [abs(float(x["delta_dex"])) for x in ledger if x["delta_dex"]]
    verdict = {"ticket": "RYA-1160", **stats,
               "nist_rows_quarantined_unphysical_fik": len(QUARANTINE),
               "cno_rows_considered": sum(1 for r in rows if r["species"] in CNO),
               "value_deltas_recorded": len(deltas),
               "max_abs_delta_dex": round(max(deltas), 4) if deltas else None,
               "median_abs_delta_dex": round(sorted(deltas)[len(deltas)//2], 4) if deltas else None,
               "log_gf_values_changed": 0,
               "policy": "evidence attached only; RYA-161 validate-don't-tune",
               "matching": "wavelength + EP, single-match-or-refuse, no argmin",
               "tier_awarded": "NIST-C+ (never LAB -- these are evaluated theory)"}
    (OUT / "adjudication_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
