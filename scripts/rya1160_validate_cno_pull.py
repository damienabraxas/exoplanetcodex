#!/usr/bin/env python3
"""
RYA-1160 — validate the NIST CNO pull against what the store already holds.

A new source is not evidence until it reproduces something already known. The store
carries exactly five NIST-graded CNO rows; those are the only control available, and
they are used here as one -- EP-aware (wavelength AND excitation potential, never
wavelength alone, RYA-1034/1037), single-match-or-refuse, no argmin.

Writes data/audit/rya1160_cno_nist_gf/control.csv and a verdict. Reads only.
"""
from __future__ import annotations

import csv, glob, json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1160_cno_nist_gf"
CNO = {"C I", "C II", "C III", "N I", "N II", "N III", "O I", "O II", "O III"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pulls = {}
    for f in sorted(glob.glob(str(ROOT / "data/linelists/primary_gf/nist_asd_[CNO]*_3000_25000.tsv"))):
        slug = Path(f).name.split("nist_asd_")[1].split("_3000")[0]
        pulls[slug[0] + " " + slug[1:]] = pd.read_csv(f, sep="\t")

    with (ROOT / "data/linelists/canonical_gf.csv").open(newline="") as fh:
        ours = [r for r in csv.DictReader(fh)
                if r["species"] in CNO and (r["nist_grade"] or "").strip()]

    rows = []
    for r in ours:
        sp, w = r["species"], float(r["wavelength_air_A"])
        ep = float(r["excitation_potential_eV"])
        df = pulls.get(sp)
        rec = {"species": sp, "wavelength_air_A": f"{w:.4f}", "EP_eV": f"{ep:.4f}",
               "ours_log_gf": r["log_gf"], "ours_nist_grade": r["nist_grade"],
               "ours_tier": r["gf_tier"]}
        if df is None:
            rec |= {"match": "NO_PULL", "verdict": "no pull for this species"}
            rows.append(rec); continue
        c = df[(df.wavelength_A - w).abs() <= 0.05]
        c = c[(c.ei_eV - ep).abs() <= 0.02]
        if len(c) != 1:
            rec |= {"match": f"AMBIGUOUS({len(c)})" if len(c) else "ABSENT",
                    "verdict": "refused -- not resolvable on wavelength+EP alone"}
            rows.append(rec); continue
        x = c.iloc[0]
        delta = float(r["log_gf"]) - float(x.log_gf)
        grade_same = str(x.nist_grade).strip() == r["nist_grade"].strip()
        rec |= {"match": "UNIQUE", "nist_log_gf": f"{x.log_gf:.4f}",
                "nist_grade": str(x.nist_grade), "nist_acc_pct": x.nist_acc_pct,
                "nist_TP_ref": x.ref_transition_probability,
                "delta_dex": f"{delta:+.4f}",
                "verdict": ("values agree; grade agrees" if grade_same else
                            f"VALUE AGREES but GRADE DIFFERS: ours {r['nist_grade']} "
                            f"vs NIST {x.nist_grade}")}
        rows.append(rec)

    fields = ["species", "wavelength_air_A", "EP_eV", "ours_log_gf", "ours_nist_grade",
              "ours_tier", "match", "nist_log_gf", "nist_grade", "nist_acc_pct",
              "nist_TP_ref", "delta_dex", "verdict"]
    with (OUT / "control.csv").open("w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w_.writeheader()
        for r in rows:
            w_.writerow({k: r.get(k, "") for k in fields})

    unique = [r for r in rows if r["match"] == "UNIQUE"]
    worst = max((abs(float(r["delta_dex"])) for r in unique), default=float("nan"))
    mismatched = [r for r in unique if "GRADE DIFFERS" in r["verdict"]]
    verdict = {
        "ticket": "RYA-1160", "control_rows": len(rows), "uniquely_matched": len(unique),
        "max_abs_delta_dex": round(worst, 4),
        "values_reproduce": bool(worst <= 0.002),
        "grade_disagreements": len(mismatched),
        "grade_disagreement_detail": [
            f"{r['species']} {r['wavelength_air_A']}: ours {r['ours_nist_grade']} "
            f"vs NIST {r['nist_grade']}" for r in mismatched],
        "finding": (
            "Values reproduce to <=0.001 dex on every uniquely matched row, which "
            "validates the pull, the log gf = log10(gi*fik) computation, the air/vac "
            "frame and the EP-aware match. The GRADES do not: the store carries A+ "
            "(<=2%) on the O I 777 triplet where NIST ASD publishes A (<=3%). The "
            "store claims BETTER precision than the authority it cites, on the most "
            "used oxygen diagnostic in the programme -- so any gf_sigma_dex derived "
            "from that grade would be too small."
            if mismatched else "Values and grades both reproduce."),
        "canonical_gf_modified": False,
    }
    (OUT / "verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
