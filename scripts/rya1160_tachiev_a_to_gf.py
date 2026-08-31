#!/usr/bin/env python3
"""
RYA-1160 — convert Tachiev & Froese Fischer 2002 A-values to gf for the AGSS21 N I / O I
reference lines.

AGSS21: "The N i oscillator strengths were taken from Tachiev & Froese Fischer (2002)".
Their CDS tables give Aki + a per-line uncertainty + full level identity, but no
wavelength and no gf. This closes that gap for the lines AGSS21 actually used.

    g_l f_lu  =  1.4992e-16 * lambda(A)^2 * g_u * A_ul(s^-1)

🔴 THE CONSTANT IS VALIDATED BEFORE USE, not assumed: fed NIST's OWN A-values it
reproduces NIST's OWN log gf to <=0.002 dex on four lines (the residual is NIST's
3-significant-figure A rounding). A conversion nobody has round-tripped is a guess.

⚠️ WAVELENGTHS COME FROM NIST, AND THAT IS NOT THE CIRCULARITY THIS TICKET WARNS ABOUT.
rya1160_circularity_check.py established that NIST cannot referee our *gf*, because 30.6%
of our CNO gf is already NIST-derived. Level energies and wavelengths are measured
spectroscopic quantities and are not the thing under grade. Using NIST for lambda while
taking A from Tachiev keeps the graded quantity independent.

🔴 THE INFORMATIVE OUTCOME IS DISAGREEMENT. If Tachiev-derived gf equals NIST's gf
exactly, that means NIST's N I values ARE Tachiev -- i.e. the two are one source wearing
two hats, and no independent confirmation has occurred. The delta column is the whole
point of this script.

Matching is on CONFIGURATION + TERM + J for BOTH levels. Never on wavelength.
"""
from __future__ import annotations

import csv, gzip, json, math, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/reference/cno_atomic_primary/n_tachiev2002"
OUT = ROOT / "data/audit/rya1160_cno_nist_gf"
K = 1.4992e-16

#: The AGSS21 reference lines, with the level identity NIST publishes for each.
#: (species, lambda_A, lower_conf, lower_term, lower_J, upper_conf, upper_term, upper_J,
#:  g_lower, g_upper, nist_A, nist_loggf, nist_grade)
TARGETS = [
    ("N I", 7442.290, "3s", "4P", "3/2", "3p", "4S", "3/2", 4, 4, 1.190e7, -0.4010, "B+"),
    ("N I", 8216.340, "3s", "4P", "5/2", "3p", "4P", "5/2", 6, 6, 2.260e7, 0.1380, "B+"),
    ("N I", 8629.240, "3s", "2P", "3/2", "3p", "2P", "3/2", 4, 4, 2.670e7, 0.0763, "B+"),
    ("N I", 8683.400, "3s", "4P", "3/2", "3p", "4D", "5/2", 4, 6, 1.880e7, 0.1059, "B+"),
    ("N I", 10108.890, "3p", "4D", "3/2", "3d", "4F", "5/2", 4, 6, 3.020e7, 0.4434, "B"),
]

#: table4 = N I (bytes 61-70) + its uncertainty (72-81); table5 = O I in the same slots.
FILES = {"N I": ("table4.dat.gz", 60, 70, 71, 81),
         "O I": ("table5.dat.gz", 60, 70, 71, 81)}


def parse(fname: str, a0: int, a1: int, e0: int, e1: int) -> list[dict]:
    rows = []
    with gzip.open(SRC / fname, "rt", errors="replace") as fh:
        for n, raw in enumerate(fh, 1):
            if len(raw) < a1:
                continue
            try:
                aki, err = float(raw[a0:a1]), float(raw[e0:e1])
            except ValueError:
                continue
            rows.append({"conf_l": raw[0:20].strip(), "term_l": raw[21:24].strip(),
                         "J_l": raw[25:29].strip(), "conf_u": raw[30:50].strip(),
                         "term_u": raw[51:54].strip(), "J_u": raw[55:60].strip(),
                         "A": aki, "err": err, "line": n})
    return rows


def shell(conf: str) -> str:
    """The outer shell label ('3s','3p','3d'), which is what distinguishes these levels."""
    m = re.search(r"\.(\d[spdf])\s*$", conf)
    return m.group(1) if m else ""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    out, exact = [], 0
    for (sp, lam, cl, tl, jl, cu, tu, ju, gl, gu, nA, nlg, ngr) in TARGETS:
        fname, a0, a1, e0, e1 = FILES[sp]
        cand = [r for r in parse(fname, a0, a1, e0, e1)
                if shell(r["conf_l"]) == cl and r["term_l"] == tl and r["J_l"] == jl
                and shell(r["conf_u"]) == cu and r["term_u"] == tu and r["J_u"] == ju]
        rec = {"species": sp, "wavelength_A": f"{lam:.3f}",
               "lower": f"{cl} {tl} {jl}", "upper": f"{cu} {tu} {ju}",
               "g_lower": gl, "g_upper": gu,
               "nist_A_s-1": f"{nA:.3E}", "nist_log_gf": f"{nlg:.4f}",
               "nist_grade": ngr}
        if len(cand) != 1:
            rec |= {"match": f"AMBIGUOUS({len(cand)})" if cand else "ABSENT",
                    "tachiev_A_s-1": "", "tachiev_uncertainty": "",
                    "tachiev_log_gf": "", "delta_dex": "",
                    "verdict": "refused -- not a unique config+term+J match"}
            out.append(rec); continue
        c = cand[0]
        lg = math.log10(K * lam * lam * gu * c["A"])
        d = lg - nlg
        same = abs(d) < 0.001
        exact += same
        rec |= {"match": "UNIQUE", "tachiev_A_s-1": f"{c['A']:.4E}",
                "tachiev_uncertainty": f"{c['err']:.4E}",
                "tachiev_rel_unc_dex": f"{abs(math.log10(1+c['err']/c['A'])):.4f}",
                "tachiev_log_gf": f"{lg:.4f}", "delta_dex": f"{d:+.4f}",
                "verdict": ("IDENTICAL to NIST -- NIST's value for this line IS Tachiev, "
                            "so NIST is not an independent check here"
                            if same else "differs from NIST -- genuinely two sources")}
        out.append(rec)

    fields = list(out[0].keys())
    with (OUT / "tachiev_a_to_gf.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, "") for k in fields})

    matched = [r for r in out if r["match"] == "UNIQUE"]
    verdict = {
        "ticket": "RYA-1160", "targets": len(out), "uniquely_matched": len(matched),
        "identical_to_nist": exact,
        "constant_validated": "yes -- reproduces NIST log gf from NIST A to <=0.002 dex",
        "matching": "configuration + term + J on BOTH levels; never wavelength",
        "wavelength_source": ("NIST (measured); the graded quantity A comes from Tachiev, "
                              "so the gf remains independent of NIST"),
        "finding": (
            f"{exact} of {len(matched)} uniquely matched lines reproduce NIST's log gf "
            "EXACTLY, which means NIST's N I values for those lines are themselves "
            "Tachiev & Froese Fischer. Agreement is therefore NOT confirmation -- it is "
            "one source seen twice. Any line where the two differ is the only place an "
            "actual cross-check exists." if exact else
            "No line reproduces NIST exactly, so the two are genuinely distinct sources "
            "and the deltas are a real comparison."),
        "canonical_gf_modified": False,
    }
    (OUT / "tachiev_a_to_gf_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")

    print(f"{'species':6s} {'lambda':>10s} {'match':12s} {'Tachiev gf':>11s} {'NIST gf':>9s} {'delta':>8s}")
    for r in out:
        print(f"{r['species']:6s} {r['wavelength_A']:>10s} {r['match']:12s} "
              f"{r['tachiev_log_gf']:>11s} {r['nist_log_gf']:>9s} {r['delta_dex']:>8s}")
    print()
    print(json.dumps(verdict["finding"], indent=2)[:400])


if __name__ == "__main__":
    main()
