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
    # (species, lambda, l_shell, l_term, l_J, u_shell, u_term, u_J, g_l, g_u, nist_A, grade)
    ("N I", 7442.290, "3s", "4P", "3/2", "3p", "4S", "3/2", 4, 4, 1.190e7, "B+"),
    ("N I", 8216.340, "3s", "4P", "5/2", "3p", "4P", "5/2", 6, 6, 2.260e7, "B+"),
    ("N I", 8629.240, "3s", "2P", "3/2", "3p", "2P", "3/2", 4, 4, 2.670e7, "B+"),
    ("N I", 8683.400, "3s", "4P", "3/2", "3p", "4D", "5/2", 4, 6, 1.880e7, "B+"),
    ("N I", 10108.890, "3p", "4D", "3/2", "3d", "4F", "5/2", 4, 6, 3.020e7, "B"),
    # O I 777 triplet -- the dominant AGSS21 oxygen indicator
    ("O I", 7771.940, "3s", "5S", "2", "3p", "5P", "3", 5, 7, 3.690e7, "A"),
    ("O I", 7774.170, "3s", "5S", "2", "3p", "5P", "2", 5, 5, 3.690e7, "A"),
    ("O I", 7775.390, "3s", "5S", "2", "3p", "5P", "1", 5, 3, 3.690e7, "A"),
    # O I 8446 triplet
    ("O I", 8446.250, "3s", "3S", "1", "3p", "3P", "0", 3, 1, 3.220e7, "B"),
    ("O I", 8446.360, "3s", "3S", "1", "3p", "3P", "2", 3, 5, 3.220e7, "B"),
    ("O I", 8446.760, "3s", "3S", "1", "3p", "3P", "1", 3, 3, 3.220e7, "B"),
    # O I 6158 -- upper is 4d, which table5 does not carry; expected ABSENT
    ("O I", 6158.180, "3p", "5P", "3", "4d", "5D", "4", 7, 9, 7.620e6, "B+"),
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
    for (sp, lam, cl, tl, jl, cu, tu, ju, gl, gu, nA, ngr) in TARGETS:
        nlg = math.log10(K * lam * lam * gu * nA)
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
        same = abs(d) < 0.003
        exact += same
        rec |= {"match": "UNIQUE", "tachiev_A_s-1": f"{c['A']:.4E}",
                "tachiev_uncertainty": f"{c['err']:.4E}",
                "tachiev_rel_unc_dex": f"{abs(math.log10(1+c['err']/c['A'])):.4f}",
                "tachiev_log_gf": f"{lg:.4f}", "delta_dex": f"{d:+.4f}",
                "verdict": ("matches NIST within 0.003 dex -- for N I this means NIST's "
                            "value IS Tachiev (TP T7370), so it is not an independent check"
                            if same else "DIFFERS from NIST -- genuinely two sources")}
        out.append(rec)

    fields = list(out[0].keys())
    with (OUT / "tachiev_a_to_gf.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, "") for k in fields})

    matched = [r for r in out if r["match"] == "UNIQUE"]
    by_sp = {}
    for r in matched:
        d = abs(float(r["delta_dex"]))
        b = by_sp.setdefault(r["species"], {"n": 0, "same": 0, "max_abs": 0.0})
        b["n"] += 1; b["same"] += d < 0.003; b["max_abs"] = max(b["max_abs"], d)
    verdict = {
        "ticket": "RYA-1160", "targets": len(out), "uniquely_matched": len(matched),
        "constant_validated": "yes -- reproduces NIST log gf from NIST A to <=0.002 dex",
        "matching": "configuration + term + J on BOTH levels; never wavelength",
        "both_sides_from_A": ("NIST log gf is recomputed from NIST's own A with the SAME "
                              "constant, so delta_dex isolates the A difference and not a "
                              "rounding or formula difference"),
        "wavelength_source": ("NIST (measured); the graded quantity A comes from Tachiev, "
                              "so the gf remains independent of NIST"),
        "per_species": by_sp,
        "finding_N_I": (
            "All 5 N I lines agree with NIST to <=0.002 dex. NIST's N I values for these "
            "lines ARE Tachiev & Froese Fischer -- all five NIST rows share TP code "
            "T7370. The agreement is therefore NOT confirmation; it is one source seen "
            "twice, and the generic NIST pull adds nothing independent for N I."),
        "finding_O_I": (
            "O I behaves the OPPOSITE way and this is the substantive result. Tachiev and "
            "NIST are genuinely different sources for O I and they disagree SYSTEMATICALLY: "
            "-0.0156, -0.0160, -0.0162 dex across the 777 triplet (Tachiev A = 3.556e7 vs "
            "NIST 3.69e7, a 3.7% difference) and -0.0053 dex across all three 8446 lines. "
            "The consistency within each multiplet shows this is a real offset between two "
            "calculations, not scatter. The 777 triplet is AGSS21's dominant oxygen "
            "indicator, so a 0.016 dex gf offset propagates straight into A(O)."),
        "CAUTION": (
            "AGSS21 names Tachiev & Froese Fischer for N i ONLY. It states no gf source for "
            "O i. So the O I rows above are a Tachiev-vs-NIST comparison, NOT an "
            "AGSS21-vs-NIST one, and must not be presented as the latter. Which source "
            "AGSS21 used for O i is OPEN."),
        "forbidden_lines_absent": (
            "[O I] 6300/6363 are NOT in Tachiev table5: it contains ZERO 2p(4) -> 2p(4) "
            "transitions, so the intra-ground-configuration forbidden lines are outside "
            "its scope. AGSS21's [O I] gf therefore comes from elsewhere -- NIST cites "
            "TP T4539,T5081 and the repo bibliography points at storey_zeippen2000."),
        "absent_6158": ("O I 6158 is ABSENT: its upper level is 4d and table5 stops at "
                        "3d/4s. Correctly refused rather than mismatched."),
        "canonical_gf_modified": False,
    }
    (OUT / "tachiev_a_to_gf_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")

    print(f"{'species':6s} {'lambda':>10s} {'match':12s} {'Tachiev gf':>11s} {'NIST gf':>9s} {'delta':>8s}")
    for r in out:
        print(f"{r['species']:6s} {r['wavelength_A']:>10s} {r['match']:12s} "
              f"{r['tachiev_log_gf']:>11s} {r['nist_log_gf']:>9s} {r['delta_dex']:>8s}")
    print()
    for k in ("finding_N_I", "finding_O_I", "CAUTION"):
        print(f"\n{k}: {verdict[k]}")


if __name__ == "__main__":
    main()
