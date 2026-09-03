#!/usr/bin/env python3
"""
RYA-1160 — ingest the gf sources AGSS21 ITSELF names for CNO.

🔴 WHY THESE TWO AND NOT NIST.  AGSS21 Sect. 4 names its gf source per indicator, and
those are the Asplund-grade gf:

  * C i (14 lines) -> "new g f-values from large-scale atomic structure calculations by
    Li et al. (2021)"  =  Li, Amarsi, Papoulia, Ekman & Jonsson 2021, MNRAS 502, 3780,
    doi 10.1093/mnras/stab214.  Supplementary tables A2-A5 cover C I, II, III and IV.
  * N i -> "The N i oscillator strengths were taken from Tachiev & Froese Fischer (2002)
    and are expected to be reliable at the 0.03 dex level based on rankings from NIST"
    =  A&A 385, 716, doi 10.1051/0004-6361:20011816, CDS J/A+A/385/716.

RYA-946 precedence is lab > NIST-C+ > Kurucz > VALD, and a generic NIST ASD pull sits at
the NIST-C+ rung AND is not independent of us (30.6% of our CNO rows are already
NIST-derived -- see rya1160_circularity_check.py).  These two are what the reference
analysis actually used, so they are the correct rung for an Asplund-grade replication.

⚠️ BOTH ARE THEORY, and that is not a defect -- it is what the field uses for C/N/O.  Li
is MCDHF/RCI; Tachiev & Froese Fischer is Breit-Pauli MCHF.  Neither may be tiered LAB.
What makes them gradeable is that BOTH SHIP A PER-LINE UNCERTAINTY, which NIST letter
grades only bucket:
  * Li table A2-A5 carry `dT`, the length-velocity gauge difference.
  * Tachiev table4/table5 carry an explicit uncertainty column PER ION.

🔴 KNOWN BLOCKER, recorded rather than papered over.  The Tachiev CDS tables give Aki +
uncertainty + full level identity (config/term/J on both sides) but NO wavelength and NO
gf.  Converting A -> gf needs level energies, which are NOT in the CDS holding (it ships
only ReadMe, table4, table5 and their .tex twins).  So N I/O I land here as A-values
pending a level join; they are NOT gf yet and must not be presented as such.

Acquisition + normalisation only.  canonical_gf is not touched.
"""
from __future__ import annotations

import csv, gzip, hashlib, json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/reference/cno_atomic_primary"
OUT = ROOT / "data/reference/cno_atomic_primary/derived"

#: Physically meaningless rows exist in the Li tables by construction: transitions between
#: near-degenerate levels give |dE| ~ 0 and therefore absurd or NEGATIVE lambda. They are
#: real table rows, not parse errors, and are KEPT with a flag rather than silently cut --
#: dropping them would misreport the source's own content.
LAMBDA_MIN_A, LAMBDA_MAX_A = 0.0, 1.0e6


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def clean_tex(s: str) -> str:
    s = re.sub(r"[$~]", "", s)
    return " ".join(s.replace("\\", "").split())


def parse_li() -> list[dict]:
    rows = []
    base = SRC / "c_li2021/CI-IV_supplementary_tables"
    for tab, species in (("tableA2_C_I.tex", "C I"), ("tableA3_C_II.tex", "C II"),
                         ("tableA4_C_III.tex", "C III"), ("tableA5_C_IV.tex", "C IV")):
        path = base / tab
        for n, line in enumerate(path.open(errors="replace"), 1):
            p = [x.strip() for x in line.split("&")]
            if len(p) < 8:
                continue
            try:
                dE, lam, S, gf, A, dT = (float(p[2]), float(p[3]), float(p[4]),
                                         float(p[5]), float(p[6]), float(p[7]))
            except ValueError:
                continue
            usable = LAMBDA_MIN_A < lam < LAMBDA_MAX_A
            rows.append({
                "source": "Li2021", "species": species,
                "upper_level": clean_tex(p[0]), "lower_level": clean_tex(p[1]),
                "delta_E_cm-1": f"{dE:.1f}", "wavelength_A": f"{lam:.4f}",
                "line_strength_S": f"{S:.4E}", "gf": f"{gf:.4E}",
                "log_gf": f"{__import__('math').log10(gf):.4f}" if gf > 0 else "",
                "A_s-1": f"{A:.4E}", "dT_gauge_difference": f"{dT:.4f}",
                "quantity": "gf", "wavelength_usable": "YES" if usable else "NO",
                "note": "" if usable else
                        "near-degenerate levels: |dE|~0 gives an unphysical lambda; "
                        "row kept to preserve the source's content, excluded from joins",
                "source_file": str(path.relative_to(ROOT)), "source_line": n,
            })
    return rows


def parse_tachiev() -> list[dict]:
    rows = []
    base = SRC / "n_tachiev2002"
    spec = {"table4.dat.gz": [("N I", 60, 70, 71, 81), ("O II", 82, 92, 93, 103)],
            "table5.dat.gz": [("O I", 60, 70, 71, 81), ("F II", 82, 92, 93, 103)]}
    for fname, ions in spec.items():
        path = base / fname
        if not path.exists():
            continue
        with gzip.open(path, "rt", errors="replace") as fh:
            for n, raw in enumerate(fh, 1):
                if len(raw) < 60:
                    continue
                cl, tl, jl = raw[0:20].strip(), raw[21:24].strip(), raw[25:29].strip()
                cu, tu, ju = raw[30:50].strip(), raw[51:54].strip(), raw[55:60].strip()
                for species, a0, a1, e0, e1 in ions:
                    try:
                        aki = float(raw[a0:a1]); err = float(raw[e0:e1])
                    except ValueError:
                        continue
                    rows.append({
                        "source": "TachievFroeseFischer2002", "species": species,
                        "upper_level": f"{cu} {tu} J={ju}",
                        "lower_level": f"{cl} {tl} J={jl}",
                        "delta_E_cm-1": "", "wavelength_A": "",
                        "line_strength_S": "", "gf": "", "log_gf": "",
                        "A_s-1": f"{aki:.4E}", "dT_gauge_difference": f"{err:.4E}",
                        "quantity": "A_only",
                        "wavelength_usable": "NO",
                        "note": ("Aki + explicit per-line uncertainty + full level "
                                 "identity, but the CDS holding ships NO level energies, "
                                 "so no wavelength and no gf. A->gf needs a level join "
                                 "before this can be graded."),
                        "source_file": str(path.relative_to(ROOT)), "source_line": n,
                    })
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    li, ta = parse_li(), parse_tachiev()
    fields = ("source", "species", "upper_level", "lower_level", "delta_E_cm-1",
              "wavelength_A", "line_strength_S", "gf", "log_gf", "A_s-1",
              "dT_gauge_difference", "quantity", "wavelength_usable", "note",
              "source_file", "source_line")
    with (OUT / "agss21_cno_atomic_gf.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader(); w.writerows(li + ta)

    def tally(rows):
        out = {}
        for r in rows:
            k = r["species"]
            d = out.setdefault(k, {"rows": 0, "usable_wavelength": 0, "with_gf": 0})
            d["rows"] += 1
            d["usable_wavelength"] += r["wavelength_usable"] == "YES"
            d["with_gf"] += bool(r["gf"])
        return out

    prov = {
        "ticket": "RYA-1160",
        "why": ("These are the gf sources AGSS21 itself names for CNO, not a generic "
                "NIST pull. RYA-946 precedence is lab > NIST-C+ > Kurucz > VALD."),
        "sources": [
            {"id": "Li2021", "citation": "Li, Amarsi, Papoulia, Ekman & Jonsson 2021, "
                                         "MNRAS 502, 3780",
             "doi": "10.1093/mnras/stab214",
             "agss21_role": "C i gf -- AGSS21: 'new g f-values from large-scale atomic "
                            "structure calculations by Li et al. (2021)'",
             "method": "MCDHF/RCI -- THEORY, not laboratory. Never tier LAB.",
             "per_line_uncertainty": "dT, the length-velocity gauge difference",
             "asset": "data/reference/cno_atomic_primary/c_li2021/"
                      "stab214_supplemental_tables.zip",
             "sha256": sha256(SRC / "c_li2021/stab214_supplemental_tables.zip"),
             "coverage": "C I, C II, C III, C IV (tables A2-A5)"},
            {"id": "TachievFroeseFischer2002",
             "citation": "Tachiev & Froese Fischer 2002, A&A 385, 716",
             "doi": "10.1051/0004-6361:20011816", "cds": "J/A+A/385/716",
             "agss21_role": "N i gf -- AGSS21: 'taken from Tachiev & Froese Fischer "
                            "(2002) and are expected to be reliable at the 0.03 dex "
                            "level based on rankings from NIST'",
             "method": "Breit-Pauli MCHF -- THEORY, not laboratory. Never tier LAB.",
             "per_line_uncertainty": "explicit uncertainty column, per ion",
             "assets": {f: sha256(SRC / "n_tachiev2002" / f) for f in
                        ("table4.dat.gz", "table5.dat.gz", "ReadMe",
                         "TachievFroeseFischer-2002-AA385-716.pdf")},
             "coverage": "table4 = N I, O II; table5 = O I, F II",
             "BLOCKER": ("Aki + uncertainty + level identity only. The CDS holding ships "
                         "NO level energies, so no wavelength and no gf. A->gf requires a "
                         "level join. These rows are NOT gf and must not be graded yet.")},
        ],
        "li2021_by_species": tally(li),
        "tachiev_by_species": tally(ta),
        "total_rows": len(li) + len(ta),
        "canonical_gf_modified": False, "adjudicated": False,
    }
    (OUT / "agss21_cno_atomic_gf.prov.json").write_text(json.dumps(prov, indent=2) + "\n")

    print("=== Li 2021 (C I-IV) — gf in hand ===")
    for k, v in tally(li).items():
        print(f"   {k:6s} rows {v['rows']:5d}   with gf {v['with_gf']:5d}   "
              f"usable lambda {v['usable_wavelength']:5d}")
    print("=== Tachiev & Froese Fischer 2002 — A-values, gf BLOCKED on a level join ===")
    for k, v in tally(ta).items():
        print(f"   {k:6s} rows {v['rows']:5d}   with gf {v['with_gf']:5d}")
    print(f"\n   total {len(li)+len(ta)} rows -> "
          f"{(OUT/'agss21_cno_atomic_gf.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
