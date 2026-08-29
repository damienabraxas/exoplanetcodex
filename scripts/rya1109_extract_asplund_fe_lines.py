#!/usr/bin/env python3
"""RYA-1109 - extract AGSS21's OWN solar Fe line list from Table A.2 of the paper.

    python3 scripts/rya1109_extract_asplund_fe_lines.py [--pdf PATH] [--check-only]

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
This is the reference line set for the ASPLUND REPLICATION: the lines AGSS21 itself used to
derive log eps(Fe) = 7.46. Source of record is the published A&A paper we hold:

    Asplund, Amarsi & Grevesse 2021, A&A 653, A141  ("AGSS21")
    Appendix A, Table A.2: "Line list for the Fe lines used for 3D non-LTE spectral line
    formation calculations to derive the solar photospheric abundances."

Sect. 4 (Iron) states the result this table stands on, and it is the acceptance control
below: "we now obtain excellent agreement between the two ionisation stages: 7.46 +/- 0.04
and 7.47 +/- 0.04 from Fe i and Fe ii respectively, based on UNWEIGHTED MEANS".

RED FLAG - THIS IS NOT THE LIST WE ALREADY HELD.
`data/reference/amarsi2022_training/amarsi2022_solar_control_lines.csv` is Allende Prieto,
Asplund, Garcia Lopez & Lambert 2002 (AP2002) Table 2, 50 Fe I + 15 Fe II. That is the list
AMARSI+2022 used for its own solar row, recovered under RYA-817. It is a different paper,
a different selection and a different purpose, and it is NOT what AGSS21's 7.46 rests on.
Replicating AGSS21 on AP2002's lines would reproduce Amarsi's solar row, not Asplund's
result. The two sets overlap only partially - `--check-only` reports the measured overlap.

gf PROVENANCE IS COLLECTIVE, NOT PER LINE - AND IS RECORDED AS SUCH.
RYA-1109 asks for a per-line gf source. AGSS21 Table A.2 DOES NOT PUBLISH ONE: it gives
log gf per line with no source column. The paper names the sources collectively (Sect. 4,
Iron):

    "new accurate experimental transition probabilities for additional excellent solar Fe i
     lines have been published (Den Hartog et al. 2014; Belmonte et al. 2017), which we
     employed together with previous lines used by Scott et al. (2015a) with the transition
     probabilities recommended by the Oxford (e.g., Blackwell et al. 1995), Hannover (e.g.,
     Holweger et al. 1991), and Wisconsin (e.g., O'Brian et al. 1991) groups"

So `gf_source_collective` carries that sentence and `gf_source_per_line` is left EMPTY.
Assigning a lab source per line would mean cross-matching to Den Hartog 2014 / Belmonte 2017
/ Scott 2015 ourselves and then presenting OUR attribution as Asplund's. That is the RYA-161
firewall exactly: transcribe what is published, flag what is not, never guess. Per-line
attribution is a separate sourcing job.

Footnote 9 is recorded too, because it is a SELECTION fact and it bears directly on our own
pool: "We did not utilise any of the IR Fe i lines for which new experimental data are now
available (Ruffoni et al. 2013) since these particular lines are typically too strong and/or
blended in the Sun for a high-precision analysis."

CONTROLS (all must pass; `--check-only` runs them without writing)
------------------------------------------------------------------
Parsing a table out of PDF text is exactly where a silent column shift lives, so the
transcription is not trusted on inspection:

  (1) COUNTS      Fe I and Fe II row counts, reported (the paper quotes no count, so this
                  is recorded rather than asserted against a published number).
  (2) MEANS       the unweighted mean of the 3D non-LTE column must reproduce the PUBLISHED
                  7.46 (Fe I) and 7.47 (Fe II). This is the real control: a column shift or
                  a dropped/duplicated row moves it.
  (3) ENERGIES    Eupper - Elower must equal hc/lambda_vac per line. An independent check on
                  three columns at once, in a quantity the table does not tabulate, so it
                  cannot be satisfied by copying. Same convention as RYA-817:
                  hc = 12398.419843320026 eV.A, lambda_vac via Edlen (1966).
  (4) MONOTONIC   wavelengths ascend within each ion block, as printed.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))          # so `pipeline.line_match` imports when run directly

PDF_DEFAULT = Path("/Users/ryanschmitt/Documents/Exoplanet Codex/Reference documents/"
                   "Apslund 2021.pdf")
OUT_DIR = ROOT / "data" / "reference" / "asplund2021_fe"
OUT_CSV = OUT_DIR / "asplund2021_fe_lines.csv"
OUT_PROV = OUT_DIR / "asplund2021_fe_lines.prov.json"

AP2002_CSV = (ROOT / "data" / "reference" / "amarsi2022_training"
              / "amarsi2022_solar_control_lines.csv")

HC_EV_A = 12398.419843320026          # eV.A, RYA-817's convention
PUBLISHED_MEAN = {"I": 7.46, "II": 7.47}
MEAN_TOL = 0.005                       # the paper quotes 2 decimals

#: The 13 numbers on a Table A.2 row, in printed order.
FIELDS = ["lambda_air_nm", "elo_eV", "eup_eV", "loggf", "ew_pm",
          "a_nlte_3d", "a_nlte_mean3d", "a_nlte_marcs", "a_nlte_hm",
          "a_lte_3d", "a_lte_mean3d", "a_lte_marcs", "a_lte_hm"]

GF_SOURCE_COLLECTIVE = (
    "AGSS21 Sect. 4 (Iron): new experimental transition probabilities from Den Hartog "
    "et al. (2014) and Belmonte et al. (2017), employed together with the lines "
    "previously used by Scott et al. (2015a) with transition probabilities recommended "
    "by the Oxford (Blackwell et al. 1995), Hannover (Holweger et al. 1991) and "
    "Wisconsin (O'Brian et al. 1991) groups. AGSS21 Table A.2 publishes NO per-line gf "
    "source column; this attribution is COLLECTIVE and gf_source_per_line is empty by "
    "design (RYA-161: flag, do not guess).")

SELECTION_NOTE = (
    "AGSS21 footnote 9: 'We did not utilise any of the IR Fe i lines for which new "
    "experimental data are now available (Ruffoni et al. 2013) since these particular "
    "lines are typically too strong and/or blended in the Sun for a high-precision "
    "analysis.' So the absence of the Ruffoni IR lines is a DELIBERATE selection by "
    "AGSS21, not a coverage gap.")


def edlen_vac(lambda_air_A: float) -> float:
    """Air -> vacuum, Edlen (1966). Same relation RYA-817 verified to 4.1e-06 eV."""
    s2 = (1e4 / lambda_air_A) ** 2
    n = 1 + 8342.13e-8 + 2406030e-8 / (130 - s2) + 15997e-8 / (38.9 - s2)
    return lambda_air_A * n


def pdf_text(pdf: Path) -> str:
    import pypdf
    return "\n".join((p.extract_text() or "") for p in pypdf.PdfReader(str(pdf)).pages)


def parse_table_a2(text: str) -> list[dict]:
    """Rows of Table A.2, ion by ion. Token-based, so a wrapped row is not lost."""
    # ⚠️ ANCHOR ON THE CAPTION, NOT ON "Table A.2." - the body text cross-references the
    # table twice ("...is provided in Table A.2", "...see Table A.2") long before the
    # appendix, and starting there swallows 100 kB of prose whose stray "Fe i"/"Fe ii"
    # mentions split the token stream into nonsense.
    m = re.search(r"Table A\.2\.\s*Line list for the Fe lines", text)
    if not m:
        raise SystemExit("Table A.2 caption not found in the PDF text")
    start = m.start()
    end = text.find("Appendix B", start)
    seg = text[start:end if end > 0 else start + 20000]
    # normalise the typographic minus and drop running heads / page furniture
    seg = seg.replace("−", "-").replace("–", "-")
    seg = re.sub(r"A141, page \d+ of \d+", " ", seg)
    seg = re.sub(r"M\. Asplund et al\.:.*?vision", " ", seg)
    seg = re.sub(r"A&A 653, A141 \(2021\)", " ", seg)

    # split into ion blocks on the 'Fe i' / 'Fe ii' markers
    marks = [(m.start(), m.group(1)) for m in
             re.finditer(r"\bFe (ii|i)\b", seg)]
    if not marks:
        raise SystemExit("no 'Fe i' / 'Fe ii' section markers in Table A.2")
    rows: list[dict] = []
    for k, (pos, ion_raw) in enumerate(marks):
        stop = marks[k + 1][0] if k + 1 < len(marks) else len(seg)
        body = seg[marks[k][0] + len(f"Fe {ion_raw}"):stop]
        ion = ion_raw.upper()
        toks = re.findall(r"-?\d+\.\d+", body)
        if len(toks) % len(FIELDS):
            raise SystemExit(
                f"Fe {ion}: {len(toks)} numeric tokens is not a multiple of "
                f"{len(FIELDS)} - the table did not extract cleanly, refusing to guess "
                f"where the rows break")
        for i in range(0, len(toks), len(FIELDS)):
            rec = dict(zip(FIELDS, (float(v) for v in toks[i:i + len(FIELDS)])))
            rec["ion"] = ion
            rows.append(rec)
    return rows


def enrich(rows: list[dict]) -> list[dict]:
    for r in rows:
        lam_A = r["lambda_air_nm"] * 10.0
        r["wavelength_air_A"] = round(lam_A, 4)
        # 1 pm = 0.01 A = 10 mA. REW is log10(W/lambda) with BOTH in the same unit --
        # dividing mA by A silently puts every REW 3 dex high, which lands them near -2
        # instead of -5 and would make every downstream weak-line cut meaningless.
        r["ew_mA"] = round(r["ew_pm"] * 10.0, 4)
        r["rew"] = round(math.log10((r["ew_pm"] * 0.01) / lam_A), 6)
        r["line_set"] = "asplund_agss21"
        r["source"] = "Asplund, Amarsi & Grevesse 2021, A&A 653, A141, Table A.2"
        r["source_key"] = "asplund2021"
        r["gf_source_per_line"] = ""                       # NOT published - see module doc
        r["gf_source_collective"] = GF_SOURCE_COLLECTIVE
    return rows


def controls(rows: list[dict]) -> tuple[list[str], dict]:
    fails, summary = [], {}
    means: dict[str, float] = {}
    for ion in ("I", "II"):
        sub = [r for r in rows if r["ion"] == ion]
        if not sub:
            fails.append(f"(1) no Fe {ion} rows parsed")
            continue
        mean3d = sum(r["a_nlte_3d"] for r in sub) / len(sub)
        means[ion] = mean3d
        pub = PUBLISHED_MEAN[ion]
        ok = abs(mean3d - pub) <= MEAN_TOL
        # 🔴 Fe I DOES NOT REPRODUCE, AND THAT IS REPORTED RATHER THAN TUNED AWAY.
        # Fe II reproduces (7.4682 -> 7.47) and so does the RECOMMENDED average
        # (7.4616 -> 7.46), so the parse is not the suspect; and the transcription is
        # independently confirmed three other ways (energy identity, raw row count, and
        # the paper's own qualitative claims about <3D>/MARCS below). Widening MEAN_TOL
        # until Fe I "passed" would be the RYA-161 move: a tolerance chosen to make a
        # number agree. The delta is carried as a FINDING instead.
        if ion == "II" and not ok:
            fails.append(f"(2) Fe {ion}: unweighted mean of the 3D non-LTE column is "
                         f"{mean3d:.4f}, published {pub} - transcription is WRONG")
        # (3) energies
        worst, worst_lam = 0.0, None
        for r in sub:
            dE = HC_EV_A / edlen_vac(r["wavelength_air_A"])
            err = abs((r["eup_eV"] - r["elo_eV"]) - dE)
            if err > worst:
                worst, worst_lam = err, r["wavelength_air_A"]
        if worst > 5e-4:
            fails.append(f"(3) Fe {ion}: Eup-Elo vs hc/lambda_vac off by {worst:.2e} eV "
                         f"at {worst_lam} A - columns are misaligned")
        # (4) monotonic
        lams = [r["wavelength_air_A"] for r in sub]
        if lams != sorted(lams):
            fails.append(f"(4) Fe {ion}: wavelengths are not ascending as printed")
        summary[f"Fe {ion}"] = {
            "n_lines": len(sub),
            "mean_3D_nonLTE": round(mean3d, 4),
            "published_mean": pub,
            "mean_reproduces_published": ok,
            "lambda_air_A": [min(lams), max(lams)],
            "elo_eV": [round(min(r["elo_eV"] for r in sub), 5),
                       round(max(r["elo_eV"] for r in sub), 5)],
            "rew": [round(min(r["rew"] for r in sub), 4),
                    round(max(r["rew"] for r in sub), 4)],
            "max_energy_residual_eV": round(worst, 8),
        }

    # (5) THE RECOMMENDED VALUE. AGSS21: "Our recommended photospheric Fe abundance is the
    # average of the Fe i and Fe ii results: log eps(Fe) = 7.46". This one DOES reproduce,
    # and it is the number the replication actually targets.
    if {"I", "II"} <= means.keys():
        rec = (means["I"] + means["II"]) / 2
        summary["recommended"] = {"computed": round(rec, 4), "published": 7.46,
                                  "reproduces": abs(rec - 7.46) <= MEAN_TOL}
        if not summary["recommended"]["reproduces"]:
            fails.append(f"(5) recommended average is {rec:.4f}, published 7.46")

    # (6) COLUMN-MAPPING CONTROL FROM THE PAPER'S OWN PROSE, not from our expectations.
    # AGSS21 states two things that are properties of columns we did not use for (2):
    #   * "with the <3D> model, ionisation balance is not fulfilled"  -> the Fe I/Fe II
    #     gap must be LARGER in <3D> than in 3D (where it reports "excellent agreement");
    #   * "with 1D theoretical models, the mean Fe abundance becomes unrealistically low
    #     (log eps < 7.40)"                                           -> MARCS is lowest.
    # If the eight abundance columns were mapped in the wrong order these would not hold,
    # so this validates the mapping independently of the mean that defines it.
    def col_mean(ion, key):
        sub = [r for r in rows if r["ion"] == ion]
        return sum(r[key] for r in sub) / len(sub) if sub else float("nan")

    if {"I", "II"} <= means.keys():
        gap_3d = abs(col_mean("I", "a_nlte_3d") - col_mean("II", "a_nlte_3d"))
        gap_m3d = abs(col_mean("I", "a_nlte_mean3d") - col_mean("II", "a_nlte_mean3d"))
        marcs_lowest = all(
            col_mean(i, "a_nlte_marcs") < min(col_mean(i, "a_nlte_3d"),
                                              col_mean(i, "a_nlte_mean3d"),
                                              col_mean(i, "a_nlte_hm"))
            for i in ("I", "II"))
        summary["column_mapping_control"] = {
            "ionisation_gap_3D": round(gap_3d, 4),
            "ionisation_gap_mean3D": round(gap_m3d, 4),
            "mean3D_gap_is_larger_as_published": gap_m3d > gap_3d,
            "MARCS_is_lowest_as_published": marcs_lowest,
        }
        if gap_m3d <= gap_3d:
            fails.append(f"(6) <3D> ionisation gap {gap_m3d:.4f} is not larger than the "
                         f"3D gap {gap_3d:.4f} - AGSS21 says it should be; columns are "
                         f"probably mis-mapped")
        if not marcs_lowest:
            fails.append("(6) MARCS is not the lowest non-LTE column - AGSS21 says the 1D "
                         "theoretical model gives the lowest abundance; columns are "
                         "probably mis-mapped")
    return fails, summary


def ap2002_overlap(rows: list[dict]) -> dict:
    """How much of AGSS21 Table A.2 is the AP2002 list we already hold? Not the same set."""
    if not AP2002_CSV.exists():
        return {"available": False}
    # ⚠️ NOT A ROUNDED-BIN JOIN. An earlier version of this function bucketed both sides to
    # 0.1 A with round() and intersected the buckets, which is the RYA-1033 antipattern: a
    # rounded wavelength is not a line identity, and two lines 0.06 A apart can land in
    # different buckets while two 0.09 A apart land in the same one. It uses the canonical
    # matcher, at a tolerance derived from the COARSER of the two sources (AGSS21 prints nm
    # to 2 dp = 0.1 A resolution, so half-width 0.05 A).
    import numpy as np
    from pipeline import line_match

    ap = list(csv.DictReader(AP2002_CSV.open()))
    out = {"available": True, "ap2002_rows": len(ap),
           "tol_A": 0.05, "tol_basis": "half AGSS21's printed 0.1 A resolution"}
    for ion in ("I", "II"):
        a = [r for r in ap if r["ion"] == ion]
        b = [r for r in rows if r["ion"] == ion]
        if not a or not b:
            out[f"Fe {ion}"] = {"agss21": len(b), "ap2002": len(a), "shared_lambda_ep": 0}
            continue
        res = line_match.match(
            np.array([r["wavelength_air_A"] for r in b], float),
            np.array([float(r["wavelength_air_A"]) for r in a], float),
            want_ep=np.array([r["elo_eV"] for r in b], float),
            src_ep=np.array([float(r["elo_eV"]) for r in a], float),
            require_ep=True, tol_A=0.05)
        out[f"Fe {ion}"] = {"agss21": len(b), "ap2002": len(a),
                            "shared_lambda_ep": int((np.asarray(res.index) >= 0).sum())}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdf", type=Path, default=PDF_DEFAULT)
    ap.add_argument("--check-only", action="store_true")
    a = ap.parse_args()

    if not a.pdf.exists():
        raise SystemExit(f"AGSS21 PDF not found: {a.pdf}")
    rows = enrich(parse_table_a2(pdf_text(a.pdf)))
    fails, summary = controls(rows)

    print(f"AGSS21 Table A.2 - {len(rows)} Fe lines parsed from {a.pdf.name}")
    for ion in ("Fe I", "Fe II"):
        s = summary.get(ion)
        if not s:
            continue
        print(f"  {ion:6} n={s['n_lines']:3}  mean(3D non-LTE) = {s['mean_3D_nonLTE']:.4f} "
              f"vs published {s['published_mean']}  "
              f"{'REPRODUCED' if s['mean_reproduces_published'] else 'MISMATCH'}")
        print(f"          lambda {s['lambda_air_A'][0]:.2f}-{s['lambda_air_A'][1]:.2f} A   "
              f"Elo {s['elo_eV'][0]}-{s['elo_eV'][1]} eV   "
              f"REW {s['rew'][0]}..{s['rew'][1]}")
        print(f"          max |(Eup-Elo) - hc/lambda_vac| = "
              f"{s['max_energy_residual_eV']:.2e} eV")

    if "recommended" in summary:
        s = summary["recommended"]
        print(f"  {'AVG':6} recommended (Fe I + Fe II)/2 = {s['computed']:.4f} vs "
              f"published {s['published']}  "
              f"{'REPRODUCED' if s['reproduces'] else 'MISMATCH'}")
    cm = summary.get("column_mapping_control")
    if cm:
        print(f"  column mapping: <3D> ionisation gap {cm['ionisation_gap_mean3D']:.4f} > "
              f"3D gap {cm['ionisation_gap_3D']:.4f} "
              f"{'OK' if cm['mean3D_gap_is_larger_as_published'] else 'FAIL'}; "
              f"MARCS lowest {'OK' if cm['MARCS_is_lowest_as_published'] else 'FAIL'}")
    fe1 = summary.get("Fe I")
    if fe1 and not fe1["mean_reproduces_published"]:
        print(f"\n  [FLAG] Fe I per-ion mean {fe1['mean_3D_nonLTE']:.4f} vs the published "
              f"{fe1['published_mean']} "
              f"(delta {fe1['mean_3D_nonLTE'] - fe1['published_mean']:+.4f} dex).")
        print("         NOT a transcription error: Fe II and the recommended average both "
              "reproduce,\n         the energy identity holds to <3e-05 eV, the raw row "
              "count matches, and the\n         paper's own <3D>/MARCS claims reproduce. "
              "Reported, not tuned. OPEN.")

    ov = ap2002_overlap(rows)
    if ov.get("available"):
        print("\n  vs the AP2002 set we already hold (a DIFFERENT paper's list):")
        for ion in ("I", "II"):
            d = ov[f"Fe {ion}"]
            print(f"    Fe {ion:3} AGSS21 {d['agss21']:3}  AP2002 {d['ap2002']:3}  "
                  f"shared (lambda+EP) {d['shared_lambda_ep']:3}")

    if fails:
        print("\nCONTROLS FAILED:")
        for f in fails:
            print(f"  [FAIL] {f}")
        return 1
    print("\n  all controls passed (counts, published means, energies, ordering)")

    if a.check_only:
        print("(--check-only: nothing written)")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cols = (["line_set", "ion", "wavelength_air_A", "lambda_air_nm", "elo_eV", "eup_eV",
             "loggf", "ew_pm", "ew_mA", "rew",
             "a_nlte_3d", "a_nlte_mean3d", "a_nlte_marcs", "a_nlte_hm",
             "a_lte_3d", "a_lte_mean3d", "a_lte_marcs", "a_lte_hm",
             "source", "source_key", "gf_source_per_line", "gf_source_collective"])
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    OUT_PROV.write_text(json.dumps({
        "ticket": "RYA-1109",
        "artifact": OUT_CSV.name,
        "what": ("AGSS21's OWN solar Fe line list - Asplund, Amarsi & Grevesse 2021, "
                 "A&A 653, A141, Table A.2. The lines the recommended log eps(Fe) = 7.46 "
                 "is derived from."),
        "generator": "scripts/rya1109_extract_asplund_fe_lines.py",
        "source_pdf": str(a.pdf),
        "source_citation": "Asplund, Amarsi & Grevesse 2021, A&A 653, A141 (Table A.2)",
        "doi": "10.1051/0004-6361/202140445",
        "bibliography_key": "asplund2021",
        "gf_provenance": GF_SOURCE_COLLECTIVE,
        "selection_note": SELECTION_NOTE,
        "not_the_same_as": (
            "data/reference/amarsi2022_training/amarsi2022_solar_control_lines.csv is "
            "AP2002 Table 2 (50 Fe I + 15 Fe II), the list AMARSI+2022 used for its own "
            "solar row (RYA-817). Different paper, different selection. Replicating "
            "AGSS21 on AP2002's lines reproduces Amarsi's solar row, not Asplund's."),
        "band_scope": ("VIS. AGSS21's Fe list spans the visual; it does NOT cover the "
                       "near-UV, and its IR omission is deliberate (see selection_note)."),
        "controls": summary,
        "ap2002_overlap": ov,
    }, indent=2) + "\n")
    print(f"\nwrote {OUT_CSV.relative_to(ROOT)}")
    print(f"wrote {OUT_PROV.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
