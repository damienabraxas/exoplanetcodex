#!/usr/bin/env python3
"""
RYA-853 scope 2 — referee the Fe I lab pool against the source papers themselves
================================================================================
The NIST pass covered ~60 rows of two hand-maintained extracts and found 70% of the grades
wrong. It did NOT touch the pool that matters: the lines graded on PRIMARY LAB sources
(Ruffoni 2014 / Den Hartog 2014 / Belmonte 2017), which is what RYA-850 promotes to the
reported value. Those grades come from the lab papers, not from NIST, so a NIST cross-check
could never have cleared them — and after a 70% failure rate next door there is no basis
for assuming they are exempt.

The papers are on the local Mac (`Reference documents/`), so this reads them off disk. No
VizieR, and no concluding a line is uncovered from a failed search (RYA-833).

    stu780.pdf                      Ruffoni et al. 2014, MNRAS 441, 3127   doi 10.1093/mnras/stu780
    Den_Hartog_2014_ApJS_215_23.pdf Den Hartog et al. 2014, ApJS 215, 23   doi 10.1088/0067-0049/215/2/23
    Belmonte_2017_ApJ_848_125.pdf   Belmonte et al. 2017, ApJ 848, 125     doi 10.3847/1538-4357/aa8cd3

💡 Those three DOIs are exactly what `pipeline/gf_grades.py` already stores — verified
against the papers' own title pages. Our citations are right.

⚠️ WHAT THE PDFs ACTUALLY CONTAIN, WHICH IS NOT EVERYTHING. Ruffoni's full line list is
Table 3, and its caption says "Only the first upper level is shown here. The full [table is
available online]" — so the complete 142-line set is supplementary, not in the paper. What
IS in the PDFs is:

    Ruffoni    Table 5   solar-synthesis subset, lambda_air / log gf / +-
    Den Hartog Table 5   solar-synthesis subset, lambda_air / log gf +- unc
    Belmonte   Table 4   full BF table, WAVELENGTH IN NANOMETRES / log gf +- unc

So this referees the OVERLAP, and reports coverage honestly rather than implying the whole
pool was cleared. A partial referee that says so is worth more than a complete one that
quietly interpolated.

⚠️ UNIT AND MEDIUM TRAPS. Belmonte's table is in NANOMETRES (358.1193 nm = 3581.193 A), and
the two Wisconsin tables are in AIR Angstroms. Mixing them silently would shift every
Belmonte line by a factor of ten and match nothing — which would read as "Belmonte does not
cover our lines" (RYA-833 again). The conversion is explicit and the match rate is reported
so a wrong convention shows up as zero coverage rather than as a clean pass.

Usage:
    python3 scripts/rya853_fe1_labpool_referee.py
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import codex_root  # noqa: E402

LAB_CSV = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"
OUT = ROOT / "data" / "results" / "rya853"

#: The papers live beside the repo, not in it — they are copyrighted PDFs.
PAPERS_DIR = Path(codex_root("local_data")).parent / "Reference documents"

PAPERS = {
    "Ruffoni2014": ("stu780.pdf", "Ruffoni et al. 2014, MNRAS 441, 3127"),
    "DenHartog2014": ("Den_Hartog_2014_ApJS_215_23.pdf",
                      "Den Hartog et al. 2014, ApJS 215, 23"),
    "Belmonte2017": ("Belmonte_2017_ApJ_848_125.pdf",
                     "Belmonte et al. 2017, ApJ 848, 125"),
}

#: Matching tolerance. These are FTS wavelengths quoted to 3-4 decimals; anything looser
#: risks a neighbouring line, and the lab table carries no EP for every source.
WAVE_TOL_A = 0.02

#: A stored value further than this from the paper's own number is a transcription defect.
VALUE_TOL_DEX = 0.02


def pdf_text(path: Path) -> str:
    """`pdftotext -layout` keeps columns aligned, which is what makes these tables
    parseable at all. pypdf's extract_text collapses them into unusable prose."""
    try:
        return subprocess.run(["pdftotext", "-layout", str(path), "-"],
                              capture_output=True, text=True, timeout=120).stdout
    except FileNotFoundError:
        pass
    try:
        from pypdf import PdfReader
        return "\n".join((p.extract_text() or "") for p in PdfReader(str(path)).pages)
    except Exception as e:
        raise SystemExit(f"cannot read {path.name}: {e}")


# 🔴 THE MINUS SIGN IS NOT AN ASCII HYPHEN. These journals typeset negatives with U+2212
# MINUS SIGN, and an ASCII `-?` silently fails to match it — the capture then starts at the
# digit and every negative log gf comes back POSITIVE. That turned 95 of 99 lines into
# "mismatches" of exactly twice the value (ours -0.31 vs "paper" +0.31) on the first run.
# A defect rate that high with a suspiciously regular signature is a parser bug, not a data
# finding; en/em dashes are accepted too because typesetters use those as well.
_MINUS = r"[-\u2212\u2013\u2014]"

# λair(Å)  ... log(gf)  ±unc     — Ruffoni Table 5 and Den Hartog Table 5
_ANGSTROM_ROW = re.compile(
    rf"^\s*(\d{{4}}\.\d{{3,4}})\s+.*?({_MINUS}?\d\.\d{{2}})\s*(?:±\s*)?(\d\.\d{{2}})\b")
# λ(nm)  ... log(gf)±unc         — Belmonte Table 4
_NM_ROW = re.compile(
    rf"^\s*(\d{{3}}\.\d{{3,4}})\s+.*?({_MINUS}?\d\.\d{{2}})\s*±\s*(\d\.\d{{2,3}})")


def _num(tok: str) -> float:
    """Normalise the unicode minus before float() — float('\u22120.31') raises."""
    return float(tok.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-"))


def parse_paper(key: str, text: str) -> pd.DataFrame:
    rows = []
    nm = key == "Belmonte2017"
    pat = _NM_ROW if nm else _ANGSTROM_ROW
    for line in text.splitlines():
        m = pat.match(line)
        if not m:
            continue
        w = _num(m.group(1))
        if nm:
            w *= 10.0                      # ⚠️ nm -> Angstrom
        if not (2000.0 <= w <= 13000.0):
            continue
        rows.append({"paper_wave_A": w, "paper_loggf": _num(m.group(2)),
                     "paper_unc_dex": _num(m.group(3)), "source": key})
    cols = ["paper_wave_A", "paper_loggf", "paper_unc_dex", "source"]
    d = pd.DataFrame(rows, columns=cols)
    # a paper can list the same line more than once across tables; keep the first
    return d.drop_duplicates(subset=["paper_wave_A"]).reset_index(drop=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not LAB_CSV.exists():
        raise SystemExit(f"lab pool absent: {LAB_CSV}")
    lab = pd.read_csv(LAB_CSV)
    print(f"[pool] {len(lab)} Fe I lab-graded lines: "
          + ", ".join(f"{k}={v}" for k, v in lab.source.value_counts().items()))

    if not PAPERS_DIR.is_dir():
        raise SystemExit(f"reference papers not found at {PAPERS_DIR}")

    parsed = {}
    for key, (fname, cite) in PAPERS.items():
        p = PAPERS_DIR / fname
        if not p.exists():
            print(f"  [{key}] MISSING: {fname}")
            continue
        d = parse_paper(key, pdf_text(p))
        parsed[key] = d
        flag = ""
        if not len(d):
            # 0 rows is a PARSE failure until shown otherwise, never "the paper lacks
            # these lines" (RYA-833). pypdf collapses the columns; pdftotext -layout
            # keeps them, and without it these tables are unreadable.
            flag = "   ⚠️ ZERO ROWS — parse failure, NOT absence; needs pdftotext -layout"
        print(f"  [{key}] {len(d):4d} per-line rows parsed from {fname}   ({cite}){flag}")

    rows = []
    for key, d in parsed.items():
        ours = lab[lab.source == key]
        for _, o in ours.iterrows():
            rec = {"source": key, "wavelength_air_A": float(o.wavelength_air_A),
                   "our_loggf": float(o.loggf), "our_unc_dex": float(o.e_loggf_dex)}
            m = d[np.abs(d.paper_wave_A - o.wavelength_air_A) <= WAVE_TOL_A]
            rec["n_paper_candidates"] = int(len(m))
            if len(m) == 1:
                rec["paper_loggf"] = float(m.iloc[0].paper_loggf)
                rec["paper_unc_dex"] = float(m.iloc[0].paper_unc_dex)
            rows.append(rec)
    t = pd.DataFrame(rows)
    t["d_loggf"] = t.get("paper_loggf") - t.our_loggf
    t["d_unc"] = t.get("paper_unc_dex") - t.our_unc_dex

    matched = t[t.get("paper_loggf").notna()] if "paper_loggf" in t else t.iloc[0:0]
    print(f"\n=== coverage ===")
    for key in parsed:
        sub = t[t.source == key]
        got = int(sub.get("paper_loggf").notna().sum()) if "paper_loggf" in sub else 0
        print(f"  {key:<15} {got:4d} of {len(sub):4d} pool lines refereed "
              f"({100*got/max(len(sub),1):5.1f}%)")
    print(f"  {'TOTAL':<15} {len(matched):4d} of {len(t):4d} "
          f"({100*len(matched)/max(len(t),1):5.1f}%)")
    print(f"  ⚠️ the remainder are in the papers' ONLINE-ONLY tables, not the PDFs — "
          f"absent from\n     this referee, NOT shown to be wrong (RYA-833)")

    bad_v = matched[matched.d_loggf.abs() > VALUE_TOL_DEX]
    bad_u = matched[matched.d_unc.abs() > 0.005]
    print(f"\n=== stored vs the paper's own table ===")
    print(f"  log gf mismatches (> {VALUE_TOL_DEX} dex): {len(bad_v)} of {len(matched)}")
    print(f"  cited-sigma mismatches:                    {len(bad_u)} of {len(matched)}")
    if len(bad_v):
        print(f"\n  🔴 value mismatches")
        for _, r in bad_v.sort_values("wavelength_air_A").iterrows():
            print(f"    {r.source:<14}{r.wavelength_air_A:10.3f}  ours {r.our_loggf:+.3f}"
                  f"  paper {r.paper_loggf:+.3f}  delta {r.d_loggf:+.3f}")
    if len(bad_u):
        print(f"\n  🔴 cited-sigma mismatches")
        for _, r in bad_u.sort_values("wavelength_air_A").iterrows():
            print(f"    {r.source:<14}{r.wavelength_air_A:10.3f}  ours {r.our_unc_dex:.3f}"
                  f"  paper {r.paper_unc_dex:.3f}")

    verdict = ("CLEAN over the refereed subset" if not len(bad_v) and not len(bad_u)
               else "DEFECTS FOUND")
    print(f"\n=== VERDICT: the Fe I lab pool is {verdict} ===")
    if not len(bad_v) and not len(bad_u):
        print(f"  Every one of the {len(matched)} refereed lines reproduces its source "
              f"paper's log gf\n  AND its cited sigma. That is the opposite of the NIST "
              f"extract result (70% wrong),\n  and it is the pool RYA-850 actually "
              f"promotes.")
    print(f"  ⚠️ Refereed subset only. The online-only remainder is UNVERIFIED, not clean.")

    t.to_csv(OUT / "rya853_fe1_labpool_referee.csv", index=False)
    (OUT / "rya853_fe1_labpool_referee.json").write_text(json.dumps({
        "ticket": "RYA-853 scope 2 (Fe I lab pool)",
        "pool_size": int(len(lab)),
        "by_source": {k: int(v) for k, v in lab.source.value_counts().items()},
        "papers": {k: {"file": v[0], "cite": v[1],
                       "rows_parsed": int(len(parsed.get(k, [])))}
                   for k, v in PAPERS.items()},
        "n_refereed": int(len(matched)), "n_pool_rows": int(len(t)),
        "coverage_frac": float(len(matched) / max(len(t), 1)),
        "n_loggf_mismatch": int(len(bad_v)),
        "n_sigma_mismatch": int(len(bad_u)),
        "loggf_mismatches": [
            {"source": r.source, "wavelength_air_A": float(r.wavelength_air_A),
             "ours": float(r.our_loggf), "paper": float(r.paper_loggf),
             "delta": float(r.d_loggf)} for _, r in bad_v.iterrows()],
        "verdict": verdict,
        "caveat": ("refereed against the tables PRESENT IN THE PDFs. Ruffoni's full Table 3 "
                   "is online-only ('Only the first upper level is shown here'), so the "
                   "unrefereed remainder is UNVERIFIED rather than clean (RYA-833)."),
        "traps": [
            "Belmonte's table is in NANOMETRES while the Wisconsin tables are in air "
            "Angstroms; mixing them matches nothing and reads as 'Belmonte does not cover "
            "our lines'",
            "pdftotext -layout is required — pypdf's extract_text collapses the columns "
            "and the tables become unparseable",
        ],
    }, indent=2))
    print(f"\n[out] {OUT}/rya853_fe1_labpool_referee.csv")
    print(f"[out] {OUT}/rya853_fe1_labpool_referee.json")


if __name__ == "__main__":
    main()
