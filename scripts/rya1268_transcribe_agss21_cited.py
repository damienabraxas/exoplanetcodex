#!/usr/bin/env python3
"""RYA-1268 -- transcribe the Asplund-cited solar line tables (Scott+2015a, Scott+2015b,
Grevesse+2015, "The elemental composition of the Sun" I-III) from the authors' arXiv LaTeX.

AGSS21 (A&A 653, A141) adopts these three papers' line sets for every element it did not
re-analyse in Table A.1/A.2 (Sect. 2.4-2.6). Parsing the LaTeX table, not the PDF text,
keeps every printed digit: the PDF extraction scrambles the multi-column rows.

Output: data/reference/agss21_cited_sets/agss21_cited_lines.csv (one row per printed line)
and agss21_cited_gf_refs.csv (each paper's numbered gf reference list). Wavelengths are the
printed air nm x 10; EW stays in pm. Nothing is matched to canonical_gf here -- that is
pipeline.line_sets' job once a set is declared in SOURCES.yaml.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "reference" / "agss21_cited_sets"
PAPERS = {
    "scott2015_NatoCa": ("tex/scott2015_NatoCa_1405.0279.tex", "10.1051/0004-6361/201424109",
                         "Scott et al. 2015 (A&A 573, A25) Table 1 -- Na to Ca"),
    "scott2015_FePeak": ("tex/scott2015_FePeak_1405.0287.tex", "10.1051/0004-6361/201424110",
                         "Scott et al. 2015 (A&A 573, A26) Table 1 -- Sc to Ni"),
    "grevesse2015_Heavies": ("tex/grevesse2015_Heavies_1405.0288.tex", "10.1051/0004-6361/201424111",
                             "Grevesse et al. 2015 (A&A 573, A27) Table 1 -- Cu to Th"),
}
ROMAN = {"i": "I", "ii": "II", "iii": "III"}
# macro = <element stem><roman>; the papers spell Sc 'scan', P 'phos', S 'sul', so the stem
# is matched against the known symbols (longest first) -- "tii" is Ti I, "nii" Ni I, "sii" Si I.
STEM = {"scan": "Sc", "phos": "P", "sul": "S"}
SYMBOLS = sorted({"na", "mg", "al", "si", "k", "ca", "ti", "v", "cr", "mn", "fe", "co", "ni",
                  "cu", "zn", "ga", "ge", "rb", "sr", "y", "zr", "nb", "mo", "ru", "rh", "pd",
                  "ag", "cd", "sn", "ba", "w", "os", "ir", "au", "pb", "th", *STEM},
                 key=len, reverse=True)


def species(macro: str) -> str:
    for stem in SYMBOLS:
        ion = macro[len(stem):]
        if macro.startswith(stem) and ion in ROMAN:
            return f"{STEM.get(stem, stem.capitalize())} {ROMAN[ion]}"
    raise ValueError(f"unknown species macro \\{macro}")


def num(s: str):
    s = re.sub(r"\\ph|\\phz|\\pst|\\phantom\{[^}]*\}|\$|\\,|\{|\}", "", s)
    s = s.replace("−", "-").replace("--", "-").strip()
    star = "*" in s
    s = s.replace("*", "").strip()
    try:
        return float(s), star
    except ValueError:
        return None, star


def table_rows(tex: str):
    block = tex.split(r"\begin{mpsupertabular}", 1)[1].split(r"\end{mpsupertabular}", 1)[0]
    block = re.sub(r"\\footnote(?:\[[^]]*\])?\{(?:[^{}]|\{[^{}]*\})*\}", "FOOTNOTE", block)
    sp = None
    for rec in block.split(r"\\"):
        rec = " ".join(rec.split())
        h = re.search(r"\\multicolumn\{\d+\}\{c\} *\{\\([a-z]+)\}", rec)
        if h:
            sp = species(h.group(1))
            continue
        cells = [c.strip() for c in rec.split("&")]
        if len(cells) < 10:
            continue
        lam, _ = num(cells[0])
        if lam is None:
            continue
        yield sp, cells, "FOOTNOTE" in rec


def parse(key: str):
    rel, doi, cite = PAPERS[key]
    tex = (BASE / rel).read_text()
    out = []
    for sp, c, fn in table_rows(tex):
        wide = len(c) >= 18  # Scott II prints the level designations (4 extra columns)
        o = 6 if wide else 2
        get = lambda i: num(c[o + i])[0]
        wt, star = num(c[o + 4])
        dn, _ = num(c[o + 10])
        out.append({
            "paper": key, "species": sp,
            "wavelength_air_A": round(num(c[0])[0] * 10, 4),
            "ep_eV": get(0), "published_loggf": get(1),
            "gf_ref": int(get(2)) if get(2) is not None else "",
            "ew_pm": get(3), "weight": wt, "gf_weight_reduced": star,
            "A_3d_lte": get(5), "A_mean3d_lte": get(6), "A_hm_lte": get(7),
            "A_marcs_lte": get(8), "A_miss_lte": get(9), "dnlte": dn,
            "A_3d_nlte_or_lte": num(c[o + 11])[0] if len(c) > o + 11 else None,
            "blend_note": "footnote (multiplet computed as one line / forbidden)" if fn else "",
            "doi": doi, "citation": cite,
        })
    return out


def refs(key: str):
    tex = (BASE / PAPERS[key][0]).read_text()
    part = tex.split(r"\textbf{References:}", 1)[1].split(r"\begin{mpsupertabular}", 1)[0]
    items = re.findall(r"\\item\s*(.+?)(?=\\item|\\end\{enumerate\})", part, flags=re.S)
    return [{"paper": key, "gf_ref": i + 1, "reference": " ".join(t.split())}
            for i, t in enumerate(items)]


def main() -> None:
    rows, rr = [], []
    for k in PAPERS:
        rows += parse(k)
        rr += refs(k)
    for name, data in (("agss21_cited_lines.csv", rows), ("agss21_cited_gf_refs.csv", rr)):
        with open(BASE / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
    from collections import Counter
    print(Counter((r["paper"], r["species"]) for r in rows))
    print(len(rows), "lines;", len(rr), "gf refs")


if __name__ == "__main__":
    main()
