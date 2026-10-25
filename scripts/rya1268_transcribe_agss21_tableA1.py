#!/usr/bin/env python3
"""RYA-1268 -- transcribe AGSS21 Table A.1 (Na, K, Ca; Mg is RYA-1250's) from the local
aa40445-21.pdf. AGSS21 re-analysed these elements in 3D non-LTE and CHANGED Scott+2015's
selection (K: 404.4/580.2 dropped, NIST gf adopted; Ca: 616.13 dropped, 616.90+616.95 fitted
as one feature), so Table A.1 -- not Scott+2015 -- is the Asplund set for them.

Printed wavelengths are nm to 2 decimals; the 4-decimal wavelength is carried from the
Scott+2015 transcription (agss21_cited_lines.csv) where the same line is there, so the
canonical join does not rest on a 0.01 nm print.
Run: python scripts/rya1268_transcribe_agss21_tableA1.py <path to aa40445-21.pdf>
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "reference" / "agss21_cited_sets"
HDR = {"Nai": "Na I", "Ki": "K I", "Cai": "Ca I", "Caii": "Ca II", "Mgi": "Mg I", "Mgii": "Mg II"}
KEEP = {"Na I", "K I", "Ca I", "Ca II"}
COLS = ["nlte_3d", "nlte_mean3d", "nlte_marcs", "nlte_hm", "lte_3d", "lte_mean3d", "lte_marcs", "lte_hm"]
NOTE_A = "doublet 616.9042 + 616.9563 nm; combined gf and EW (AGSS21 Table A.1 note a)"


def main(pdf: str) -> None:
    import fitz
    text = "\n".join(p.get_text() for p in fitz.open(pdf))
    tab = text.split("Table A.1. Line list", 1)[1].split("Table A.2. Line list", 1)[0]
    toks = [t.strip() for t in tab.split("\n") if t.strip()]
    scott = {}
    with open(BASE / "agss21_cited_lines.csv") as fh:
        for r in csv.DictReader(fh):
            if r["paper"] == "scott2015_NatoCa":
                scott.setdefault(r["species"], []).append(float(r["wavelength_air_A"]))
    rows, sp, k = [], None, 0
    while k < len(toks):
        t = toks[k]
        if t in HDR:
            sp, k = HDR[t], k + 1
            continue
        m = re.fullmatch(r"(\d{3,4}\.\d+)( \(a\))?", t)
        if sp and m and k + 12 < len(toks):
            try:
                f = [float(v.replace("−", "-")) for v in toks[k + 1:k + 13]]
            except ValueError:
                k += 1
                continue
            k += 13
            if sp not in KEEP:
                continue
            wl = float(m.group(1)) * 10
            near = [w for w in scott.get(sp, []) if abs(w - wl) <= 0.1]
            rows.append({
                "species": sp, "wavelength_printed_A": round(wl, 2),
                "wavelength_air_A": near[0] if len(near) == 1 and not m.group(2) else round(wl, 2),
                "ep_eV": f[0], "upper_eV": f[1], "published_loggf": f[2], "ew_pm": f[3],
                **{c: f[4 + i] for i, c in enumerate(COLS)},
                "note": NOTE_A if m.group(2) else "",
                "doi": "10.1051/0004-6361/202140445",
                "citation": "Asplund, Amarsi & Grevesse 2021 (A&A 653, A141) Table A.1",
            })
            continue
        k += 1
    with open(BASE / "agss21_tableA1_NaKCa.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(Counter(r["species"] for r in rows))


if __name__ == "__main__":
    main(sys.argv[1])
