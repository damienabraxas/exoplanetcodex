#!/usr/bin/env python3
"""RYA-470 / RYA-733 -- Holmes et al. 2016 (ApJS 224, 35) V I lab gf into the lab table.

Holmes+2016 Table 4 is an independent V I lab set (Imperial College FTS branching fractions on
TR-LIF lifetimes): 181 lines, 3043-20231 A, nine of them past 1.6 um -- the only lab V I gf in the
H/K bands. CDS (cdsarc) was unreachable on 2026-11-17 from the Mac and from Sirius, so the table is
vendored as the CfA VizieR mirror's ASU-TSV export (vizier/J_ApJS_224_35/asu_all_tables.tsv)
instead of ReadMe + table4.dat, and parsed here into transcribed/holmes2016_v1.csv, which
rya1268_transcribe_lab_gf.py appends to lab_gf_lines.csv.

Checks (the script stops if one fails):
  * every wavelength is AIR: it must equal the air wavelength (Peck & Reeder 1972) computed from
    the row's own level energies within AIR_TOL_A. A row that fails is listed and LEFT OUT -- a
    wavelength that disagrees with its own levels is a transcription or assignment error, never
    repaired here.
  * the row count and the Holmes log(gf) count match the catalogue (208 rows, 181 with gf).

Columns: log(gf) and its uncertainty are Holmes's own (BF error and lifetime error combined,
Sect. 2). A and dA are derived from them, A = 10^loggf / (1.4992e-16 lambda^2 g_up), so the
table's A-based fields stay consistent; sigma_dex is Holmes's stated value, not re-derived.
Holmes's flag `b` (BF found by fitting the HFS of a blend with another V I line, Sect. 2.3) is
kept in `note`. Lawler, Whaling and Doerr values that Holmes reprints are NOT transcribed: a
reprint is not the primary source.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "reference" / "lab_gf_tables" / "vizier" / "J_ApJS_224_35" / "asu_all_tables.tsv"
OUT = ROOT / "data" / "reference" / "lab_gf_tables" / "transcribed" / "holmes2016_v1.csv"
CM_PER_EV = 8065.544
AIR_TOL_A = 0.01
N_ROWS, N_GF = 208, 181
DOI = "10.3847/0067-0049/224/2/35"


def air(vac_a: float) -> float:
    """Peck & Reeder (1972) standard air, the index Holmes's own notes cite."""
    s2 = (1e4 / vac_a) ** 2
    n = 1 + (8060.51 + 2480990 / (132.274 - s2) + 17455.7 / (39.32957 - s2)) * 1e-8
    return vac_a / n


def table4() -> list[dict]:
    lines = SRC.read_text().split("\n")
    start = next(i for i, ln in enumerate(lines) if ln.startswith("#Name: J/ApJS/224/35/table4"))
    hdr, rows = None, []
    for ln in lines[start:]:
        if ln.startswith("#") or not ln.strip():
            if rows:
                break
            continue
        f = [x.strip() for x in ln.split("\t")]
        if hdr is None:
            hdr = f
            continue
        try:
            float(f[1])
        except ValueError:          # the units and dashes rows under the header
            continue
        rows.append(dict(zip(hdr, f)))
    return rows


def main() -> None:
    rows = table4()
    with_gf = [r for r in rows if r["TS-log(gf)"]]
    if len(rows) != N_ROWS or len(with_gf) != N_GF:
        raise SystemExit(f"expected {N_ROWS} rows / {N_GF} with gf, found {len(rows)} / {len(with_gf)}")
    out, bad = [], []
    for r in with_gf:
        w, eu, el = float(r["lambda"]), float(r["Eup"]), float(r["Elow"])
        if abs(w - air(1e8 / (eu - el))) > AIR_TOL_A:
            bad.append((w, eu, el, round(w - air(1e8 / (eu - el)), 3)))
            continue
        gf, sig = float(r["TS-log(gf)"]), float(r["e_TS-log(gf)"])
        g_up = 2 * float(r["Jup"]) + 1
        a = 10 ** gf / (1.4992e-16 * w * w * g_up) / 1e6
        out.append({"paper": "holmes2016_v1", "species": "V I", "wavelength_air_A": round(w, 4),
                    "ep_eV": round(el / CM_PER_EV, 4), "j_lower": float(r["Jlo"]),
                    "ep_upper_eV": round(eu / CM_PER_EV, 4), "A_1e6_s": round(a, 6),
                    "dA_1e6_s": round(a * (10 ** sig - 1), 6), "loggf": gf, "loggf_sigma_dex": sig,
                    "doi": DOI, "vizier": "J/ApJS/224/35",
                    "note": "BF from HFS fit of a blended V I feature (Holmes flag b)" if r["f_TS-log(gf)"] == "b" else ""})
    for b in bad:
        print("LEFT OUT (wavelength disagrees with its own levels):", b)
    print(f"{len(out)} Holmes V I lines; {sum(1 for o in out if o['wavelength_air_A'] > 16000)} beyond 1.6 um")
    if "--write" not in sys.argv:
        return
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print("wrote", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
