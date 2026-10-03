#!/usr/bin/env python3
"""RYA-1233 -- the Si Reference lines ARE Asplund's: adopt Amarsi & Asplund 2017 Table 1.

Ryan, 2026-10-03: graded lines are the best lines, as close to Asplund's as possible.
For Si that is the AGSS21 line set (Asplund 2021 -> Amarsi & Asplund 2017, MNRAS 464,
264, doi 10.1093/mnras/stw2445 -> the Scott 2015 lines), with the gf THEY used: Si I =
Garz 1973 laboratory values + 0.097 dex (renormalised to the O'Brian & Lawler 1991 lab
lifetimes); Si II 6371 = mean of Schulz-Gulde 1969, Blanco 1995, Matheron 2001
(sigma 0.02). RYA-1169 already transcribed the table and matched each row to its
canonical line id: data/audit/rya1169_si_intake/si_agss21_reference_lines.csv.

Every row Amarsi & Asplund USED becomes gf_tier LAB with their log gf. The two violet
lines they explicitly rejected (3905, 4102) are left as they are (Den Hartog 2023 LAB).
Si I per-line Garz uncertainties were never transcribed: gf_sigma_dex stays empty and the
row says so. Idempotent.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
TABLE = ROOT / "data" / "audit" / "rya1169_si_intake" / "si_agss21_reference_lines.csv"
DOI = "10.1093/mnras/stw2445"           # bibliography: amarsi2017_si (Crossref-verified)


def main() -> None:
    cg = pd.read_csv(CANONICAL, low_memory=False)
    with TABLE.open(newline="", encoding="utf-8") as fh:
        table = [r for r in csv.DictReader(fh) if r["reference_status"] == "used"]
    done = []
    for r in table:
        hit = cg.index[cg.line_id.astype(str) == r["canonical_line_id"]]
        if len(hit) != 1:
            raise SystemExit(f"{r['species']} {r['wavelength_air_A']}: canonical line id "
                             f"{r['canonical_line_id']!r} matched {len(hit)} rows -- refusing")
        i = hit[0]
        if abs(float(cg.at[i, "wavelength_air_A"]) - float(r["wavelength_air_A"])) > 0.02:
            raise SystemExit(f"{r['canonical_line_id']}: wavelength disagrees with the table")
        old = cg.at[i, "log_gf"]
        cg.at[i, "log_gf"] = float(r["published_loggf"])
        cg.at[i, "gf_tier"] = "LAB"
        cg.at[i, "nist_grade"] = float("nan")
        cg.at[i, "loggf_reference"] = (f"Amarsi & Asplund 2017 Table 1 (AGSS21 Si set): "
                                       f"{r['published_gf_source']}")
        cg.at[i, "gf_source_doi"] = DOI
        cg.at[i, "lab_source_tag"] = "AGSS21_Si_Amarsi2017"
        sig = r.get("published_gf_sigma_dex") or ""
        cg.at[i, "gf_sigma_dex"] = float(sig) if sig else float("nan")
        cg.at[i, "adjudication_status"] = (
            "RYA-1233: Reference line = Asplund's; "
            + ("" if sig else "per-line Garz 1973 sigma not yet transcribed"))
        done.append((r["species"], r["wavelength_air_A"], old, r["published_loggf"]))
    cg.to_csv(CANONICAL, index=False)
    for d in done:
        print(f"  {d[0]:<6} {d[1]:>9}  log gf {d[2]} -> {d[3]}  LAB")
    print(f"{len(done)} Asplund Si lines adopted as LAB")


if __name__ == "__main__":
    main()
