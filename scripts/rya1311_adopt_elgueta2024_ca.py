#!/usr/bin/env python3
"""RYA-1311 -- grade the Y/J-band Ca lines with Elgueta et al. 2024 (MNRAS 532, 3694; Table 8),
the only published Ca gf with a per-line uncertainty beyond Den Hartog et al. 2021's reach.

    python3 scripts/rya1311_adopt_elgueta2024_ca.py [--check]

Why this source: no laboratory gf exists for any Ca line from 0.9 to 2.5 um. NIST ASD has no
Ca I transition probability from 9000 to 25000 A, and Ca II has four, all accuracy D
(theoretical). Den Hartog et al. 2021 stops at 10000 A, and its 9099-9416 A lines are themselves
theory. The APOGEE list (Shetrone 2015; Smith 2021) is tuned to the solar spectrum, so it is
circular for a solar measurement; Elgueta et al. 2026 refuse it for the same reason. Elgueta 2024
calibrate log gf on WINERED Cepheid spectra against Luck (2018) abundances (not the Sun) and
print per line the SD of the per-spectrum offsets: a published gf uncertainty, i.e. Codex Grade
(Ryan 2026-10-09, `pipeline.gf_grades.is_gf_graded`).

Per canonical Ca row, the calibration of the VALD gf is preferred (our row IS the VALD value),
else the calibration of the Melendez & Barbuy 1999 gf. sigma_dex = the paper's SD. gf_tier
E2024-ASTRO: priced, not LAB, so it never enters the Reference pool. A row is matched on air
wavelength, unique within 0.03 A with no second candidate within 0.10 A; anything else is
reported, never guessed. Every replaced value is written to
data/audit/rya1311_ca_elgueta2024/adoption.csv.

CAVEAT (measured, see the source README): on the three Fe I lines of their Table 7 that we hold
with a laboratory gf, their calibrated values sit 0.14 dex below the lab values.
"""
from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/linelists/primary_gf/elgueta2024/ca_table8.csv"
GF = ROOT / "data/linelists/canonical_gf.csv"
AUDIT = ROOT / "data/audit/rya1311_ca_elgueta2024/adoption.csv"
DOI = "10.1093/mnras/stae1674"
REF = "Elgueta et al. 2024, MNRAS 532, 3694, Table 8 (astrophysical, Cepheid-calibrated)"
TIER, TAG = "E2024-ASTRO", "E2024_astro"
TOL_A, CLEAR_A = 0.03, 0.10
STATUS = "adjudicated_rya1311_elgueta2024"


def read_table8() -> pd.DataFrame:
    d = pd.read_csv(SRC)
    vald = d.vald_calibrated_loggf.notna()
    d["loggf"] = d.vald_calibrated_loggf.where(vald, d.mb99_calibrated_loggf)
    d["sigma"] = d.vald_sd_dex.where(vald, d.mb99_sd_dex)
    d["calibrated_from"] = ["VALD" if v else "MB99" for v in vald]
    return d


def match(cg_sp: pd.DataFrame, w: float):
    dw = (cg_sp.wavelength_air_A - w).abs()
    near = dw[dw <= CLEAR_A].sort_values()
    if len(near) == 0 or near.iloc[0] > TOL_A:
        return None, "no canonical row within 0.03 A (not in our line list)"
    if len(near) > 1:
        return None, f"{len(near)} canonical rows within 0.10 A -- ambiguous"
    return near.index[0], ""


def _splice(edits: dict) -> None:
    """Rewrite ONLY the edited rows, byte for byte elsewhere: a pandas round-trip reformats
    floats on unrelated rows (-0.890 -> -0.89), which is a diff nobody reviewed."""
    lines = GF.read_text().splitlines(keepends=True)
    header = next(csv.reader([lines[0]]))
    done = 0
    for k, line in enumerate(lines[1:], 1):
        lid = line.split(",", 1)[0]
        if lid not in edits:
            continue
        row = dict(zip(header, next(csv.reader([line]))))
        row.update(edits[lid])
        buf = io.StringIO()
        csv.writer(buf, lineterminator=line[len(line.rstrip("\r\n")):]).writerow(
            [row[h] for h in header])
        lines[k] = buf.getvalue()
        done += 1
    if done != len(edits):
        raise RuntimeError(f"spliced {done} of {len(edits)} rows; refusing to write")
    GF.write_text("".join(lines))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="report only, write nothing")
    a = ap.parse_args(argv)
    cg = pd.read_csv(GF, low_memory=False)
    audit, edits = [], {}
    for _, r in read_table8().iterrows():
        rec = {"species": r.species, "air_A": r.wavelength_air_A, "ep_eV": r.ep_eV,
               "e24_loggf": r.loggf, "e24_sigma_dex": r.sigma, "calibrated_from": r.calibrated_from}
        i, why = match(cg[cg.species == r.species], float(r.wavelength_air_A))
        if i is None:
            audit.append({**rec, "action": "unmatched", "why": why})
            continue
        old = cg.loc[i]
        if str(old.gf_tier) == "LAB" or pd.notna(old.gf_sigma_dex) and str(old.lab_source_tag) != TAG:
            audit.append({**rec, "action": "kept", "why": f"already graded ({old.gf_tier})"})
            continue
        audit.append({**rec, "action": "adopted", "why": "", "canonical_A": old.wavelength_air_A,
                      "old_loggf": old.log_gf, "old_tier": old.gf_tier,
                      "old_reference": old.loggf_reference,
                      "d_loggf": round(float(r.loggf) - float(old.log_gf), 3)})
        if a.check:
            continue
        edits[old.line_id] = {"log_gf": f"{float(r.loggf):.3f}", "gf_sigma_dex": f"{float(r.sigma):.3f}",
                              "loggf_reference": f"{REF}, calibrated from {r.calibrated_from}",
                              "gf_source_doi": DOI, "gf_tier": TIER, "lab_source_tag": TAG,
                              "adjudication_status": STATUS}
    au = pd.DataFrame(audit)
    print(au[["species", "air_A", "action", "e24_loggf", "e24_sigma_dex", "why"]].to_string(index=False))
    if not a.check:
        _splice(edits)
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        au.to_csv(AUDIT, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
