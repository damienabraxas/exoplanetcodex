"""Sirius-only diagnostic for two identified NIST N I levels missing in VALD."""
import csv
from dataclasses import replace
import hashlib
from pathlib import Path
import runpy
import sys

import numpy as np

from config.synth_bands import SYNTH_BANDS
from pipeline import nearuv_linelist as deck


def main():
    base_path, components, out = map(Path, sys.argv[1:])
    out = out.resolve()
    out.mkdir()
    deck._ensure_ispec_on_path()
    import ispec

    base = ispec.read_atomic_linelist(str(base_path))
    mask = (base["element"] == "N 1") & (abs(base["wave_A"] - 10539.575) < .002)
    if mask.sum() != 1:
        raise RuntimeError("The existing same-multiplet VALD damping template is required")
    with components.open() as source:
        records = list(csv.DictReader(source))
    rows = np.repeat(base[mask], len(records))
    input_csv = out / "input.csv"
    with input_csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["species", "wavelength_air_A", "ep_eV", "match_tol_A"])
        for row, record in zip(rows, records, strict=True):
            w, elo, eup, jl, ju, gf = (float(record[k]) for k in
                ("wavelength_air_A", "lower_cm-1", "upper_cm-1", "lower_J", "upper_J", "loggf"))
            if np.any((base["element"] == "N 1") & (abs(base["wave_A"] - w) < .015)):
                raise RuntimeError(f"Duplicate nitrogen transition {w}")
            for field, value in {
                "wave_A":w, "wave_nm":w/10, "lower_state_cm1":elo, "upper_state_cm1":eup,
                "lower_state_eV":elo/8065.544, "upper_state_eV":eup/8065.544,
                "lower_j":jl, "upper_j":ju, "upper_g":2*ju+1, "loggf":gf,
                "reference_code":"NIST T7370",
            }.items():
                row[field] = value
            writer.writerow(["N I", f"{w:.3f}", f"{elo/8065.544:.5f}", ".005"])

    modified = out / "diagnostic_atomic_lines.tsv"
    merged = np.concatenate([base, rows])
    deck.write(merged[np.argsort(merged["wave_A"], kind="stable")], modified)
    if len(ispec.read_atomic_linelist(str(modified))) != len(merged):
        raise RuntimeError("Synth deck row count changed on round trip")
    with (out/"sha256.txt").open("w") as stream:
        for source in (base_path, components, modified, input_csv):
            stream.write(f"{hashlib.sha256(source.read_bytes()).hexdigest()}  {source}\n")
    print("Line-specific broadening unavailable. Use same-multiplet VALD template as a diagnostic-only model.",
          flush=True)
    SYNTH_BANDS["NIR"] = replace(SYNTH_BANDS["NIR"], linelist_spec=str(modified))
    runner = "scripts/derive_band_products.py"
    sys.path.insert(0, str(Path(runner).resolve().parent))
    sys.argv = [runner, "--element", "N", "--ion", "I", "--lo", "10501", "--hi", "10527",
                "--instrument", "crires_plus", "--holding", "solar_crires_plus_y_wide_rya1054",
                "--force-synthesis", "--skip-engine-b", "--lines-from-set",
                f"NI_CAFFAU_RYA1397_NIST_Y_DIAGNOSTIC={input_csv}", "--out", str(out)]
    print(sys.argv, flush=True)
    runpy.run_path(runner, run_name="__main__")


if __name__ == "__main__":
    main()
