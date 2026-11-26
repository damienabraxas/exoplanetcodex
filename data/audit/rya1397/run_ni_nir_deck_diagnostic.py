"""Build a scratch-only N I deck from a held VALD extract and run the diagnostic."""
import argparse
import csv
from dataclasses import replace
import hashlib
from pathlib import Path
import runpy
import sys

import numpy as np

from config.synth_bands import SYNTH_BANDS
from pipeline import nearuv_linelist as deck

TARGETS = ((10105.130, 11.750), (10112.480, 11.758),
           (10114.640, 11.764), (10539.570, 11.844))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vald_extract", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(exist_ok=True)
    original = SYNTH_BANDS["NIR"].linelist
    modified = out / "diagnostic_atomic_lines.tsv"
    if modified.exists():
        raise SystemExit(f"Choose a new versioned output directory: {modified}")
    records, report = deck.read_band(args.vald_extract, 10100.0, 10545.0)
    selected = []
    for wave, ep in TARGETS:
        match = [r for r in records if r["element"] == "N" and r["ion"] == "I"
                 and abs(r["wavelength"] - wave) <= .015
                 and abs(r["e_low_eV"] - ep) <= .01]
        if len(match) != 1:
            raise SystemExit(f"N I {wave}: {len(match)} agreeing VALD tuples")
        selected.append(match[0])
    print(f"VALD parser: {report}", flush=True)

    deck._ensure_ispec_on_path()
    import ispec

    base = ispec.read_atomic_linelist(str(original))
    add = deck.to_ispec_array(selected).astype(base.dtype)
    for row in add:
        nearby = (base["element"] == "N 1") & (abs(base["wave_A"] - row["wave_A"]) <= .015)
        if np.any(nearby):
            raise SystemExit(f"Unexpected existing N I deck row at {row['wave_A']}")
    merged = np.concatenate((base, add))
    deck.write(merged[np.argsort(merged["wave_A"], kind="stable")], modified)
    check = ispec.read_atomic_linelist(str(modified))
    assert len(check) == len(base) + len(TARGETS)
    with (out / "sha256.txt").open("w") as stream:
        for path in [args.vald_extract, original, modified]:
            stream.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path}\n")
    source = out / "input.csv"
    with source.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["species", "wavelength_air_A", "ep_eV", "match_tol_A"])
        for row in add:
            writer.writerow(["N I", f"{row['wave_A']:.3f}", f"{row['lower_state_eV']:.4f}", .005])

    # The configuration override is scoped to this process. Shared decks stay intact.
    SYNTH_BANDS["NIR"] = replace(SYNTH_BANDS["NIR"], linelist_spec=str(modified))
    runner = "scripts/derive_band_products.py"
    sys.path.insert(0, str(Path(runner).resolve().parent))
    sys.argv = [runner, "--element", "N", "--ion", "I", "--lo", "10100", "--hi", "10545",
                "--instrument", "crires_plus", "--holding", "solar_crires_plus_y_wide_rya1054",
                "--force-synthesis", "--skip-engine-b",
                "--lines-from-set", f"NI_CAFFAU_RYA1397_Y_DIAGNOSTIC={source}",
                "--out", str(out)]
    print(f"Diagnostic command: {sys.argv}", flush=True)
    runpy.run_path(runner, run_name="__main__")


if __name__ == "__main__":
    main()
