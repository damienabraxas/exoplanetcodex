"""Compare the canonical overlay with distinct Sobeck transitions in scratch."""
import argparse
import csv
from pathlib import Path
import runpy
import sys

import numpy as np

from config.synth_bands import SYNTH_BANDS
from pipeline import abundances_derive as synthesis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--holding", default="solar_kpno_molecfit_corrected",
                        choices=("solar_kpno_molecfit_corrected", "solar_kpno_kurucz2005_corrected"))
    parser.add_argument("--nominal-only", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(exist_ok=False)
    include = out / "input.csv"
    include.write_text("species,wavelength_air_A,ep_eV,match_tol_A\nCr I,4829.309,2.545,0.005\n")
    deck = SYNTH_BANDS["VIS"].linelist
    load = synthesis._load_synth_resources
    terms = ((4829.309, 4.0, -1.420, .0629), (4829.372, 5.0, -.510, .058))
    mode = "nominal"

    def load_pair(linelist_file=None, *args, **kwargs):
        raw, isotope, chem = load(linelist_file, *args, **kwargs)
        if Path(linelist_file or deck).resolve() != deck.resolve():
            return raw, isotope, chem
        explicit = raw.copy()
        records = []
        for wave, lower_j, gf, sigma in terms:
            indices = np.flatnonzero((raw["element"] == "Cr 1")
                                     & (abs(raw["wave_A"] - wave) <= .002)
                                     & (raw["lower_j"] == lower_j)
                                     & (raw["upper_j"] == 4)
                                     & (abs(raw["upper_state_eV"] - 5.1112) <= .002))
            if len(indices) != 1:
                raise RuntimeError(f"Expected one matching Cr I component at {wave}: {indices}")
            idx, = indices
            shift = {"both_minus": -sigma, "both_plus": sigma}.get(mode, 0.0)
            if mode != "control":
                explicit["loggf"][idx] = gf + shift
            records.append([mode, wave, lower_j, raw["loggf"][idx], explicit["loggf"][idx], sigma])
        with (out / "oscillator_strengths.csv").open("a", newline="") as handle:
            csv.writer(handle).writerows(records)
        return explicit, isotope, chem

    (out / "oscillator_strengths.csv").write_text(
        "mode,air_A,lower_J,canonical_overlay_loggf,explicit_lab_loggf,source_sigma_dex\n")
    synthesis._load_synth_resources = load_pair
    runner = "scripts/derive_band_products.py"
    sys.path.insert(0, str(Path(runner).resolve().parent))
    modes = ("nominal",) if args.nominal_only else ("control", "nominal", "both_minus", "both_plus")
    for mode in modes:
        stage = out / mode
        stage.mkdir()
        sys.argv = [runner, "--element", "Cr", "--ion", "I", "--lo", "4823", "--hi", "4836",
                    "--instrument", "kpno_solar_atlas", "--holding", args.holding,
                    "--force-synthesis", "--skip-engine-b",
                    "--lines-from-set", f"CRI_RYA1397_LAB={include}", "--out", str(stage)]
        print(f"mode={mode} argv={sys.argv}", flush=True)
        runpy.run_path(runner, run_name="__main__")


if __name__ == "__main__":
    main()
