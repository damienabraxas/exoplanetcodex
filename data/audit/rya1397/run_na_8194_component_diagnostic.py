"""Measure the Na I 8194 Å fine-structure strength response on Kitt Peak."""
import csv
import math
from pathlib import Path
import runpy
import sys

import numpy as np

from config.synth_bands import SYNTH_BANDS
from pipeline import abundances_derive as synthesis


def main():
    out = Path(sys.argv[1]).resolve()
    out.mkdir()
    deck = SYNTH_BANDS["VIS"].linelist.resolve()
    load = synthesis._load_synth_resources
    selector = out / "input.csv"
    selector.write_text("species,wavelength_air_A,ep_eV,match_tol_A\nNa I,8194.824,2.104,.005\n")
    strengths = out / "oscillator_strengths.csv"
    strengths.write_text("mode,air_A,upper_J,overlay_loggf,used_loggf\n")

    def load_pair(linelist_file=None, *args, **kwargs):
        rows, isotope, chemistry = load(linelist_file, *args, **kwargs)
        if Path(linelist_file or deck).resolve() != deck:
            return rows, isotope, chemistry
        rows = rows.copy()
        records = []
        # Lind et al. oscillator strengths times the lower-level weight 2J+1 = 4.
        for wavelength, upper_j, oscillator_strength in (
            (8194.790, 1.5, 0.0863), (8194.824, 2.5, 0.777),
        ):
            match = (rows["element"] == "Na 1") & (abs(rows["wave_A"] - wavelength) < .002)
            match &= (rows["lower_j"] == 1.5) & (rows["upper_j"] == upper_j)
            match &= abs(rows["lower_state_eV"] - 2.104) < .002
            indices = np.flatnonzero(match)
            if len(indices) != 1:
                raise RuntimeError(f"Expected one Na I component at {wavelength}, found {indices}")
            index, = indices
            overlay = float(rows["loggf"][index])
            if mode == "lind":
                rows["loggf"][index] = math.log10(4 * oscillator_strength)
            records.append((mode, wavelength, upper_j, overlay, rows["loggf"][index]))
        with strengths.open("a", newline="") as stream:
            csv.writer(stream).writerows(records)
        return rows, isotope, chemistry

    synthesis._load_synth_resources = load_pair
    runner = Path("scripts/derive_band_products.py").resolve()
    sys.path.insert(0, str(runner.parent))
    for mode in ("control", "lind"):
        stage = out / mode
        stage.mkdir()
        sys.argv = [str(runner), "--element", "Na", "--ion", "I", "--lo", "8190",
                    "--hi", "8199", "--instrument", "kpno_solar_atlas", "--holding",
                    "solar_kpno_kurucz2005_corrected", "--force-synthesis", "--skip-engine-b",
                    "--lines-from-set", f"NAI_RYA1397_8194_DIAGNOSTIC={selector}", "--out", str(stage)]
        print(f"{mode}: {sys.argv}", flush=True)
        runpy.run_path(str(runner), run_name="__main__")


if __name__ == "__main__":
    main()
