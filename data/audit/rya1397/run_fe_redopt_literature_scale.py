"""Test the AGSS21 strength scale on the exact four measured iron transitions."""
import csv
from pathlib import Path
import runpy
import sys

import numpy as np

from config.synth_bands import SYNTH_BANDS
from pipeline import abundances_derive as synthesis

# deck wavelength, E_low, J_low, E_up, J_up, canonical log gf, literature log gf
ITEMS = {
    "I": [(8293.513, 3.301, 2.0, 4.796, 2.0, -2.164, -2.203)],
    "II": [(7222.391, 3.889, 1.5, 5.605, .5, -3.36, -3.26),
           (7224.478, 3.889, .5, 5.605, .5, -3.24, -3.20),
           (7515.830, 3.903, 3.5, 5.553, 2.5, -3.46, -3.39)],
}


def main():
    root = Path(sys.argv[1]).resolve()
    root.mkdir()
    deck = SYNTH_BANDS["red-optical"].linelist
    load = synthesis._load_synth_resources
    runner = "scripts/derive_band_products.py"
    sys.path.insert(0, str(Path(runner).resolve().parent))
    ion = "I"

    def inject(linelist_file=None, *args, **kwargs):
        rows, isotope, chem = load(linelist_file, *args, **kwargs)
        if Path(linelist_file or deck).resolve() != deck.resolve():
            return rows, isotope, chem
        output = rows.copy()
        records = []
        token = "Fe 1" if ion == "I" else "Fe 2"
        for wave, ep, jl, upper, ju, canonical, literature in ITEMS[ion]:
            hits = np.flatnonzero((rows["element"] == token) & (abs(rows["wave_A"] - wave) <= .002)
                                  & (abs(rows["lower_state_eV"] - ep) <= .002)
                                  & (abs(rows["upper_state_eV"] - upper) <= .002)
                                  & (rows["lower_j"] == jl) & (rows["upper_j"] == ju))
            if len(hits) != 1 or not np.isclose(rows["loggf"][hits[0]], canonical, atol=.002):
                raise RuntimeError(f"Unexpected canonical payload for Fe {ion} {wave}: {hits}")
            output["loggf"][hits[0]] = literature
            records.append([ion, wave, ep, jl, upper, ju, rows["loggf"][hits[0]], literature])
        with (root / "strength-override.csv").open("a", newline="") as stream:
            csv.writer(stream).writerows(records)
        return output, isotope, chem

    (root / "strength-override.csv").write_text(
        "ion,wavelength_air_A,lower_eV,lower_J,upper_eV,upper_J,canonical_loggf,literature_loggf\n")
    synthesis._load_synth_resources = inject
    for ion, lo, hi, holding in [
        ("I", 8288, 8299, "solar_kpno_molecfit_corrected"),
        ("II", 7217, 7521, "solar_kpno_kurucz2005_corrected"),
    ]:
        selector = root / f"Fe{ion}_input.csv"
        with selector.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["species", "wavelength_air_A", "ep_eV", "match_tol_A"])
            for w, ep, *_ in ITEMS[ion]:
                writer.writerow([f"Fe {ion}", w, ep, .005])
        out = root / f"Fe{ion}"
        out.mkdir()
        sys.argv = [runner, "--element", "Fe", "--ion", ion, "--lo", str(lo), "--hi", str(hi),
                    "--instrument", "kpno_solar_atlas", "--holding", holding,
                    "--force-synthesis", "--skip-engine-b", "--lines-from-set",
                    f"FE_RYA1397_AGSS21_DIAGNOSTIC={selector}", "--out", str(out)]
        print(sys.argv, flush=True)
        runpy.run_path(runner, run_name="__main__")


if __name__ == "__main__":
    main()
