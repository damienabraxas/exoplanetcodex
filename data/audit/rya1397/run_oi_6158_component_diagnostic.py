"""Compare the standard O I 6158 Å deck with all three NIST components."""
import csv
from pathlib import Path
import runpy
import sys

import numpy as np

from config.synth_bands import SYNTH_BANDS
from pipeline import abundances_derive as synthesis


def main():
    component_path, out = map(Path, sys.argv[1:])
    out = out.resolve()
    out.mkdir()
    deck = SYNTH_BANDS["VIS"].linelist
    load = synthesis._load_synth_resources
    runner = "scripts/derive_band_products.py"
    with component_path.open() as stream:
        components = list(csv.DictReader(stream))
    mode = "control"

    def resolve_oxygen(linelist_file=None, *args, **kwargs):
        rows, isotope, chem = load(linelist_file, *args, **kwargs)
        if Path(linelist_file or deck).resolve() != deck.resolve() or mode == "control":
            return rows, isotope, chem
        mask = (rows["element"] == "O 1") & (abs(rows["wave_A"] - 6158.18) < .07)
        if mask.sum() != 1:
            raise RuntimeError("Expected exactly one original deck representative for the O I feature")
        added = np.repeat(rows[mask], len(components))
        for row, comp in zip(added, components, strict=True):
            wave, elo, eup, jl, ju, gf = (float(comp[key]) for key in (
                "wavelength_Ritz_air_A", "lower_cm-1", "upper_cm-1", "J_lower", "J_upper", "NIST_loggf"))
            for name, value in {
                "wave_A": wave, "wave_nm": wave/10, "lower_state_cm1": elo, "upper_state_cm1": eup,
                "lower_state_eV": elo/8065.544, "upper_state_eV": eup/8065.544, "lower_j": jl,
                "upper_j": ju, "upper_g": 2*ju+1, "loggf": gf, "reference_code": "NIST T5803 LS",
            }.items():
                row[name] = value
        new = np.concatenate([rows[~mask], added])
        print("Diagnostic assumption: use the existing multiplet representative's broadening "
              "for all three components, with distinct Ritz wavelengths, J and gf.", flush=True)
        return new[np.argsort(new["wave_A"], kind="stable")], isotope, chem

    synthesis._load_synth_resources = resolve_oxygen
    sys.path.insert(0, str(Path(runner).resolve().parent))
    for mode, wave in (("control", 6158.186), ("split", 6158.187)):
        stage = out / mode
        stage.mkdir()
        selector = stage / "input.csv"
        selector.write_text(f"species,wavelength_air_A,ep_eV,match_tol_A\nO I,{wave},10.741,.005\n")
        sys.argv = [runner, "--element", "O", "--ion", "I", "--lo", "6152", "--hi", "6164",
                    "--instrument", "kpno_solar_atlas", "--holding", "solar_kpno_molecfit_corrected",
                    "--force-synthesis", "--skip-engine-b", "--lines-from-set",
                    f"OI_RYA1397_6158_DIAGNOSTIC={selector}", "--out", str(stage)]
        print(f"{mode}: {sys.argv}", flush=True)
        runpy.run_path(runner, run_name="__main__")


if __name__ == "__main__":
    main()
