"""Run the V I HARPS component insertion test on Sirius with no publication."""

import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

from pipeline.line_patch import LinePatch, LinePatchError
from pipeline.run_descriptor import RunDescriptor, resolve
from pipeline.run_matrix import verify_numpy_ceiling

WINDOW = (4200.0, 6908.0)
STEM = "VI_4200_6908_harps_solar_harps_molecfit_corrected_SYNTH_SET-V_AGSS21"
PRIOR = Path("data/results/orchestrator/budget/solar_V/stage").resolve()
SOURCE = Path("data/reference/line_sets/v_agss21_VI_graded.csv")
# Lawler Table 3 identities, with the corresponding GES energy pair.
FEATURES = {
    4875.486: (322.622, 20825.235, 3.5, 2.5, 5),
    5731.249: (8581.739, 26019.446, 2.5, 3.5, 1),
    6111.651: (8412.362, 24769.286, 0.5, 0.5, 2),
}


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    interpreter = os.environ["CODEX_SYNTH_PYTHON"]
    ok, reason = verify_numpy_ceiling(interpreter)
    if not ok:
        raise SystemExit(reason)
    holding = "solar_harps_molecfit_corrected"
    descriptor = RunDescriptor("V", "I", "harps", holding, *WINDOW,
                               "ts-lte", "synthesis", "reference")
    plan = resolve(descriptor, interpreter=interpreter, ispec_dir=os.environ["ISPEC_DIR"])
    if plan.blocked_reason:
        raise SystemExit(plan.blocked_reason)
    step, = plan.steps

    prior_tables = {t: PRIOR / f"{STEM}_{t}_lines.csv" for t in ("1D-LTE", "ENGINE-B")}
    prior_waves = {
        t: {round(float(row["wavelength_air_A"]), 3) for row in read(path)}
        for t, path in prior_tables.items()
    }
    if len(prior_waves["1D-LTE"]) != 27 or prior_waves["1D-LTE"] != prior_waves["ENGINE-B"]:
        raise LinePatchError("Unexpected baseline pool membership")

    sys.path.insert(0, str(Path("scripts").resolve()))
    import derive_band_products as derive

    select = derive._cand_from_ew_artifact
    bind = LinePatch.bind
    load_prior = LinePatch.prior
    admitted = set()
    overlay = out / "v_agss21_harps_components.csv"
    selection = out / "refit_components.csv"

    def select_with_components(linelist, path, **kw):
        source_rows = read(path)
        component_records = []
        refits = []
        waves = linelist["wave_A"] if "wave_A" in linelist.dtype.names else linelist["wave_nm"] * 10
        for row in source_rows:
            center = round(float(row["wavelength_air_A"]), 3)
            if center not in FEATURES:
                continue
            low, up, low_j, up_j, count = FEATURES[center]
            mask = ((linelist["element"] == "V 1") & (np.abs(waves - center) < 0.5)
                    & np.isclose(linelist["lower_state_cm1"], low, rtol=0, atol=0.1)
                    & np.isclose(linelist["upper_state_cm1"], up, rtol=0, atol=0.1)
                    & (linelist["lower_j"] == low_j) & (linelist["upper_j"] == up_j))
            group = np.flatnonzero(mask)
            if len(group) != count:
                raise LinePatchError(f"{center}: expected {count} deck components, found {len(group)}")
            loggf = linelist["loggf"][group]
            strongest = group[np.argmax(loggf)]
            representative = round(float(waves[strongest]), 3)
            total_loggf = math.log10(float(np.sum(10.0 ** loggf)))
            if abs(total_loggf - float(row["loggf_published"])) > 0.02:
                raise LinePatchError(
                    f"{center}: component sum {total_loggf:.5f} differs from published {row['loggf_published']}"
                )
            admitted.add(representative)
            row["wavelength_air_A"] = f"{representative:.3f}"
            refits.append({"species": "V I", "wavelength_air_A": representative,
                           "ep_eV": float(linelist["lower_state_eV"][strongest])})
            for i in group:
                component_records.append({
                    "laboratory_air_A": center, "fit_representative_air_A": representative,
                    "component_air_A": float(waves[i]), "component_loggf": float(linelist["loggf"][i]),
                    "total_loggf": total_loggf, "lower_cm1": low, "upper_cm1": up,
                    "lower_J": low_j, "upper_J": up_j,
                })
        if len(admitted) != len(FEATURES) or admitted & prior_waves["1D-LTE"]:
            raise LinePatchError("Expected three new, distinct features")
        write(overlay, source_rows[0], source_rows)
        write(selection, refits[0], refits)
        write(out / "component_mapping.csv", component_records[0], component_records)
        return select(linelist, overlay, **kw)

    def strict_bind(self, waves, energies, *args):
        current = {round(float(w), 3) for w in waves}
        if self.stem != STEM or current != prior_waves["1D-LTE"] | admitted:
            raise LinePatchError(f"Pool changed unexpectedly: {self.stem}, new {current - prior_waves['1D-LTE']}")
        bind(self, waves, energies, *args)
        if self.refit_wl != admitted or self.not_in_pool:
            raise LinePatchError(f"New feature binding differs: {self.refit_wl}, {self.not_in_pool}")

    def strict_prior(self, treatment, wave, ep):
        row = load_prior(self, treatment, wave, ep)
        if row is None and round(float(wave), 3) not in admitted:
            raise LinePatchError(f"Stored fit unavailable for {treatment} {wave:.4f}")
        return row

    derive._cand_from_ew_artifact = select_with_components
    LinePatch.bind = strict_bind
    LinePatch.prior = strict_prior
    args = list(step["args"])
    args[args.index("--out") + 1] = str(out)
    i = args.index("--lines-tier")
    args[i:i + 2] = ["--lines-from-set", f"V_AGSS21={SOURCE}"]
    sys.argv = [step["script"], *args, "--refit-lines", str(selection),
                "--reuse-fits-from", str(PRIOR)]
    (out / "inputs.json").write_text(json.dumps({
        "command": sys.argv, "source_sha256": checksum(SOURCE),
        "prior_sha256": {t: checksum(p) for t, p in prior_tables.items()},
    }, indent=2) + "\n")
    derive.main()


if __name__ == "__main__":
    main()
