"""Inventory the held Bride environment/atom/snapshot without invoking RT or SE.

Run with the held evaluation interpreter, from the checkout. Snapshot identity
is explicit; paths resolve through the path register. JSON is audit evidence,
not an executable science-run manifest or product admission.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from importlib.metadata import version
import inspect
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def fingerprint(path: Path) -> dict:
    return {"name": path.name, "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def audit(stem: str) -> dict:
    import numpy as np
    import lightweaver as lw
    from lightweaver import rh_atoms
    from config.constants import get_star_params, require_codex_path
    from pipeline.bride_foundation import validate_snapshot
    from rya1118_multi3d_read import read_atmos_multi3d

    if Path(stem).name != stem or not stem or stem in (".", ".."):
        raise ValueError("snapshot-stem must be an explicit filename stem")
    held = require_codex_path("engines.bride_eval")
    paths = {kind: held / f"{kind}.{stem}" for kind in ("mesh", "atm3d", "info")}
    snapshot = read_atmos_multi3d(paths["mesh"], paths["atm3d"])
    shape = validate_snapshot(snapshot)
    atoms = []
    for name in ("O_atom", "OI_fine_atom", "OI_fine_crd_atom", "OI_ord_atom", "O_no_777_atom"):
        atom = getattr(rh_atoms, name)()
        lines = []
        # This band is a diagnostic inventory of the named O I 777 multiplet,
        # not a wavelength-matching/adoption tolerance or a synthesis window.
        for line in atom.lines:
            if 776 < float(line.lambda0) < 779:
                lower, upper = atom.levels[line.i], atom.levels[line.j]
                lines.append({"i": line.i, "j": line.j,
                              "vacuum_nm": float(line.lambda0), "f": float(line.f),
                              "lower_label": lower.label, "upper_label": upper.label,
                              "lower_g": lower.g, "upper_g": upper.g,
                              "lower_stage": lower.stage, "upper_stage": upper.stage,
                              "broadening": str(line.broadening)})
        atoms.append({"factory": name, "levels": len(atom.levels), "lines": len(atom.lines),
                      "continua": len(atom.continua), "collisions": len(atom.collisions),
                      "collision_classes": dict(Counter(type(c).__name__ for c in atom.collisions)),
                      "multiplet_region_transitions": lines,
                      "verdict": "HOLD: no demonstrated resolved-triplet and benchmark-rate parity"})
    canonical = []
    with (ROOT / "data/linelists/canonical_gf.csv").open() as f:
        for row in csv.DictReader(f):
            if row["species"] == "O I" and 7770 < float(row["wavelength_air_A"]) < 7777:
                canonical.append({k: row[k] for k in ("line_id", "physical_id", "wavelength_air_A",
                    "excitation_potential_eV", "log_gf", "loggf_reference", "nist_grade",
                    "gf_sigma_dex", "gf_tier")})
    params = get_star_params("solar")
    return {
        "schema": "codex.bride_foundation_audit/1", "science_executed": False,
        "code_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "generator": fingerprint(Path(__file__)),
        "input_validator": fingerprint(ROOT / "pipeline/bride_foundation.py"),
        "canonical_gf_source": fingerprint(ROOT / "data/linelists/canonical_gf.csv"),
        "star_params_source": fingerprint(ROOT / "config/stars.yaml"),
        "held_path_key": "engines.bride_eval", "snapshot_stem": stem,
        "files": {k: fingerprint(p) for k, p in paths.items()},
        "info_header": paths["info"].read_text().splitlines(),
        "shape_zyx": list(shape),
        "axis_endpoints_m": {k: [float(snapshot[k][0]), float(snapshot[k][-1])] for k in ("x", "y", "z")},
        "mean_temperature_endpoints_K": snapshot["temperature"].mean(axis=(1, 2))[[0, -1]].tolist(),
        "environment": {p: version(p) for p in ("lightweaver", "numpy", "scipy", "numba")},
        "python": sys.version.split()[0],
        "atom_source": fingerprint(Path(inspect.getsourcefile(rh_atoms))),
        "atmosphere_source": fingerprint(Path(inspect.getsourcefile(lw.Atmosphere))),
        "make_1d_signature": str(inspect.signature(lw.Atmosphere.make_1d)),
        "gravity": {"source": "config/stars.yaml via get_star_params(solar)",
                    "canonical_logg_cgs": params["logg"],
                    "lightweaver_logG_si": params["logg"] - float(np.log10(100.))},
        "atoms": atoms, "canonical_triplet": canonical,
        "solver_timing": None,
        "timing_reason": "No O NLTE column run authorized; Ca/Mg proxy times are not O timings",
        "first_light_verdict": "HOLD: reference atom and independent baseline gates remain open",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot-stem", required=True)
    parser.add_argument("--output", type=Path, help="New audit file; existing files are never overwritten")
    args = parser.parse_args()
    doc = audit(args.snapshot_stem)
    data = json.dumps(doc, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as f:
            f.write(data)
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
