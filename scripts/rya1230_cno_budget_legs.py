#!/usr/bin/env python3
"""RYA-1230 -- the RYA-587 perturbation legs for every re-run C/N/O band product.

    python3 scripts/rya1230_cno_budget_legs.py --units UNITS.txt --out DIR [--jobs 8]

Ryan, 2026-09-26: a re-run owes its uncertainty, not just its value. Each UNIT is one
`derive_band_products.py` invocation exactly as the nominal re-run made it (one line of
UNITS.txt, `--out` stripped). Every leg re-runs that SAME invocation with ONE input varied,
into its own directory, so the per-line pools pair (RYA-1120: dA/dp is a paired
per-line differential, never a difference of aggregates):

  xi_minus / xi_plus   CODEX_XI_OVERRIDE = 0.90 / 1.10 km/s, stamped with
                       pipeline.xi_pairing.write_stamp (RYA-1178 A) -> stellar.xi
  core                 --half-width-A 0.25 (RYA-1220's core window) -> profile_ew
  nominal              no override: the product itself, run by the same code as its legs
  q80 / q97            CODEX_CONT_Q = 80 / 97: the model-guided continuum's pixel-selection
                       quantile either side of the nominal 90 -> continuum
  marcs                CODEX_MODEL_GRID = MARCS.GES -> model_atmosphere
  tmask                CODEX_TELLURIC_MASK = 1: pixels a MEASURED sky absorbs are blanked
                       -> telluric
  cscale               CODEX_CONT_SCALE = 1.001: dA/df per line, to convert a measured
                       telluric flux residual into dex
  c_minus / c_plus     CODEX_ABUND_OFFSET = C:-/+0.10 (N I units: the CN blended into
                       the profiles follows A(C)) -> blends

A leg is complete when its directory holds DONE; a re-invocation skips completed legs.
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

XI_NOMINAL, XI_STEP = 1.0, 0.10

LEGS = {
    "xi_minus": ({"CODEX_XI_OVERRIDE": f"{XI_NOMINAL - XI_STEP:.2f}"}, []),
    "xi_plus": ({"CODEX_XI_OVERRIDE": f"{XI_NOMINAL + XI_STEP:.2f}"}, []),
    "core": ({}, ["--half-width-A", "0.25"]),
    "nominal": ({}, []),
    "q80": ({"CODEX_CONT_Q": "80"}, []),
    "q97": ({"CODEX_CONT_Q": "97"}, []),
    # RYA-1232: continuum REFERENCE leg -- the IAG telluric-free atlas instead of the
    # synthesis (Amarsi+2021's two-atlas spread); equals nominal outside 5001-11086 A
    "contref": ({"CODEX_CONT_REF": "iag"}, []),
    "marcs": ({"CODEX_MODEL_GRID": "MARCS.GES"}, []),
    "tmask": ({"CODEX_TELLURIC_MASK": "1"}, []),
    "cscale": ({"CODEX_CONT_SCALE": "1.001"}, []),
    # N I only: the CN inside the line profiles follows A(C) -> blends
    "c_minus": ({"CODEX_ABUND_OFFSET": "C:-0.10"}, []),
    "c_plus": ({"CODEX_ABUND_OFFSET": "C:+0.10"}, []),
    # O I 844.6 only (RYA-1232): the blending Fe I 8446.575 at Ruffoni+2014's lab gf instead
    # of the adopted VALD value -> the blend term
    "blendgf": ({"CODEX_GF_OVERRIDE": "pk_b3d47c33d602=-1.44"}, []),
}
DEFAULT_LEGS = ("nominal", "xi_minus", "xi_plus", "core", "q80", "q97", "contref", "marcs", "cscale")


def parse_units(path: Path) -> list[list[str]]:
    """One argv per line: the nominal invocation's derive_band_products arguments."""
    units = []
    for line in path.read_text().splitlines():
        line = line.split(">")[0].strip()
        if not line:
            continue
        argv = shlex.split(line)
        i = argv.index("scripts/derive_band_products.py")
        args = argv[i + 1:]
        if "--out" in args:
            j = args.index("--out")
            del args[j:j + 2]
        units.append(args)
    return units


def run_leg(unit: int, args: list[str], leg: str, out: Path, py: str) -> str:
    env_extra, argv_extra = LEGS[leg]
    d = out / f"unit{unit:02d}" / leg
    if (d / "DONE").exists():
        return f"skip unit{unit:02d}/{leg}"
    d.mkdir(parents=True, exist_ok=True)
    if leg.startswith("xi_"):
        from pipeline.xi_pairing import write_stamp
        write_stamp(d, xi_kms=float(env_extra["CODEX_XI_OVERRIDE"]),
                    leg="minus" if leg == "xi_minus" else "plus",
                    xi_nominal=XI_NOMINAL, step_kms=XI_STEP, unit=f"unit{unit:02d}")
    env = {**os.environ, **env_extra}
    cmd = [py, "scripts/derive_band_products.py", *args, *argv_extra, "--out", str(d)]
    (d / "argv.txt").write_text(shlex.join(cmd) + "\n" + repr(env_extra) + "\n")
    with open(d / "run.log", "w") as log:
        rc = subprocess.call(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    if rc == 0:
        (d / "DONE").write_text("0\n")
    return f"rc={rc} unit{unit:02d}/{leg}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--units", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--legs", default=",".join(DEFAULT_LEGS))
    ap.add_argument("--only-units", default=None, help="comma list of unit indices")
    a = ap.parse_args()
    units = parse_units(a.units)
    a.out.mkdir(parents=True, exist_ok=True)
    if not (a.out / "units.txt").exists():
        (a.out / "units.txt").write_text("\n".join(shlex.join(u) for u in units) + "\n")
    keep = None if a.only_units is None else {int(x) for x in a.only_units.split(",")}
    tasks = [(i, u, leg) for i, u in enumerate(units) for leg in a.legs.split(",")
             if keep is None or i in keep]
    print(f"{len(units)} units x {len(a.legs.split(','))} legs = {len(tasks)} runs", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for msg in ex.map(lambda t: run_leg(t[0], t[1], t[2], a.out, sys.executable), tasks):
            print(msg, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
