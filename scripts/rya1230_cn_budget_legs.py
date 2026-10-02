#!/usr/bin/env python3
"""RYA-1230 -- RYA-587 perturbation legs for the molecular CN A-X products (cno_synthesis).

    python3 scripts/rya1230_cn_budget_legs.py --regions nir_cn_iag nir_cn_kp j_cn_crires \
        --out DIR [--jobs 9]

The CN band fit is ONE indicator, so every component is a paired single-value response.
Each leg re-runs `python -m pipeline.cno_synthesis --region R --pin C=... --pin O=...
--no-systematics` with ONE input varied, into its own directory:

  nominal              A(C) 8.488 (our CH G-band carbon, the RYA-1214 pin), A(O) 8.690
  xi_minus / xi_plus   CODEX_XI_OVERRIDE 0.90 / 1.10 (stamped, pipeline.xi_pairing) -> stellar.xi
  q80 / q97            CODEX_CONT_Q, model-guided continuum quantile               -> continuum
  marcs                CODEX_MODEL_GRID=MARCS.GES                                  -> model_atmosphere
  cscale               CODEX_CONT_SCALE=1.001 (dA/df)                              -> telluric
  win60                CODEX_CNO_WINDOW_SCALE=0.6 (fit windows, continuum unchanged) -> profile_ew
  c_minus / c_plus     --pin C = 8.388 / 8.588                                      -> molecular_coupling
  o_minus / o_plus     --pin O = 8.590 / 8.790                                      -> molecular_coupling
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

C0, O0, STEP = 8.488, 8.690, 0.10
XI_NOMINAL, XI_STEP = 1.0, 0.10
LEGS = {
    "nominal": ({}, C0, O0),
    "xi_minus": ({"CODEX_XI_OVERRIDE": "0.90"}, C0, O0),
    "xi_plus": ({"CODEX_XI_OVERRIDE": "1.10"}, C0, O0),
    "q80": ({"CODEX_CONT_Q": "80"}, C0, O0),
    "q97": ({"CODEX_CONT_Q": "97"}, C0, O0),
    "marcs": ({"CODEX_MODEL_GRID": "MARCS.GES"}, C0, O0),
    "cscale": ({"CODEX_CONT_SCALE": "1.001"}, C0, O0),
    "win60": ({"CODEX_CNO_WINDOW_SCALE": "0.6"}, C0, O0),
    "c_minus": ({}, round(C0 - STEP, 3), O0),
    "c_plus": ({}, round(C0 + STEP, 3), O0),
    "o_minus": ({}, C0, round(O0 - STEP, 3)),
    "o_plus": ({}, C0, round(O0 + STEP, 3)),
}


def run_leg(region: str, leg: str, out: Path) -> str:
    env_extra, c, o = LEGS[leg]
    d = out / region / leg
    if (d / "DONE").exists():
        return f"skip {region}/{leg}"
    d.mkdir(parents=True, exist_ok=True)
    if leg.startswith("xi_"):
        from pipeline.xi_pairing import write_stamp
        write_stamp(d, xi_kms=float(env_extra["CODEX_XI_OVERRIDE"]),
                    leg="minus" if leg == "xi_minus" else "plus",
                    xi_nominal=XI_NOMINAL, step_kms=XI_STEP, unit=region)
    cmd = [sys.executable, "-m", "pipeline.cno_synthesis", "--star", "solar", "--region", region,
           "--pin", f"C={c}", "--pin", f"O={o}", "--no-systematics", "--out", str(d)]
    (d / "argv.txt").write_text(shlex.join(cmd) + "\n" + repr(env_extra) + "\n")
    with open(d / "run.log", "w") as log:
        rc = subprocess.call(cmd, cwd=ROOT, env={**os.environ, **env_extra},
                             stdout=log, stderr=subprocess.STDOUT)
    if rc == 0:
        (d / "DONE").write_text("0\n")
    return f"rc={rc} {region}/{leg}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--jobs", type=int, default=9)
    ap.add_argument("--legs", default=",".join(LEGS))
    a = ap.parse_args()
    tasks = [(r, leg) for r in a.regions for leg in a.legs.split(",")]
    print(f"{len(a.regions)} regions x {len(a.legs.split(','))} legs = {len(tasks)} runs", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for msg in ex.map(lambda t: run_leg(t[0], t[1], a.out), tasks):
            print(msg, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
