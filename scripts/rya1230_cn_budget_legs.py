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

STEP = 0.10
XI_NOMINAL, XI_STEP = 1.0, 0.10
#: base pins per region: the CN regions pin C and O (RYA-1214); `vis` fits C, N and O
#: jointly and pins nothing. A coupling leg pins ONE other element at its solar value
#: +/- STEP (the value printed in each nominal run's "seed A" line is the reference).
BASE_PINS = {"nir_cn_iag": {"C": 8.488, "O": 8.690}, "nir_cn_kp": {"C": 8.488, "O": 8.690},
             "j_cn_crires": {"C": 8.488, "O": 8.690}, "vis": {},
             "vis_kp_k05": {}, "vis_kp_mf": {}, "vis_iag": {}}
#: the pin reference for a `vis` coupling leg: the nominal joint solution is read from the
#: nominal run's per-band output at launch time (see main)
VIS_REF: dict = {}
ENV_LEGS = {
    "nominal": {}, "xi_minus": {"CODEX_XI_OVERRIDE": "0.90"}, "xi_plus": {"CODEX_XI_OVERRIDE": "1.10"},
    "q80": {"CODEX_CONT_Q": "80"}, "q97": {"CODEX_CONT_Q": "97"},
    "contref": {"CODEX_CONT_REF": "iag"},
    "marcs": {"CODEX_MODEL_GRID": "MARCS.GES"}, "cscale": {"CODEX_CONT_SCALE": "1.001"},
    "win60": {"CODEX_CNO_WINDOW_SCALE": "0.6"},
}
PIN_LEGS = {f"{el.lower()}_{side}": (el, sgn) for el in ("C", "N", "O", "Ni")
            for side, sgn in (("minus", -1), ("plus", +1))}
LEGS = {**{k: None for k in ENV_LEGS}, **{k: None for k in PIN_LEGS}}


def leg_spec(region: str, leg: str):
    pins = dict(BASE_PINS[region])
    env = dict(ENV_LEGS.get(leg, {}))
    if leg in PIN_LEGS:
        el, sgn = PIN_LEGS[leg]
        ref = pins.get(el, VIS_REF.get(el))
        if ref is None:
            raise SystemExit(f"{region}/{leg}: no reference value for A({el})")
        pins[el] = round(ref + sgn * STEP, 3)
    return env, pins


def run_leg(region: str, leg: str, out: Path) -> str:
    env_extra, pins = leg_spec(region, leg)
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
           *[x for el, v in pins.items() for x in ("--pin", f"{el}={v}")],
           "--no-systematics", "--out", str(d)]
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
    if "vis" in a.regions:
        import pandas as pd
        nom = a.out / "vis" / "nominal" / "solar_vis_cno_per_band.csv"
        if nom.exists():
            d = pd.read_csv(nom)
            for el in ("C", "N", "O"):
                prim = d[(d.element == el) & (d.role == "primary")]
                if len(prim):
                    VIS_REF[el] = float(prim.A_X.iloc[0])
        import re
        seed = (a.out / "vis" / "nominal" / "run.log")
        if seed.exists():
            m = re.search(r"Ni=([0-9.]+)", seed.read_text())
            if m:
                VIS_REF["Ni"] = float(m.group(1))
    tasks = [(r, leg) for r in a.regions for leg in a.legs.split(",")]
    print(f"{len(a.regions)} regions x {len(a.legs.split(','))} legs = {len(tasks)} runs", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for msg in ex.map(lambda t: run_leg(t[0], t[1], a.out), tasks):
            print(msg, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
