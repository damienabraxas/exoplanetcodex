"""Diagnostic, never production (RYA-1311): one 1D-LTE synthesis cell, flux vs mu=1.

Runbook 7b. Defaults to Ca I, set:CA_AGSS21, HARPS VIS.
    ISPEC_DIR=... python scripts/diag_flux_vs_mu1.py flux
    CODEX_SYNTH_CACHE=0 ISPEC_DIR=... python scripts/diag_flux_vs_mu1.py mu1 [vmac]
Env DIAG_ELEMENT / DIAG_ION / DIAG_SET override the cell; DIAG_OUT the output root.

mode=flux : the production cell (HARPS, flux) -- the CONTROL; must reproduce the published lines.
mode=mu1  : same code on the Liege disk-centre intensity atlas, Turbospectrum intensity at mu=1
            (RYA-441 injection, wrapping the fit AND the local-continuum synthesis), vsini=0
            (no rotational broadening at one point on the disk), R told 500000.
Optional 2nd arg overrides vmac (km/s) for a sensitivity leg.
"""
import gzip
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path.cwd()
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
mode = sys.argv[1]
if mode == "mu1":
    assert os.environ.get("CODEX_SYNTH_CACHE") == "0", "synth cache key ignores flux/intensity: run mu1 with CODEX_SYNTH_CACHE=0"
vmac = float(sys.argv[2]) if len(sys.argv) > 2 else None
out = Path(os.environ.get("DIAG_OUT", Path.home() / "scratch/flux_vs_mu1")) / (mode + (f"_vmac{vmac:g}" if vmac is not None else ""))
out.mkdir(parents=True, exist_ok=True)

from pipeline.run_descriptor import RunDescriptor, resolve  # noqa: E402

d = RunDescriptor(element=os.environ.get("DIAG_ELEMENT", "Ca"),
                  ion=os.environ.get("DIAG_ION", "I"), instrument="harps",
                  holding="solar_harps_molecfit_corrected", lo_A=3782.6, hi_A=6910.0,
                  engine_deck="ts-lte", method="synthesis", pool="set:" + os.environ.get("DIAG_SET", "CA_AGSS21"))
step = [s for s in resolve(d, interpreter=sys.executable,
                           ispec_dir=os.environ["ISPEC_DIR"]).steps
        if s["name"] == "derive_products"][0]
args = list(step["args"])
assert "--lines-from-set" in args, args
args[args.index("--out") + 1] = str(out)
print("[mu1-diag] args:", " ".join(args))

import pipeline.nearuv_synth as ns  # noqa: E402

_bsc = ns.build_solar_context


def bsc(*a, **k):
    ctx = _bsc(*a, **k)
    if mode == "mu1":
        ctx["vsini"] = 0.0
    if vmac is not None:
        ctx["macroturbulence"] = vmac
    print(f"[mu1-diag] ctx vsini={ctx['vsini']} vmac={ctx['macroturbulence']} "
          f"R={ctx['resolving_power']}")
    return ctx


ns.build_solar_context = bsc

import derive_band_products as dbp  # noqa: E402

sys.argv = ["derive_band_products.py"] + args
if mode == "flux":
    dbp.main()
    sys.exit(0)

import measure_band_ew as mbe  # noqa: E402
import pandas as pd  # noqa: E402

with gzip.open("/mnt/codex-data/solar_reference/delbouille_bass2000/"
               "delbouille_bass2000_visible.csv.gz", "rt") as fh:
    L = pd.read_csv(fh).to_numpy(float)
LW, LF = L[:, 0], L[:, 1]
spec = mbe.HoldingSpec(holding_id="solar_delbouille_liege", reader="liege_diag",
                       pre_normalised=True, span_A=(3000.0, 10000.0))


def sel(*a, **k):
    return spec


def lwx(instrument, centre, pad, segs=None, allow_uncorrected=False, holding=None):
    m = (LW >= centre - pad) & (LW <= centre + pad)
    return mbe.Window(LW[m], LF[m], "holding=solar_delbouille_liege BASS2000 mu=1 (DIAG)", spec)


mbe.select_holding, mbe.load_window_ex = sel, lwx
_rc = pd.read_csv


def rc(p, *a, **k):
    df = _rc(p, *a, **k)
    if str(p).endswith("instrument_catalog.csv"):
        df.loc[df.iloc[:, 0].astype(str) == "harps", "resolving_power_max"] = 500000
    return df


pd.read_csv = rc
from ace_co_mu1_bracket import turbospectrum_intensity_at_mu  # noqa: E402

with turbospectrum_intensity_at_mu(1.0) as st:
    try:
        dbp.main()
    finally:
        print(f"[mu1-diag] intensity injections fired={st['fired']} remapped={st['remapped']}")
assert st["remapped"] > 0, "no intensity remap -> this would be flux"
