#!/usr/bin/env python3
"""RYA-1214 — the K synthesis band: its line list, and the numbers its BandPolicy records.

K did not exist as a synthesis band. `band_policy` resolved 19500-24000 A to "NIR", whose
line list is RYA-762's curated 9203-12976 A set, so a K run would synthesise against a list
that contains no K line at all — the exact blocker RYA-1094 hit for H.

LINE LIST — the same builder and the same kind of input H used: `pipeline.nearuv_linelist.
build()` on our solar hfs-ON VALD extract `vald_solar_ir_17000_25000_hfson_raw.txt`
(17000.57-24997.99 A), cut to the K product's MEASURED extent (the RYA-1214 rest-frame J/K
conditioning wrote 19452.4-24845.6 A). Molecules are excluded by that builder by design and
reported; CO enters as a separate .bsyn list (RYA-1207), not here.

BAND NUMBERS — measured, as RYA-1094 measured H's, and on the same two kinds of material:
the line list (lines per A, median gap) and the product (continuum median and p95 over
pixels the co-add wrote). ⚠️ The K product is RYA-1219's per-segment normalisation after
molecfit, so these continuum figures describe that product, not the sky.

Must run where iSpec is importable (Sirius).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import nearuv_linelist as nl  # noqa: E402

RAW = ROOT / "data/linelists/vald_solar_ir_17000_25000_hfson_raw.txt"
PRODUCT = ROOT / "data/results/rya1214_crires_jk/solar_crires_plus_k_rya1219_rest.csv"
OUT = ROOT / "data/linelists/ispec_k_19452_24846/atomic_lines.tsv"
AUDIT = ROOT / "data/audit/rya1214_crires_jk/k_band_build.json"


def main() -> int:
    prod = pd.read_csv(PRODUCT)
    lo = float(prod.wavelength_air_A.min())
    hi = float(prod.wavelength_air_A.max())
    path, rep = nl.build(RAW, lo, hi, out_path=OUT)
    import ispec
    arr = ispec.read_atomic_linelist(str(OUT))
    w = np.sort(np.asarray(arr["wave_A"], float))
    flux = prod.flux_normalized.to_numpy(float)
    stats = {
        "lo_A": round(lo, 2), "hi_A": round(hi, 2),
        "n_lines": int(len(w)), "lines_per_A": round(len(w) / (hi - lo), 3),
        "median_gap_A": round(float(np.median(np.diff(w))), 3),
        "continuum_median": round(float(np.median(flux)), 3),
        "continuum_p95": round(float(np.percentile(flux, 95)), 3),
        "frac_points_below_0p5": round(float(np.mean(flux < 0.5)), 4),
        "product_pixels": int(len(prod)),
        "anchor_A": round(0.5 * (lo + hi), 1),
        # the H convention (synth_bands.yaml): 20.51 sigma_D at the anchor, sigma_D for
        # the same 1.7 km/s (thermal + micro) the other bands use
        "half_width_A": round(20.51 * 0.5 * (lo + hi) * 1.7 / 299792.458, 2),
    }
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT.write_text(json.dumps({"ticket": "RYA-1214", "source": str(RAW.relative_to(ROOT)),
                                 "product": str(PRODUCT.relative_to(ROOT)),
                                 "linelist": str(OUT.relative_to(ROOT)),
                                 "builder_report": {k: v for k, v in rep.items()
                                                    if isinstance(v, (int, float, str, dict))},
                                 "band": stats}, indent=2, default=str) + "\n")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
