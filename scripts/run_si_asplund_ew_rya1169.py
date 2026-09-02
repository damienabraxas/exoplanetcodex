#!/usr/bin/env python3
"""Measure the exact RYA-1169 Asplund Si set on staged Solar atlases.

This is the observational run: source membership and published gf are inputs,
while EW and observed depth come from each holding.  No generic Si candidate is
substituted when a source line is unserved.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))

DEFAULT_KP = Path("/Users/ryanschmitt/codex/spectra/Kitt Peak Flux Atlas")
DEFAULT_K2005 = Path("/Users/ryanschmitt/codex/spectra/KITTS/irradrelwl.dat")
OUT = ROOT / "data/results/rya1169/asplund_ew_run"

# Scott et al. (2015), electronic Table 2, W_lambda in pm.  These are carried
# only as an external comparison; the EWs emitted by this run come from our
# atlas integration and are never replaced by these values.
SCOTT_EW_MA = {
    5645.6128: 35.0, 5684.4840: 63.7, 5690.4250: 52.6,
    5701.1040: 39.5, 5772.1460: 56.0, 5793.0727: 45.8,
    6741.6280: 16.3, 7034.9006: 74.0, 7226.2079: 38.7,
    6371.3714: 36.6,
}


def iag_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    from astropy.io import fits
    from pipeline.uv_conditioning import vac_to_air
    with fits.open(path) as hdul:
        d = hdul[1].data
        w = vac_to_air(1e8 / np.asarray(d["v"], float))
        f = np.asarray(d["s"], float)
    order = np.argsort(w)
    return w[order], f[order]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kp-dir", type=Path, default=DEFAULT_KP)
    ap.add_argument("--kurucz2005", type=Path, default=DEFAULT_K2005)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    os.environ["CODEX_KP_ATLAS"] = str(a.kp_dir)

    from measure_band_ew import (kp_segments, load_kp1984_composite_window,
                                 _read_kurucz2005_residual, window_half_width)
    from pipeline.band_products import equivalent_width

    src = pd.read_csv(ROOT / "data/results/rya1169/si_asplund_line_test.csv")
    src = src[src.asplund_grade == "asplund"].copy()
    ll = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)
    allw = ll.wavelength_air_A.astype(float).to_numpy()
    segs = kp_segments()
    rows = []
    for line in src.itertuples():
        centre = float(line.depth_source_wavelength_A)
        hw = window_half_width(allw, centre)
        for holding in ("kpno_1984_molecfit_corrected", "kpno_kurucz2005_residual"):
            status = "served"; reason = ""; ew = depth = np.nan; prov = ""
            try:
                if holding == "kpno_1984_molecfit_corrected":
                    w, f, prov = load_kp1984_composite_window(centre, hw * 3, segs)
                else:
                    w, f = _read_kurucz2005_residual(a.kurucz2005,
                                                     centre - hw * 3, centre + hw * 3)
                    if w.size < 12:
                        raise LookupError(f"Kurucz 2005 has only {w.size} pixels in window")
                    prov = a.kurucz2005.name
                ew, method, concern = equivalent_width(w, f, centre, hw, pre_normalised=True)
                core = np.abs(w - centre) <= min(hw, 0.10)
                depth = float(1 - np.nanmin(f[core]))
                if concern: reason = concern
            except Exception as exc:
                status = "unserved"; reason = f"{type(exc).__name__}: {exc}"
                method = ""
            rows.append({
                "holding": holding, "species": line.species,
                "wavelength_air_A": centre, "ep_eV": float(line.ep_eV),
                "published_loggf": float(line.published_loggf),
                "source_band": line.source_band_recomputed,
                "predicted_feature_depth": float(line.solar_feature_depth),
                "observed_core_depth": round(depth, 5) if np.isfinite(depth) else "",
                "ew_mA": round(ew, 4) if np.isfinite(ew) else "",
                "scott2015_published_ew_mA": SCOTT_EW_MA[centre],
                "ew_ratio_ours_vs_scott": round(ew / SCOTT_EW_MA[centre], 5) if np.isfinite(ew) else "",
                "half_width_A": round(hw, 4), "status": status,
                "reason_or_concern": reason, "method": method, "spectrum_source": prov,
            })
    a.out.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame(rows)
    out.to_csv(a.out / "si_asplund_ew_per_line.csv", index=False)
    holdings = []
    for name, g in out.groupby("holding", sort=False):
        served = g[g.status == "served"]
        holdings.append({"holding": name, "n_source_lines": len(g), "n_served": len(served),
                         "n_unserved": int((g.status != "served").sum()),
                         "median_ew_mA": round(float(pd.to_numeric(served.ew_mA).median()), 4) if len(served) else None,
                         "median_ew_ratio_vs_scott": round(float(pd.to_numeric(served.ew_ratio_ours_vs_scott).median()), 5) if len(served) else None,
                         "median_observed_depth": round(float(pd.to_numeric(served.observed_core_depth).median()), 5) if len(served) else None,
                         "unserved": g.loc[g.status != "served", ["wavelength_air_A", "reason_or_concern"]].to_dict("records")})
    doc = {"ticket":"RYA-1169", "line_set":"asplund", "run_type":"observed EW/depth",
           "n_lines":len(src), "holdings":holdings,
           "raw_1984_policy":"QUARANTINED; never loaded by this runner",
           "abundance_inversion":"PENDING: requires source-gf override through synthesis/EW engine",
           "firewall":"no missing line substituted; distance from A(Si)=7.51 was not used"}
    (a.out / "si_asplund_ew_summary.json").write_text(json.dumps(doc, indent=2)+"\n")
    print(json.dumps(doc, indent=2))


if __name__ == "__main__": main()
