#!/usr/bin/env python3
"""Measure the RYA-1169 Si source lines on every registered solar holding.

This is an observation/holding matrix only. It does not infer abundances or
select an engine. Every row names the requested holding and records a refusal
when that holding cannot serve the exact window.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

HOLDINGS = (
    ("harps", "solar_harps"),
    ("harps", "solar_harps_molecfit_corrected"),
    ("iag_fts_solar_atlas", "solar_iag"),
    ("iag_fts_solar_atlas", "solar_iag_reiners2016"),
    ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
    ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
    ("crires_plus", "solar_crires_plus_y_rya794"),
    ("crires_plus", "solar_crires_plus_y_wide_rya1054"),
    ("crires_plus", "solar_crires_plus_h_rya1094"),
    ("crires_plus", "solar_vesta_crires_plus_idp"),
)


def main() -> None:
    from measure_band_ew import load_window_ex, window_half_width
    from pipeline.band_products import equivalent_width

    src = pd.read_csv(ROOT / "data/results/rya1169/si_asplund_line_test.csv")
    src = src[src.asplund_grade == "asplund"].copy()
    ll = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)
    allw = ll.wavelength_air_A.astype(float).to_numpy()
    rows = []
    for line in src.itertuples():
        centre = float(line.depth_source_wavelength_A)
        hw = window_half_width(allw, centre)
        for instrument, holding in HOLDINGS:
            row = {
                "instrument": instrument, "holding": holding,
                "species": line.species, "wavelength_air_A": centre,
                "ep_eV": float(line.ep_eV), "published_loggf": float(line.published_loggf),
                "source_band": line.source_band_recomputed, "half_width_A": round(hw, 4),
                "status": "unserved", "ew_mA": "", "observed_core_depth": "",
                "method": "", "spectrum_source": "", "reason_or_concern": "",
            }
            try:
                win = load_window_ex(instrument, centre, hw * 3.0, holding=holding)
                ew, method, concern = equivalent_width(
                    win.wave, win.flux, centre, hw, pre_normalised=win.pre_normalised)
                core = np.abs(win.wave - centre) <= min(hw, 0.10)
                row.update({
                    "status": "served", "ew_mA": round(float(ew), 4),
                    "observed_core_depth": round(float(1 - np.nanmin(win.flux[core])), 5),
                    "method": method, "spectrum_source": win.provenance,
                    "reason_or_concern": concern or "",
                })
            except Exception as exc:
                row["reason_or_concern"] = f"{type(exc).__name__}: {exc}"
            rows.append(row)
    out = ROOT / "data/results/rya1218/si_multiholding_ew"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "si_multiholding_ew_per_line.csv", index=False)
    summary = []
    for (instrument, holding), g in df.groupby(["instrument", "holding"], sort=False):
        served = g[g.status == "served"]
        summary.append({"instrument": instrument, "holding": holding,
                        "n_source_lines": len(g), "n_served": len(served),
                        "n_unserved": int((g.status != "served").sum()),
                        "served_wavelengths": [round(float(x), 4) for x in served.wavelength_air_A],
                        "unserved": g.loc[g.status != "served",
                                           ["wavelength_air_A", "reason_or_concern"]].to_dict("records")})
    (out / "si_multiholding_ew_summary.json").write_text(json.dumps({
        "ticket": "RYA-1218", "line_set": "RYA-1169 Asplund source lines",
        "run_type": "observed EW/depth; no abundance or engine selection", "holdings": summary
    }, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
