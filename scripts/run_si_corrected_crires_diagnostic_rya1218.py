#!/usr/bin/env python3
"""Observed Si EW/depth diagnostics on the corrected CRIRES Y and H products.

The canonical Si census is still HOLD for laboratory/identity adjudication, so
these rows are measurements of the corrected spectra only; they are not an
abundance or grade pool.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HOLDINGS = (
    ("Y", "solar_crires_plus_y_wide_rya1054", 9800.0, 10796.0),
    ("H", "solar_crires_plus_h_rya1094", 15007.0, 17494.0),
)


def main() -> None:
    from measure_band_ew import load_window_ex, window_half_width
    from pipeline.band_products import equivalent_width

    census = pd.read_csv(ROOT / "data/audit/rya1218_si_protocol/canonical_si_census.csv",
                         low_memory=False)
    linelist = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)
    allw = linelist.wavelength_air_A.astype(float).to_numpy()
    rows = []
    for band, holding, lo, hi in HOLDINGS:
        lines = census[(census.wavelength_air_A >= lo) &
                       (census.wavelength_air_A <= hi)]
        for line in lines.itertuples():
            centre = float(line.wavelength_air_A)
            hw = window_half_width(allw, centre)
            row = {
                "band": band, "holding": holding, "species": line.species,
                "wavelength_air_A": centre,
                "ep_eV": float(line.excitation_potential_eV),
                "loggf": float(line.log_gf),
                "canonical_line_id": line.line_id,
                "measurement_eligibility": line.measurement_eligibility,
                "status": "HOLD", "ew_mA": "", "core_depth": "",
                "method": "", "provenance": "", "reason": "",
            }
            try:
                win = load_window_ex("crires_plus", centre, hw * 3.0,
                                     holding=holding)
                ew, method, concern = equivalent_width(
                    win.wave, win.flux, centre, hw, pre_normalised=win.pre_normalised)
                core = np.abs(win.wave - centre) <= min(hw, 0.10)
                row.update({"status": "DIAGNOSTIC_SERVED", "ew_mA": round(float(ew), 4),
                            "core_depth": round(float(1 - np.nanmin(win.flux[core])), 5),
                            "method": method, "provenance": win.provenance,
                            "reason": concern or "canonical line remains HOLD for identity/gf adjudication"})
            except Exception as exc:
                row["reason"] = f"{type(exc).__name__}: {exc}"
            rows.append(row)
    out = ROOT / "data/results/rya1218/si_crires_corrected_diagnostic"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "si_corrected_crires_per_line.csv", index=False)
    summary = []
    for (band, holding), g in df.groupby(["band", "holding"], sort=False):
        summary.append({"band": band, "holding": holding, "n_canonical_lines": len(g),
                        "n_served": int((g.status == "DIAGNOSTIC_SERVED").sum()),
                        "n_unserved": int((g.status != "DIAGNOSTIC_SERVED").sum()),
                        "species": sorted(g.species.unique().tolist()),
                        "abundance_status": "HOLD"})
    (out / "si_corrected_crires_summary.json").write_text(json.dumps({
        "ticket": "RYA-1218", "purpose": "corrected-spectrum Si diagnostic",
        "abundance_or_grade_pool": False, "holdings": summary}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
