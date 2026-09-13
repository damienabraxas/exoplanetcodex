#!/usr/bin/env python3
"""Audit published Si J-band lines against every registered CRIRES holding.

This deliberately records coverage and conditioning status only.  Raw Vesta
IDPs are never promoted to measurements, and the existing corrected Y/H
products cannot answer a J-band line request.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HOLDINGS = (
    "solar_crires_plus_y_rya794",
    "solar_crires_plus_y_wide_rya1054",
    "solar_crires_plus_h_rya1094",
    "solar_vesta_crires_plus_idp",
)


def main() -> None:
    from measure_band_ew import load_window_ex

    pool = pd.read_csv(ROOT / "data/audit/rya1218_si_protocol/si_nir_line_pool.csv")
    rows = []
    for line in pool.itertuples():
        for holding in HOLDINGS:
            row = {
                "holding": holding,
                "source": line.source,
                "wavelength_air_A": float(line.wavelength_air_A),
                "ep_eV": float(line.ep_eV),
                "loggf": float(line.loggf),
                "band": line.band,
                "status": "HOLD",
                "conditioning": "",
                "reason": "",
                "provenance": "",
            }
            try:
                win = load_window_ex("crires_plus", float(line.wavelength_air_A), 1.0,
                                     holding=holding)
                row.update({
                    "status": "REACHABLE",
                    "conditioning": "corrected_normalized" if win.pre_normalised else "uncorrected_or_unnormalized",
                    "provenance": win.provenance,
                    "reason": "coverage reached; EW still withheld pending a validated J-band corrected product"
                    if holding != "solar_vesta_crires_plus_idp" else
                    "raw IDP coverage only; telluric correction and rest-frame conditioning are not applied",
                })
            except Exception as exc:
                row["reason"] = f"{type(exc).__name__}: {exc}"
            rows.append(row)
    out = ROOT / "data/results/rya1218/si_crires_nir"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "si_j_band_crires_coverage.csv", index=False)
    summary = []
    for holding, group in df.groupby("holding", sort=False):
        summary.append({
            "holding": holding,
            "n_lines": len(group),
            "n_reachable": int((group.status == "REACHABLE").sum()),
            "n_hold": int((group.status == "HOLD").sum()),
            "disposition": "HOLD_MEASUREMENT" if holding == "solar_vesta_crires_plus_idp" or
            not (group.status == "REACHABLE").all() else "HOLD_NO_J_COVERAGE",
            "reasons": sorted(set(group.reason)),
        })
    payload = {
        "ticket": "RYA-1218",
        "line_pool": "Bergemann et al. 2013 Table 1 J-band Si I lines",
        "measurement_policy": "No EW or abundance is reported from raw/topocentric IDPs or without a corrected J holding.",
        "holdings": summary,
    }
    (out / "si_j_band_crires_coverage.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
