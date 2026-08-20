#!/usr/bin/env python3
"""Compare corrected HARPS EWs to the frozen website baseline (RYA-927).

The baseline is selected before correction: all seven website-published lines inside
the O2 B band.  This tool never promotes or overwrites either EW product.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "data/audit/rya927_harps_telluric/website_ew_baseline.csv"
KEY = ["element", "ion", "wavelength_air_A"]


def compare(corrected_path: Path) -> pd.DataFrame:
    old = pd.read_csv(BASELINE, comment="#")
    new = pd.read_csv(corrected_path, comment="#")
    required = set(KEY + ["ew_mA", "ew_err_mA"])
    missing = required - set(new.columns)
    if missing:
        raise ValueError(f"corrected EW product lacks columns: {sorted(missing)}")
    joined = old.merge(new, on=KEY, how="left", suffixes=("_uncorrected", "_corrected"),
                       validate="one_to_one", indicator=True)
    if not (joined["_merge"] == "both").all():
        absent = joined.loc[joined["_merge"] != "both", KEY]
        raise ValueError(f"corrected product omitted predeclared O2-B lines:\n{absent}")
    joined["delta_ew_mA"] = joined["ew_mA_corrected"] - joined["ew_mA_uncorrected"]
    joined["ratio_corrected_to_uncorrected"] = (
        joined["ew_mA_corrected"] / joined["ew_mA_uncorrected"]
    )
    return joined.drop(columns="_merge")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("corrected_ews", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    result = compare(args.corrected_ews)
    print(result[KEY + ["ew_mA_uncorrected", "ew_mA_corrected", "delta_ew_mA"]]
          .to_string(index=False))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
