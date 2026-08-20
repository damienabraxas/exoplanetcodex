#!/usr/bin/env python3
"""Extract and validate the 2023-08-02T15 La Silla GDAS1 profile for RYA-927."""
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

from pipeline.telluric.arl_gdas import (extract_profile, interpolate_at_height,
                                        write_molecfit_ascii)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("archive", type=Path)
    ap.add_argument("output_dir", type=Path)
    args = ap.parse_args()
    profile = extract_profile(args.archive, datetime(2023, 8, 2, 15), -29.26, -70.73)
    at_site = interpolate_at_height(profile, 2400.0)
    observed = {"pressure_hpa": 773.2, "temperature_c": 20.9,
                "relative_humidity_pct": 8.0}
    delta = {key: at_site[key] - observed[key] for key in observed}
    # Wide, predeclared identity tolerances: detect a wrong site/slot, not tune GDAS.
    limits = {"pressure_hpa": 15.0, "temperature_c": 5.0,
              "relative_humidity_pct": 10.0}
    if any(abs(delta[key]) > limits[key] for key in limits):
        raise SystemExit(f"GDAS-vs-ESO telemetry validation failed: {delta}")
    profile["observatory_elevation_validation"] = {
        "gdas_interpolated": at_site, "eso_harps_observed": observed,
        "delta_gdas_minus_observed": delta, "limits": limits, "verdict": "PASS",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_molecfit_ascii(profile, args.output_dir / "noaa_C-70.7-29.3_2023-08-02T15.txt")
    (args.output_dir / "noaa_C-70.7-29.3_2023-08-02T15.json").write_text(
        json.dumps(profile, indent=2) + "\n"
    )
    print(json.dumps({k: v for k, v in profile.items() if k != "levels"}, indent=2))


if __name__ == "__main__":
    main()
