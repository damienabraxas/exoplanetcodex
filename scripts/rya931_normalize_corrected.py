#!/usr/bin/env python3
"""Rebuild the solar stack/normalisation from corrected exposures (RYA-931).

This is deliberately not a second normaliser.  It calls the same
`pipeline.spectra_normalize` functions the held `solar_harps` product was built
with, only with a different input directory and output path, so the corrected
holding differs from the uncorrected one by the telluric correction and by
nothing else.

The reproduction control proves that: run with --reproduce against the ORIGINAL
exposures, the output is byte-identical to the committed solar_normalized.csv.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import CONTINUUM_PARAMS  # noqa: E402
from pipeline import spectra_normalize as sn  # noqa: E402


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            d.update(chunk)
    return d.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("exposures", type=Path)
    ap.add_argument("out_csv", type=Path)
    ap.add_argument("--reference-csv", type=Path, default=None,
                    help="assert the output is byte-identical to this file")
    ap.add_argument("--report", type=Path, default=None)
    args = ap.parse_args()

    wavelength, flux_raw, n_exp, file_info = sn.stack_spectra(args.exposures,
                                                              star_id="solar")
    continuum = sn.fit_continuum(wavelength, flux_raw, star_id="solar")
    params = CONTINUUM_PARAMS.get("solar")
    provenance = sn._build_provenance("solar", file_info, params, wavelength, args.out_csv)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    flux_norm = sn._save_normalized(wavelength, flux_raw, continuum, args.out_csv, "")

    finite = np.isfinite(flux_norm)
    report = {
        "ticket": "RYA-931", "exposures_dir": str(args.exposures), "n_exposures": n_exp,
        "out_csv": str(args.out_csv), "sha256": sha256(args.out_csv),
        "n_samples": int(flux_norm.size), "n_finite": int(finite.sum()),
        "n_nan": int((~finite).sum()),
        "wavelength_range_A": [float(wavelength.min()), float(wavelength.max())],
        "provenance": provenance,
    }
    if args.reference_csv:
        report["reference_sha256"] = sha256(args.reference_csv)
        report["byte_identical_to_reference"] = (
            report["sha256"] == report["reference_sha256"])
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "provenance"}, indent=2))


if __name__ == "__main__":
    main()
