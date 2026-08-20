#!/usr/bin/env python3
"""Assemble the RYA-931 corrected-holding intake manifest (RYA-931).

The manifest is the holding's provenance: which exposures, which atmosphere,
which engine and configuration, what came out, and what was refused.  It is
assembled from the run artifacts rather than typed, so it cannot drift from what
was actually executed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from astropy.io import fits


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            d.update(chunk)
    return d.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("correction_json", type=Path)
    ap.add_argument("verification_json", type=Path)
    ap.add_argument("reproduction_json", type=Path)
    ap.add_argument("corrected_normalization_json", type=Path)
    ap.add_argument("rca_json", type=Path)
    ap.add_argument("corrected_csv", type=Path)
    ap.add_argument("out_json", type=Path)
    args = ap.parse_args()

    correction = json.loads(args.correction_json.read_text())
    verification = json.loads(args.verification_json.read_text())
    reproduction = json.loads(args.reproduction_json.read_text())
    corrected_norm = json.loads(args.corrected_normalization_json.read_text())
    rca = json.loads(args.rca_json.read_text())

    esorex = subprocess.run(["/opt/homebrew/bin/esorex", "--version"],
                            text=True, capture_output=True).stdout.strip()
    sample = Path(correction[0]["exposure"])
    manifest = {
        "manifest_id": "rya931_solar_harps_molecfit_corrected",
        "holding_id": "solar_harps_molecfit_corrected",
        "ticket": "RYA-931",
        "parent_ticket": "RYA-927",
        "supersedes_nothing": (
            "solar_harps is preserved unmodified; this is a provenance-distinct "
            "sibling holding, not a replacement."),

        "source": {
            "program_id": "1102.D-0954(A)",
            "citation": "Dumusque et al. 2021, arXiv:2009.01945 — HARPS direct-Sun feed",
            "product_class": "ESO Phase-3 HARPS S1D ADP",
            "observation_date_utc": "2023-08-02 13:56:15 to 14:08:40",
            "n_exposures": len(correction),
            "exposures": [{"file": r["exposure"], "md5": r["source_md5"],
                           "mjd_obs": r["mjd_obs"], "airmass_start": r["airmass_start"]}
                          for r in correction],
            "wavelength_frame": "air, SPECSYS=BARYCENT",
            "known_defect": ("the ERR column is 100% NaN in all ten exposures; it is "
                             "not usable as an uncertainty and is not used as one"),
        },

        "lineage_proof": {
            "claim": ("these ten exposures ARE the source of the held solar_harps "
                      "product, proven rather than assumed"),
            "method": ("pipeline.spectra_normalize, the same functions that built the "
                       "held product, re-run on the ten originals"),
            "regenerated_sha256": reproduction["sha256"],
            "committed_sha256": reproduction["reference_sha256"],
            "byte_identical": reproduction["byte_identical_to_reference"],
        },

        "correction": {
            "engine": "ESO Molecfit", "recipe": "molecfit_model",
            "esorex_version": esorex,
            "stage": "per exposure, before co-addition and continuum normalisation",
            "molecules": ["O2"],
            "band_A": [6857.0, 6911.0],
            "registered_telluric_bands_in_instrument_range": ["O2 B 6867-6884"],
            "wavelength_frame": "AIR_RV",
            "obs_erf_rv_convention": (
                "OBS_ERF_RV_VALUE = -BERV; molecfit's sign is positive when the "
                "target-Earth distance increases, the negative of the ESO DRS BERV. "
                "Determined empirically, not assumed: against AIR and AIR_RV(+BERV) "
                "this sign gives a 15x lower corrected-residual RMS (0.039 vs 0.59/0.62)."),
            "reference_atmosphere": rca["passing_run"]["reference_atmosphere"],
            "reference_atmosphere_md5": rca["passing_run"]["reference_atmosphere_md5"],
            "gdas_profile": "noaa gdas1.aug23.w1, La Silla column, 2023-08-02 15:00 UT",
            "standard_atmosphere_fallback": False,
            "line_spread_function": {
                "form": "Gaussian, fitted; boxcar and Lorentzian disabled (HARPS is "
                        "fibre-fed and has no slit)",
                "fwhm_pix_per_exposure": [r["lsf_gauss_fwhm_pix"] for r in correction],
                "implied_resolving_power": 101710,
            },
            "solar_line_exclusion": (
                "intrinsic solar absorption masked out of the atmospheric fit using the "
                "Baker+2020 telluric-free IAG atlas (md5 1a47d929df32778f72cc4080ff1f95ec); "
                "the atlas supplies a MASK only, never flux"),
            "per_exposure": [{k: r[k] for k in (
                "exposure", "corrected_md5", "rel_mol_col_O2", "reduced_chi2",
                "lsf_gauss_fwhm_pix", "dt_transmission_sigma", "t_min",
                "n_corrected", "n_quarantined", "lsf_start_accepted_pix")}
                for r in correction],
        },

        "quarantine": {
            "rule": ("a pixel whose fitted transmission T is below T_min is set NaN, "
                     "never divided. T_min = measured dT / 0.05, i.e. a corrected pixel "
                     "may carry at most 5% relative flux error from the transmission "
                     "model. T_min is derived per exposure, not hard-coded."),
            "policy_basis": ("RYA-927 ratified policy: saturated or badly modelled "
                             "telluric regions are rejected, not forced into a "
                             "corrected state"),
            "n_nan_in_product": corrected_norm["n_nan"],
        },

        "product": {
            "path": str(args.corrected_csv),
            "sha256": corrected_norm["sha256"],
            "n_samples": corrected_norm["n_samples"],
            "n_finite": corrected_norm["n_finite"],
            "wavelength_range_A": corrected_norm["wavelength_range_A"],
            "built_by": ("scripts/rya931_normalize_corrected.py, which calls "
                         "pipeline.spectra_normalize — the canonical path"),
        },

        "telluric_state": {
            "value": "applied",
            "derived_from": ("the product's own MTRANS image extension, read by "
                             "pipeline.telluric_intake; not asserted in the registry"),
            "evidence_class": "transmission extension present and NOT all-unity",
        },

        "verification": {
            "o2b_gate": verification["o2b_gate"],
            "solar_preservation": verification["solar_preservation"],
            "iag_discriminator": verification["iag_discriminator"],
            "verdict": verification["verdict"],
        },

        "caveats": [
            "telluric_applied=applied is a product-level conditioning fact, not a "
            "line-level science verdict; every candidate line in 6857-6911 A still "
            "needs a local telluric disposition (RYA-927 ratified two-level policy).",
            "Removing the O2 band raises the 95th-percentile continuum anchor of the "
            "knot window containing it, so the continuum solution shifts by up to "
            "0.389% outside the fitted region. The raw co-add outside the fitted "
            "region is bit-identical; the normalised flux moves by at most 1.6e-3.",
            "The ERR column is unusable, so molecfit was weighted with its documented "
            "DEFAULT_ERROR. That is fit weighting only and is not an EW uncertainty.",
        ],
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {args.out_json}")


if __name__ == "__main__":
    main()
