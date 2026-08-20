#!/usr/bin/env python3
"""Assemble the RYA-925 Kitt Peak Al evidence without reading frozen gold."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.treatment_axes import axes_for
from pipeline.validate_element import BandProduct, validate_element

OUT = ROOT / "data/results/rya925"


def _product_inputs():
    return sorted((OUT / "vis").glob("*products.csv")) + sorted(
        (OUT / "red_optical").glob("*products.csv")) + sorted(
        OUT.glob("*3000_3800*products.csv"))


def product_matrix() -> pd.DataFrame:
    frames = [pd.read_csv(p) for p in _product_inputs()]
    d = pd.concat(frames, ignore_index=True)
    # RYA-925 follow-up: the sole near-UV candidate does not measure Al. Across
    # A(Al)=1.43..6.43 the peak synthetic response is 2.7e-5 while the window RMS
    # residual is 0.046; its apparent minimum is numerical structure in a blended,
    # model-inadequate window. Keep the attempted line in the appendix, never emit a
    # one-line abundance product from it.
    d = d[~d["band"].eq("near-UV")].copy()
    d["display_name"] = [
        axes_for(r.treatment, handler=r.handler, model=r.model).display
        for r in d.itertuples()
    ]
    d["holding"] = "solar_kpno"
    d["instrument"] = "Kitt Peak Solar Flux Atlas (NSO/KPNO FTS)"
    d["coverage_A"] = "2960-13000"
    d["delta_same_line_LTE"] = 0.0
    # Amarsi coverage is a subset; compare only with the LTE values for those same lines.
    same = {("VIS", "ENGINE-A"): -0.0223,
            ("red-optical", "ENGINE-A"): -0.0282,
            ("VIS", "ENGINE-B"): -0.2598,
            ("red-optical", "ENGINE-B"): -0.0827}
    for i, r in d.iterrows():
        d.loc[i, "delta_same_line_LTE"] = same.get((r.band, r.treatment), 0.0)
    d["literature_comparison"] = np.where(
        d.band.eq("VIS"), "validate against RYA-716 6.43 [6.39,6.47]", "report-and-note")
    d["domain_status"] = "in-domain"
    cols = ["holding", "instrument", "coverage_A", "band", "route", "scale", "model", "atmos", "gf", "display_name",
            "n_lines", "A", "stat_dex", "syst_dex", "delta_same_line_LTE",
            "literature_comparison", "domain_status"]
    return d[cols].rename(columns={"n_lines": "n", "stat_dex": "sigma_stat",
                                   "syst_dex": "sigma_syst"})


def calibration_matrix(products: pd.DataFrame) -> pd.DataFrame:
    """Replication gate: agreement with the external value within combined 1-sigma.

    Keep the two-part internal bars in product_matrix; this projection combines them
    only for the explicitly named calibration comparison.
    """
    d = products.copy()
    d["literature_value"] = 6.43
    d["literature_sigma"] = 0.04
    d["literature_scale"] = "3D-NLTE"
    d["delta_dex"] = d["A"] - d["literature_value"]
    d["sigma_combined"] = np.sqrt(
        d["sigma_stat"] ** 2 + d["sigma_syst"] ** 2 + d["literature_sigma"] ** 2)
    d["z_combined"] = d["delta_dex"].abs() / d["sigma_combined"]
    d["calibration_verdict"] = np.where(
        d["z_combined"] <= 1.0, "REPLICATED-WITHIN-1SIGMA", "NOT-REPLICATED")
    d["comparison_caveat"] = np.where(
        d["scale"].eq("3D-NLTE"), "scale-matched",
        "diagnostic cross-scale comparison; 3D-NLTE route not yet executed")
    return d[["holding", "instrument", "coverage_A", "band", "display_name", "scale",
              "n", "A", "sigma_stat", "sigma_syst", "literature_value",
              "literature_sigma", "literature_scale", "delta_dex", "sigma_combined",
              "z_combined", "calibration_verdict", "comparison_caveat"]]


def per_line_matrix() -> pd.DataFrame:
    rows = []
    for folder in (OUT / "vis", OUT / "red_optical"):
        for p in sorted(folder.glob("*lines.csv")):
            d = pd.read_csv(p)
            for r in d.itertuples():
                rows.append(dict(
                    band="VIS" if folder.name == "vis" else "red-optical",
                    wavelength_air_A=r.wavelength_air_A, treatment=r.treatment,
                    route="synth" if r.treatment == "ENGINE-B" else "ew",
                    abundance=r.abundance, in_aggregate=r.in_aggregate,
                    excluded_reason=getattr(r, "excluded_reason", ""),
                    nlte_delta_dex=getattr(r, "nlte_delta_dex", np.nan)))
    p = next(OUT.glob("*3000_3800*1D-LTE_lines.csv"))
    for r in pd.read_csv(p).itertuples():
        reason = r.excluded_reason
        in_aggregate = r.in_aggregate
        if abs(float(r.wavelength_air_A) - 3057.144) < 0.01:
            in_aggregate = False
            reason = ("ABUNDANCE-INSENSITIVE / BLENDED MODEL FAILURE: 146 catalogue "
                      "transitions in +/-0.72 A; changing A(Al)=1.43..6.43 moves peak "
                      "synthetic flux only 2.7e-5 while best-window RMS residual is "
                      "0.046. Numerical minimum does not measure Al; appendix only.")
        rows.append(dict(band="near-UV", wavelength_air_A=r.wavelength_air_A,
                         treatment=r.treatment, route="synth", abundance=r.abundance,
                         in_aggregate=in_aggregate, excluded_reason=reason,
                         nlte_delta_dex=0.0))
    return pd.DataFrame(rows)


def disposition() -> pd.DataFrame:
    rows = []
    for band, p in (
        ("VIS", ROOT / "data/measured/band_ew/AlI_3800_6910_kpno_solar_atlas_PROFILEFIT_ew.csv"),
        ("red-optical", ROOT / "data/measured/band_ew/AlI_6910_9199_kpno_solar_atlas_PROFILEFIT_ew.csv")):
        for r in pd.read_csv(p).itertuples():
            rows.append(dict(lambda_A=r.wavelength_air_A, band=band, holding="solar_kpno",
                gf_source_grade="canonical/Kurucz; UNGRADED except RYA-835 line-specific evidence",
                permitted_routes="ew|synth", model_domain="LTE; Amarsi per-line gated",
                disposition="measured" if r.in_aggregate else "quarantined",
                reason="accepted by current profile/feature gates" if r.in_aggregate else r.excluded_reason))
    rows += [
        dict(lambda_A=3057.144, band="near-UV", holding="solar_kpno",
             gf_source_grade="VALD3; ungraded", permitted_routes="synth only",
             model_domain="1D-LTE attempted; abundance-insensitive",
             disposition="rejected",
             reason="146 transitions in +/-0.72 A; Al response 2.7e-5 vs 0.046 RMS residual; numerical minimum does not measure Al; diagnostic plot banked"),
        dict(lambda_A=3944.006, band="VIS", holding="solar_kpno",
             gf_source_grade="Nordlander-Lind/NIST context", permitted_routes="synth",
             model_domain="model line exists", disposition="skipped",
             reason="known CH blend; RYA-716 says do not attempt"),
        dict(lambda_A=3961.520, band="VIS", holding="solar_kpno",
             gf_source_grade="Nordlander-Lind", permitted_routes="synth",
             model_domain="Gerber/mean3D level exists; production Al path not registered",
             disposition="skipped",
             reason="resonance line requires dedicated blend-aware optical synthesis; generic VIS intake is EW-keyed"),
        dict(lambda_A=10872.973, band="NIR", holding="solar_kpno",
             gf_source_grade="VALD3; theoretical depth 0.061", permitted_routes="synth",
             model_domain="LTE line present", disposition="skipped",
             reason="below predeclared synthesis theoretical-depth floor 0.15"),
        dict(lambda_A=11254.926, band="NIR", holding="solar_kpno",
             gf_source_grade="VALD3", permitted_routes="synth",
             model_domain="not evaluated", disposition="quarantined",
             reason="inside registered H2O 11120-11560 A telluric interval"),
    ]
    return pd.DataFrame(rows).sort_values("lambda_A")


def validation(products: pd.DataFrame) -> dict:
    inputs = [BandProduct(r.band, r.display_name, r.A, r.sigma_stat, r.sigma_syst,
                          r.n, r.scale) for r in products.itertuples()]
    legacy = validate_element("Al", inputs, asplund=6.43).to_dict()
    cal = calibration_matrix(products)
    science = cal[cal.band.isin(["VIS", "red-optical"])]
    legacy["classifier_scope"] = (
        "legacy central-value interval diagnostic; not the calibration verdict")
    legacy["calibration_verdict"] = (
        "REPLICATED-WITHIN-1SIGMA" if science.calibration_verdict.eq(
            "REPLICATED-WITHIN-1SIGMA").all() else "NOT-REPLICATED")
    legacy["calibration_gate"] = (
        "abs(Codex-literature) / sqrt(stat^2+syst^2+literature_sigma^2) <= 1")
    legacy["calibration_products_passed"] = int(
        science.calibration_verdict.eq("REPLICATED-WITHIN-1SIGMA").sum())
    legacy["calibration_products_total"] = int(len(science))
    return legacy


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    prod = product_matrix()
    prod.to_csv(OUT / "product_matrix.csv", index=False)
    cal = calibration_matrix(prod)
    cal.to_csv(OUT / "calibration_matrix.csv", index=False)
    tracker_rows = cal.where(pd.notna(cal), None).to_dict(orient="records")
    (OUT / "tracker_data.json").write_text(json.dumps({
        "updated": "2026-08-19",
        "literature": {"value": 6.43, "sigma": 0.04, "scale": "3D-NLTE"},
        "planned_instruments": [
            {"id": "solar_kpno", "label": "Kitt Peak FTS", "coverage": "2960–13000 Å",
             "readiness": "measured", "source": "RYA-925"},
            {"id": "solar_harps", "label": "HARPS solar", "coverage": "3780–6910 Å (full VIS)",
             "readiness": "loader ready; Al run owed", "source": "RYA-911 / PR #313"},
            {"id": "solar_iag", "label": "IAG FTS", "coverage": "4047–10650 Å",
             "readiness": "holding ready; Al run owed", "source": "RYA-708"},
            {"id": "solar_crires_plus", "label": "CRIRES+ solar/Vesta", "coverage": "held Y: 10280–10680 Å",
             "readiness": "holding ready; Al line coverage gap", "source": "RYA-794/805"},
        ],
        "bands": ["near-UV", "VIS", "red-optical", "NIR"],
        "products": tracker_rows,
        "open_models": [
            "Gerber 1D-NLTE — Al abundance-axis adapter not wired",
            "Nordlander–Lind mean-3D NLTE — verified grid not held/registered",
            "3D-NLTE — scale-matched production route owed",
        ],
    }, indent=2) + "\n")
    per_line_matrix().to_csv(OUT / "per_line_matrix.csv", index=False)
    disposition().to_csv(OUT / "line_disposition.csv", index=False)
    (OUT / "validation.json").write_text(json.dumps(validation(prod), indent=2) + "\n")
    print(prod.to_string(index=False))


if __name__ == "__main__":
    main()
