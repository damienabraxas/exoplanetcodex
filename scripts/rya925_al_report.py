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
    d["display_name"] = [
        axes_for(r.treatment, handler=r.handler, model=r.model).display
        for r in d.itertuples()
    ]
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
    d.loc[d.band.eq("near-UV"), "domain_status"] = (
        "in-domain LTE; n=1 constrained of 3; no registered NLTE synthesis deck")
    cols = ["band", "route", "scale", "model", "atmos", "gf", "display_name",
            "n_lines", "A", "stat_dex", "syst_dex", "delta_same_line_LTE",
            "literature_comparison", "domain_status"]
    return d[cols].rename(columns={"n_lines": "n", "stat_dex": "sigma_stat",
                                   "syst_dex": "sigma_syst"})


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
        rows.append(dict(band="near-UV", wavelength_air_A=r.wavelength_air_A,
                         treatment=r.treatment, route="synth", abundance=r.abundance,
                         in_aggregate=r.in_aggregate, excluded_reason=r.excluded_reason,
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
    return validate_element("Al", inputs, asplund=6.43).to_dict()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    prod = product_matrix()
    prod.to_csv(OUT / "product_matrix.csv", index=False)
    per_line_matrix().to_csv(OUT / "per_line_matrix.csv", index=False)
    disposition().to_csv(OUT / "line_disposition.csv", index=False)
    (OUT / "validation.json").write_text(json.dumps(validation(prod), indent=2) + "\n")
    print(prod.to_string(index=False))


if __name__ == "__main__":
    main()
