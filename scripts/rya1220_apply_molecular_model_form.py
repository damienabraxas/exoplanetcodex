#!/usr/bin/env python3
"""
RYA-1220 -- apply the measured 1D->3D model-form term to every molecular CNO route.

Closes `3d_nlte_model: HOLD_MOLECULAR`, which the gate carried on every route while the
repo already held the measurement (Amarsi et al. 2021 Table 2, five model atmospheres).

The term is a BOUND on the systematic we carry by synthesizing in 1D, not a correction:
no abundance moves (RYA-161). It is admitted only where the route's window lies inside
the transitions the reference analysis actually used -- borrowing a CN A-X (0-0) NIR
term for the CN red system at 6125-6200 A would be a real measurement of the wrong
population.

This does NOT promote any route. ⚠️ And it does not judge the fits on red_chi2: that
guard was tried and refuted on this data (fit_validity.py -- 1262 of 1366 in-aggregate
lines flagged, because the HARPS ERR column is all NaN and the absolute chi2 is
meaningless). red_chi2 is reported here only because the route files carry it.
See cno_hold_readjudication.json for what actually discriminates.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from pipeline.molecular_identity import identity_for
from pipeline.molecular_model_form import model_form_term

#: route directory -> diagnostic key, instrument, band
ROUTES = [
    ("cn_crires_j_continuum_v1", "CN_AX_J", "CRIRES+", "J"),
    ("cn_iag_3d_anchor_v1", "CN_AX_IR", "IAG", "NIR"),
    ("cn_kp_continuum_v8", "CN_AX_IR", "Kitt Peak", "NIR"),
    ("cn_harps_continuum_v2", "CN_red", "HARPS", "VIS"),
    ("cn_nearuv_continuum_v1", "NH_AX", "Kitt Peak", "near-UV"),
    ("oh_crires_h_continuum_v1", "OH_H", "CRIRES+", "H"),
    ("co_crires_k_continuum_v1", "CO_K", "CRIRES+", "K"),
]


def red_chi2_values(node, out=None):
    out = [] if out is None else out
    if isinstance(node, dict):
        if isinstance(node.get("red_chi2"), (int, float)):
            out.append(node["red_chi2"])
        for v in node.values():
            red_chi2_values(v, out)
    elif isinstance(node, list):
        for v in node:
            red_chi2_values(v, out)
    return out


def build(root: pathlib.Path) -> dict:
    rows = []
    for route, key, instrument, band in ROUTES:
        d = root / route
        if not (d / "provenance.json").exists():
            continue
        prov = json.loads((d / "provenance.json").read_text())
        fits = json.loads((d / "fits.json").read_text())
        windows = prov["pool"]["windows_air_A"]
        lo, hi = min(w[0] for w in windows), max(w[1] for w in windows)
        ident = identity_for(key)
        term = model_form_term(ident.molecule, lo, hi)
        chis = red_chi2_values(fits)
        rows.append({
            "route": route, "diagnostic": key, "instrument": instrument, "band": band,
            "molecule": ident.molecule, "system": ident.system, "vib_band": ident.band,
            "element": ident.element, "window_air_A": [lo, hi],
            "red_chi2_max_UNCALIBRATED": round(max(chis), 2) if chis else None,
            "in_reference_band": term["state"] == "MEASURED",
            "model_atmosphere_state": term["state"],
            "model_atmosphere_sigma_dex": term["sigma_dex"],
            "model_atmosphere_reason": term["reason"],
        })
    return {
        "schema": "rya1220.molecular_route_model_form.v1", "ticket": "RYA-1220",
        "source": ("Amarsi et al. 2021 Table 2 via "
                   "data/reference/molecular_cno_literature_rya1220/"),
        "note": ("model_atmosphere is a BOUND on the systematic carried by synthesizing "
                 "in 1D, not a correction; no abundance is moved (RYA-161). Closing this "
                 "component does not promote a route -- the remaining RYA-587 components "
                 "and the diagnostic's own standing still apply."),
        "routes": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="data/output/rya1220")
    ap.add_argument("--out", default="data/output/rya1220/molecular_route_model_form.json")
    args = ap.parse_args()
    doc = build(pathlib.Path(args.root))
    if not doc["routes"]:
        print("no molecular routes found", file=sys.stderr)
        return 2
    pathlib.Path(args.out).write_text(json.dumps(doc, indent=2) + "\n")

    print("RYA-1220 molecular model-form application")
    print(f"  {'diagnostic':10s} {'instrument':10s} {'mol':7s} {'sys':4s} "
          f"{'chi2*':>8s} {'in-band':>8s}  {'model_atm':9s} {'sigma':>6s}")
    for r in doc["routes"]:
        s = (f"{r['model_atmosphere_sigma_dex']:.3f}"
             if r["model_atmosphere_sigma_dex"] is not None else "  -")
        print(f"  {r['diagnostic']:10s} {r['instrument']:10s} {r['molecule']:7s} "
              f"{r['system']:4s} {r['red_chi2_max_UNCALIBRATED']:8.1f} "
              f"{'YES' if r['in_reference_band'] else 'NO':>8s}  "
              f"{r['model_atmosphere_state']:9s} {s:>6s}")
    n = sum(1 for r in doc["routes"] if r["model_atmosphere_state"] == "MEASURED")
    print(f"\n  model_atmosphere MEASURED: {n}/{len(doc['routes'])} (was 0 -- all HOLD_MOLECULAR)")
    print("  * red_chi2 is UNCALIBRATED on this data and decides nothing.")
    print(f"  -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
