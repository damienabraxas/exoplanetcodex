#!/usr/bin/env python3
"""
RYA-1220 -- assemble the RYA-587 budget for every molecular CNO route, on every
instrument, from the perturbation legs each route already ran.

Ryan, 2026-09-25: "if there is more uncertainty in IR that is fine. Especially if the
literature ran into the same problem ... we do the best we can and post results."

So this prices what the legs support and BOUNDS what they do not, rather than holding a
component because one input is unmeasured. An unmeasured correlation is not an unpriceable
term: rho is unknown in [-1, +1], so the term is evaluated at the WORST rho. That is the
RYA-1226 rule -- bound it, and say what was bounded.

WHAT IS PRICED PER ROUTE
  measurement        the fit's own sigma_fit
  continuum          half-range over the continuum legs; ONE-SIDED where only one leg ran,
                     reported as a one-sided bound and never halved into a fake symmetry
  molecular_coupling dA/dC and dA/dO from the C/O legs, combined with sigma_C and sigma_O
                     at the worst rho in [-1, +1]
  stellar.*          from the teff / vturb legs where they ran
  model_atmosphere   the measured Amarsi 2021 1D->3D bound (molecular_model_form)
  holding_instrument the spread across holdings that ran the SAME diagnostic

A component with no leg stays HOLD and names the leg that would close it. Nothing is
combined in quadrature across components here -- that is the contract's job, not this
script's; this produces the component table the contract consumes.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from pipeline.molecular_identity import identity_for
from pipeline.molecular_model_form import model_form_term

ROUTES = {
    "cn_crires_j_continuum_v1": ("CN_AX_J", "CRIRES+", "J", "solar_crires_plus_j_rya1219"),
    "cn_iag_3d_anchor_v1": ("CN_AX_IR", "IAG", "NIR", "solar_iag"),
    "cn_iag_coupling_v2": ("CN_AX_IR", "IAG", "NIR", "solar_iag"),
    "cn_iag_stellar_v2": ("CN_AX_IR", "IAG", "NIR", "solar_iag"),
    "cn_kp_continuum_v8": ("CN_AX_IR", "Kitt Peak", "NIR", "solar_kpno_molecfit_corrected"),
    "cn_harps_continuum_v2": ("CN_red", "HARPS", "VIS", "solar_harps_molecfit_corrected"),
    "cn_nearuv_continuum_v1": ("NH_AX", "Kitt Peak", "near-UV", "solar_kpno_kurucz2005_corrected"),
    "oh_crires_h_continuum_v1": ("OH_H", "CRIRES+", "H", "solar_crires_plus_h"),
    "co_crires_k_continuum_v1": ("CO_K", "CRIRES+", "K", "solar_crires_plus_k"),
}

#: Accepted anchor marginals for the CNO coupling (c_o_covariance_assessment.json).
SIGMA_C, SIGMA_O = 0.034, 0.421
ANCHOR_SOURCE = "data/output/rya1220/c_o_covariance_assessment.json (solar VIS CNO product)"


def legs(recs) -> dict:
    return {r["name"]: r["fit"]["A_X"] for r in recs
            if r.get("name") and isinstance(r.get("fit", {}).get("A_X"), (int, float))}


def half_range(vals: dict, minus: str, plus: str, nominal: float):
    """Symmetric half-range when both legs ran; a one-sided BOUND when only one did."""
    lo, hi = vals.get(minus), vals.get(plus)
    if lo is not None and hi is not None:
        return abs(hi - lo) / 2.0, "MEASURED", f"half-range over {minus}/{plus}"
    only = lo if lo is not None else hi
    if only is None:
        return None, "HOLD", f"neither {minus} nor {plus} ran"
    return abs(only - nominal), "MEASURED_ONE_SIDED", (
        f"ONE-SIDED: only {minus if lo is not None else plus} ran; the full excursion is "
        f"carried as the bound rather than halved into a symmetry that was not measured")


def coupling(vals: dict, nominal: float):
    """dA/dC and dA/dO from the C/O legs, at the worst rho in [-1, +1]."""
    if not {"C_-1", "C_+1", "O_-1", "O_+1"} <= set(vals):
        return None, "HOLD", "C/O perturbation legs did not run on this route", None
    # legs are +/-1 sigma of the anchor, so the response per sigma is the half-difference
    dC = (vals["C_+1"] - vals["C_-1"]) / 2.0
    dO = (vals["O_+1"] - vals["O_-1"]) / 2.0
    worst = math.sqrt(dC ** 2 + dO ** 2 + 2 * abs(dC * dO))   # rho at its worst sign
    return worst, "MEASURED_BOUNDED", (
        f"dA/dC={dC:+.3f}, dA/dO={dO:+.3f} per anchor sigma; rho is UNMEASURED so the "
        f"term is evaluated at the worst rho in [-1,+1]. Bounded, not assumed."), \
        {"dA_dC_per_sigma": round(dC, 4), "dA_dO_per_sigma": round(dO, 4),
         "sigma_C": SIGMA_C, "sigma_O": SIGMA_O, "anchor_source": ANCHOR_SOURCE}


def build(root: pathlib.Path) -> dict:
    per_diagnostic: dict[str, list] = {}
    rows = []
    for route, (key, instrument, band, holding) in ROUTES.items():
        path = root / route / "fits.json"
        if not path.exists():
            continue
        recs = json.loads(path.read_text())
        vals = legs(recs)
        nominal = vals.get("nominal")
        if nominal is None:
            continue
        nom_fit = [r for r in recs if r.get("name") == "nominal"][0]["fit"]
        ident = identity_for(key)
        comp: dict[str, dict] = {}

        comp["measurement"] = {"sigma_dex": nom_fit.get("sigma_fit"), "state": "MEASURED",
                               "note": "the fit's own sigma_fit"}

        for label, (m, p) in {"continuum": ("continuum_-1", "continuum_+1"),
                              "continuum_level": ("continuum_level_-1", "continuum_level_+1"),
                              "continuum_slope": ("continuum_slope_-1", "continuum_slope_+1")}.items():
            if m in vals or p in vals:
                s, st, why = half_range(vals, m, p, nominal)
                comp.setdefault("continuum", {"sigma_dex": 0.0, "state": "HOLD", "note": ""})
                if s is not None and s >= comp["continuum"].get("sigma_dex", 0):
                    comp["continuum"] = {"sigma_dex": round(s, 4), "state": st,
                                         "note": f"{label}: {why}"}

        comp.setdefault("continuum", {"sigma_dex": None, "state": "HOLD",
                                      "note": "no continuum leg ran"})

        s, st, why, detail = coupling(vals, nominal)
        comp["molecular_coupling"] = {"sigma_dex": round(s, 4) if s else None,
                                      "state": st, "note": why, "detail": detail}

        for label, (m, p) in {"stellar.teff": ("teff_K_-1", "teff_K_+1"),
                              "stellar.xi": ("vturb_kms_-1", "vturb_kms_+1")}.items():
            if m in vals or p in vals:
                s, st, why = half_range(vals, m, p, nominal)
                comp[label] = {"sigma_dex": round(s, 4) if s is not None else None,
                               "state": st, "note": why}

        mf = model_form_term(ident.molecule, *(lambda w: (min(x[0] for x in w), max(x[1] for x in w)))(
            json.loads((root / route / "provenance.json").read_text())["pool"]["windows_air_A"]))
        comp["model_atmosphere"] = {"sigma_dex": mf["sigma_dex"], "state": mf["state"],
                                    "note": mf["reason"]}

        row = {"route": route, "diagnostic": key, "instrument": instrument, "band": band,
               "holding": holding, "molecule": ident.molecule, "element": ident.element,
               "A_X": nominal, "components": comp}
        rows.append(row)
        per_diagnostic.setdefault(key, []).append((instrument, nominal))

    # Separate WITHIN-holding run scatter from ACROSS-holding spread. Collapsing them
    # would charge a convergence problem to the instrument.
    for row in rows:
        same_diag = per_diagnostic[row["diagnostic"]]
        by_inst: dict[str, list] = {}
        for inst, a in same_diag:
            by_inst.setdefault(inst, []).append(a)

        mine = by_inst.get(row["instrument"], [])
        if len(mine) >= 2:
            start = max(mine) - min(mine)
            row["fit_start_dependence_dex"] = round(start, 4)
            row["components"]["measurement"] = {
                "sigma_dex": round(max(row["components"]["measurement"]["sigma_dex"] or 0.0,
                                       start), 4),
                "state": "MEASURED_BOUNDED",
                "note": (f"sigma_fit is {row['components']['measurement']['sigma_dex']}, but "
                         f"{len(mine)} runs of this SAME diagnostic on this SAME holding land "
                         f"{start:.3f} dex apart ({', '.join(f'{v:.3f}' for v in sorted(mine))}). "
                         f"The chi2 surface is shallow, so the optimum depends on where the fit "
                         f"starts. The run spread bounds the measurement term; sigma_fit alone "
                         f"understates it.")}

        across = {k: (sum(v) / len(v)) for k, v in by_inst.items()}
        if len(across) >= 2:
            spread = max(across.values()) - min(across.values())
            row["components"]["holding_instrument"] = {
                "sigma_dex": round(spread, 4), "state": "MEASURED",
                "note": ("spread across holdings running the same diagnostic (per-holding "
                         "mean, so within-holding run scatter is not double counted): "
                         + ", ".join(f"{k} {v:.3f}" for k, v in sorted(across.items())))}
        else:
            row["components"]["holding_instrument"] = {
                "sigma_dex": None, "state": "HOLD",
                "note": (f"only one holding ({row['instrument']}) runs {row['diagnostic']}; "
                         f"a second is needed to measure repeatability")}
    return {"schema": "rya1220.molecular_route_budgets.v1", "ticket": "RYA-1220",
            "policy": ("price what the legs support, BOUND what they do not. An unmeasured "
                       "correlation is bounded at its worst value, never held."),
            "routes": rows}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="data/output/rya1220")
    ap.add_argument("--out", default="data/output/rya1220/molecular_route_budgets.json")
    args = ap.parse_args()
    doc = build(pathlib.Path(args.root))
    pathlib.Path(args.out).write_text(json.dumps(doc, indent=2) + "\n")

    order = ["measurement", "continuum", "molecular_coupling", "stellar.teff",
             "stellar.xi", "model_atmosphere", "holding_instrument"]
    print("RYA-1220 molecular route budgets -- all instruments\n")
    hdr = f"  {'diagnostic':10s} {'instrument':10s} {'A_X':>6s} " + " ".join(f"{c.split('.')[-1][:9]:>9s}" for c in order)
    print(hdr)
    for r in doc["routes"]:
        cells = []
        for c in order:
            g = r["components"].get(c)
            if not g or g.get("sigma_dex") is None:
                cells.append(f"{'HOLD':>9s}")
            else:
                mark = "*" if "ONE_SIDED" in g["state"] else ("b" if "BOUNDED" in g["state"] else " ")
                cells.append(f"{g['sigma_dex']:8.3f}{mark}")
        print(f"  {r['diagnostic']:10s} {r['instrument']:10s} {r['A_X']:6.3f} " + " ".join(cells))
    print("\n  * one-sided bound   b bounded over unmeasured rho")
    n_hold = sum(1 for r in doc["routes"] for c in order
                 if (r["components"].get(c) or {}).get("sigma_dex") is None)
    print(f"  components still HOLD: {n_hold} of {len(doc['routes'])*len(order)}")
    print(f"  -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
