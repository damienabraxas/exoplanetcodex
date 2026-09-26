#!/usr/bin/env python3
"""
RYA-1220 -- the SOLAR stellar-parameter Jacobians for the CNO diagnostics.

Closes `stellar.teff` and `stellar.xi` on the RYA-587 budget, which the gate matrix
carried as HOLD on every CNO route, by running the perturbation legs the harness already
supports and differencing them.

WHICH ROUTE EACH ELEMENT IS MEASURED ON, AND WHY IT IS NOT THE SAME REGION FOR ALL THREE:

  C   region `vis`         CH G-band, the primary solar carbon indicator
  O   region `vis`         [O I] 6300, the forbidden blend
  N   region `nir_cn_iag`  CN A-X, the RATIFIED PRIMARY

🔴 N IS DELIBERATELY NOT MEASURED ON `vis`. That region's only N diagnostic is CN_red,
and its ksi leg returned NON-MINIMUM -- `frac_rise_weaker = -5.6e-06`, i.e. chi2 at the
bracket end was not above the reported minimum, so the optimiser did not find a minimum
and the fit determined nothing. There is no slope to take on a flat surface. Running it
again would reproduce a non-result on a diagnostic the ratified policy already rejects,
so the Jacobian is taken on the route the policy actually publishes.

A LEG PAIR THAT STRADDLES NOMINAL IS A SLOPE. ONE THAT DOES NOT IS CURVATURE.
Where both perturbed legs land on the SAME side of nominal, the response is not
monotonic across the step and the central difference is meaningless at this step size --
it is reported as UNRESOLVED with a bound, never as a small slope. Solar C's Teff legs
do exactly that: 8.485 and 8.483 against a nominal 8.477.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

RUNS = {
    "C": ("stellar_jacobian_vis_C", "vis", "CH_Gband"),
    "O": ("stellar_jacobian_vis_O", "vis", "OI_6300"),
    "N": ("stellar_jacobian_nir_cn_N", "nir_cn_iag", "CN_AX_IR"),
}
PARAMS = {"stellar.teff": ("teff_K_-1", "teff_K_+1"),
          "stellar.xi": ("vturb_kms_-1", "vturb_kms_+1")}


def legs(recs) -> dict:
    return {r["name"]: r["fit"] for r in recs if r.get("name")}


def jacobian(fits: dict, minus: str, plus: str) -> dict:
    nom = fits.get("nominal", {}).get("A_X")
    lo = fits.get(minus, {}).get("A_X")
    hi = fits.get(plus, {}).get("A_X")
    if nom is None or lo is None or hi is None:
        return {"state": "HOLD", "sigma_dex": None,
                "why": f"missing leg: {minus if lo is None else plus}"}
    unconstrained = [n for n in (minus, plus, "nominal")
                     if fits.get(n, {}).get("constrained") is False]
    if unconstrained:
        return {"state": "HOLD", "sigma_dex": None, "nominal_A": nom,
                "why": (f"leg(s) {unconstrained} did not converge; a non-minimum fit has "
                        f"no slope to take")}
    central = (hi - lo) / 2.0
    straddles = (lo - nom) * (hi - nom) < 0
    out = {"nominal_A": round(nom, 4), "minus_A": round(lo, 4), "plus_A": round(hi, 4),
           "central_difference_dex_per_sigma": round(central, 5),
           "same_side_of_nominal": not straddles}
    if straddles:
        out.update({"state": "MEASURED", "sigma_dex": abs(round(central, 5)),
                    "why": "both legs straddle nominal; the central difference is a slope"})
    else:
        bound = max(abs(lo - nom), abs(hi - nom))
        out.update({"state": "UNRESOLVED", "sigma_dex": round(bound, 5),
                    "why": (f"both legs sit on the same side of nominal ({lo:.4f}, "
                            f"{hi:.4f} vs {nom:.4f}), so the response is not monotonic "
                            f"across this step and the central difference "
                            f"({central:+.5f}) is curvature, not a slope. The larger "
                            f"single-leg excursion is carried as a BOUND.")})
    return out


def build(root: pathlib.Path) -> dict:
    rows = []
    for element, (folder, region, diagnostic) in RUNS.items():
        path = root / folder / "fits.json"
        if not path.exists():
            rows.append({"element": element, "region": region, "diagnostic": diagnostic,
                         "state": "NOT_RUN", "components": {}})
            continue
        fits = legs(json.loads(path.read_text()))
        comps = {name: jacobian(fits, m, p) for name, (m, p) in PARAMS.items()}
        rows.append({"element": element, "region": region, "diagnostic": diagnostic,
                     "n_legs": len(fits),
                     "sigma_fit": fits.get("nominal", {}).get("sigma_fit"),
                     "components": comps})
    return {"schema": "rya1220.solar_stellar_jacobians.v1", "ticket": "RYA-1220",
            "star": "solar",
            "note": ("stellar.teff and stellar.xi for the CNO diagnostics, from the "
                     "harness's own perturbation legs. A pair that straddles nominal "
                     "gives a slope; one that does not is reported UNRESOLVED with a "
                     "bound rather than as a small slope."),
            "n_is_on_the_ratified_route": (
                "N is measured on nir_cn_iag (CN A-X), not vis. The vis region's only N "
                "diagnostic is CN_red, whose ksi leg returned NON-MINIMUM "
                "(frac_rise_weaker = -5.6e-06): no minimum, so no slope."),
            "elements": rows}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="data/output/rya1220")
    ap.add_argument("--out", default="data/output/rya1220/solar_stellar_jacobians.json")
    a = ap.parse_args()
    doc = build(pathlib.Path(a.root))
    pathlib.Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
    print("RYA-1220 solar stellar-parameter Jacobians\n")
    for row in doc["elements"]:
        if row.get("state") == "NOT_RUN":
            print(f"  {row['element']}  {row['diagnostic']:10s}  NOT RUN")
            continue
        print(f"  {row['element']}  {row['diagnostic']:10s} ({row['region']})  "
              f"sigma_fit={row['sigma_fit']}")
        for name, g in row["components"].items():
            s = f"{g['sigma_dex']:.5f}" if g.get("sigma_dex") is not None else "   -"
            print(f"      {name:14s} {g['state']:10s} {s}")
    print(f"\n  -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
