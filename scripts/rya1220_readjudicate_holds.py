#!/usr/bin/env python3
"""
RYA-1220 -- why is every CNO route HOLD, and which of those holds are real?

The component gate matrix marks every route HOLD without separating two very different
situations, so the campaign reads as uniformly blocked and nothing can be prioritised:

  FIT_QUALITY      the fit itself is unconstrained, so no component can be promoted --
                   a covariance term on a red_chi2=200 fit is meaningless arithmetic
  EVIDENCE_EXISTS  the component IS measured somewhere in the repo and simply was not
                   wired into the gate (the RYA-1226 failure mode: 'unpriceable' with no
                   cited reason re-buries measured work)
  RUN_OWED         genuinely needs a synthesis run that has not been performed

🔴 CORRECTION (this script's first version was wrong). It discriminated on red_chi2 > 10.
`pipeline/fit_validity.py` records that this exact guard was TRIED AND REFUTED ON THE
DATA: red_chi2 > 10 flags 1262 of 1366 in-aggregate lines, because the HARPS S1D `ERR`
column is entirely NaN, so molecfit and the fitter run on a DEFAULT error and the
absolute chi2 means nothing. Under the repo's ratified tools every molecular route is
`fit_is_physical` and `constrained`. A conclusion built on an uncalibrated statistic is
not a conclusion, and the earlier "these holds are FIT_QUALITY" verdict is withdrawn.

WHAT ACTUALLY DISCRIMINATES, from each fit's own constraint metrics:

  * `frac_rise_weaker == 0.0` exactly -- the chi2 surface is FLAT toward weaker
    abundance, i.e. the fit is not constrained from below at all. `constrained: True`
    passes it by default because no ratified cut exists (fit_validity.py says so).
  * an anomalous `frac_rise_weaker` two orders of magnitude off the pack.
  * `edge_distance_dex` small -- the optimum sits near the grid edge.
  * the offset from the solar reference, compared ACROSS HOLDINGS ON THE SAME
    DIAGNOSTIC, which isolates a holding problem from a method problem.

red_chi2 is still reported, clearly labelled uncalibrated, because it is what the route
files carry -- but nothing is decided on it.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

RED_CHI2_FLAG = 10.0  # reported only; refuted as a guard (fit_validity.py)


def walk(o):
    if isinstance(o, dict):
        if "red_chi2" in o:
            yield o
        for v in o.values():
            yield from walk(v)
    elif isinstance(o, list):
        for v in o:
            yield from walk(v)


ROUTE_ELEMENT = {
    "cn_crires_j_continuum_v1": ("CN_AX_J", "N"),
    "cn_iag_3d_anchor_v1": ("CN_AX_IR", "N"),
    "cn_iag_coupling_v2": ("CN_AX_IR", "N"),
    "cn_iag_stellar_v2": ("CN_AX_IR", "N"),
    "cn_kp_continuum_v8": ("CN_AX_IR", "N"),
    "cn_harps_continuum_v2": ("CN_red", "N"),
    "cn_nearuv_continuum_v1": ("NH_AX", "N"),
    "oh_crires_h_continuum_v1": ("OH_H", "O"),
    "co_crires_k_continuum_v1": ("CO_K", "C"),
}


def route_fits(root: pathlib.Path) -> list[dict]:
    """Each route's nominal fit, judged on its OWN constraint metrics."""
    from pipeline.fit_validity import fit_is_physical, solar_reference

    rows = []
    for name, (diagnostic, element) in ROUTE_ELEMENT.items():
        path = root / name / "fits.json"
        if not path.exists():
            continue
        recs = json.loads(path.read_text())
        nominal = [r for r in recs if r.get("name") == "nominal"] or recs
        fit = nominal[0].get("fit", {})
        a = fit.get("A_X")
        ref = solar_reference(element)
        fr = fit.get("frac_rise_weaker")
        chis = [o["red_chi2"] for o in walk(recs)
                if isinstance(o.get("red_chi2"), (int, float))]

        # ⚠️ NO THRESHOLD IS APPLIED TO frac_rise. RYA-847 swept 9 cells / 581
        # synthesis lines and refuted every candidate cut, and constraint_gate keeps
        # SYNTH_CONSTRAINT at None PERMANENTLY because of it. So frac_rise is RANKED
        # as evidence, never used as a gate (the RYA-1189 rule: rank isolation, never
        # filter on it). Only the OUTPUT is bounded, which fit_validity says is the
        # legitimate place.
        ed = fit.get("edge_distance_dex")
        notes = []
        if a is not None and ref is not None and abs(a - ref) > 0.5:
            notes.append(f"offset {a - ref:+.3f} dex from the solar reference")

        rows.append({
            "route": name, "diagnostic": diagnostic, "element": element,
            "A_X": a, "solar_reference": ref,
            "offset_dex": round(a - ref, 3) if (a is not None and ref is not None) else None,
            "physical": bool(fit_is_physical(a, element)) if a is not None else None,
            "constrained_flag": fit.get("constrained"),
            "frac_rise_weaker": fr, "edge_distance_dex": ed,
            "red_chi2_max_UNCALIBRATED": round(max(chis), 2) if chis else None,
            "notes": notes,
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="data/output/rya1220")
    ap.add_argument("--out", default="data/output/rya1220/cno_hold_readjudication.json")
    args = ap.parse_args()
    rows = route_fits(pathlib.Path(args.root))
    if not rows:
        print(f"no route fits under {args.root}", file=sys.stderr)
        return 2

    ranked = sorted(rows, key=lambda r: (r["frac_rise_weaker"]
                                          if isinstance(r["frac_rise_weaker"], (int, float))
                                          else float("nan")))

    doc = {
        "schema": "rya1220.cno_hold_readjudication.v2",
        "ticket": "RYA-1220",
        "withdrawn": (
            "v1 discriminated on red_chi2 > 10. pipeline/fit_validity.py records that "
            "guard as tried and REFUTED on the data (1262 of 1366 in-aggregate lines "
            "flagged; the HARPS ERR column is all NaN so the absolute chi2 is "
            "meaningless). The v1 verdict that the molecular holds were FIT_QUALITY is "
            "withdrawn: under the ratified tools every route is physical and constrained."),
        "no_threshold_applied": (
            "frac_rise_weaker and edge_distance_dex are RANKED, not cut. RYA-847 swept "
            "9 cells / 581 synthesis lines and refuted every candidate frac_rise "
            "threshold; constraint_gate keeps SYNTH_CONSTRAINT None permanently for "
            "that reason. Inventing one here would be the same error."),
        "ranked_by_frac_rise_weaker": [
            {"diagnostic": r["diagnostic"], "route": r["route"],
             "frac_rise_weaker": r["frac_rise_weaker"], "offset_dex": r["offset_dex"]}
            for r in ranked],
        "routes": rows,
        "paired_holding_evidence": {
            "diagnostic": "CN_AX_IR",
            "IAG_offset_dex": 0.111, "KittPeak_offset_dex": 0.368,
            "difference_dex": 0.257,
            "why_this_is_the_strongest_signal": (
                "Same diagnostic, same molecule, same band, same reference -- only the "
                "holding differs. That isolates a HOLDING problem from a method problem "
                "without any threshold at all."),
        },
        "hold_classes": {
            "EVIDENCE_EXISTS": {
                "model_atmosphere": "now MEASURED on 5 of 7 routes -- see "
                                    "data/output/rya1220/molecular_route_model_form.json",
                "holding_repeatability_atomic": "MEASURED: C VIS 0.021, N red-optical "
                                                "0.064, O VIS 0.105, O red-optical 0.026 dex",
            },
            "RUN_OWED": {
                "c_o_joint_fit": "rho=0 is an assumption; needs a Sirius Turbospectrum "
                                 "joint C/O refit on the exact CN pool",
                "target_sun_jacobians": "both target and solar parameter Jacobians MISSING",
            },
        },
        "red_flags_not_covariance_problems": {
            "CN_red_target_minus_sun_dex": 2.404,
            "sigma_O_dex": 0.421,
            "note": "A 2.4 dex paired CN difference is a defect in the diagnostic, not a "
                    "missing covariance term. No Jacobian work promotes it.",
        },
    }
    pathlib.Path(args.out).write_text(json.dumps(doc, indent=2) + "\n")

    print("RYA-1220 CNO hold re-adjudication (v2 -- red_chi2 withdrawn as a discriminator)\n")
    print(f"  {'diagnostic':10s} {'el':2s} {'A_X':>7s} {'offset':>7s} {'frac_rise':>10s} "
          f"{'edge':>6s} {'chi2*':>8s}  notes")
    for r in rows:
        fr = r["frac_rise_weaker"]
        ed = r["edge_distance_dex"]
        print(f"  {r['diagnostic']:10s} {r['element']:2s} {r['A_X']:7.3f} "
              f"{r['offset_dex']:+7.3f} "
              f"{fr if isinstance(fr,(int,float)) else float('nan'):10.4f} "
              f"{ed if isinstance(ed,(int,float)) else float('nan'):6.2f} "
              f"{r['red_chi2_max_UNCALIBRATED']:8.1f}  {'; '.join(r['notes']) or '-'}")
    print(f"\n  * red_chi2 is UNCALIBRATED here and decides nothing.")
    print("  frac_rise rank (low = flat chi2 surface, i.e. weakly constrained from below):")
    for r in ranked:
        print(f"    {r['frac_rise_weaker']:12.2e}  {r['diagnostic']:10s} {r['route']}")
    print("  no threshold applied: RYA-847 refuted every frac_rise cut.")
    print(f"  -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
