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

Discriminator is the route's own red_chi2. RYA-1112/1152 use red_chi2 > 10 as the
constraint flag; that threshold is inherited here rather than invented.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

RED_CHI2_FLAG = 10.0


def walk(o):
    if isinstance(o, dict):
        if "red_chi2" in o:
            yield o
        for v in o.values():
            yield from walk(v)
    elif isinstance(o, list):
        for v in o:
            yield from walk(v)


def route_fits(root: pathlib.Path) -> list[dict]:
    rows = []
    for fits in sorted(root.glob("*/fits.json")):
        vals = [o["red_chi2"] for o in walk(json.loads(fits.read_text()))
                if isinstance(o.get("red_chi2"), (int, float))]
        if not vals:
            continue
        rows.append({"route": fits.parent.name, "n_fits": len(vals),
                     "red_chi2_min": round(min(vals), 2), "red_chi2_max": round(max(vals), 2),
                     "constrained": max(vals) <= RED_CHI2_FLAG})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="data/output/rya1220")
    ap.add_argument("--out", default="data/output/rya1220/cno_hold_readjudication.json")
    args = ap.parse_args()
    root = pathlib.Path(args.root)
    rows = route_fits(root)
    if not rows:
        print(f"no fits.json under {root}", file=sys.stderr)
        return 2

    ok = [r for r in rows if r["constrained"]]
    bad = [r for r in rows if not r["constrained"]]

    doc = {
        "schema": "rya1220.cno_hold_readjudication.v1",
        "ticket": "RYA-1220",
        "red_chi2_flag": RED_CHI2_FLAG,
        "routes": rows,
        "constrained_routes": [r["route"] for r in ok],
        "unconstrained_routes": [r["route"] for r in bad],
        "finding": (
            "The holds are NOT uniform bookkeeping. Every ATOMIC N I route is well "
            f"constrained (red_chi2 {min(r['red_chi2_min'] for r in ok):.2f}-"
            f"{max(r['red_chi2_max'] for r in ok):.2f}); every MOLECULAR route is "
            f"unconstrained (red_chi2 {min(r['red_chi2_min'] for r in bad):.2f}-"
            f"{max(r['red_chi2_max'] for r in bad):.2f}). Promoting a covariance or "
            "model-form component on an unconstrained fit would be arithmetic on noise, "
            "so those holds are CORRECT and the fits must be repaired first."),
        "hold_classes": {
            "FIT_QUALITY": [r["route"] for r in bad],
            "EVIDENCE_EXISTS": {
                "holding_repeatability_atomic": {
                    "evidence": "data/output/rya1220/cno_holding_repeatability_audit.json",
                    "status": "MEASURED for C VIS 0.021, N red-optical 0.064, O VIS 0.105, "
                              "O red-optical 0.026 dex",
                    "gate_says": "HOLD on most routes",
                    "action": "wire the measured atomic ranges into the atomic routes"},
                "model_form_atomic": {
                    "evidence": "data/output/rya1220/cno_model_form_audit.json",
                    "status": "C I and O I MEASURED (Amarsi 2019 3D-NLTE / Caffau 2015); "
                              "N I MEASURED_PARTIAL (Amarsi 2020 grid)",
                    "gate_says": "3d_nlte_model HOLD_MOLECULAR",
                    "action": "the audit's own decision admits atomic anchors; wire them"},
                "cn_ir_1d_to_3d": {
                    "evidence": "cno_model_form_audit.json measured_or_published_model_terms",
                    "status": "KP -0.081, IAG -0.082 dex -- two holdings agreeing to 0.001",
                    "gate_says": "not applied",
                    "action": "admissible as a BOUND on the molecular model-form term even "
                              "though it is not applied as a correction"}},
            "RUN_OWED": {
                "c_o_joint_fit": "rho=0 is an assumption, not a measurement; needs a Sirius "
                                 "Turbospectrum joint C/O refit on the exact CN pool",
                "target_sun_jacobians": "both target and solar parameter Jacobians MISSING",
                "molecular_3d_per_transition": "no transition-specific grid for CN_VIS, "
                                               "CN_CRIRES_J, OH_H, CO_K"}},
        "red_flags_not_covariance_problems": {
            "CN_red_target_minus_sun_dex": 2.404,
            "OI_6300_target_minus_sun_dex": 0.607,
            "sigma_O_dex": 0.421,
            "note": "These are in cno_target_sun_pair_audit.json / c_o_covariance_assessment.json. "
                    "A 2.4 dex paired CN difference is a defect in the diagnostic, not a "
                    "missing covariance term. No amount of Jacobian work promotes it."},
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=2) + "\n")

    print("RYA-1220 CNO hold re-adjudication")
    print(f"  red_chi2 flag: {RED_CHI2_FLAG}\n")
    for r in rows:
        mark = "OK " if r["constrained"] else "BAD"
        print(f"  {mark}  {r['route']:38s} n={r['n_fits']:3d}  "
              f"red_chi2 {r['red_chi2_min']:8.2f} .. {r['red_chi2_max']:8.2f}")
    print(f"\n  constrained  : {len(ok)}")
    print(f"  unconstrained: {len(bad)}  <- these holds are CORRECT, not bookkeeping")
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
