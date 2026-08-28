#!/usr/bin/env python3
"""RYA-1089 — what the solar gate margin actually depends on.

The ticket asks a decision question ("which delta_xi") and then a factual one ("state
sigma_reported against the 0.05 gate, honestly"). This script answers the SECOND for
every candidate delta_xi on the table, so the decision is made against the consequence
rather than against a single pre-selected number.

🔴 THE FINDING THIS SCRIPT EXISTS TO RECORD. The recorded `solar.e_xi` (0.0588, RYA-311's
formal Delta-chi2 = 1 error) is SMALLER THAN EVERY OTHER DEFENSIBLE ESTIMATE of how far
the ADOPTED xi = 1.0 could sit from the truth -- including the spread of our own solve's
answers about that same 1.0. A formal error describes how tightly a fit pins its OWN root;
it is not a statement about a central value that fit did not return. Our solve returns
0.7088 at the production ceiling, 0.291 away from the adopted 1.0 -- five times the error
being carried. Those two numbers cannot both be right, and the budget currently carries
the smaller one.

Nothing here is chosen to land inside the gate (RYA-161): the script reports every
candidate, and the honest span straddles the gate rather than clearing it.

Inputs are read, never retyped:
  * |dA/dxi| and the per-element budget -> data/audit/uncertainty/solar_uncertainty_rya158.json
  * sigma_stat of each <3D> member      -> data/products/solar/Fe.json
  * the xi determinations               -> data/results/rya311/rya311_xi_ceiling_sensitivity_all.json
  * the recorded delta_xi candidates    -> config/stars.yaml (solar.e_xi*)
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config.constants import get_star_params           # noqa: E402

SOLAR_GATE = 0.05      # dex — the solar reported-uncertainty gate (RYA-282)
ADOPTED_XI = 1.0       # km/s — PINNED (RYA-196). NOT changed here; out of scope.
OUT = ROOT / "data/results/rya1089/rya1089_delta_xi_gate_margin.json"


def _load():
    budget = json.loads((ROOT / "data/audit/uncertainty/solar_uncertainty_rya158.json").read_text())
    fe_I = next(r for r in budget["per_element"] if r["element"] == "Fe" and r["ion"] == "I")
    dAdxi = abs(float(fe_I["dA_dvmic_per_kms"]))
    prods = json.loads((ROOT / "data/products/solar/Fe.json").read_text())["products"]
    mean3d = [p for p in prods if "mean3D" in str(p.get("treatment", ""))]
    sens = json.loads((ROOT / "data/results/rya311/rya311_xi_ceiling_sensitivity_all.json").read_text())
    return dAdxi, mean3d, sens


def _candidates(sens):
    """Every delta_xi with a defensible claim, with what it MEASURES (not its size)."""
    sp = get_star_params("solar")
    legs = [l for l in sens["legs"] if l["root_found"]]
    fx = [l["xi_measured"] for l in legs]
    ols = [l["xi_ols_root"] for l in legs]
    allx = fx + ols
    rms = lambda v: math.sqrt(sum((x - ADOPTED_XI) ** 2 for x in v) / len(v))
    return [
        (0.05, "RETIRED: uncited RYA-158 spec line 'vturb = 0.9 +/- 0.05'",
         "nothing -- no citation, and its central value (0.9) is not the adopted xi (1.0)"),
        (float(sp["e_xi"]), "RECORDED NOW: RYA-311 FITEXY formal, Delta-chi2 = 1",
         "how tightly the both-axis fit pins ITS OWN root (0.7088) given the EW errors"),
        (0.18, "RETIRED for the Sun (Ryan, 2026-08-27): Jofre+2014 sigma_vmic",
         "dispersion of vmic across seven independent GBS pipelines -- not our solve"),
        (rms(ols), "our OLS roots vs adopted 1.0 (n=3 ceilings)",
         "spread of the estimator the production bisection actually fits, about the pin"),
        (rms(allx), "our 6 determinations vs adopted 1.0 (3 ceilings x 2 estimators)",
         "spread of every defensible analysis of OUR data, about the pin"),
        (float(sp["e_xi_selection"]), "RYA-311 ceiling-selection spread (recorded, unused)",
         "how far xi moves across equally defensible 80/100/120 mA saturation cuts"),
        (abs(fx[1] - ADOPTED_XI) if len(fx) > 1 else float("nan"),
         "|production-leg xi - adopted xi| = |0.7088 - 1.0|",
         "the distance the budget must cover for the PIN to be defensible at all"),
        (float(sp["e_xi_chi2_scaled"]), "RYA-311 sqrt(chi2_red)-scaled (recorded, unused)",
         "the formal error inflated for a fit whose chi2_red is 27.9, i.e. a poor fit"),
        (max(allx) - min(allx), "full range of our 6 determinations",
         "the outer bound: best-to-worst across our own analysis choices"),
    ]


def main() -> int:
    dAdxi, mean3d, sens = _load()
    rows = []
    for dxi, source, measures in sorted(_candidates(sens), key=lambda r: r[0]):
        sigma_params = dAdxi * dxi
        per_member = {}
        for m in mean3d:
            st = float(m["sigma_stat"])
            sr = math.sqrt(st ** 2 + sigma_params ** 2)     # RYA-282 convention
            per_member[m["treatment"]] = {
                "sigma_stat": round(st, 4),
                "sigma_reported": round(sr, 4),
                "gate_ratio": round(sr / SOLAR_GATE, 3),
                "verdict": "INSIDE" if sr < SOLAR_GATE else "OVER",
            }
        rows.append({
            "delta_xi": round(dxi, 4),
            "source": source,
            "what_it_measures": measures,
            "sigma_params": round(sigma_params, 4),
            "members": per_member,
        })

    verdicts = {v["verdict"] for r in rows for v in r["members"].values()}
    out = {
        "ticket": "RYA-1089",
        "what": ("solar sigma_reported vs the 0.05 gate for every candidate delta_xi -- "
                 "the decision's consequence, computed rather than asserted"),
        "solar_gate_dex": SOLAR_GATE,
        "adopted_xi_kms": ADOPTED_XI,
        "adopted_xi_note": "PINNED (RYA-196); NOT changed by this ticket",
        "dA_dxi_per_kms": round(dAdxi, 4),
        "dA_dxi_source": "data/audit/uncertainty/solar_uncertainty_rya158.json (Fe I dA_dvmic_per_kms)",
        "recorded_e_xi": round(float(get_star_params("solar")["e_xi"]), 4),
        "candidates": rows,
        "finding": (
            "The gate verdict is NOT settled by the data: candidates span "
            f"{min(r['delta_xi'] for r in rows):.4f}-{max(r['delta_xi'] for r in rows):.4f} km/s and the "
            f"resulting sigma_reported straddles the 0.05 gate ({sorted(verdicts)}). The "
            "currently recorded e_xi is the SMALLEST candidate on the table and the only "
            "class that ignores where the adopted value sits relative to our own solve."
        ),
        "rya1104_dependency": (
            "⚠️ Ryan's RYA-311 sequencing caveat (2026-08-28): a reduced-EW/FITEXY solve is "
            "TILTED by gf systematics -- if gf quality correlates with line strength, the "
            "slope that determines xi is biased, so 0.7088 may carry the gf zero-point "
            "RYA-1104 is investigating. That gates the PIN-MOVE decision, and it qualifies "
            "the |0.7088 - 1.0| row here: the gap is real, but its CAUSE is unsettled. "
            "It does NOT rescue the formal error -- 0.0588 propagates EW errors only and "
            "knows nothing about gf, so a gf-contaminated solve makes it a WORSE "
            "underestimate of delta_p, not a better one. Either xi is genuinely ~0.71 (the "
            "pin is wrong) or our ability to determine xi is gf-limited (the uncertainty is "
            "larger than the fit reports). Neither reading supports 0.0588."
        ),
        "open_decision": (
            "Which delta_xi enters Type B is Ryan's call (RYA-311 spec 6 / RYA-1089 spec 4). "
            "This artifact states the consequence of each; it selects none."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n")

    w = max(len(r["source"]) for r in rows)
    print(f"|dA/dxi| = {dAdxi} dex per km/s     solar gate = {SOLAR_GATE} dex")
    print(f"adopted xi = {ADOPTED_XI} (pinned)  recorded e_xi = {out['recorded_e_xi']}\n")
    print(f"{'delta_xi':>8}  {'source':<{w}}  {'sig_par':>7} {'sig_rep':>7}  {'x gate':>6}  verdict")
    for r in rows:
        m = r["members"]["synth-mean3D-NLTE-gerber-stagger"]
        mark = "  <-- RECORDED" if abs(r["delta_xi"] - out["recorded_e_xi"]) < 1e-9 else ""
        print(f"{r['delta_xi']:8.4f}  {r['source']:<{w}}  {r['sigma_params']:7.4f} "
              f"{m['sigma_reported']:7.4f}  {m['gate_ratio']:6.2f}  {m['verdict']}{mark}")
    print(f"\n{out['finding']}")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
