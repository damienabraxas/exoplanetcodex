#!/usr/bin/env python3
"""RYA-1088: record the solar sigma_params on the ⟨3D⟩ product — as a MEASURED value.

🔴 WHY A ZERO HAS TO BE WRITTEN DOWN. For the Sun the Type B term is ~0 by construction:
logg and [Fe/H] are exact by definition (delta_p = 0, so `uncertainty_stack` correctly
skips the derivative runs) and Teff is known to +/-1 K. The ⟨3D⟩ product then carried NO
parameter term at all — and an absent term reads as an omission, not as a measured zero.

That is not hypothetical: a readiness checklist reading "sigma_params never measured" on
the SOLAR product manufactured a false alarm exactly this way, and the term had been built
since RYA-158. RYA-907 already ratified the principle for the other direction — UNMEASURED
IS NOT ZERO — and this is its mirror: A MEASURED ZERO IS NOT AN ABSENCE. Both states have
to be sayable.

MEASURED (pipeline.uncertainty_stack --star solar, Fe I, n=62):

    sigma_B Teff  0.0007   dA/dTeff = 0.0665 dex/100K, delta_p = 1 K
    sigma_B logg  0.0      delta_p = 0 -- exact by definition, no run
    sigma_B vmic  0.0120   dA/dvmic = -0.24 dex/(km/s), delta_p = 0.05
    sigma_B FeH   0.0      delta_p = 0 -- the Sun DEFINES the zero point, no run
    ------------------------------------------------------------------
    sigma_params  0.0120   quadrature; vmic-dominated

⚠️ THE DERIVATIVES COME FROM THE EW ROUTE, AND THE PRODUCT IS A SYNTHESIS ROUTE PRODUCT.
`uncertainty_stack._derive_A` re-derives through `abundances_derive.run` (EW->abundance),
so dA/dp is measured there while the ⟨3D⟩ product is a flux fit. RYA-1071 measured the
route systematic at +0.054 median with sd 0.176, so the two routes are NOT interchangeable
and this transfer is an approximation. It is recorded rather than hidden.

It is also a small approximation HERE and the numbers say why: sigma_B(Teff) is 0.0007 at
delta_p = 1 K and would stay negligible under any plausible dA/dTeff, and sigma_B(vmic) is
0.0120 against a sigma_stat of 0.0218 — so even a large error in dA/dvmic moves
sigma_reported in the third decimal. That argument does NOT transfer to a real star, where
delta_p are 20-80x larger.

⚠️ ADDS FIELDS ONLY. The abundance and the per-band sigma_syst are untouched — that budget
(`error_budget.py`) is a different quantity and is correct.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRODUCT = ROOT / "data" / "products" / "solar" / "Fe.json"
TYPE_B = ROOT / "data" / "audit" / "uncertainty" / "solar_uncertainty_rya158.json"
TREATMENTS = ("synth-mean3D-NLTE-gerber-stagger", "synth-mean3D-LTE-gerber-stagger")

REASON = ("MEASURED, not omitted: solar logg and [Fe/H] are exact by definition so "
          "their delta_p = 0 and no derivative was run; Teff is +/-1 K (sigma_B 0.0007); "
          "the term is vmic-DOMINATED. RYA-1093 re-stamps delta_xi from RYA-158's "
          "inherited 0.05 km/s to the 0.2912 km/s our own pool actually leaves "
          "unresolved, so sigma_B_vmic goes 0.0120 -> ~0.0699 and sigma_params is no "
          "longer small against sigma_stat. dA/dp measured on the EW route "
          "(uncertainty_stack), transferred to this synthesis product -- an approximation "
          "that does NOT transfer to a real star (RYA-1071).")

#: 🔴 RYA-1093 2D — THE PROVENANCE-TAGGED RE-STAMP ALLOWANCE.
#: RYA-1080's guard exists to stop a published number DRIFTING silently while artifacts are
#: reconciled. A recompute that carries an explicit reason and a source ticket is not drift,
#: and the distinction has to be recorded IN THE PRODUCT or the next reader cannot tell the
#: two apart from the file alone.
#:
#: ⚠️ AND THE GUARD DID NOT NEED WIDENING. RYA-1088's note proposed either (a) make the
#: guard precise -- refuse CHANGED values, allow ADDED keys -- or (b) record sigma_params
#: elsewhere. Measured on this branch: `published_value_edits` ALREADY iterates the BASELINE
#: fields and skips anything absent there ("a key ADDED later is another ticket's
#: business"), so option (a) has already happened and the stamp is unblocked. Nothing in
#: RYA-1093 relaxes a guard -- the tests that were xfailed simply pass now, and widening the
#: guard to make my own change fit would have been the antipattern RYA-1088 refused.
RESTAMP_PROVENANCE = {
    "source_ticket": "RYA-1093",
    "reason": ("delta_xi re-derived from RYA-158's inherited 0.05 km/s to 0.2912 km/s, the "
               "method+selection spread the 58-line Fe I pool cannot resolve. This is a "
               "RECOMPUTE WITH A STATED REASON, not drift: the audit that produced it is "
               "data/results/rya1093/xi_robustness_audit.json and the value is reproducible "
               "from it."),
    "audit_artifact": "data/results/rya1093/xi_robustness_audit.json",
    "supersedes": {"delta_xi_kms": 0.05, "source_ticket": "RYA-158"},
    "not_drift_because": ("a silent overwrite carries no reason and no ticket. This one "
                          "carries both, and the previous value is recorded in "
                          "`supersedes` so the change is legible from the product alone "
                          "(RYA-1080's intent: reconcile the artifacts, not the numbers)."),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--product", default=str(PRODUCT))
    ap.add_argument("--type-b", default=str(TYPE_B))
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    tb = json.loads(Path(a.type_b).read_text())
    row = next((r for r in tb["per_element"]
                if str(r.get("element")) == a.element and str(r.get("ion")) == "I"), None)
    if row is None:
        raise SystemExit(f"no {a.element} I row in {a.type_b} -- refusing to invent one")

    sig_b = {k: float(row.get(f"sigma_B_{k}", 0.0))
             for k in ("Teff", "logg", "vmic", "FeH")}
    sigma_params = round(math.sqrt(sum(v * v for v in sig_b.values())), 4)

    doc = json.loads(Path(a.product).read_text())
    touched = []
    for p in doc["products"]:
        if p.get("treatment") not in TREATMENTS:
            continue
        before_A, before_syst = p.get("A"), p.get("sigma_syst")
        p["sigma_params"] = sigma_params
        p["sigma_params_terms"] = {k: round(v, 4) for k, v in sig_b.items()}
        p["sigma_params_reason"] = REASON
        # RYA-1093 2D — the allowance travels WITH the value it permits.
        p["sigma_params_restamp"] = RESTAMP_PROVENANCE
        # sigma_reported = sqrt(sigma_stat^2 + sigma_params^2) -- the RYA-282 field
        # convention. sigma_syst is the SEPARATE per-band budget and is NOT folded in here.
        st = float(p.get("sigma_stat") or 0.0)
        p["sigma_reported"] = round(math.sqrt(st * st + sigma_params ** 2), 4)
        assert p.get("A") == before_A and p.get("sigma_syst") == before_syst, \
            "this script adds fields only -- it must never move A or sigma_syst"
        touched.append((p["treatment"], st, p["sigma_reported"]))

    if not touched:
        raise SystemExit(f"no ⟨3D⟩ product found in {a.product} -- nothing stamped")
    print(f"sigma_params = {sigma_params}  terms {sig_b}")
    # ⚠️ RYA-1093 2E — THE GATE IS REPORTED, NEVER MET BY CHOOSING A SMALLER BAR.
    _GATE = 0.05
    print(f"  vmic term alone: {sig_b['vmic']:.4f} dex  vs the {_GATE} dex solar gate "
          f"-> {'FAILS' if sig_b['vmic'] > _GATE else 'inside'}")
    if sig_b["vmic"] > _GATE:
        print(f"  ⚠️ the xi term ALONE exceeds the solar gate, before any other "
              f"systematic. That is a finding about the GATE (the field's own solar Fe "
              f"floor is 0.04-0.06 dex), not a number to tune: RYA-1093 forbids picking "
              f"the 0.0588 formal error because it would pass (RYA-161). The "
              f"adopt-or-relax call is Ryan's.")
    for t, st, rep in touched:
        print(f"  {t}: sigma_stat {st} -> sigma_reported {rep}")
    if a.dry_run:
        print("(dry run -- not written)")
        return 0
    Path(a.product).write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {a.product}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
