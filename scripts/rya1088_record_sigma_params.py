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

REASON = ("MEASURED ~0, not omitted: solar logg and [Fe/H] are exact by definition so "
          "their delta_p = 0 and no derivative was run; Teff is +/-1 K (sigma_B 0.0007); "
          "the term is vmic-dominated (sigma_B 0.0120 at delta_p 0.05 km/s, Jofre+2014 "
          "scale). dA/dp measured on the EW route (uncertainty_stack), transferred to this "
          "synthesis product -- an approximation that is small here because sigma_params "
          "is half sigma_stat, and that does NOT transfer to a real star (RYA-1071).")


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
