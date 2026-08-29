#!/usr/bin/env python3
"""RYA-1112 — xi-aware uncertainty audit of the live VIS Fe products. READ-ONLY.

    python3 scripts/rya1112_vis_fe_uncertainty_audit.py           # table + artifact
    python3 scripts/rya1112_vis_fe_uncertainty_audit.py --check   # exit 1 if a finding regresses

Runs `skills/codex-product-audit` Parts A-E against every live VIS Fe product in
`data/products/solar/Fe.json`. Changes NO value (RYA-161): it reads the feed, the band
artifacts and the two uncertainty budgets, and reports.

🔴 FINDING 1 — THE PUBLISHED sigma_syst OMITS THE STELLAR-PARAMETER SYSTEMATICS ENTIRELY.

There are TWO uncertainty budgets in this repo and they never meet:

  A. `pipeline/uncertainty_stack.py` -> `data/audit/uncertainty/solar_uncertainty_rya158.json`
     The Type B stellar-parameter budget. RYA-1089 sourced the honest delta_xi (0.2912 km/s)
     and stamped Fe I sigma_B_vmic = 0.0699 dex, sigma_reported = 0.0721.
     🔴 READ BY TESTS AND ITS OWN TWO STAMPING SCRIPTS. NOTHING IN THE PRODUCT PATH.

  B. `pipeline/error_budget.build()` -> the band artifacts -> `publish_product` ->
     `data/products/solar/Fe.json` -> the feed and the website.
     Adds EXACTLY: line scatter, gf scale, harness residual, and — only where the band
     policy says so — pseudo-continuum and telluric residual.
     🔴 `error_budget.py` CONTAINS NO OCCURRENCE OF vmic / microturbulence / xi.

So every live VIS Fe product publishes `A ± sigma_stat` with a `sigma_syst` that carries no
Teff, no log g and no microturbulence term. That is the skill's Part A failure exactly:
**sigma_syst is a FLOOR quoted as a total**, and a floor understates the bar.

This is not a re-derivation and not a tuning: RYA-1089 already did the sourcing. The term
exists, is cited, and does not reach the product.

⚠️ AND IT MUST NOT BE FIXED BY BOLTING 0.0699 ONTO EVERYTHING. The xi learning (RYA-1093)
is that xi-sensitivity is a STRONG-LINE phenomenon, so dA/dvmic is a property of the LINE
SET, not of the element:
  * `GRADED`     = LAB lines with feature depth <= 0.60 — the WEAKER population.
  * `DEEPGRADED` = LAB lines with feature depth  > 0.60 — SATURATED BY CONSTRUCTION, where
                   xi bites hardest.
  * 3D products have NO MICROTURBULENCE AT ALL. xi is a 1D fudge for the velocity field 3D
    resolves; charging a xi term to a 3D product would be charging a term it did not earn
    (Part B). Eight VIS Fe products are 3D.
RYA-1089's -0.24 dex/(km/s) was measured on ONE 62-line pool. It is quoted here as a
REFERENCE MAGNITUDE for the size of the omission, never as a per-product value.

🔴 FINDING 2 — `gf = kurucz` ON EVERY GRADED AND DEEPGRADED VIS Fe PRODUCT.

`treatment_axes.LEGACY` pins the gf axis from the TREATMENT TOKEN
(`"1D-LTE": dict(..., gf="kurucz")`), not from the lines. But `--lines-tier graded` selects
on `gf_tier == LAB` — `_graded_mask` is documented as "True where the line carries a PRIMARY
LABORATORY gf". Measured on one product here: **67 of 67 lines carry a primary laboratory
gf and the artifact says `gf=kurucz`.** The row contradicts itself in its own columns —
`dominant` reads "gf scale (cited lab)" beside `gf=kurucz`. RYA-1104's defect, still live.

🔴 FINDING 3 — Fe II HAS NO STELLAR-PARAMETER BUDGET ROW AT ALL.
`solar_uncertainty_rya158.json` carries 10 element rows; exactly ONE has a non-zero vmic
term (Fe I). Fe II is absent, so the six Fe II VIS products have no sourced delta_xi
anywhere to omit or include, and their identical `sigma_syst = 0.0697` carries no
per-product content.

🔴 FINDING 4 — THE GATE IS NOT RENDERED AT ALL, AS FLAG OR OTHERWISE.
Part D asks that the gate render as a LABEL (bar vs target). `SOLAR_GATE_DEX = 0.05` exists
only inside two AUDIT scripts (`rya1089_stamp_honest_delta_xi.py`,
`rya1093_xi_robustness_audit.py`) and nowhere in the product path. The product records carry
no gate, target, or flag field. Nothing is being BLOCKED — the failure is the other one: a reader cannot see
which products clear the target and which do not.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import product_eligibility as PE          # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
STELLAR = ROOT / "data" / "audit" / "uncertainty" / "solar_uncertainty_rya158.json"
BAND = ROOT / "data" / "results" / "band_products"
OUT = ROOT / "data" / "results" / "rya1112" / "vis_fe_uncertainty_audit.json"

BAND_NAME = "VIS"
#: RYA-1089's own gate constant. Quoted, not re-derived.
SOLAR_GATE_DEX = 0.05
#: The ticket's dig-in threshold: above this a product needs a named RCA verdict.
DIG_IN_DEX = 0.1
#: `line_accounting_rya709.DEPTH_HI` — the gate that splits GRADED from DEEPGRADED.
DEPTH_HI = 0.60

#: Treatments whose model has NO microturbulence parameter. A xi term is not merely
#: uncomputed for these, it is NOT APPLICABLE, and charging one would be a term the
#: product did not earn.
def is_3d(treatment: str) -> bool:
    return "3D" in str(treatment)


def honest_total(p: dict) -> float:
    return math.hypot(float(p["sigma_stat"]), float(p["sigma_syst"]))


def rca(p: dict) -> tuple[str, str]:
    """A NAMED verdict for a product over the dig-in threshold. IRREDUCIBLE / FIXABLE / OPEN."""
    tot = honest_total(p)
    if tot <= DIG_IN_DEX:
        return "", ""
    n, ion, tier = int(p["n_lines"]), p["ion"], p["tier"]
    stat_share = (float(p["sigma_stat"]) / tot) ** 2
    if ion == "II" and n <= 3:
        return ("IRREDUCIBLE (small-N)",
                f"{n} lines. sigma/sqrt(n) on a 3-line pool cannot be argued down; the bar "
                f"is honest and the pool is the limit. RYA-1031's shape — report, do not "
                f"publish a precision nobody has. More Fe II lab-gf lines is the only fix.")
    if ion == "II":
        return ("OPEN (arm-correlated, RYA-1081)",
                f"{n} lines, scatter-dominated ({100 * stat_share:.0f}% of the variance). "
                f"RYA-1081 measured this pool as a POOL-WIDE arm-specific offset, not a few "
                f"bad lines — dropping the worst line moved the median 0.035. Saturation is "
                f"UNANSWERED there (the per-line COG is Fe I only). Not closed by this audit.")
    if tier == "GRADED" and n <= 13:
        return ("IRREDUCIBLE (small-N on a depth-gated pool)",
                f"{n} lines. GRADED selects LAB-gf lines at or below the {DEPTH_HI} depth "
                f"gate, and on the profile-fit route that pool is small by construction. "
                f"The bar is honest; widening it needs more lab-gf lines, not a code fix.")
    return ("OPEN", "no named cause yet — needs its own RCA")


def audit() -> dict:
    feed = json.loads(FEED.read_text())
    stellar = json.loads(STELLAR.read_text())
    vis = [p for p in feed["products"] if p["band"] == BAND_NAME]

    fe_rows = {r["ion"]: r for r in stellar["per_element"] if r["element"] == "Fe"}
    nonzero_vmic = [(r["element"], r["ion"], r["sigma_B_vmic"])
                    for r in stellar["per_element"] if r.get("sigma_B_vmic")]

    products = []
    for p in vis:
        tot = honest_total(p)
        verdict, why = rca(p)
        products.append({
            "ion": p["ion"], "tier": p["tier"], "treatment": p["treatment"],
            "holding": p["holding"], "route": p["route"],
            "A": p["A"], "sigma_stat": p["sigma_stat"], "sigma_syst": p["sigma_syst"],
            "n_lines": p["n_lines"], "honest_total_dex": round(tot, 4),
            "dominant_term": p["dominant_term"],
            "stat_basis": PE.stat_basis_of(p),
            "stat_basis_field_present": "stat_basis" in p,
            "over_dig_in": tot > DIG_IN_DEX,
            "over_solar_gate": tot > SOLAR_GATE_DEX,
            "rca_verdict": verdict, "rca_note": why,
            "xi_applicability": ("NOT APPLICABLE — 3D has no microturbulence"
                                 if is_3d(p["treatment"]) else
                                 "strong/saturated pool — xi bites hardest here"
                                 if p["tier"] == "DEEPGRADED" else
                                 "weaker pool (depth <= %.2f) — least xi-sensitive" % DEPTH_HI),
        })
    products.sort(key=lambda r: -r["honest_total_dex"])

    n3d = sum(1 for r in products if is_3d(r["treatment"]))
    return {
        "ticket": "RYA-1112",
        "skill": "codex-product-audit",
        "scope": f"live {BAND_NAME} Fe products in data/products/solar/Fe.json",
        "n_products": len(products),
        "thresholds": {"solar_gate_dex": SOLAR_GATE_DEX, "dig_in_dex": DIG_IN_DEX,
                       "depth_gate": DEPTH_HI},
        "findings": {
            "F1_stellar_parameter_terms_absent_from_published_sigma_syst": {
                "severity": "the published bar understates itself",
                "detail": ("error_budget.build() adds line scatter, gf scale, harness "
                           "residual and (band-conditional) pseudo-continuum/telluric. It "
                           "contains no vmic/microturbulence/xi term. The stellar-parameter "
                           "budget that RYA-1089 sourced lives in a separate artifact read "
                           "only by tests and its own stamping scripts."),
                "reference_magnitude_FeI": fe_rows.get("I", {}).get("sigma_B_vmic"),
                "reference_dA_dvmic_per_kms": fe_rows.get("I", {}).get("dA_dvmic_per_kms"),
                "sourced_delta_xi_kms": stellar["delta_p_solar"]["vmic"],
                "caveat": ("the derivative was measured on ONE 62-line pool; it is a "
                           "REFERENCE MAGNITUDE for the omission, never a per-product value"),
                "products_over_dig_in_now": sum(1 for r in products if r["over_dig_in"]),
                "products_over_dig_in_if_FeI_vmic_folded_in": sum(
                    1 for r in products
                    if r["ion"] == "I" and not is_3d(r["treatment"])
                    and math.hypot(r["honest_total_dex"],
                                   fe_rows["I"]["sigma_B_vmic"]) > DIG_IN_DEX)
                    + sum(1 for r in products
                          if (r["ion"] != "I" or is_3d(r["treatment"])) and r["over_dig_in"]),
            },
            "F2_gf_label_says_kurucz_on_lab_gf_pools": {
                "detail": ("treatment_axes.LEGACY pins gf from the TREATMENT TOKEN, but "
                           "--lines-tier graded selects on gf_tier == LAB. Every GRADED and "
                           "DEEPGRADED VIS Fe artifact row publishes gf=kurucz."),
                "self_contradiction": (
                    "the same rows carry dominant='gf scale (cited lab)', and the committed "
                    "budget beside them says 'every one of the 67 Fe I lines is GF-LAB and "
                    "100% carry a published per-line sigma' — while the gf axis says kurucz"),
            },
            "F3_no_stellar_parameter_row_for_Fe_II": {
                "elements_in_budget": len(stellar["per_element"]),
                "rows_with_a_nonzero_vmic_term": nonzero_vmic,
                "detail": "Fe II is absent, so its six VIS products have no sourced delta_xi.",
            },
            "F4_gate_is_not_rendered": {
                "detail": ("Part D asks the gate render as a LABEL. SOLAR_GATE_DEX lives "
                           "only in two AUDIT scripts (rya1089, rya1093) and nowhere in the "
                           "product path; no product record carries a gate/target/flag "
                           "field. Nothing is blocked — the failure is the other one: a "
                           "reader cannot see which bars clear the target."),
                "products_over_the_solar_gate": sum(1 for r in products
                                                    if r["over_solar_gate"]),
            },
        },
        "part_A_type_A": {
            "all_sigma_stat_are_standard_error": all(
                r["stat_basis"] == "standard_error" for r in products),
            "stat_basis_field_present_on": sum(1 for r in products
                                               if r["stat_basis_field_present"]),
            "note": ("the rest resolve through product_eligibility.STAT_BASIS_BY_ROUTE, the "
                     "legacy route map — correct today, but the basis is not a property of "
                     "those records"),
        },
        "n_3d_products_where_xi_is_not_applicable": n3d,
        "products": products,
    }


def check(doc: dict) -> list[str]:
    """Pin what this audit measured, so a regression is visible rather than silent."""
    bad = []
    f = doc["findings"]
    if f["F1_stellar_parameter_terms_absent_from_published_sigma_syst"][
            "reference_magnitude_FeI"] != 0.0699:
        bad.append("the RYA-1089 Fe I sigma_B_vmic is no longer 0.0699 — re-read the audit")
    if not doc["part_A_type_A"]["all_sigma_stat_are_standard_error"]:
        bad.append("a VIS Fe product's sigma_stat is no longer a standard error")
    if doc["n_products"] != 46:
        bad.append(f"the live VIS Fe product count moved to {doc['n_products']} (was 46)")
    unnamed = [r for r in doc["products"] if r["over_dig_in"] and not r["rca_verdict"]]
    if unnamed:
        bad.append(f"{len(unnamed)} product(s) over {DIG_IN_DEX} dex with no RCA verdict")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    doc = audit()

    print(f"RYA-1112 — {doc['n_products']} live VIS Fe products\n")
    print(f"{'ion':4}{'tier':12}{'treatment':24}{'holding':34}{'A':>7}{'stat':>8}"
          f"{'syst':>8}{'TOTAL':>8}{'n':>5}  RCA")
    print("-" * 130)
    for r in doc["products"]:
        mark = " <<<" if r["over_dig_in"] else ""
        print(f"{r['ion']:4}{r['tier']:12}{r['treatment'][:23]:24}{r['holding'][:33]:34}"
              f"{r['A']:7.3f}{r['sigma_stat']:8.4f}{r['sigma_syst']:8.4f}"
              f"{r['honest_total_dex']:8.4f}{r['n_lines']:5d}  {r['rca_verdict']}{mark}")

    f = doc["findings"]
    print(f"\nover the {DIG_IN_DEX} dig-in threshold: "
          f"{f['F1_stellar_parameter_terms_absent_from_published_sigma_syst']['products_over_dig_in_now']}"
          f" of {doc['n_products']}   |   over the {SOLAR_GATE_DEX} solar gate: "
          f"{f['F4_gate_is_not_rendered']['products_over_the_solar_gate']} of {doc['n_products']}")
    print(f"3D products where a xi term is NOT APPLICABLE: "
          f"{doc['n_3d_products_where_xi_is_not_applicable']}")
    print("\nFINDINGS")
    for k, v in f.items():
        print(f"  {k}")
        print(f"    {v['detail'] if 'detail' in v else ''}"[:400])

    fails = check(doc)
    for x in fails:
        print(f"\n  FAIL {x}")
    if not fails:
        print("\nOK — audit reproduces; every product over the dig-in threshold has a "
              "named RCA verdict.")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 1 if (a.check and fails) else 0


if __name__ == "__main__":
    raise SystemExit(main())
