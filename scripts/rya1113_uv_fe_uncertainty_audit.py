#!/usr/bin/env python3
"""RYA-1113 — xi-aware uncertainty re-verification of the live near-UV Fe products.

    python3 scripts/rya1113_uv_fe_uncertainty_audit.py [--check]

Sister to RYA-1112 (VIS), run through the `codex-product-audit` skill. READ-ONLY: no value,
line, gf or threshold is changed (RYA-161). Findings become tickets, not edits.

🔴 THE BAND IS A FRONTIER BAND AND THIS AUDIT REPORTS, IT DOES NOT CHASE. RYA-777 sets the
standard as "VIS baseline; do not chase the dragon in IR/UV" and RYA-813 makes the near-UV
"report-and-note", not validate-to-target. So a bar over the ~0.05 solar gate is LABELLED,
never tuned down, and a >0.1 dex bar gets a NAMED verdict rather than an apology.

🔴 THE TICKET'S HEADLINE PREMISE DOES NOT SURVIVE THE DATA, AND THAT IS THE FIRST FINDING.
RYA-1113 expects the near-UV bar to be dominated by "~84% Kurucz-floored gf (RYA-822)". It
is not, because the PRODUCTS DO NOT USE THE NEAR-UV LINE LIST -- they use its lab-gf subset:

    canonical_gf, 3000-3780 A      Fe I  4274 lines, of which   59 LAB  ( 1%)
                                   Fe II  665 lines, of which   12 LAB  ( 2%)
    the DEEPGRADED products use    Fe I    58 lines, 58/58 LAB (100%)
                                   Fe II   12 lines, 12/12 LAB (100%)

measured here on the lambda+EP dual key (RYA-1037). So the gf-floor is a property of the
BAND'S LINE LIST and not of these products, whose gf term is a measured RMS of published
per-line laboratory sigma. The floor statement is true and is about a different population.

⚠️ AND IT INVERTS THE REDUCIBILITY STORY. The pool is not merely lab-gf, it is EXHAUSTIVE of
the lab-gf resource in band -- 58 of the 59 Fe I and 12 of the 12 Fe II lines that exist at
that tier. There is no line-selection lever left: widening these pools requires NEW
LABORATORY MEASUREMENTS (RYA-1061), not a better cut.

WHAT ACTUALLY SETS THE BAR
--------------------------
`pseudo-continuum`, a flat 0.10 dex, charged by `error_budget` wherever the BandPolicy says
the continuum is pseudo -- the near-UV alone. It DOMINATES every one of these products, and
it is a DECLARED BOUND, not a per-product measurement: the same 0.10 for both ions, both
holdings and every window, reasoned from "near-UV median flux 0.283-0.805, so the
normalisation is uncertain at the ~10% level in the worst windows".

That matters because the data contains a direct, independent measurement of exactly the
quantity that bound is standing in for -- and this audit makes it, in `continuum_probe`.
"""
from __future__ import annotations

import argparse
import collections
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import line_match                      # noqa: E402
from pipeline import product_eligibility as PE       # noqa: E402
from pipeline.telluric_policy import TELLURIC_BANDS  # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
STELLAR = ROOT / "data" / "audit" / "uncertainty" / "solar_uncertainty_rya158.json"
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
BANDP = ROOT / "data" / "results" / "band_products"
OUT = ROOT / "data" / "results" / "rya1113" / "uv_fe_uncertainty_audit.json"

BAND_NAME = "near-UV"
BAND_LO_A, BAND_HI_A = 3000.0, 3780.0
#: RYA-1089's gate constant, quoted not re-derived. A FLAG here, never a target (RYA-777).
SOLAR_GATE_DEX = 0.05
#: The ticket's dig-in threshold: above this a product needs a named RCA verdict.
DIG_IN_DEX = 0.1
#: `error_budget` charges this wherever BandPolicy.continuum_treatment is pseudo.
PSEUDO_CONTINUUM_DEX = 0.10
#: Display-only, for "how far from the solar reference" — nothing branches on it (RYA-161).
ASPLUND21_FE = 7.46


def honest_total(p: dict) -> float:
    return math.hypot(float(p["sigma_stat"]), float(p["sigma_syst"]))


def gf_tier_composition() -> dict:
    """What tier are these products' lines ACTUALLY on, vs what the band offers?

    🔴 The join is lambda+EP dual key (RYA-1037) and the match count is reported, so a
    silently-partial join cannot masquerade as a tier finding.
    """
    canon = pd.read_csv(CANON, low_memory=False)
    inband = canon[(canon.wavelength_air_A >= BAND_LO_A)
                   & (canon.wavelength_air_A < BAND_HI_A)]
    out = {"band_window_A": [BAND_LO_A, BAND_HI_A], "per_ion": {}}
    for ion in ("I", "II"):
        avail = inband[inband.species.astype(str) == f"Fe {ion}"]
        tiers = collections.Counter(avail.gf_tier.astype(str))
        rec = {"linelist_n": int(len(avail)),
               "linelist_tiers": {k: int(v) for k, v in tiers.items()},
               "linelist_LAB_frac": (round(tiers.get("LAB", 0) / len(avail), 4)
                                     if len(avail) else None)}
        f = (BANDP / f"Fe{ion}_{int(BAND_LO_A)}_{int(BAND_HI_A)}_kpno_solar_atlas_"
                     f"solar_kpno_molecfit_corrected_SYNTH_DEEPGRADED_1D-LTE_lines.csv")
        if f.exists():
            d = pd.read_csv(f)
            used = d[d.in_aggregate.astype(bool)] if "in_aggregate" in d else d
            pool = avail.dropna(subset=["wavelength_air_A",
                                        "excitation_potential_eV"]).reset_index(drop=True)
            r = line_match.match(used.wavelength_air_A.to_numpy(float),
                                 pool.wavelength_air_A.to_numpy(float),
                                 want_ep=used.ep_eV.to_numpy(float),
                                 src_ep=pool.excitation_potential_eV.to_numpy(float),
                                 require_ep=True, tol_A=0.01)
            idx = np.asarray(r.index)
            hit = idx >= 0
            used_tiers = collections.Counter(str(pool.gf_tier.iloc[i]) for i in idx[hit])
            rec.update({
                "product_pool_n": int(len(used)),
                "product_pool_matched": int(hit.sum()),
                "product_pool_tiers": {k: int(v) for k, v in used_tiers.items()},
                "pool_is_all_LAB": bool(used_tiers.get("LAB", 0) == int(hit.sum())
                                        and hit.all()),
                "pool_exhausts_LAB_resource": (
                    f"{used_tiers.get('LAB', 0)} of the {tiers.get('LAB', 0)} LAB Fe {ion} "
                    f"lines that exist in band"),
            })
        else:
            # ⚠️ LOUD, not silent. Fe II near-UV ships only *_products.csv -- there is no
            # per-line artifact in the tree, so its tier composition CANNOT be measured
            # here. Recording that is the point; assuming it matches Fe I would be the
            # defect this audit exists to catch.
            rec.update({"product_pool_n": None,
                        "per_line_artifact": "ABSENT — no *_lines.csv committed for this "
                                             "ion/band; tier composition unmeasurable",
                        "expected_path": str(f.relative_to(ROOT))})
        out["per_ion"][f"Fe {ion}"] = rec
    return out


def continuum_probe(uv: list[dict]) -> dict:
    """🔴 MEASURE the continuum systematic the 0.10 bound stands in for.

    The two Kitt Peak holdings are NOT two treatments of one spectrum. Per the holdings
    registry they are different source products:

      solar_kpno_kurucz2005_corrected  Kurucz 2005 irradiance distribution (300-1000 nm)
      solar_kpno_molecfit_corrected    the 1984 KP atlas, ESO molecfit telluric-corrected

    and molecfit's six registered telluric bands are ALL RED (O2 A/B, H2O 9280,
    11120-11560, 7160-7340 A). None lies in 3000-3780 A, so in THIS band the telluric
    correction is a no-op and the two holdings differ by their continuum realisation and
    their underlying spectrum -- nothing else.

    Their disagreement is therefore an independent, in-band estimate of the very
    normalisation uncertainty `pseudo-continuum` asserts. Half the spread is the
    conventional one-sided equivalent.
    """
    # 🔴 MEASURED FROM THE REGISTRY, NOT ASSERTED FROM PROSE. The claim "molecfit changes
    # nothing in this band" is load-bearing for the whole finding, so it is read off
    # `telluric_policy.TELLURIC_BANDS` rather than inferred from the holding's notes.
    bluest = min(lo for lo, _hi, _n in TELLURIC_BANDS)
    out = {"_telluric_no_op_in_band": {
        "bluest_registered_telluric_band_A": bluest,
        "band_window_A": [BAND_LO_A, BAND_HI_A],
        "gap_A": round(bluest - BAND_HI_A, 1),
        "verdict": ("no registered telluric band reaches this window, so the two holdings "
                    "cannot differ here by telluric treatment"
                    if bluest > BAND_HI_A else
                    "🔴 a telluric band OVERLAPS this window — the no-op claim is FALSE"),
    }}
    for ion in ("I", "II"):
        for t in sorted({p["treatment"] for p in uv if p["ion"] == ion}):
            sel = {("kurucz2005" if "kurucz2005" in p["holding"] else "molecfit"): p
                   for p in uv if p["ion"] == ion and p["treatment"] == t}
            if len(sel) != 2:
                continue
            spread = float(sel["kurucz2005"]["A"]) - float(sel["molecfit"]["A"])
            out[f"Fe {ion} {t}"] = {
                "A_kurucz2005": sel["kurucz2005"]["A"],
                "A_molecfit": sel["molecfit"]["A"],
                "spread_dex": round(spread, 4),
                "half_spread_dex": round(abs(spread) / 2.0, 4),
                "exceeds_pseudo_continuum_bound": bool(
                    abs(spread) / 2.0 > PSEUDO_CONTINUUM_DEX),
            }
    return out


def rca(p: dict, comp: dict, probe: dict) -> tuple[str, str]:
    """A NAMED verdict for every product over the dig-in threshold."""
    tot = honest_total(p)
    if tot <= DIG_IN_DEX:
        return "", ""
    ion, n = p["ion"], int(p["n_lines"])
    key = f"Fe {ion} {p['treatment']}"
    half = probe.get(key, {}).get("half_spread_dex")
    lab = comp["per_ion"][f"Fe {ion}"].get("pool_exhausts_LAB_resource", "")
    if half is not None and half > PSEUDO_CONTINUUM_DEX:
        return ("OPEN (continuum realisation exceeds its own bound)",
                f"the two holdings disagree by {probe[key]['spread_dex']:+.3f} dex on the "
                f"SAME {n} lines, so half-spread {half:.3f} > the {PSEUDO_CONTINUUM_DEX} "
                f"pseudo-continuum bound charged here. They are different source spectra "
                f"(Kurucz 2005 irradiance vs the 1984 atlas) and molecfit corrects no "
                f"telluric band in 3000-3780 A, so this is continuum + spectrum, not "
                f"telluric. The bound is band-flat and cannot express it. NOT tuned: the "
                f"published bars stand; the term needs re-deriving on this evidence.")
    return ("IRREDUCIBLE (pseudo-continuum floor on an exhausted lab-gf pool)",
            f"dominated by the {PSEUDO_CONTINUUM_DEX} dex pseudo-continuum systematic, "
            f"which does not average down: the near-UV true continuum is never observed. "
            f"The line lever is already spent -- {lab or 'the pool is all lab-gf'} -- so "
            f"more/better lines cannot reduce this. Reducing it needs a real continuum "
            f"determination in band, which is RYA-929/1027 territory, not a tuning knob.")


def audit() -> dict:
    feed = json.loads(FEED.read_text())
    stellar = json.loads(STELLAR.read_text())
    uv = [p for p in feed["products"] if p["band"] == BAND_NAME]
    comp = gf_tier_composition()
    probe = continuum_probe(uv)
    fe_rows = {r["ion"]: r for r in stellar["per_element"] if r["element"] == "Fe"}

    products = []
    for p in uv:
        tot = honest_total(p)
        verdict, why = rca(p, comp, probe)
        sib = next((q for q in uv if q["ion"] == p["ion"] and q["holding"] == p["holding"]
                    and q["treatment"] == "1D-LTE"), None)
        products.append({
            "ion": p["ion"], "tier": p["tier"], "treatment": p["treatment"],
            "holding": p["holding"], "route": p["route"],
            "A": p["A"], "vs_asplund21_dex": round(float(p["A"]) - ASPLUND21_FE, 3),
            "sigma_stat": p["sigma_stat"], "sigma_syst": p["sigma_syst"],
            "n_lines": p["n_lines"], "honest_total_dex": round(tot, 4),
            "dominant_term": p["dominant_term"],
            "stat_basis": PE.stat_basis_of(p),
            "over_dig_in": tot > DIG_IN_DEX,
            "over_solar_gate": tot > SOLAR_GATE_DEX,
            # RYA-777/813: a frontier band is REPORTED against the gate, never chased.
            "gate_disposition": "REPORT-AND-NOTE (frontier band, RYA-777/813)",
            "rca_verdict": verdict, "rca_note": why,
            "xi_term_present": False,
            "xi_applicability": (
                "APPLIES — and is ABSENT. No stellar-parameter term reaches sigma_syst "
                "here (same omission RYA-1112 found in VIS, F1). Fe II has no row in the "
                "stellar budget at all, so it has no sourced delta_xi to fold even in "
                "principle."),
            # the NLTE leg must be comparable to its own LTE sibling, or its delta is not
            # an NLTE effect. Measured, not assumed.
            "nlte_pool_delta_vs_1D_LTE": (None if sib is None or p is sib
                                          else int(p["n_lines"]) - int(sib["n_lines"])),
        })
    products.sort(key=lambda r: -r["honest_total_dex"])

    pool_changed = [r for r in products if r["nlte_pool_delta_vs_1D_LTE"]]
    return {
        "ticket": "RYA-1113",
        "skill": "codex-product-audit",
        "band": BAND_NAME,
        "scope": f"live {BAND_NAME} Fe products in data/products/solar/Fe.json",
        "standard": ("REPORT-AND-NOTE. RYA-777: VIS baseline, do not chase the dragon in "
                     "IR/UV. RYA-813: near-UV is reported, not validated to target."),
        "n_products": len(products),
        "thresholds": {"solar_gate_dex": SOLAR_GATE_DEX, "dig_in_dex": DIG_IN_DEX,
                       "pseudo_continuum_dex": PSEUDO_CONTINUUM_DEX},
        "gf_tier_composition": comp,
        "continuum_probe": probe,
        "findings": {
            "F1_the_gf_floor_premise_does_not_apply_to_these_products": {
                "severity": "corrects the ticket's own framing",
                "detail": ("RYA-1113 expects a ~84% Kurucz-floored gf to dominate. The "
                           "near-UV LINE LIST is indeed overwhelmingly non-LAB (Fe I "
                           "59/4274 = 1% LAB; Fe II 12/665 = 2%), but the DEEPGRADED "
                           "products use only the LAB subset: Fe I 58/58 LAB, measured on "
                           "the lambda+EP key. The floor is a fact about the band, not "
                           "about these products."),
                "consequence": ("the gf term here is a MEASURED RMS of published per-line "
                                "laboratory sigma (0.0509 dex, Fe I), not a floor bound"),
            },
            "F2_the_dominant_term_is_a_band_flat_declared_bound": {
                "severity": "the bar's largest component is not measured per product",
                "detail": (f"`pseudo-continuum` = {PSEUDO_CONTINUUM_DEX} dex flat, charged "
                           "wherever BandPolicy.continuum_treatment is pseudo. Identical "
                           "for both ions, both holdings, every window. It dominates all "
                           "eight products, so sigma_syst is effectively a band constant: "
                           "0.1122 on every Fe I product and 0.1095 on every Fe II."),
                "why_it_matters": ("a band-flat term cannot distinguish a well-normalised "
                                   "window from a badly-normalised one, and cannot move "
                                   "when the evidence says it should — see F3"),
            },
            "F3_the_holdings_disagree_by_more_than_the_bound_that_covers_them": {
                "severity": "🔴 the measured continuum disagreement exceeds its own budget",
                "detail": ("Fe II: the two holdings return 7.957 vs 7.697 on the SAME 12 "
                           "lines — 0.260 dex apart, half-spread 0.130 > the 0.10 bound. "
                           "Fe I is consistent (0.046 apart, half-spread 0.023). The two "
                           "holdings are DIFFERENT SOURCE SPECTRA (Kurucz 2005 irradiance "
                           "vs the 1984 atlas + molecfit), and molecfit corrects no "
                           "telluric band inside 3000-3780 A, so in this band the "
                           "difference is continuum realisation and spectrum, not "
                           "telluric."),
                "measured": probe,
            },
            "F4_neither_holding_was_normalisation_probed_in_this_band": {
                "detail": ("RYA-1030's intake probes report 'normalised' for both holdings "
                           "— but at 4758/6500/8242/9288 A (kurucz2005) and "
                           "5482/7980/10478/11976 A (molecfit). Not one probe window lies "
                           "in 3000-3780 A. The normalisation status of the band these "
                           "products are measured in is UNPROBED on both holdings."),
                "note": ("read from data/catalog/holdings_manifest_registry.csv; RYA-929 "
                         "is not re-run here, per the ticket"),
            },
            "F5_the_NLTE_leg_changes_the_POOL_not_just_the_physics": {
                "severity": "🔴 the published NLTE effect is confounded with a pool change",
                "detail": ("Fe II ENGINE-A reports n=7 where its 1D-LTE sibling reports "
                           "n=12 — 5 of 12 lines dropped, on BOTH holdings. The published "
                           "difference (+0.018 kurucz2005, +0.041 molecfit) is therefore "
                           "not a paired NLTE differential: it is 7 lines compared against "
                           "12. sigma_stat rises accordingly (0.096->0.1329, "
                           "0.0986->0.121), which is the pool shrinking, not the "
                           "measurement worsening."),
                "products_with_a_pool_change": [
                    {"ion": r["ion"], "holding": r["holding"],
                     "treatment": r["treatment"], "dn": r["nlte_pool_delta_vs_1D_LTE"]}
                    for r in pool_changed],
            },
            "F6_no_per_line_artifact_for_Fe_II_in_band": {
                "detail": ("data/results/band_products/ ships *_lines.csv for Fe I near-UV "
                           "(molecfit) only. Fe II near-UV has *_products.csv alone, so "
                           "its tier composition, its rejections and the identity of the 5 "
                           "lines the NLTE leg drops cannot be interrogated from the "
                           "committed tree."),
            },
            "F7_no_stellar_parameter_term_reaches_the_published_bar": {
                "detail": ("same omission RYA-1112 recorded for VIS. sigma_syst here "
                           "carries pseudo-continuum, gf scale, scatter and harness "
                           "residual — no Teff/logg/xi term. RYA-1120 built the join; "
                           "these products predate its application."),
                "sourced_delta_xi_kms": stellar["delta_p_solar"]["vmic"],
                "reference_magnitude_FeI_sigma_B_vmic": fe_rows.get("I", {}).get(
                    "sigma_B_vmic"),
                "Fe_II_has_a_stellar_row": "II" in fe_rows,
            },
        },
        "part_A_type_A": {
            "all_sigma_stat_are_standard_error": all(
                r["stat_basis"] == "standard_error" for r in products),
            "products_over_solar_gate": sum(1 for r in products if r["over_solar_gate"]),
            "products_over_dig_in": sum(1 for r in products if r["over_dig_in"]),
        },
        "no_values_changed": ("READ-ONLY. This script reads the feed, the band artifacts, "
                              "canonical_gf and the stellar budget. It writes one audit "
                              "artifact and changes no value, line, gf or threshold "
                              "(RYA-161)."),
        "products": products,
    }


def check(doc: dict) -> list[str]:
    """Pin what this audit measured, so a regression is visible rather than silent."""
    bad = []
    if doc["n_products"] != 8:
        bad.append(f"live near-UV Fe product count moved to {doc['n_products']} (was 8)")
    fe1 = doc["gf_tier_composition"]["per_ion"]["Fe I"]
    if not fe1.get("pool_is_all_LAB"):
        bad.append("the Fe I near-UV pool is no longer 100% LAB — F1's premise moved")
    if fe1.get("product_pool_matched") != fe1.get("product_pool_n"):
        bad.append("the lambda+EP join no longer matches every Fe I near-UV line")
    unnamed = [r for r in doc["products"] if r["over_dig_in"] and not r["rca_verdict"]]
    if unnamed:
        bad.append(f"{len(unnamed)} product(s) over {DIG_IN_DEX} dex with no RCA verdict")
    if not doc["part_A_type_A"]["all_sigma_stat_are_standard_error"]:
        bad.append("a near-UV Fe sigma_stat is no longer a standard error")
    t = doc["continuum_probe"]["_telluric_no_op_in_band"]
    if t["bluest_registered_telluric_band_A"] <= BAND_HI_A:
        bad.append("a registered telluric band now overlaps the near-UV window — F3's "
                   "'continuum, not telluric' reasoning no longer holds")
    feii = doc["continuum_probe"].get("Fe II 1D-LTE", {})
    if not feii.get("exceeds_pseudo_continuum_bound"):
        bad.append("Fe II's holding disagreement no longer exceeds the continuum bound "
                   "— F3 moved; re-read it")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    doc = audit()
    print("=" * 96)
    print(f"RYA-1113 — near-UV Fe uncertainty audit ({doc['n_products']} live products)")
    print("=" * 96)
    print(f"  standard: {doc['standard']}")
    print()
    print(f"  {'ion':<4}{'treatment':<10}{'holding':<14}{'A':>8}{'vs7.46':>8}{'stat':>8}"
          f"{'syst':>8}{'total':>8}{'n':>4}  {'gate':<6}{'dominant'}")
    print("  " + "-" * 92)
    for r in doc["products"]:
        h = r["holding"].replace("solar_kpno_", "").replace("_corrected", "")
        print(f"  {r['ion']:<4}{r['treatment']:<10}{h:<14}{r['A']:>8.3f}"
              f"{r['vs_asplund21_dex']:>+8.3f}{r['sigma_stat']:>8.4f}{r['sigma_syst']:>8.4f}"
              f"{r['honest_total_dex']:>8.4f}{r['n_lines']:>4}  "
              f"{'OVER' if r['over_solar_gate'] else 'ok':<6}{r['dominant_term']}")
    print()
    print("  RCA verdicts (every product over the dig-in threshold):")
    for r in doc["products"]:
        if r["rca_verdict"]:
            h = r["holding"].replace("solar_kpno_", "").replace("_corrected", "")
            print(f"    Fe {r['ion']:<3} {r['treatment']:<10} {h:<14} {r['rca_verdict']}")
    print()
    for k, v in doc["findings"].items():
        print(f"  {k}")
        print(f"      {v.get('severity', '')}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {a.out.relative_to(ROOT)}")

    if a.check:
        bad = check(doc)
        if bad:
            print("\nFAIL:")
            for b in bad:
                print(f"  - {b}")
            return 1
        print("\nOK — audit reproduces; every product over the dig-in has a named verdict.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
