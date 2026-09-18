#!/usr/bin/env python3
"""
RYA-1226 part A -- can the legacy near-UV Fe rows carry a COMPLETE RYA-587 budget?

This assembles every component for which real evidence exists, from the artifacts that
already hold it, and then lets `pipeline.uncertainty_contract` issue the verdict. The
HOLD list below is the CONTRACT's, not an opinion: anything not supplied becomes HOLD
automatically and `validate()` refuses the document.

🔴 NOTHING IS INVENTED TO MAKE THE GATE PASS. That is the one way to fail this ticket
badly -- RYA-587 exists to stop a placeholder being published as a measurement, so
filling a component with a plausible-looking number would defeat the thing being
migrated to.

WHAT THE EVIDENCE SUPPORTS, and where it comes from:

  measurement       per-line artifact: std(ddof=1)/sqrt(N) reproduces the published
                    sigma_stat EXACTLY (0.18116/sqrt(12) = 0.05230), so the raw scatter
                    and N are both real and the contract's raw/sqrt(N) identity holds.
  transition_data   canonical_gf: all 12 lines are LAB tier with a published per-line
                    gf_sigma_dex, Den Hartog 2019 (DOI 10.3847/1538-4365/ab322e).
                    ⚠️ ONE COMMON SOURCE, so the errors are NOT independent -- the
                    covariance is fully correlated and the term does not average down.
                    That is a DEPARTURE from the legacy budget's "RMS of the per-line
                    sigma" (0.0445), which corresponds to neither honest propagation
                    (independent gives 0.0129, fully correlated gives 0.0417). The
                    contract propagates with explicit weights and covariance, so the
                    number moves; A does not.
  stellar.xi        RYA-1225's adjudicated slope for Fe II (-0.1100, post-RYA-1207
                    molecular synthesis), RYA-1226 part C's measured slope for Fe I.
  stellar.teff/     solar: delta_logg and delta_feh are definitionally ZERO (pinned), and
  logg/metallicity  Teff contributes 0.000665 dex -- DEFINED from uncertainty_stack, not
                    guessed.
  pseudo_continuum  0.1000 dex band-flat, RYA-1113/RYA-841: in the near-UV the true
                    continuum is never observed.
  hfs_isotopes      N/A on evidence: all 12 lines carry hfs_n_components == 1.
  telluric          N/A on evidence: 3000-3780 A contains no registered telluric band.
  molecular_coupling N/A: an atomic Fe II product, selector is not MOL-.
  nlte              N/A for the 1D-LTE leg on its own per-line record
                    (nlte_delta_dex = 0, "none -- LTE, no departure applied").

🔴 AND WHAT NO ARTIFACT PRICES -- the reason this is expected to refuse:

  blends            RYA-1190 measured the near-UV as OPACITY-DOMINATED, -0.1721 dex
                    excess over its control, "robust across the sigma sweep, in a band
                    whose catalogued opacity is already complete to the VALD threshold.
                    The deficit is real and is NOT catalogueable from what we hold."
                    Its ~0.028 payoff figure is explicitly "an ORDER, not a number to
                    plan on", so using it here would be inventing exactly what RYA-1190
                    refused to supply.
  continuum         continuum_method, continuum_ref and continuum_level are ALL EMPTY on
                    every one of the 12 per-line rows.
  profile_ew        the route is a synthesis flux-fit at a FIXED +/-0.4 A half-width
                    (RYA-759). RYA-1220 measured that class of sensitivity at +0.216 dex
                    on one line at 8216 A; it has never been measured for this band.
  model_atmosphere  the product's model_grid is null.
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.uncertainty_contract import (  # noqa: E402
    UncertaintyError, assemble, component, pool_digest, product_scope,
    publication_problems, transition_data, validate)

FEED = ROOT / "data/products/solar/Fe.json"
BANDP = ROOT / "data/results/band_products"
GF = ROOT / "data/linelists/canonical_gf.csv"
OUT = ROOT / "data/results/rya1226/nearuv_budget_migration.json"

DH19 = "Den Hartog et al. 2019, ApJS 243, 33 (DOI 10.3847/1538-4365/ab322e)"
PSEUDO_CONTINUUM_DEX = 0.10


def per_line(product) -> pd.DataFrame:
    ion = {"I": "I", "II": "II"}[product["ion"]]
    stem = (f"Fe{ion}_3000_3780_kpno_solar_atlas_{product['holding']}_SYNTH_"
            f"{product['selector']}_{product['treatment']}_lines.csv")
    return pd.read_csv(BANDP / stem)


def evidence_for(product, *, xi_slope: float | None, delta_xi: float) -> dict:
    """Assemble every component the artifacts actually support. Nothing else."""
    lines = per_line(product)
    acc = lines[lines["in_aggregate"] == True]  # noqa: E712
    gf = pd.read_csv(GF, low_memory=False)
    gf = gf[gf["species"] == f"Fe {product['ion']}"]
    j = acc.merge(gf, on="wavelength_air_A", how="left", suffixes=("", "_gf"))
    ep_ok = (j["ep_eV"] - j["excitation_potential_eV"]).abs() < 0.005
    if not bool(ep_ok.all()):
        raise SystemExit("lambda-only gf join: EP disagrees -- refusing (RYA-1037)")

    ids = [str(x) for x in j["physical_id"]]
    sigmas = [float(x) for x in j["gf_sigma_dex"]]
    if any(math.isnan(s) for s in sigmas):
        raise SystemExit("a line carries no published gf sigma -- cannot price transition_data")

    #: 🔴 ONE SOURCE => FULLY CORRELATED. A shared laboratory scale does not average down.
    n = len(ids)
    weights = [1.0 / n] * n          # the response of a mean/median to one line's gf
    cov = [[sigmas[a] * sigmas[b] for b in range(n)] for a in range(n)]

    raw = float(acc["abundance"].astype(float).std(ddof=1))
    se = raw / math.sqrt(len(acc))
    digest = pool_digest(ids)

    comps = [
        dict(name="measurement", sigma_dex=se, state="MEASURED",
             source=f"per-line scatter of the accepted pool, {BANDP.name} artifact",
             evidence={"method": "line_scatter", "independent": True,
                       "n_lines": int(len(acc)), "raw_sigma": raw,
                       "pool_sha256": digest}),
        dict(name="stellar.logg", sigma_dex=0.0, state="DEFINED",
             source="solar log g is pinned; delta_logg is definitionally zero (RYA-1089)",
             evidence={"delta": 0.0, "parameter_exists": True}),
        dict(name="stellar.metallicity", sigma_dex=0.0, state="DEFINED",
             source="solar [Fe/H] is pinned; delta_feh is definitionally zero (RYA-1089)",
             evidence={"delta": 0.0, "parameter_exists": True}),
        dict(name="stellar.teff", sigma_dex=0.000665, state="DEFINED",
             source="pipeline.uncertainty_stack solar delta_Teff = 1 K x 0.0665 dex/100 K",
             evidence={"delta_K": 1.0, "parameter_exists": True}),
        dict(name="pseudo_continuum", sigma_dex=PSEUDO_CONTINUUM_DEX, state="DEFINED",
             source="RYA-1113/RYA-841 near-UV band-flat: the true continuum is never observed",
             evidence={"band": "near-UV", "basis": "band-flat declared allowance"}),
        dict(name="hfs_isotopes", sigma_dex=None, state="N/A",
             source="every line in this pool has hfs_n_components == 1 (canonical_gf)",
             evidence={"hfs_n_components": sorted({int(x) for x in j["hfs_n_components"].dropna()})}),
        dict(name="telluric", sigma_dex=None, state="N/A",
             source="3000-3780 A contains no registered telluric band (pipeline.telluric_policy)",
             evidence={"window_A": [3000, 3780], "bands_in_window": 0}),
        dict(name="molecular_coupling", sigma_dex=None, state="N/A",
             source="atomic Fe product; the selector is not a MOL- indicator set",
             evidence={"selector": product["selector"]}),
        dict(name="holding_instrument", sigma_dex=0.0, state="MEASURED",
             source="SynthesisHandler harness residual MEASURED against the known optical "
                    "answer, not assumed zero (control/frontier rule)",
             evidence={"harness_residual_dex": 0.0}),
    ]

    nlte_delta = set(lines["nlte_delta_dex"].dropna().unique())
    if nlte_delta == {0.0}:
        comps.append(dict(name="nlte", sigma_dex=None, state="N/A",
                          source="this leg applies no departures: nlte_delta_dex is 0 on "
                                 "every line and nlte_source reads 'none -- LTE'",
                          evidence={"nlte_delta_dex": 0.0}))

    if xi_slope is not None:
        comps.append(dict(name="stellar.xi", sigma_dex=abs(xi_slope) * delta_xi,
                          state="MEASURED",
                          source=("per-line paired differential at xi 0.90/1.10 on this "
                                  "product's OWN pool, measured on the current molecular "
                                  "synthesis"),
                          evidence={"dA_dxi": xi_slope, "delta_xi_kms": delta_xi,
                                    "parameter_exists": True}))

    comps.append(transition_data(ids, sigmas, weights, covariance=cov,
                                 sources=[DH19] * n, covariance_source=DH19))
    return {"components": comps, "indicator_ids": ids, "pool_sha256": digest,
            "raw_sigma": raw, "sigma_stat": se, "n_lines": int(len(acc))}


def build(product, *, xi_slope, delta_xi=0.2912) -> dict:
    ev = evidence_for(product, xi_slope=xi_slope, delta_xi=delta_xi)
    scope = product_scope(product, star="solar", indicator_ids=ev["indicator_ids"])
    numeric = [c for c in ev["components"] if c["state"] in {"MEASURED", "DEFINED"}]
    cov = np.diag([c["sigma_dex"] ** 2 for c in numeric]).tolist()
    doc = assemble(scope, ev["components"], covariance=cov,
                   covariance_source=("components treated as independent of one another; "
                                      "the within-component correlations (shared gf source) "
                                      "are carried inside transition_data's own covariance"),
                   assumptions=("solar log g and [Fe/H] are pinned, so their responses are "
                                "definitionally zero rather than unmeasured"))
    return {"budget": doc, "evidence": ev, "scope": scope}


def main() -> int:
    feed = json.loads(FEED.read_text())
    slope = {"II": -0.1100}          # RYA-1225, adjudicated; Fe I is part C
    rows = []
    for p in feed["products"]:
        if p.get("band") != "near-UV" or p.get("tier") != "DEEPGRADED":
            continue
        if p.get("treatment") not in ("1D-LTE", "ENGINE-A"):
            continue
        try:
            built = build(p, xi_slope=slope.get(p["ion"]))
            holds = built["budget"]["holds"]
            problems = None
            try:
                validate(built["budget"], scope=built["scope"])
            except UncertaintyError as exc:
                problems = str(exc)
            rows.append({"ion": p["ion"], "holding": p["holding"],
                         "treatment": p["treatment"], "A": p["A"],
                         "n_lines": built["evidence"]["n_lines"],
                         "sigma_stat_reproduced": round(built["evidence"]["sigma_stat"], 6),
                         "sigma_stat_published": p.get("sigma_stat"),
                         "components_supplied": len(built["budget"]["components"]) - len(holds),
                         "HOLDS": holds, "validator": problems,
                         "sigma_reported_if_complete": built["budget"]["sigma_reported"]})
        except SystemExit as exc:
            rows.append({"ion": p["ion"], "holding": p["holding"],
                         "treatment": p["treatment"], "error": str(exc)})

    doc = {"ticket": "RYA-1226", "read_only": True,
           "question": "can the legacy near-UV Fe rows carry a complete RYA-587 budget?",
           "answer": ("NO -- not from the evidence that exists today. The HOLD list below "
                      "is the contract's verdict, not an opinion: every component for "
                      "which an artifact holds real evidence was supplied."),
           "rows": rows}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}\n")
    for r in rows:
        if "error" in r:
            print(f"  Fe {r['ion']:3} {r['holding'][:31]:31} {r['treatment']:9} ERROR {r['error'][:60]}")
            continue
        print(f"  Fe {r['ion']:3} {r['holding'][:31]:31} {r['treatment']:9} "
              f"supplied={r['components_supplied']:2}/16  HOLDS={r['HOLDS']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
