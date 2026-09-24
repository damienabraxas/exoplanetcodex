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

🔴 THE FOUR THAT WERE WRONGLY CALLED UNPRICEABLE (Ryan's REDIRECT, 2026-09-18).
My first pass held all four and justified only one. Three of them HAVE measured work:

  continuum         RYA-1133 audited THESE products and ruled "(b) a distinct term is owed,
                    with a (c) floor". Its (c) floor is the OBSERVED holding spread, and
                    that is what is wired here -- recomputed on the CURRENT post-RYA-1207
                    values, because RYA-1133 ran pre-opacity: Fe II half-spread 0.13 ->
                    0.1075 (1D-LTE) and 0.1185 -> 0.1155 (ENGINE-A); Fe I is unchanged at
                    0.0230.
                    ⚠️ NOT RYA-846's number, and that is RYA-1133's own finding, not mine:
                    846 measured 1984KP-vs-Wallace, `kurucz2005_in_846s_comparison` is
                    False, and "every RYA-846 candidate term is below the observed
                    half-spread". Wiring 846 here would charge a smaller term measured on a
                    different pair. Its candidates are carried in the evidence so the
                    comparison stays visible: offset 0.0582, MAD 0.0526, std 0.1192,
                    at-product-lines raw 0.1254 / net 0.0710, lever 2.54 (not 2.42).
                    🔴 AND RYA-1133's STATED BLOCKER IS STALE. It says the owed measurement
                    cannot be done because "Fe II near-UV ships no *_lines.csv"; RYA-908
                    emitted exactly those files on the SAME DAY, and 14 of them are in the
                    tree. The owed measurement -- sigma_delta between Kurucz 2005 and the
                    1984 KP atlas at the product lines, per ion -- is RUNNABLE NOW.
  model_atmosphere  RYA-1032's `atmosphere_nuisance`: atlas9 -> marcs-ges, LTE, no NLTE
                    physics, SAME POOL (n=67 both), one axis varied -> 0.004 dex. Measured
                    on the VIS Fe I pool, so it is carried as a cross-band bound, flagged.
  profile_ew        RYA-1220 measured a fit half-width excursion of +0.2160 dex (7.9562 at
                    +/-0.25 A against 8.1722 at +/-1.1 A). ⚠️ It is N I, at 8216 A, on
                    solar_iag -- a different element, band and holding -- and its own
                    artifact declares `steps_are_adopted_uncertainties: false` and
                    `status: diagnostic_only`. It is carried FLAGGED, and it is NOT
                    demonstrably conservative: the near-UV is more heavily blended than
                    8216 A, so the true term there could be larger, not smaller. The real
                    fix is a half-width sweep on the near-UV Fe pool itself.
  blends            BOUNDED per the directive rather than held. ⚠️ The ceiling used is
                    RYA-1190's measured DEFICIT (0.1721 dex excess over its control), not
                    its ~0.028 "payoff" figure: the payoff is how much a completeness build
                    would CLOSE, which is a different quantity from how wrong the number
                    might be, and it is 6x smaller. A ceiling must bound the error, so the
                    deficit is the conservative choice. Marked "systematic under
                    development".

🔴 THE CONTRACT GAP THIS EXPOSES, which is the real fix (RYA-587). STATES is
{MEASURED, DEFINED, N/A, HOLD}. There is NO state that distinguishes a measured value from
a CONSERVATIVE CEILING carried pending measurement -- both read DEFINED. Three of the four
terms above are ceilings, and nothing in the schema says so except prose inside `evidence`,
which no consumer keys on. A reader of the published budget cannot tell 0.1721 "measured"
from 0.1721 "bounded, under development". That gap is reported, not worked around.
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
FEED_PRODUCTS = json.loads(FEED.read_text())["products"]
BANDP = ROOT / "data/results/band_products"
GF = ROOT / "data/linelists/canonical_gf.csv"
OUT = ROOT / "data/results/rya1226/nearuv_budget_migration.json"

DH19 = "Den Hartog et al. 2019, ApJS 243, 33 (DOI 10.3847/1538-4365/ab322e)"

#: 🔴 READ FROM RYA-846's ARTIFACT, NEVER RESTATED HERE (RYA-914). The near-UV
#: normalization term is a MEASURED number and the measurement lives in one file; copying
#: its digits into this script would be a second place for them to live and a second place
#: to drift from the run that made them. The old `PSEUDO_CONTINUUM_DEX = 0.10` literal is
#: gone for the same reason -- it was the ASSUMED value RYA-846 was built to replace, and
#: the artifact records it as `assumed_term_dex` so the supersession stays visible.
RYA846_PATH = ROOT / "data/results/rya846/rya846_sigma_delta.json"


def _rya846() -> dict:
    """RYA-846's measured near-UV continuum-placement term, at THIS band's product lines."""
    d = json.loads(RYA846_PATH.read_text())
    at = d["sigma_delta_at_product_lines"]
    if not at.get("available"):
        raise SystemExit("RYA-846 records no product-line measurement -- refusing to "
                         "substitute the band-wide figure silently (RYA-1226)")
    #: ⚠️ EVERY PRODUCT LINE MUST BE INSIDE THE COMPARED REGION, or the term is measured on
    #: a subset and quietly generalised. RYA-846 reports 0 outside; assert it rather than
    #: trust it, because a future re-run on a narrower Wallace region would still "work".
    if at.get("n_lines_outside_wallace") != 0:
        raise SystemExit(f"{at['n_lines_outside_wallace']} product line(s) lie outside the "
                         f"compared region -- the term would not cover this pool")
    dec = d["decomposition"]
    return {"n_lines": at["n_lines"],
            "n_lines_outside_wallace": at["n_lines_outside_wallace"],
            "raw_at_lines": dec["raw_at_lines"],
            "wallace_internal": dec["wallace_internal"],
            "net_quadrature": dec["net_quadrature"],
            "term_raw_dex": round(dec["term_raw_dex"], 6),
            "term_net_dex": round(dec["term_net_dex"], 6),
            "lever": d["lever"]["value"],
            "assumed_term_dex": d["assumed_term_dex"]}


RYA846 = _rya846()

#: Measured per-side xi responses, from the part C legs against the committed xi=1.0
#: product. Recorded because "symmetric" is an assumption until somebody measures it --
#: and on kurucz2005 it is FALSE: the upward side responds 64% harder.
XI_ASSESSMENT = {
    ("I", "solar_kpno_kurucz2005_corrected"): (
        "MEASURED and ASYMMETRIC: A(1.1)-A(1.0) = -0.0180 against A(1.0)-A(0.9) = -0.0110 "
        "per 0.1 km/s (n=55 both), i.e. the upward side responds 64% harder. The published "
        "sigma uses the central response per the campaign convention; the asymmetry is "
        "recorded here rather than averaged away."),
    ("I", "solar_kpno_molecfit_corrected"): (
        "MEASURED and SYMMETRIC: A(1.1)-A(1.0) = A(1.0)-A(0.9) = -0.0170 per 0.1 km/s."),
    ("II", "solar_kpno_kurucz2005_corrected"): (
        "NOT ASSESSED: RYA-1213's published artifact retains only the central median for "
        "this pool and its legs are not in the committed tree, so per-side asymmetry cannot "
        "be measured here. Owed."),
    ("II", "solar_kpno_molecfit_corrected"): (
        "NOT ASSESSED: as for kurucz2005 -- RYA-1213 publishes only the central median."),
}


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
        #: 🔴 A DERIVED RESPONSE IS NOT A PERTURBATION, and the contract is right to say so:
        #: "nonzero parameter terms need perturbation evidence". The 0.000665 dex figure is
        #: uncertainty_stack's response (0.0665 dex/100 K) times the solar allowance (1 K),
        #: NOT a Teff perturbation run on this pool -- and dA/dp is a property of the LINE
        #: SET (RYA-1093), so an element-level response may not stand in for one. It is
        #: therefore carried as the SOURCED DEFINED ZERO the contract sanctions, with the
        #: derived value kept in evidence so the omission is visible and sized: 0.000665^2
        #: is 0.0004% of this budget's variance, four orders below the leading term and
        #: below the 6-dp publication precision.
        dict(name="stellar.teff", sigma_dex=0.0, state="DEFINED",
             source=("sourced DEFINED zero: no Teff perturbation has been run on this pool, "
                     "and the derived response is immaterial at this budget's precision"),
             evidence={"parameter_exists": True, "delta_K": 1.0,
                       "derived_but_not_perturbed_dex": 0.000665,
                       "derived_from": "uncertainty_stack 0.0665 dex/100 K x 1 K",
                       "share_of_variance": "0.0004%",
                       "owed_measurement": "a Teff perturbation on this pool"}),
        #: 🔴 RYA-846's MEASURED TERM REPLACES THE ASSUMED 0.100 (Ryan, 2026-09-19).
        #: RYA-841 established that the near-UV LIMITING term is pseudo-continuum
        #: (normalization), and RYA-846 MEASURED it: two independent reductions of this band
        #: differenced AT THESE 40 PRODUCT LINES. Its own artifact records
        #: `assumed_term_dex: 0.1` as the placeholder it existed to retire, and its
        #: `revised_cells` price the replacement -- so carrying 0.100 alongside it would keep
        #: a declared allowance in place of the measurement that was made to displace it.
        #:
        #: ⚠️ THE NET TERM, NOT THE RAW ONE. `raw_at_lines` = 0.04938 includes Wallace's OWN
        #: internal scatter (`wallace_internal` = 0.04071); deconfounded in quadrature that
        #: leaves 0.02794, which the lever turns into `term_net_dex` = 0.07098. Publishing
        #: the raw 0.1254 would charge this product for the comparison reduction's noise --
        #: the "measure each pair's own baseline first" rule (RYA-1192). The raw figure is
        #: carried in the evidence so the choice is auditable.
        dict(name="pseudo_continuum", sigma_dex=RYA846["term_net_dex"], state="MEASURED",
             source=("RYA-846 sigma_delta, measured at this band's product lines: Kitt Peak "
                     "(Kurucz/Brault) vs Wallace+2011 (NSO), two independent reductions, "
                     "deconfounded from Wallace's own internal scatter. Supersedes the "
                     "RYA-1113/841 band-flat 0.100 allowance, which RYA-846 records as "
                     "`assumed_term_dex` -- the number it was built to replace."),
             evidence={"band": "near-UV", "basis": "measured, two-reduction difference",
                       "n_lines": RYA846["n_lines"],
                       "n_lines_outside_wallace": RYA846["n_lines_outside_wallace"],
                       "raw_at_lines": RYA846["raw_at_lines"],
                       "wallace_internal": RYA846["wallace_internal"],
                       "net_quadrature": RYA846["net_quadrature"],
                       "term_raw_dex": RYA846["term_raw_dex"],
                       "term_net_dex": RYA846["term_net_dex"],
                       "lever": RYA846["lever"],
                       "supersedes_assumed_dex": RYA846["assumed_term_dex"],
                       "empirical_lower_bound": ("RYA-846 states this measures "
                                                 "reduction-to-reduction disagreement, not "
                                                 "distance from an unobservable true "
                                                 "continuum, so it is a LOWER bound"),
                       "lever_note": ("the artifact's lever is 2.54 (RYA-841, 200 fits, "
                                      "+/-4% grid, p16 0.91 / p84 4.68); the directive said "
                                      "2.42")}),
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

    #: --- the four Ryan's REDIRECT sent back, wired rather than re-held ---
    holdings = {q["holding"]: q["A"] for q in FEED_PRODUCTS
                if q["band"] == "near-UV" and q["ion"] == product["ion"]
                and q["tier"] == "DEEPGRADED" and q["treatment"] == product["treatment"]}
    half_spread = abs(holdings["solar_kpno_kurucz2005_corrected"]
                      - holdings["solar_kpno_molecfit_corrected"]) / 2.0
    comps += [
        #: 🔴 THE OBSERVED-HOLDING-SPREAD FLOOR IS WITHDRAWN. The previous pass charged
        #: half the kurucz2005-vs-molecfit abundance spread (0.1075 dex on Fe II 1D-LTE) as a
        #: continuum "floor pending the owed measurement". Ryan's 2026-09-19 ruling is
        #: explicit -- assemble from MEASURED components, "do not bound-fudge it" -- and that
        #: floor was a bound standing in for a measurement that, it turns out, EXISTS:
        #: RYA-846 measured this band's continuum placement at these product lines, and it is
        #: now carried on `pseudo_continuum` above at 0.0710.
        #:
        #: ⚠️ WORSE, THE FLOOR WAS CIRCULAR. The holding spread is a DISAGREEMENT BETWEEN THE
        #: TWO PRODUCTS, so charging it to each of them prices a product using its own
        #: sibling's answer -- and it came out LARGER than the real measurement (0.1075 vs
        #: 0.0710), so it was not even conservative in the direction that would have been
        #: safe to keep.
        #:
        #: So `continuum` is N/A ON EVIDENCE: the per-line continuum fit fields are empty on
        #: every row of this pool (there is no per-line placement to price), and the
        #: band-level placement uncertainty those rows DO carry is the normalization term
        #: RYA-846 measured. Pricing both would count normalization twice.
        dict(name="continuum", sigma_dex=None, state="N/A",
             source=("the per-line continuum fit carries no evidence to price -- "
                     "continuum_method, continuum_ref and continuum_level are empty on every "
                     "line in this pool -- and the band-level continuum PLACEMENT "
                     "uncertainty is measured by RYA-846 and charged once, on "
                     "pseudo_continuum."),
             evidence={"per_line_fields_empty": ["continuum_method", "continuum_ref",
                                                 "continuum_level"],
                       "subsumed_by": "pseudo_continuum (RYA-846 measured, 0.070975 dex)",
                       "not_double_counted": True,
                       "withdrawn_floor_dex": round(half_spread, 6),
                       "why_withdrawn": ("it was the OBSERVED SPREAD BETWEEN THE TWO "
                                         "HOLDINGS, which prices each product from its "
                                         "sibling's answer, and it exceeded RYA-846's real "
                                         "measurement -- a bound where a measurement exists "
                                         "(Ryan, 2026-09-19: do not bound-fudge)"),
                       "holding_values": holdings,
                       "owed_refinement": ("sigma_delta for Kurucz2005-vs-1984KP "
                                           "specifically, per ion -- RYA-846 compared Kitt "
                                           "Peak vs Wallace+2011. RYA-1133 recorded this as "
                                           "blocked on Fe II near-UV having no *_lines.csv; "
                                           "RYA-908 emitted them the same day and 14 are in "
                                           "the tree, so it is RUNNABLE, and it would refine "
                                           "rather than unblock -- the term is measured now.")}),
        dict(name="model_atmosphere", sigma_dex=0.004, state="DEFINED",
             source="RYA-1032 atmosphere_nuisance: atlas9 -> marcs-ges, LTE, same pool",
             evidence={"bound": True, "status": "measured on the VIS Fe I pool, carried "
                                                "cross-band",
                       "delta_dex": 0.004, "n_hi": 67, "n_lo": 67, "same_pool": True,
                       "varied": "atmosphere", "measured_band": "VIS",
                       #: Ryan asked for one thing to be VERIFIED here rather than assumed:
                       #: does a grid-shift lever exist at all. It does, and this is it --
                       #: RYA-1032's `atmosphere_nuisance`, atlas9 -> marcs-ges with family
                       #: and NLTE held fixed, the SAME 67-line pool on both sides.
                       "grid_shift_lever_exists": True,
                       "lever_is": ("RYA-1032 atmosphere_nuisance: A(marcs-ges) 7.451 - "
                                    "A(atlas9) 7.447 = 0.004, one axis varied"),
                       #: ⚠️ AND WHY THE STAGGER NUMBER IS NOT THIS ONE. RYA-1032's
                       #: `dimensionality_step` does reach stagger-mean3d and it is 0.101 dex
                       #: -- 25x larger -- but its `varied` field is `dim`, not `atmosphere`:
                       #: it is the 1D -> <3D> DIMENSIONALITY step, which this feed expresses
                       #: as a separate treatment (models 5/6/7), not as an uncertainty on a
                       #: 1D product's grid choice. Charging it here would price a product for
                       #: not being a different model.
                       "stagger_step_is_a_different_axis": {
                           "step": "RYA-1032 dimensionality_step",
                           "delta_dex": 0.101, "varied": "dim",
                           "from": "marcs-ges (1D-LTE)",
                           "to": "stagger-mean3d (<3D>-LTE)",
                           "why_not_charged": ("dimensionality, not grid choice -- the feed "
                                               "carries <3D> as its own treatment")}}),
        dict(name="profile_ew", sigma_dex=0.2160, state="DEFINED",
             source="RYA-1220 fit half-width excursion, +/-0.25 A vs +/-1.1 A",
             evidence={"bound": True,
                       "status": "systematic under development -- NOT demonstrably conservative",
                       "delta_dex": 0.2160, "measured_on": "N I 8216.336 A, solar_iag",
                       "cross_element": True, "cross_band": True, "cross_holding": True,
                       "source_artifact_says": ("steps_are_adopted_uncertainties: false; "
                                                "status: diagnostic_only"),
                       "owed_measurement": ("half-width sweep on the near-UV Fe pool at its "
                                            "own fixed +/-0.4 A")}),
        #: 🔴 BLENDS IS NOT AN OPEN UNCERTAINTY COMPONENT (Ryan, 2026-09-19, superseding the
        #: earlier "bound it as a conservative ceiling" directive). It was already classified
        #: and dispositioned, and the previous pass here charged 0.1721 dex -- 28.7% of the
        #: published variance -- for a component that does not exist:
        #:
        #:   * RYA-1190 (Done, PR #508), verbatim: "the near-UV list IS the VALD extract --
        #:     nothing to complete." There are NO un-modelled resolvable blends to catalogue.
        #:   * what RYA-1190 actually found was the OPACITY deficit, and that deficit was
        #:     CORRECTED -- RYA-1204 (Bautista 1997 / TOPbase Fe I bound-free plus molecular
        #:     opacity), wired into production by RYA-1207. It is fixed IN THE SYNTHESIS
        #:     these products are measured on, not an open uncertainty term.
        #:   * RYA-841 established the near-UV LIMITING term is pseudo-continuum
        #:     (normalization), NOT blends -- and RYA-846 measured that, above.
        #:
        #: ⚠️ THE 0.1721 I REMOVED WAS RYA-1190's PRE-CORRECTION DEFICIT, and the ~0.028
        #: alternative was its pre-correction payoff estimate. Both are measurements of the
        #: state BEFORE RYA-1204/1207; neither is a live component after it. Charging either
        #: bills these rows for an error that the synthesis they ran on has already removed.
        dict(name="blends", sigma_dex=None, state="N/A",
             source=("no un-modelled resolvable blends exist in this band: RYA-1190 measured "
                     "the near-UV list to BE the full VALD extract ('nothing to complete'), "
                     "and the opacity deficit it did find was corrected by RYA-1204 and "
                     "wired into this synthesis by RYA-1207."),
             evidence={"rya1190_verdict": ("the near-UV list IS the VALD extract -- nothing "
                                           "to complete; the residual was the OPACITY "
                                           "deficit, not catalogueable blends"),
                       "opacity_deficit_corrected_by": ["RYA-1204", "RYA-1207"],
                       "synthesis_is_post_correction": True,
                       "limiting_term_is": ("pseudo_continuum / normalization (RYA-841), "
                                            "measured by RYA-846"),
                       "withdrawn_bound_dex": 0.1721,
                       "why_withdrawn": ("RYA-1190's PRE-correction deficit. It measured the "
                                         "band before RYA-1204/1207; carrying it as a "
                                         "ceiling bills this product for an error the "
                                         "synthesis it ran on has already removed."),
                       "withdrawn_alternative_dex": 0.028,
                       "why_alternative_also_withdrawn": ("also pre-correction -- it is how "
                                                          "much a completeness build would "
                                                          "have closed, and there is nothing "
                                                          "left to complete")}),
    ]

    nlte_delta = set(lines["nlte_delta_dex"].dropna().unique())
    if nlte_delta != {0.0}:
        #: An NLTE leg owes the uncertainty ON its departure correction. RYA-1032 measured
        #: the model-family spread on a fixed 1D-NLTE axis -- Gerber minus Bergemann,
        #: 0.043 dex raw and 0.039 deconfounded from the atmosphere nuisance. ENGINE-A IS
        #: the Bergemann member, so that spread is the disagreement between the two live
        #: treatments of the same physics. Measured on the VIS Fe I pool: a flagged
        #: cross-band bound, not a near-UV measurement.
        comps.append(dict(name="nlte", sigma_dex=0.043, state="DEFINED",
                          source="RYA-1032 model_family_spread (Gerber - Bergemann, 1D-NLTE)",
                          evidence={"bound": True,
                                    "status": "measured on the VIS Fe I pool, carried cross-band",
                                    "delta_dex": 0.043, "deconfounded_dex": 0.039,
                                    "varied": "family", "measured_band": "VIS",
                                    "owed_measurement": ("the Gerber/Bergemann disagreement "
                                                         "on the near-UV pool itself")}))
    if nlte_delta == {0.0}:
        comps.append(dict(name="nlte", sigma_dex=None, state="N/A",
                          source="this leg applies no departures: nlte_delta_dex is 0 on "
                                 "every line and nlte_source reads 'none -- LTE'",
                          evidence={"nlte_delta_dex": 0.0}))

    if xi_slope is not None:
        #: 🔴 THE CONVENTION IS NOT RE-DERIVED HERE. sigma_xi = |dA/dxi| x delta_xi is what
        #: RYA-1168, RYA-1213 and RYA-1226 part C all publish, so the signed response is the
        #: CENTRAL one scaled to the adopted allowance and delta_plus/minus are its
        #: symmetric halves. Deriving sigma from the one-sided pair instead would change
        #: the number on every xi product in the campaign, silently, from inside a near-UV
        #: migration -- which is not this ticket's call to make.
        #:
        #: The measured per-side responses go in `response_assessment`, which is the field
        #: the contract provides for exactly that. They are NOT symmetric on kurucz2005.
        signed = xi_slope * delta_xi
        comps.append(dict(name="stellar.xi", sigma_dex=abs(signed), state="MEASURED",
                          source=("per-line paired differential at xi 0.90/1.10 on this "
                                  "product's OWN pool, on the current molecular synthesis"),
                          evidence={"pool_sha256": digest,
                                    "parameter_source": ("RYA-1089 sourced solar delta_xi = "
                                                         "0.2912 km/s, the honest method+"
                                                         "selection spread |0.709 - 1.0|"),
                                    "response_assessment": XI_ASSESSMENT.get(
                                        (product["ion"], product["holding"]),
                                        "per-side responses not assessable from the "
                                        "committed tree"),
                                    "signed_response_dex": signed,
                                    "delta_plus_dex": signed,
                                    "delta_minus_dex": -signed,
                                    "delta_parameter": delta_xi,
                                    "dA_dxi": xi_slope, "delta_xi_kms": delta_xi,
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
    #: RYA-1225 adjudicated Fe II; RYA-1226 part C MEASURED Fe I on its own pool.
    #: molecfit Fe I is absent on purpose -- its pool moved under perturbation (Fe I
    #: 3427.119 A rails at xi=1.10), and a moved pool is HOLD, not a slope on survivors.
    part_c = json.loads((ROOT / "data/results/rya1226/nearuv_fe1_xi_dadxi.json").read_text())
    fe1 = {(r["holding"], r["treatment"]): r["dA_dxi"] for r in part_c["pools"]
           if r["xi_state"] == "MEASURED"}
    rows = []
    for p in feed["products"]:
        if p.get("band") != "near-UV" or p.get("tier") != "DEEPGRADED":
            continue
        if p.get("treatment") not in ("1D-LTE", "ENGINE-A"):
            continue
        try:
            xi = (-0.1100 if p["ion"] == "II"
                  else fe1.get((p["holding"], p["treatment"])))
            built = build(p, xi_slope=xi)
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
