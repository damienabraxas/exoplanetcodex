"""RYA-515 — per-line provenance for the live Fe pool: the audit, pinned.

AUDIT. RYA-515 asks a provenance question, and the answer is only worth anything if it
stays answerable. These pin the three things that would otherwise rot silently:

  * the resolver that ties a product to its per-line evidence -- it must disambiguate two
    same-field artifacts on VALUE, not give up on them, because the two products this
    ticket exists to audit (Fe II VIS 7.966 and 7.617) are exactly the ambiguous pair;
  * the finding that twelve live products CONTRADICT their own committed evidence -- a
    guard that goes quiet the moment someone re-derives them, which is the point;
  * the shape of the Fe I anchor pool -- well-behaved under a median, not under a mean.

Nothing here asserts a product VALUE is right. RYA-161: validate, don't tune.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
AUD = ROOT / "data/audit/rya515_fe_perline"
PROV = AUD / "perline_provenance.csv"
REACH = AUD / "reachability.csv"
FE2 = AUD / "fe2_blend_census.csv"
FE1 = AUD / "fe1_control.csv"
BP = ROOT / "data/results/band_products"


def _load(p):
    if not p.exists():
        pytest.skip(f"{p.name} absent")
    return pd.read_csv(p)


@pytest.fixture(scope="module")
def prov():
    return _load(PROV)


@pytest.fixture(scope="module")
def reach():
    return _load(REACH)


# --- the resolver -----------------------------------------------------------------

def test_resolver_disambiguates_on_n_lines_rather_than_giving_up(reach):
    """The HARPS Fe II VIS pair is the whole reason the gate moved.

    Two artifacts share every field of these products (RYA-1191 committed a
    re-measurement at stem 4200_6908 beside the shipped 4200_6910). Bailing out on
    len(cands) > 1 made 7.966 and 7.617 -- the two numbers this ticket names -- unreachable.
    """
    harps2 = reach[(reach.ion == "II") & (reach.band == "VIS")
                   & reach.holding.str.contains("harps", na=False)]
    assert len(harps2) == 2, "expected exactly the two HARPS Fe II VIS products"
    assert harps2.perline_reachable.all(), "the ambiguous pair must resolve"
    assert harps2.resolution.str.startswith("resolved on n_lines").all()


def test_resolver_never_resolves_two_artifacts_to_one_product(reach):
    """A product resolving is only meaningful if the file it resolved to is unique to it."""
    got = reach[reach.perline_reachable & (reach.perline_artifact.astype(str) != "")]
    dupes = got.perline_artifact.value_counts()
    assert (dupes == 1).all(), f"artifact stamped onto >1 product: {dupes[dupes > 1].to_dict()}"


def test_resolved_products_agree_with_their_artifact_on_n_lines(reach):
    """The gate is value-equality; assert it actually held, rather than trusting the label."""
    for _, r in reach[reach.perline_reachable].iterrows():
        t = pd.read_csv(BP / f"{r.perline_artifact}")
        n = int(t["in_aggregate"].sum()) if "in_aggregate" in t.columns else len(t)
        assert n == int(r.n_lines), f"{r.perline_artifact}: {n} vs product {r.n_lines}"


# --- the finding ------------------------------------------------------------------

@pytest.mark.xfail(
    strict=True,
    reason="RYA-515 6: all six live Fe II VIS products still publish the pre-curation "
           "9-line/3-line pool. RYA-1203 regenerates Fe.json from the post-1191 "
           "band_products; when it lands this XPASSes and the marker must come off.",
)
def test_no_live_fe2_vis_product_stands_on_a_curated_out_line(reach):
    """6, written as the invariant rather than as the defect.

    Asserting the DEFECT (len(stale) == 6) would have to be deleted the moment someone
    fixed it, and would go red on the branch that did the fixing -- punishing the fix.
    The invariant is the durable statement: a live product's line count must match the
    evidence committed beside it. Strict xfail records that it does not hold today and
    turns the fix into a loud XPASS rather than a silent one.
    """
    vis2 = reach[(reach.ion == "II") & (reach.band == "VIS")]
    assert len(vis2) == 6, "expected 6 live Fe II VIS products"
    disagree = vis2[vis2.resolution.str.startswith("n_lines MISMATCH")
                    | vis2.resolution.str.startswith("resolved on n_lines")]
    assert disagree.empty, (
        "these products publish a line pool their own committed evidence does not contain: "
        + ", ".join(f"{r.holding}/{r.treatment} (A={r.A})" for _, r in disagree.iterrows())
    )


def test_only_the_two_harps_products_still_aggregate_4303(prov):
    """The curated line must not be aggregated anywhere else in the 1464-row table."""
    hit = prov[(abs(prov.wavelength_air_A - 4303.170) < 0.02) & prov.in_aggregate.astype(bool)]
    assert set(hit.treatment) == {"1D-LTE", "ENGINE-A"}
    assert (hit.holding == "solar_harps_molecfit_corrected").all()
    assert len(hit) == 2, f"4303.170 aggregated in {len(hit)} products, expected 2"
    assert (hit.line_A > 9.9).all(), "the artifact is a +2.5 dex outlier or it is not this line"


def test_curation_registry_actually_carries_the_ruling():
    """§3e claims the ruling is live in code. Assert that, not the prose."""
    p = ROOT / "data/catalog/line_curation_exclusions.csv"
    if not p.exists():
        pytest.skip("curation registry absent")
    reg = pd.read_csv(p, comment="#")
    row = reg[(reg.species == "Fe II") & (abs(reg.wavelength_air_A - 4303.170) < 0.02)]
    assert len(row) == 1, "Fe II 4303.170 must be in the curation registry"
    assert str(row.iloc[0].ticket).strip(), "a curated exclusion with no ticket is not a ruling"


# --- the control that makes 3a a finding -------------------------------------------

def test_red_chi2_does_not_discriminate_fe2_from_the_fe1_control():
    """3a. Reporting Fe II's red_chi2 as a defect is only refutable against the Fe I pool."""
    f2, f1 = _load(FE2), _load(FE1)
    assert (f2.red_chi2 > 5).mean() > 0.95 and (f1.red_chi2 > 5).mean() > 0.95, (
        "if the Fe I control ever stops sharing Fe II's red_chi2 distribution, 3a is wrong"
    )


def test_cog_regime_separates_the_pools_where_chi2_did_not():
    """3c. The discriminator that survived: every Fe II line is saturated, few Fe I are."""
    f2, f1 = _load(FE2), _load(FE1)
    assert (f2.obs_depth > 0.70).mean() == 1.0
    assert (f1.obs_depth > 0.70).mean() < 0.20
    assert f2.obs_rew.median() - f1.obs_rew.median() > 0.5, "REW separation collapsed"


# --- the anchor -------------------------------------------------------------------

def test_fe1_anchor_is_well_behaved_under_a_median_but_not_a_mean(prov):
    """5. The ticket's question, answered conditionally -- the tail is real."""
    a = prov[(prov.ion == "I") & prov.in_aggregate.astype(bool)].line_A.dropna()
    assert len(a) > 1000
    assert a.skew() > 2.0, "the right tail is the whole point of the conditional answer"
    core = a[(a > 7.0) & (a < 8.0)]
    assert abs(core.median() - a.median()) < 0.01, "median must be resistant to the tail"
    assert a.mean() - a.median() > 0.02, "and the mean must not be"


# --- provenance is stamped on EVERY line, which is the deliverable ------------------

def test_every_row_carries_engine_route_and_correction(prov):
    for c in ("route", "treatment", "nlte_source", "holding", "band", "ion"):
        assert prov[c].notna().all(), f"{c} unstamped on some rows"
    # A missing delta is only legitimate where the correction was NOT SERVED -- that is a
    # stamped fact, not a gap. Anywhere else an unstamped delta is an unprovenanced line.
    blank = prov[prov.nlte_delta_dex.isna()]
    assert blank.excluded_reason.astype(str).str.startswith("ENGINE-A-NOT-SERVED").all()
    assert not blank.in_aggregate.astype(bool).any(), "a NOT-SERVED line entered an aggregate"


def test_selector_separates_the_two_harps_fe1_vis_products(reach):
    """The one artifact on disk is the DEEPGRADED product's (median 7.535), not LOCALRENORM's.

    Resolving on tier alone gave both products the same 109 lines. They publish 7.339 and
    7.535 -- 0.196 dex apart -- so one of the two was being shown another product's evidence.
    """
    two = reach[(reach.ion == "I") & (reach.band == "VIS")
                & (reach.holding == "solar_harps_molecfit_corrected")
                & (reach.tier == "DEEPGRADED") & (reach.treatment == "1D-LTE")]
    assert len(two) == 2, "expected the DEEPGRADED / DEEPGRADED-LOCALRENORM pair"
    assert two.perline_reachable.sum() == 1, "exactly one of the pair has evidence on disk"
    got = two[two.perline_reachable].iloc[0]
    assert abs(got.A - 7.535) < 1e-6, "the resolved one must be the 7.535 product"



def test_engine_a_is_the_only_treatment_carrying_a_real_nlte_delta(prov):
    """4. The stamp defect, pinned: NLTE-labelled treatments record no departure."""
    ag = prov[prov.in_aggregate.astype(bool)]
    ea = ag[ag.treatment == "ENGINE-A"]
    assert (ea.nlte_delta_dex.abs() > 0).all(), "ENGINE-A must carry real per-line deltas"
    for t in ("ENGINE-B-NLTE", "synth-mean3D-NLTE-gerber-stagger"):
        rows = ag[ag.treatment == t]
        if rows.empty:
            continue
        assert (rows.nlte_delta_dex == 0).all(), f"{t} gained a delta -- §4 needs revisiting"
