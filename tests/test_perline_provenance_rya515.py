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


def test_the_harps_stem_drift_is_gone_and_nothing_resolves_by_tiebreak(reach):
    """RYA-1203 retired the duplicate stem, so the tiebreak should no longer be REACHED.

    RYA-515 had to disambiguate two artifacts of one product on n_lines because HARPS
    carried both 4200_6908 and 4200_6910. That was a workaround for a defect, not a
    feature: the fix is one artifact per product. The tiebreak code stays (it is what
    would CATCH a recurrence) but nothing may need it, and the HARPS Fe II VIS products
    -- the ones this audit exists to judge -- must resolve plainly.
    """
    harps2 = reach[(reach.ion == "II") & (reach.band == "VIS")
                   & reach.holding.str.contains("harps", na=False)]
    assert len(harps2) >= 2, "the HARPS Fe II VIS products must still be live"
    assert harps2.perline_reachable.all(), "every HARPS Fe II VIS product needs evidence"
    assert (harps2.resolution == "resolved").all(), harps2.resolution.tolist()
    assert not reach.resolution.str.startswith("resolved on n_lines").any(), (
        "a product needed the n_lines tiebreak again -- two artifacts share its key, "
        "which is the RYA-515 6 stem drift returning")


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

def test_no_live_fe2_vis_product_stands_on_a_curated_out_line(reach):
    """6, and it now HOLDS -- the marker came off in RYA-1203.

    RYA-515 could only record this as a strict xfail: all six live Fe II VIS products
    published a pre-curation pool their own committed evidence no longer contained.
    RYA-1203 re-derived them from the post-RYA-1191 artifacts, so the invariant is true
    and is asserted plainly. The pool is nine now, not six, because the Gerber LTE leg was
    added on each of the three holdings.
    """
    vis2 = reach[(reach.ion == "II") & (reach.band == "VIS")]
    assert len(vis2) == 9, f"expected 9 live Fe II VIS products, got {len(vis2)}"
    disagree = vis2[vis2.resolution.str.startswith("n_lines MISMATCH")
                    | vis2.resolution.str.startswith("resolved on n_lines")]
    assert disagree.empty, (
        "these products publish a line pool their own committed evidence does not contain: "
        + ", ".join(f"{r.holding}/{r.treatment} (A={r.A})" for _, r in disagree.iterrows()))
    assert (vis2.n_lines == 8).sum() == 6, "the six 8-line legs are 1D-LTE + Gerber-LTE"



def test_no_live_product_aggregates_the_curated_line(prov):
    """The inverse of what RYA-515 could assert. 4303.170 is curated out
    (data/catalog/line_curation_exclusions.csv, ruled RYA-1191); RYA-515 6 found the
    HARPS pair still aggregating it because the products predated the ruling. RYA-1203
    re-derived them, so the registry and the feed finally agree."""
    hit = prov[(abs(prov.wavelength_air_A - 4303.170) < 0.02) & prov.in_aggregate.astype(bool)]
    assert hit.empty, (
        "these products still aggregate a line the repo has ruled invalid: "
        + ", ".join(f"{r.holding}/{r.treatment}" for _, r in hit.iterrows()))
    #: 🔴 AND IT LEAVES NO ROW AT ALL. The curation drops the line before the artifact is
    #: written, so the per-line table does not record that it was considered and rejected --
    #: `data/catalog/line_curation_exclusions.csv` is the only trace. That is a real gap in
    #: per-line provenance (an exclusion you cannot see from the evidence), recorded here
    #: rather than asserted away. Pinned so that if the pipeline ever starts emitting
    #: curated rows as excluded-but-present, this fails and the note gets updated.
    seen = prov[abs(prov.wavelength_air_A - 4303.170) < 0.02]
    assert seen.empty, ("4303.170 now appears in the per-line table -- if it is emitted as "
                        "excluded-but-present that is an improvement; update this note")


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
    #: the gap narrowed 0.029 -> 0.019 as RYA-1203 widened the pool from 1336 to 2139
    #: lines; the point is that it is NON-ZERO and one-signed, not its exact size.
    assert a.mean() - a.median() > 0.010, "and the mean must not be"


# --- provenance is stamped on EVERY line, which is the deliverable ------------------


def test_every_row_carries_engine_route_and_correction(prov):
    for c in ("route", "treatment", "nlte_source", "holding", "band", "ion"):
        assert prov[c].notna().all(), f"{c} unstamped on some rows"
    #: A BLANK DELTA IS NOT THE SAME AS A ZERO ONE, and RYA-1203 made the difference
    #: visible: the Gerber legs leave `nlte_delta_dex` NaN where the older kurucz2005 runs
    #: wrote 0.0. Both mean "no ADDITIVE per-line correction" -- the Gerber departures are
    #: applied inside the radiative transfer and the per-line schema has no field for them
    #: (RYA-515 4) -- so the row is legitimately blank. What is NOT legitimate is a blank
    #: on a leg that corrects additively, which is ENGINE-A and only ENGINE-A.
    blank = prov[prov.nlte_delta_dex.isna()]
    bad = blank[(blank.treatment == "ENGINE-A") & blank.in_aggregate.astype(bool)]
    assert bad.empty, f"ENGINE-A rows in an aggregate with no delta: {len(bad)}"


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
    """4. ENGINE-A corrects additively and says so per line; nothing else does.

    Stated as an EXCLUSIVE claim rather than a per-treatment list, so a new treatment
    cannot quietly start carrying departures without this failing -- RYA-1203 added three
    (ENGINE-B-NLTE on five new cells, synth-1D-LTE-gerber on eight) and the claim held.
    """
    ag = prov[prov.in_aggregate.astype(bool)]
    ea = ag[ag.treatment == "ENGINE-A"]
    assert len(ea) > 100 and (ea.nlte_delta_dex.abs() > 0).all(), (
        "ENGINE-A must carry a real per-line delta on every aggregated row")
    others = ag[ag.treatment != "ENGINE-A"]
    carrying = others[others.nlte_delta_dex.abs() > 0]
    assert carrying.empty, (
        "a non-ENGINE-A treatment now carries a per-line departure: "
        + ", ".join(sorted(set(carrying.treatment))))