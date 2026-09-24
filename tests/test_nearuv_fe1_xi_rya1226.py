"""RYA-1226 part C -- the near-UV Fe I slope, measured rather than adopted.

🔴 WHAT THIS PROTECTS. Two things that each look like a smaller version of the other:
a pool that MOVED under perturbation, and a pool that was merged carelessly. The first is
a HOLD; the second is an arithmetic error that IMITATES one. Both bit during this ticket.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

pd = pytest.importorskip("pandas")

from pipeline.paired_differential import paired_differential  # noqa: E402
from pipeline.xi_pairing import assert_pair                    # noqa: E402

LEGS = ROOT / "data/results/rya1226/legs"
ART = ROOT / "data/results/rya1226/nearuv_fe1_xi_dadxi.json"
KUR = "solar_kpno_kurucz2005_corrected"
MOL = "solar_kpno_molecfit_corrected"


@pytest.fixture(scope="module")
def art():
    return json.loads(ART.read_text())


def _row(art, holding, treat):
    return next(p for p in art["pools"]
                if p["holding"] == holding and p["treatment"] == treat)


def _legs(tag, holding, treat):
    stem = f"FeI_3000_3780_kpno_solar_atlas_{holding}_SYNTH_DEEPGRADED_{treat}_lines.csv"
    return (pd.read_csv(LEGS / f"nuvI_{tag}_dg_xi0.9000" / stem),
            pd.read_csv(LEGS / f"nuvI_{tag}_dg_xi1.1000" / stem))


@pytest.mark.parametrize("tag", ["kur", "kpmf"])
def test_the_pairing_is_proven_from_the_legs_own_stamps(tag):
    """Worktree isolation is not a proof of pairing (RYA-1178 A)."""
    assert_pair(LEGS / f"nuvI_{tag}_dg_xi0.9000", LEGS / f"nuvI_{tag}_dg_xi1.1000")


@pytest.mark.parametrize("treat", ["1D-LTE", "ENGINE-A"])
def test_kurucz2005_pairs_the_products_own_pool_exactly_and_is_measured(art, treat):
    r = _row(art, KUR, treat)
    assert r["n_paired"] == r["product_n_lines"] == 55
    assert r["pool_moved"] is False
    assert r["xi_state"] == "MEASURED"
    assert r["dA_dxi"] == -0.17


@pytest.mark.parametrize("treat", ["1D-LTE", "ENGINE-A"])
def test_molecfit_pool_MOVED_so_it_is_held_not_restricted_to_survivors(art, treat):
    """🔴 A moved pool is HOLD at the caller, never a slope computed on what survived."""
    r = _row(art, MOL, treat)
    assert r["n_paired"] == 53 and r["product_n_lines"] == 54
    assert r["pool_moved"] is True
    assert r["xi_state"] == "UNMEASURED"
    assert r["dA_dxi"] is None and r["sigma_xi"] is None
    assert r["dA_dxi_unpublishable_float"] is not None, "carry the float, do not publish it"
    assert "3427.119" in r["xi_note"]


def test_the_line_that_leaves_the_molecfit_pool_is_named_correctly():
    """The disposition rests on a specific line and a specific leg -- verify both."""
    lo, hi = _legs("kpmf", MOL, "1D-LTE")
    agg = lambda d: {round(w, 3) for w, a in zip(d["wavelength_air_A"], d["in_aggregate"]) if a}
    assert 3427.119 in agg(lo), "the minus leg keeps it"
    assert 3427.119 not in agg(hi), "the plus leg is the one that drops it"
    row = hi[hi["wavelength_air_A"].round(3) == 3427.119]
    assert "NON-MINIMUM" in str(row["excluded_reason"].iloc[0])


def test_a_raw_merge_imitates_a_moved_pool_and_that_is_why_the_helper_is_used():
    """🔴 CONTROL for my own first-pass error.

    A merge that does not filter to `in_aggregate` pairs excluded rows too -- they still
    carry abundances -- and reports 57 against a 55-line product, which reads exactly like
    a moved pool. paired_differential filters; a raw merge does not.
    """
    lo, hi = _legs("kur", KUR, "1D-LTE")
    raw = lo.merge(hi, on=["wavelength_air_A", "ep_eV"], suffixes=("_lo", "_hi"))
    naive_n = int((raw["abundance_hi"] - raw["abundance_lo"]).notna().sum())
    assert naive_n == 57, naive_n
    assert paired_differential(hi, lo).n_paired == 55


def test_the_slope_was_measured_not_adopted_from_the_reference_tier(art):
    """The ticket forbids same-pool adoption for Fe I; prove the number is not REFERENCE's."""
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    ref = {(p["holding"], p["treatment"]): p.get("dA_dxi_dex_per_kms")
           for p in feed["products"]
           if p["band"] == "near-UV" and p["ion"] == "I" and p["tier"] == "REFERENCE"}
    for treat in ("1D-LTE", "ENGINE-A"):
        assert _row(art, KUR, treat)["dA_dxi"] != ref[(KUR, treat)]


def test_the_run_is_recorded_as_post_molecular_with_its_provenance(art):
    p = art["provenance"]
    assert p["commit"] == "97cde328"
    assert "use_molecules" in p["post_molecular"]
    assert "none touches pipeline/" in p["uncommitted_at_run_time"]
    assert art["xi_span_kms"] == 0.20 and art["delta_xi_kms"] == 0.2912


def test_no_abundance_moved(art):
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    live = {(p["ion"], p["holding"], p["tier"], p["treatment"], p["band"]): p["A"]
            for p in feed["products"]}
    for r in art["pools"]:
        assert live[("I", r["holding"], "DEEPGRADED", r["treatment"], "near-UV")] == r["A"]
