"""RYA-1227 — the xi campaign is COMPLETE: no 1D engine may read NOT_IN_CAMPAIGN.

`NOT_IN_CAMPAIGN` is not a disposition, it is an unfinished campaign. A product carrying it
publishes a `sigma_reported` with NO xi term at all -- not zero, absent -- and the audit
found them on 1D engines where microturbulence plainly applies, while the SAME engines were
MEASURED in other bands.

🔴 THE ONLY LEGITIMATE NON-MEASURED xi STATES ARE:
  * UNMEASURED    -- a genuine sub-floor hold (RYA-1224's min_paired floor), and
  * NOT_APPLICABLE -- models 5 and 6 (<3D> mean) and model 7 (full 3D) ONLY.
A 1D engine (models 1-4) reading either NOT_IN_CAMPAIGN or NOT_APPLICABLE is a defect.

⚠️ THE 1D/3D SPLIT IS DERIVED FROM `model_registry.csv` BY `model_id`, NOT FROM THE
ATMOSPHERE AND NOT FROM A LIST OF NAMES HERE. Model 7 (ENGINE-A-3DNLTE) has atmosphere
`atlas9`, exactly like models 1 and 3, so an atmosphere test would call full 3D a 1D engine
and let a real defect through. The registry is the SSOT; a name list here is a fourth place
the split can drift (RYA-1170).
"""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data/products/solar/Fe.json"
ARTIFACT = ROOT / "data/results/rya1227/codex_deep_xi_dadxi.json"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from rya1178_emit_fe_schema import (XI_BAND_RUNS, load_sources,  # noqa: E402
                                    xi_band_index, xi_min_paired)

#: Models whose published value legitimately carries no NET xi dependence.
XI_EXEMPT_MODEL_IDS = (5, 6, 7)


@pytest.fixture(scope="module")
def models():
    m = pd.read_csv(ROOT / "data/catalog/model_registry.csv", comment="#")
    return m.dropna(subset=["stored_token"])


@pytest.fixture(scope="module")
def one_d_treatments(models):
    """The 1D engines, by `model_id` from the registry."""
    out = {str(r.stored_token) for r in models.itertuples()
           if int(r.model_id) not in XI_EXEMPT_MODEL_IDS}
    assert out, "no 1D treatments resolved from the registry — the guard would be vacuous"
    return out


@pytest.fixture(scope="module")
def feed():
    """🔴 THE GENERATOR'S VERDICT, NOT THE PUBLISHED FILE (RYA-1224's lesson: when code and
    artifact disagree, the stale side is the one that passes the gate)."""
    from rya1178_emit_fe_schema import enrich
    hold, inst, models_df, xi_doc = load_sources()
    out, _ = enrich(json.loads(FEED.read_text()), hold, inst, models_df, xi_doc)
    return out


@pytest.fixture(scope="module")
def published():
    return json.loads(FEED.read_text())


@pytest.fixture(scope="module")
def artifact():
    return json.loads(ARTIFACT.read_text())


# ── the gate ──────────────────────────────────────────────────────────────────
def test_no_product_anywhere_reads_not_in_campaign(feed):
    """THE TICKET'S GATE. Feed-wide, not just on the 11 this ticket measured."""
    left = [f"{p['ion']}/{p['band']}/{p['tier']}/{p['treatment']}"
            for p in feed["products"] if p.get("xi_state") == "NOT_IN_CAMPAIGN"]
    assert not left, f"{len(left)} product(s) still read NOT_IN_CAMPAIGN: {left}"


def test_a_1d_engine_is_never_exempt_from_xi(feed, one_d_treatments):
    """🔴 NOT_APPLICABLE ON A 1D ENGINE IS A DEFECT, not a disposition. xi is a 1D fudge
    for the velocity field; a 1D model does not resolve that field, so the term applies."""
    bad = [f"{p['ion']}/{p['band']}/{p['tier']}/{p['treatment']}={p['xi_state']}"
           for p in feed["products"]
           if p["treatment"] in one_d_treatments
           and p.get("xi_state") in ("NOT_APPLICABLE", "NOT_IN_CAMPAIGN")]
    assert not bad, bad


def test_the_only_states_a_1d_engine_may_carry(feed, one_d_treatments):
    """Stated positively, so a NEW state cannot slip in unexamined."""
    seen = {p["xi_state"] for p in feed["products"]
            if p["treatment"] in one_d_treatments}
    assert seen, "no 1D-engine products in the feed — this guard would be vacuous"
    assert seen <= {"MEASURED", "UNMEASURED"}, seen


def test_every_exempt_product_really_is_a_3d_model(feed, one_d_treatments):
    """The converse control: NOT_APPLICABLE must belong to models 5/6/7 and nothing else."""
    na = {p["treatment"] for p in feed["products"]
          if p.get("xi_state") == "NOT_APPLICABLE"}
    assert na, "nothing is NOT_APPLICABLE — the converse guard would be vacuous"
    assert not (na & one_d_treatments), na


def test_the_campaign_moved_no_published_abundance(feed, published):
    """🔴 validate-don't-tune (RYA-161). The xi term may change; A may not."""
    def ident(q):
        return tuple(str(q.get(k) or "") for k in
                     ("element", "ion", "band", "instrument", "holding", "tier",
                      "selector", "route", "treatment"))
    was = {ident(p): p["A"] for p in published["products"]}
    assert len(was) == len(published["products"]), "identity key is not unique"
    moved = [f"{ident(p)} {was[ident(p)]} -> {p['A']}"
             for p in feed["products"] if p["A"] != was[ident(p)]]
    assert not moved, moved


# ── the artifact this ticket adds ─────────────────────────────────────────────
def test_the_artifact_emits_only_each_decks_own_treatment(artifact, models):
    """🔴 THE DELIBERATE NARROWING, ASSERTED. The same legs measured 1D-LTE and ENGINE-A on
    these pools, but those cells are ALREADY MEASURED and a second derivative would move a
    published sigma out of scope. Emitting one must fail here, not be noticed later."""
    own = set(artifact["emits_own_treatment_only"].values())
    assert own == {"synth-1D-LTE-gerber", "ENGINE-B-NLTE"}
    got = {p["treatment"] for p in artifact["pools"]}
    assert got <= own, f"the artifact emits a treatment it should not: {got - own}"
    #: and the base-treatment derivatives it declined to emit are RECORDED, so the
    #: narrowing is auditable rather than merely claimed.
    skipped = {r["treatment"] for r in artifact["base_treatments_measured_but_not_emitted"]}
    assert skipped <= {"1D-LTE", "ENGINE-A"}
    assert skipped, "no declined derivative recorded — the audit trail is missing"


def test_the_artifact_floor_is_the_derived_floor(artifact):
    """One floor, derived from the runs, never a literal (RYA-1224)."""
    assert artifact["min_paired"] == xi_min_paired()


def test_every_pool_is_measured_on_its_own_graded_or_deep_pool(artifact):
    """No Reference pool may answer for a Codex/Deep cell -- that is the borrowing
    RYA-1224 removed, and `~/scratch/rya1213_xi` is full of reference-pool legs whose
    names differ only in a tier segment."""
    for p in artifact["pools"]:
        assert p["tier"] in ("GRADED", "DEEPGRADED"), p
        assert p["tier"] in p["pool"], f"the pool string must name its tier: {p['pool']}"


def test_a_measured_pool_clears_the_floor_and_a_held_one_says_why(artifact):
    floor = artifact["min_paired"]
    for p in artifact["pools"]:
        if p["xi_state"] == "MEASURED":
            assert p["n_paired"] >= floor, p
            assert p["sigma_xi"] > 0, p
        else:
            assert p["xi_state"] == "UNMEASURED", p
            assert f"below min_paired={floor}" in p["xi_note"], p


def test_the_pairing_was_proven_not_assumed(artifact):
    """Worktree isolation is not a proof of pairing (RYA-1120). Every pool's legs were
    verified from their own `xi_run.json` stamps."""
    for p in artifact["pools"]:
        assert p["xi_pairing_verified"] is True, p
        assert p["span_kms"] == pytest.approx(0.2, abs=1e-6), p


def test_the_paired_route_carries_its_own_justification(artifact):
    """🔴 The paired differential exists because differencing aggregates flipped a SIGN on
    the first cell RYA-1213 measured. Every pool records both, so a future sign flip is
    visible rather than silently averaged away."""
    for p in artifact["pools"]:
        assert "difference_of_aggregates" in p
        assert isinstance(p["sign_disagreement_with_aggregates"], bool)
        if p["sign_disagreement_with_aggregates"]:
            assert p["dA_dxi"] * p["difference_of_aggregates"] < 0, p


# ── composition: five band-run files, one key space ───────────────────────────
def test_the_band_key_is_unique_across_every_band_run():
    """🔴 THE GUARD THAT MAKES THE COMPOSITION SAFE. Five files now feed one index; if two
    claimed the same (ion, holding, tier, treatment, band) a product would silently take
    one. `xi_band_index` refuses -- run it rather than trust it."""
    assert ARTIFACT in XI_BAND_RUNS, "the RYA-1227 artifact is not wired into XI_BAND_RUNS"
    idx = xi_band_index()          # raises SystemExit on a duplicate key
    assert idx, "the band index is empty"


def test_every_artifact_pool_reaches_a_real_product(artifact, feed):
    """No orphan derivative: a measured pool with no product is a pool measured for
    nothing, and usually means the selector did not match the published cell."""
    live = {(p["ion"], p["holding"], p["tier"], p["treatment"], p["band"])
            for p in feed["products"]}
    for p in artifact["pools"]:
        key = (p["ion"], p["holding"], p["tier"], p["treatment"], p["band"])
        assert key in live, f"the artifact measures a pool with no product: {key}"


def test_each_artifact_pool_actually_serves_its_product(artifact, feed):
    """And the reverse direction: the product must end up MEASURED from this run, with the
    artifact's own sigma. A derivative that reaches the index but loses a precedence
    contest would leave the cell exactly as unfixed as before."""
    by_key = {(p["ion"], p["holding"], p["tier"], p["treatment"], p["band"]): p
              for p in feed["products"]}
    for pool in artifact["pools"]:
        key = (pool["ion"], pool["holding"], pool["tier"], pool["treatment"], pool["band"])
        prod = by_key[key]
        if pool["xi_state"] != "MEASURED":
            continue
        assert prod["xi_state"] == "MEASURED", (key, prod["xi_state"])
        assert prod["sigma_xi"] == pytest.approx(pool["sigma_xi"], abs=1e-9), key
        assert prod["dA_dxi_dex_per_kms"] == pool["dA_dxi"], key
