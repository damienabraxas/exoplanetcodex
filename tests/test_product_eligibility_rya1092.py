"""
RYA-1092 — the live-product eligibility gate, and its enforcement.

🔴 THIS FILE IS THE ENFORCEMENT, NOT A DESCRIPTION OF IT. The pools have existed since
RYA-711 and the discipline was ratified long ago; what was missing was anything that
CHECKED `products[]` against it, so ineligible products stayed live and the public feed
rendered them. `test_every_live_product_passes_the_gate` is the check. Everything else
here exists so that check cannot pass vacuously.

Every number asserted below was MEASURED on the committed feed before the gate was
written, and the threshold test asserts an INVARIANCE WINDOW rather than a value -- if new
data ever narrows that window, this file fails and forces a re-derivation instead of
letting a stale constant quietly decide (RYA-161).
"""
from __future__ import annotations

import json
import math
import statistics
import subprocess
from pathlib import Path

import pytest

from pipeline import product_eligibility as pe

ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data" / "products" / "solar" / "Fe.json"


@pytest.fixture(scope="module")
def doc():
    return json.loads(FEED.read_text())


# ── THE GATE ─────────────────────────────────────────────────────────────────────

def test_every_live_product_passes_the_gate(doc):
    """🔴 THE CI GATE. An ineligible product cannot be merged into the live list."""
    bad = {k: [f"{r.code}: {r.detail}" for r in v]
           for k, v in pe.evaluate_feed(doc).items() if v}
    assert not bad, (
        f"{len(bad)} product(s) in products[] fail the RYA-1092 eligibility gate. They "
        f"must be MOVED to quarantine[] with their reason (never deleted) -- run "
        f"`python3 scripts/rya1092_apply_eligibility_gate.py --element Fe --apply`:\n"
        + "\n".join(f"  {k}\n    " + "\n    ".join(v) for k, v in sorted(bad.items())))


def test_no_live_product_has_a_null_sigma_syst(doc):
    """The ticket's own smoke test, as a test. Stated separately from the gate because
    this is the single condition that let the incomplete Amarsi products render an error
    bar nobody computed."""
    assert [pe.key_of(p) for p in doc["products"] if p.get("sigma_syst") is None] == []


def test_no_live_product_is_outside_the_ratified_species_set(doc):
    assert [pe.key_of(p) for p in doc["products"] if not pe.is_real_species(p)] == []


def test_no_duplicate_live_cells(doc):
    """Two live products with the same identity means the feed cannot say which one IS
    the product for that cell, and a consumer picks arbitrarily."""
    assert pe.duplicate_live_cells(doc) == {}


# ── POSITIVE CONTROLS: the gate must actually bite, and must not over-bite ────────

def _solid(doc):
    """A product the gate must NOT touch: 1D-LTE, corrected holding, real syst."""
    for p in doc["products"]:
        if p["treatment"] == "1D-LTE" and p.get("sigma_syst") is not None:
            return p
    pytest.skip("no solid 1D-LTE product in the feed")


def test_staged_ineligible_product_FAILS_the_gate(doc):
    """POSITIVE CONTROL — the ticket's CRITICAL condition. A gate that never fires is
    indistinguishable from no gate at all."""
    staged = dict(_solid(doc))
    staged["sigma_syst"] = None
    reasons = pe.evaluate(staged)
    assert "SYST_INCOMPLETE" in {r.code for r in reasons}


def test_a_solid_1D_LTE_product_is_NOT_quarantined(doc):
    """POSITIVE CONTROL the other way, and a CRITICAL condition of its own: being 1D or
    LTE is not a disqualification. These are the good products."""
    assert pe.evaluate(_solid(doc), peers=[]) == ()


def test_an_uncorrected_holding_FAILS(doc):
    """`solar_kpno` is the telluric-RETAINING 1984 atlas, whitelisted as the RYA-1026
    CONTROL. A control in products[] is a science claim by placement."""
    staged = dict(_solid(doc), holding="solar_kpno")
    assert "UNCORRECTED_HOLDING" in {r.code for r in pe.evaluate(staged)}


def test_a_pre_continuum_fix_artifact_FAILS(doc):
    staged = dict(_solid(doc))
    staged["provenance"] = dict(staged["provenance"], artifact_mtime="2026-08-01T00:00:00Z")
    assert "PRE_CONTINUUM_FIX" in {r.code for r in pe.evaluate(staged)}


def test_a_missing_artifact_mtime_FAILS_rather_than_passing(doc):
    """UNKNOWN IS NOT 'AFTER' (RYA-907). A record that cannot show it postdates the fixes
    must not be admitted because the field happens to be absent."""
    staged = dict(_solid(doc))
    staged["provenance"] = {k: v for k, v in staged["provenance"].items()
                            if k != "artifact_mtime"}
    assert "PRE_CONTINUUM_FIX" in {r.code for r in pe.evaluate(staged)}


def test_a_product_with_no_holding_RAISES_rather_than_passing(doc):
    staged = {k: v for k, v in _solid(doc).items() if k != "holding"}
    with pytest.raises(pe.EligibilityError):
        pe.evaluate(staged)


def test_a_non_Fe_non_Al_species_FAILS(doc):
    staged = dict(_solid(doc), element="Ba", ion="II")
    assert "NOT_YET_DEFENSIBLE" in {r.code for r in pe.evaluate(staged)}


def test_Al_passes_in_any_ionisation_stage(doc):
    """Ryan's set is {Fe I, Fe II, Al} -- Al is named without an ion, so both stages are
    in. Writing it as ("Al", None) rather than guessing "Al I" is the difference between
    honouring the ruling and narrowing it."""
    assert pe.is_real_species(dict(_solid(doc), element="Al", ion="I"))
    assert pe.is_real_species(dict(_solid(doc), element="Al", ion="II"))


# ── THE AMARSI CONSTRAINT (Ryan's 2026-08-28 correction) ─────────────────────────

def test_SUPERSEDED_is_never_reachable_from_an_engine_name(doc):
    """🔴 Ryan's correction, pinned so it cannot be undone by a later edit.

    RYA-1092 was written believing Amarsi 3D-NLTE was superseded by the mean-⟨3D⟩ STAGGER
    redo. It is not: Amarsi is the authoritative reference STAGGER is validated AGAINST,
    the redo was a PREREQUISITE for building STAGGER, and the two coexist. The Amarsi
    products are held out because they are INCOMPLETE, never because they are dead.

    So a complete Amarsi product must pass. If supersession were ever inferred from the
    engine's name, this fails -- and the wrong verdict would be baked into the enforcement
    mechanism, which is the worst possible place for it.
    """
    complete_amarsi = dict(_solid(doc),
                           treatment="ENGINE-A-3DNLTE",
                           route="EW-3D",
                           display="EW · 3D-NLTE · Amarsi")
    codes = {r.code for r in pe.evaluate(complete_amarsi, peers=[])}
    assert "SUPERSEDED" not in codes
    assert "ANOMALOUS_SCATTER" not in codes
    # The one thing left is a SCHEMA defect the whole EW-3D route shares -- its sigma_stat
    # is a raw scatter where the rest of the feed publishes a standard error -- and it is
    # fixed by carrying `stat_basis` into the feed, not by re-measuring Amarsi.
    assert codes == {"STAT_BASIS_MISMATCH"}, codes


def test_SUPERSEDED_comes_from_the_record(doc):
    """The code still exists for genuine supersession -- it is just read, not inferred."""
    staged = dict(_solid(doc), superseded_by="some later re-run")
    assert "SUPERSEDED" in {r.code for r in pe.evaluate(staged)}


def test_the_Amarsi_products_were_COMPLETED_not_written_off(doc):
    """🔴 THE WHOLE POINT OF HOLDING THEM OUT, AND IT HAS NOW HAPPENED.

    RYA-1092 quarantined the four EW·3D-NLTE·Amarsi products as INCOMPLETE -- never as
    superseded -- because Ryan's correction is that Amarsi is the 3D reference the STAGGER
    deck is validated against, not a retired engine. RYA-1095 then completed the leg: the
    log gf was single-sourced to canonical (it had been VALD, differing on all 50 lines),
    a real `sigma_syst` was computed, and the statistical bar was put on the same footing
    as every other route. They are LIVE.

    This test asserts the resolution rather than the old state. Deleting it would erase the
    requirement; leaving it asserting `quarantine` would assert that the fix did not happen.
    """
    live = [p for p in doc["products"] if p.get("treatment") == "ENGINE-A-3DNLTE"]
    assert len(live) == 4, f"expected the 4 Amarsi products live, found {len(live)}"
    for p in live:
        assert p.get("sigma_syst") is not None, f"{pe.key_of(p)} still has no systematic"
        assert pe.evaluate(p, peers=[]) == (), pe.evaluate(p, peers=[])


def test_the_ARCHIVED_Amarsi_records_still_say_why_they_were_held_out(doc):
    """The evidence for the withdrawal survives the fix. A completed product must not
    erase the record of what was wrong with its predecessor (RYA-711)."""
    old = [p for p in (doc.get("archive") or []) + (doc.get("quarantine") or [])
           if p.get("treatment") == "ENGINE-A-3DNLTE"
           and "RYA-1092" in str(p.get("quarantine_reason", ""))]
    assert len(old) == 4, f"the 4 withdrawn Amarsi records are gone, found {len(old)}"
    for p in old:
        codes = set(p.get("quarantine_codes") or [])
        assert "SYST_INCOMPLETE" in codes
        assert "STAT_BASIS_MISMATCH" in codes
        # ⚠️ NOT anomalous, and never superseded. The first gate flagged them at 5-7x by
        # treating a raw scatter as a standard error; their real line scatter is ordinary.
        assert "ANOMALOUS_SCATTER" not in codes
        assert "SUPERSEDED" not in codes


# ── THE ANOMALY THRESHOLD: measured null, and an invariance window ───────────────

def _pre_sweep_population(doc):
    """The population as it was BEFORE the sweep: the live products plus the ones RYA-1092
    withdrew. Judging the threshold on the post-sweep list would be circular -- the
    flagged products would already be gone."""
    return list(doc["products"]) + [p for p in doc.get("quarantine", [])
                                    if "RYA-1092" in str(p.get("quarantine_reason", ""))]


def _ratios(doc):
    """Every comparable product's leave-one-out raw-scatter ratio, pre-sweep."""
    pop = _pre_sweep_population(doc)
    groups: dict = {}
    for p in pop:
        groups.setdefault(pe.peer_group_of(p), []).append(p)
    out = []
    for p in pop:
        peers = [q for q in groups[pe.peer_group_of(p)] if q is not p]
        r = pe.scatter_ratio(p, peers)
        if r is not None:
            out.append((r, p))
    return out


def test_the_two_routes_publish_DIFFERENT_statistics(doc):
    """🔴 THE DEFECT THAT MADE FOUR SOUND POOLS LOOK LIKE 7x OUTLIERS.

    `sigma_stat` is a standard error on the band route (`error_budget.py:609`,
    scatter/sqrt(n)) and a raw standard deviation on the EW-3D route
    (`band_products.py:506`, `np.std(vals, ddof=1)`). Nothing in the feed records which.
    Multiplying every product's `sigma_stat` by sqrt(n) is therefore right for 59 products
    and wrong for 4, and the wrong 4 are exactly the ones the first version of this gate
    called anomalous.

    Verified against the artifacts, not asserted: a band product's published `sigma_stat`
    equals std/sqrt(n) of its own per-line file, and the committed RYA-817 EW-3D artifact's
    published sigma equals the plain std of its 114 lines.
    """
    band = ROOT / ("data/results/band_products/FeI_4200_6910_kpno_solar_atlas_"
                   "solar_kpno_molecfit_corrected_SYNTH_GRADED_1D-LTE_lines.csv")
    ew3d = ROOT / "data/results/rya817/rya817_3dnlte_per_line.csv"
    if not (band.exists() and ew3d.exists()):
        pytest.skip("per-line artifacts not present")
    import pandas as pd
    a = pd.read_csv(band)
    a = a[a.in_aggregate == True]["abundance"].dropna()          # noqa: E712
    se = float(a.std(ddof=1)) / math.sqrt(len(a))
    # Across pools: this record was withdrawn as PRE_CONTINUUM_FIX, and what sigma_stat
    # MEANS is a property of the record, not of which pool it currently sits in.
    rec = [p for p in _pre_sweep_population(doc)
           if p["route"] == "SYNTH" and p["treatment"] == "1D-LTE" and p["band"] == "VIS"
           and p["holding"] == "solar_kpno_molecfit_corrected" and p["tier"] == "GRADED"]
    assert rec and rec[0]["sigma_stat"] == pytest.approx(se, abs=5e-4), (
        "the band route's sigma_stat must be a STANDARD ERROR")

    b = pd.read_csv(ew3d)
    b = b[(b.element == "Fe") & (b.ion == "I") & (b.band == "VIS")
          & (b.in_domain == True)]["a_3dnlte"].dropna()          # noqa: E712
    assert float(b.std(ddof=1)) == pytest.approx(0.3418, abs=1e-3), (
        "the EW-3D route's published sigma must be the raw per-line STD")

    assert pe.STAT_BASIS_BY_ROUTE["SYNTH"] == "standard_error"
    assert pe.STAT_BASIS_BY_ROUTE["EW-3D"] == "line_scatter"


def test_raw_scatter_dispatches_on_the_route(doc):
    """The conversion is route-dependent; assuming one basis is the bug this fixes."""
    band = dict(_solid(doc), route="SYNTH", sigma_stat=0.02, n_lines=100)
    ew3d = dict(_solid(doc), route="EW-3D", sigma_stat=0.20, n_lines=100)
    assert pe.raw_scatter(band) == pytest.approx(0.20)   # 0.02 * sqrt(100)
    assert pe.raw_scatter(ew3d) == pytest.approx(0.20)   # already the scatter
    # An unrecognised route ABSTAINS rather than being assumed into a basis (RYA-907).
    assert pe.raw_scatter(dict(band, route="SOMETHING-NEW")) is None


def test_the_amarsi_pools_are_NOT_anomalous_once_the_basis_is_right(doc):
    """The retraction, pinned. Their line-to-line scatter is ordinary.

    Read from the LIVE products now that RYA-1095 completed them -- the point is unchanged
    and it is the reason `ANOMALOUS_SCATTER` was the wrong code: on the corrected runs the
    scatter is 0.156-0.167 against 0.16-0.21 for the 1D pools on the same holdings.
    """
    # ⚠️ RYA-1106 — SELECT ON THE TREATMENT, NOT THE ROUTE. This read `route == "EW-3D"`,
    # which selected these products by the one field about them that was WRONG: RYA-1104
    # traced that token to a stranded ProfileFitHandler and RYA-1106 corrected it to
    # SYNTH. Keying a test on the mislabel means the test passes only while the bug
    # survives, and goes silent -- `assert amarsi` on an empty list -- the moment it is
    # fixed. The treatment is the stable identity (RYA-874: never rewritten).
    amarsi = [p for p in doc["products"] if p.get("treatment") == "ENGINE-A-3DNLTE"]
    assert amarsi, "the Amarsi products should be live"
    peers_1d = [p for p in doc["products"]
                if p["band"] == "VIS" and p["ion"] == "I" and p["route"] == "SYNTH"]
    ref = statistics.median([pe.raw_scatter(p) for p in peers_1d
                             if pe.raw_scatter(p) is not None])
    for p in amarsi:
        r = pe.raw_scatter(p) / ref
        assert 0.4 < r < 2.0, (
            f"{pe.key_of(p)} raw scatter ratio {r:.2f} -- these were reported at ~7x by "
            f"the sqrt(n) mistake and are ordinary once the basis is right")


def test_the_null_alone_justifies_the_threshold(doc):
    """🔴 A TOLERANCE NEEDS A MEASURED NULL — and here the null is ALL there is.

    After the basis correction NO live or withdrawn product is anomalous, so there is no
    flagged population to separate from and the constant cannot be justified by a gap.
    What CAN be asserted is the null: the largest scatter ratio anywhere in the feed, and
    that `ANOMALY_RATIO` sits comfortably above it. If any product ever approaches the
    threshold this fails, which forces a re-derivation instead of letting a stale constant
    decide (RYA-161).

    Saying "the criterion currently has no instance" is the honest report. A criterion
    that has never fired on real data is not evidence that it works.
    """
    rs = _ratios(doc)
    assert len(rs) >= 25, f"only {len(rs)} comparable products -- the null is too thin"
    worst = max(r for r, _ in rs)
    assert worst < pe.ANOMALY_RATIO, (
        f"a product now sits at {worst:.2f}x, at or past ANOMALY_RATIO="
        f"{pe.ANOMALY_RATIO} -- re-derive the threshold against the new population")
    assert pe.ANOMALY_RATIO / worst > 1.5, (
        f"the margin over the measured null has narrowed to "
        f"{pe.ANOMALY_RATIO / worst:.2f}x -- re-derive")
    assert not [k for k, v in pe.evaluate_feed(doc).items()
                if any(r.code == "ANOMALOUS_SCATTER" for r in v)]


def test_raw_scatter_is_sigma_stat_times_sqrt_n_on_the_BAND_route(doc):
    p = _solid(doc)
    assert pe.stat_basis_of(p) == "standard_error"
    assert pe.raw_scatter(p) == pytest.approx(p["sigma_stat"] * math.sqrt(p["n_lines"]))


def test_the_scatter_check_ABSTAINS_on_a_thin_group(doc):
    """One product must not be able to declare another anomalous. Below the peer floor the
    check returns None -- an abstention, not a pass and not a failure."""
    p = _solid(doc)
    assert pe.scatter_ratio(p, [dict(p), dict(p)]) is None


# ── THE GATE'S OWN INPUTS ────────────────────────────────────────────────────────

def test_the_identity_is_the_SAME_one_publish_product_uses():
    """Two notions of 'which product is this' would drift, and the drift would be silent
    (RYA-845). The gate does not get its own key."""
    src = (ROOT / "scripts" / "publish_product.py").read_text()
    import ast
    tree = ast.parse(src)
    found = None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
                getattr(t, "id", None) == "KEY_FIELDS" for t in node.targets):
            found = tuple(ast.literal_eval(node.value))
    assert found is not None, "publish_product.KEY_FIELDS not found"
    assert found == pe.KEY_FIELDS


def test_the_continuum_fix_cutoff_is_DERIVED_from_the_named_commits():
    """The cutoff is an instant, but it is the instant of NAMED COMMITS. Re-resolving them
    is what stops it becoming a number nobody can check -- the failure shape RYA-311 hit
    with an uncited 0.05."""
    times = []
    for sha in pe.CONTINUUM_FIX_COMMITS:
        r = subprocess.run(["git", "show", "-s", "--format=%cI", sha],
                           cwd=ROOT, capture_output=True, text=True)
        if r.returncode != 0:
            pytest.skip(f"commit {sha[:8]} not reachable in this clone")
        times.append(r.stdout.strip())
    import datetime as dt
    latest = max(dt.datetime.fromisoformat(t).astimezone(dt.timezone.utc) for t in times)
    assert latest.strftime("%Y-%m-%dT%H:%M:%SZ") == pe.CONTINUUM_FIX_CUTOFF_UTC


def test_the_telluric_fact_is_read_through_the_one_source():
    """`applied_state` is the single source for `telluric_applied`. A hand-written clean
    set is the RYA-845 shape, and the first draft of telluric_display_policy got two of
    five entries wrong by trying it."""
    src = (ROOT / "pipeline" / "product_eligibility.py").read_text()
    assert "telluric_display_policy" in src
    for name in ("solar_kpno_molecfit_corrected", "solar_harps_molecfit_corrected",
                 "solar_iag", "solar_kpno"):
        assert name not in src, (
            f"{name} is hardcoded in the gate -- the holding's state must be read from "
            f"the registry, never listed here")


# ── QUARANTINE IS A MOVE, NOT A DELETE ───────────────────────────────────────────

def test_a_product_may_leave_the_live_list_ONLY_into_a_pool_that_says_WHY(doc):
    """🔴 CRITICAL: quarantine must be loud and reversible.

    RYA-1080's value guard deliberately does not police pool membership -- "a product
    removed later is a different question", and this is that question. A product may leave
    `products[]` (that is what withdrawal IS) but only into a pool that records the reason.
    Vanishing, or landing somewhere with no reason attached, is a silent removal, and a
    removed product that cannot say why it was removed is exactly the evidence RYA-711
    refuses to destroy.
    """
    r = subprocess.run(["git", "show", "origin/main:data/products/solar/Fe.json"],
                       cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip("origin/main not available")
    before = {pe.key_of(p) for p in json.loads(r.stdout)["products"]}
    live_now = {pe.key_of(p) for p in doc["products"]}
    withdrawn = {}
    for pool in ("quarantine", "superseded", "archive"):
        for p in doc.get(pool) or []:
            withdrawn[pe.key_of(p)] = p
    for k in before - live_now:
        rec = withdrawn.get(k)
        assert rec is not None, f"product left the live list and is in NO pool: {k}"
        assert rec.get("quarantine_reason") or rec.get("superseded_reason"), (
            f"product was withdrawn with no reason recorded: {k}")


def test_every_moved_product_carries_a_reason_and_a_timestamp(doc):
    """⚠️ SEARCHED ACROSS THE WITHDRAWAL POOLS. A record the gate moved can be moved AGAIN
    -- publishing a completed replacement sends it quarantine -> archive -- and reading
    only `quarantine[]` made this assert "the sweep moved nothing" the moment a follow-up
    ticket succeeded. Same one-pool mistake this ticket had to fix in two other guards."""
    moved = [p for pool in ("quarantine", "archive", "superseded")
             for p in doc.get(pool) or []
             if "RYA-1092" in str(p.get("quarantine_reason", ""))]
    assert moved, "the sweep moved nothing -- has it been run?"
    for p in moved:
        assert p.get("quarantined_at")
        assert p.get("quarantine_codes")
        assert set(p["quarantine_codes"]) <= {
            "UNCORRECTED_HOLDING", "PRE_CONTINUUM_FIX", "SYST_INCOMPLETE",
            "NOT_YET_DEFENSIBLE", "ANOMALOUS_SCATTER", "STAT_BASIS_MISMATCH",
            "SUPERSEDED"}
