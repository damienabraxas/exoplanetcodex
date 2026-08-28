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
    assert pe.evaluate(complete_amarsi, peers=[]) == ()


def test_SUPERSEDED_comes_from_the_record(doc):
    """The code still exists for genuine supersession -- it is just read, not inferred."""
    staged = dict(_solid(doc), superseded_by="some later re-run")
    assert "SUPERSEDED" in {r.code for r in pe.evaluate(staged)}


def test_the_quarantined_Amarsi_products_carry_INCOMPLETE_not_superseded(doc):
    """The four Amarsi products are out on SYST_INCOMPLETE (and anomalous scatter), and
    the reason text must say so -- a wrong reason is a wrong record even when the
    placement is right."""
    q = [p for p in doc.get("quarantine", [])
         if p.get("treatment") == "ENGINE-A-3DNLTE"
         and "RYA-1092" in str(p.get("quarantine_reason", ""))]
    assert len(q) == 4, f"expected the 4 Amarsi EW·3D-NLTE products, found {len(q)}"
    for p in q:
        codes = set(p.get("quarantine_codes") or [])
        assert "SYST_INCOMPLETE" in codes
        assert "SUPERSEDED" not in codes


# ── THE ANOMALY THRESHOLD: measured null, and an invariance window ───────────────

def _pre_sweep_population(doc):
    """The population as it was BEFORE the sweep: the live products plus the ones RYA-1092
    withdrew. Judging the threshold on the post-sweep list would be circular -- the
    flagged products would already be gone."""
    return list(doc["products"]) + [p for p in doc.get("quarantine", [])
                                    if "RYA-1092" in str(p.get("quarantine_reason", ""))]


def _ratio_populations(doc, statistic):
    """(null, flagged) leave-one-out ratios under an arbitrary per-product statistic.

    "Flagged" is defined by the INDEPENDENT criterion -- sigma_syst is null -- not by the
    scatter rule itself, so this cannot become a test of the threshold against its own
    output.
    """
    pop = _pre_sweep_population(doc)
    groups: dict = {}
    for p in pop:
        groups.setdefault(pe.peer_group_of(p), []).append(p)
    null, flagged = [], []
    for p in pop:
        mine = statistic(p)
        if not mine:
            continue
        peers = [q for q in groups[pe.peer_group_of(p)] if q is not p]
        others = [statistic(q) for q in peers]
        others = [x for x in others if x]
        if len(others) < pe.ANOMALY_MIN_PEERS:
            continue
        r = mine / statistics.median(others)
        (flagged if p.get("sigma_syst") is None else null).append(r)
    return null, flagged


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


def test_the_chosen_constant_would_have_MISFIRED_on_sigma_stat(doc):
    """Why the statistic is raw scatter and not the standard error, measured.

    ⚠️ NOT "they overlap on sigma_stat" -- an earlier draft asserted that and it is FALSE
    on this feed. What is true, and what is asserted here, is sharper: the separation on
    `sigma_stat` is 1.8x narrower, and `ANOMALY_RATIO` falls BELOW the sigma_stat null's
    maximum. The same constant on the standard error would have flagged a legitimate
    product -- which is the concrete cost of using a statistic that carries 1/sqrt(n).
    """
    null_se, flag_se = _ratio_populations(doc, lambda p: p.get("sigma_stat"))
    null_raw, flag_raw = _ratio_populations(doc, pe.raw_scatter)
    assert null_se and flag_se and null_raw and flag_raw
    sep_se = min(flag_se) / max(null_se)
    sep_raw = min(flag_raw) / max(null_raw)
    assert sep_raw > sep_se, (
        f"raw scatter must separate the populations more cleanly than the standard "
        f"error (raw {sep_raw:.2f}x vs SE {sep_se:.2f}x)")
    assert max(null_se) > pe.ANOMALY_RATIO, (
        f"the sigma_stat null tops out at {max(null_se):.3f}, which must exceed "
        f"ANOMALY_RATIO={pe.ANOMALY_RATIO} -- that is the measured cost of the wrong "
        f"statistic, and if it stops being true the docstring's argument has changed")
    assert max(null_raw) < pe.ANOMALY_RATIO


def test_the_null_and_the_invariance_window(doc):
    """🔴 A TOLERANCE NEEDS A MEASURED NULL.

    The null is every comparable product that nobody has challenged; the flagged set is
    the products the gate calls anomalous. What is asserted is not the threshold but the
    GAP: any value strictly between the null's maximum and the flagged set's minimum gives
    the identical verdict, so no number in that range can be tuned. `ANOMALY_RATIO` must
    sit inside it. If the window ever closes, this test fails.
    """
    rs = _ratios(doc)
    # 37 measured: only the four peer groups with >= 3 OTHER pools take part. The
    # thirteen 3-product groups abstain by design -- two peers is not a population.
    assert len(rs) >= 30, f"only {len(rs)} comparable products -- the null is too thin"
    flagged = [r for r, p in rs if r > pe.ANOMALY_RATIO]
    null = [r for r, p in rs if r <= pe.ANOMALY_RATIO]
    assert len(flagged) == 4, f"expected the 4 Amarsi pools, got {len(flagged)}"
    assert max(null) < min(flagged), "null and flagged populations must not overlap"
    assert max(null) < pe.ANOMALY_RATIO < min(flagged)
    # The window is wide, not marginal: a factor of two with nothing inside it.
    assert min(flagged) / max(null) > 2.0, (   # measured 2.69x
        f"the separation has narrowed to {min(flagged) / max(null):.2f}x "
        f"(null max {max(null):.2f}, flagged min {min(flagged):.2f}) -- re-derive the "
        f"threshold rather than keeping this constant")


def test_raw_scatter_is_sigma_stat_times_sqrt_n(doc):
    p = _solid(doc)
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

def test_nothing_was_deleted_by_the_sweep(doc):
    """🔴 CRITICAL: quarantine must be loud and reversible. Every product origin/main had
    live is still SOMEWHERE in this feed."""
    r = subprocess.run(["git", "show", "origin/main:data/products/solar/Fe.json"],
                       cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip("origin/main not available")
    before = {pe.key_of(p) for p in json.loads(r.stdout)["products"]}
    after = set()
    for pool in ("products", "quarantine", "superseded", "archive"):
        after |= {pe.key_of(p) for p in doc.get(pool) or []}
    assert not (before - after), f"products vanished entirely: {sorted(before - after)}"


def test_every_moved_product_carries_a_reason_and_a_timestamp(doc):
    moved = [p for p in doc.get("quarantine", [])
             if "RYA-1092" in str(p.get("quarantine_reason", ""))]
    assert moved, "the sweep moved nothing -- has it been run?"
    for p in moved:
        assert p.get("quarantined_at")
        assert p.get("quarantine_codes")
        assert set(p["quarantine_codes"]) <= {
            "UNCORRECTED_HOLDING", "PRE_CONTINUUM_FIX", "SYST_INCOMPLETE",
            "NOT_YET_DEFENSIBLE", "ANOMALOUS_SCATTER", "SUPERSEDED"}
