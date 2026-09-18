"""RYA-1224 — the min_paired floor is a property of the RULE, so it binds every tier.

The RYA-1221 audit found the DEEPGRADED tier carrying six ALIASED Fe II VIS products: a
BORROWED xi derivative published as though it were the product's own. Three of them sat on
two-line pools. Reference had already declined those same pools as UNMEASURED -- so one
physical pool carried two honesty standards depending only on which tier label was read.

🔴 THE CAUSE WAS NOT FE II AND NOT THE TIER. `xi_terms` had two routes to a derivative. The
band-keyed route read each run's own verdict, and every band-keyed run declares and applies
`min_paired`, so that route was right by inheritance. The RYA-1120 campaign declares no
floor and the campaign route never checked one -- so the floor existed in the ARTIFACTS and
not in the RULE, and whichever route answered decided how honest the answer was. These tests
pin the floor to the rule, and pin the identity gate on the route that reads across tiers.
"""
import copy
import json
import sys
from collections import defaultdict
from pathlib import Path

import pytest

#: 🔴 ABSOLUTE, NEVER RELATIVE — other modules in this suite chdir (RYA-1178's note).
ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data/products/solar/Fe.json"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from rya1178_emit_fe_schema import (XI_BAND_RUNS, load_sources, xi_index,  # noqa: E402
                                    xi_min_paired, xi_min_paired_hold,
                                    xi_same_artifact_entry, xi_terms)


@pytest.fixture(scope="module")
def feed():
    """🔴 THE GENERATOR'S VERDICT, NOT THE PUBLISHED FILE.

    RYA-1116's lesson, the hard way: when code and artifact disagree, the STALE side is the
    one that passes the gate. `data/products/solar/Fe.json` cannot currently be rewritten by
    its own emitter -- RYA-587's `assert_publication_feed` refuses all 160 rows, because
    `enrich` restamps `generated_at`/`code_commit` on every run and a restamped row is no
    longer `in legacy`, so it loses the legacy exemption and is then asked for an
    uncertainty budget it has never carried (reported on RYA-1224; NOT this ticket's fix).
    Reading the published feed here would therefore assert the rule against a file written
    before the rule existed, and every guard below would pass while proving nothing.
    """
    from rya1178_emit_fe_schema import enrich
    hold, inst, models, xi_doc = load_sources()
    out, _ = enrich(json.loads(FEED.read_text()), hold, inst, models, xi_doc)
    return out


@pytest.fixture(scope="module")
def published():
    return json.loads(FEED.read_text())


@pytest.fixture(scope="module")
def idx(feed):
    return xi_index(load_sources()[3], feed)


def test_the_rule_moves_no_published_abundance(feed, published):
    """🔴 THE RYA-161 FIREWALL, AND THE TICKET'S FIRST CRITICAL. The xi classification and
    the uncertainty may change; A may not. Compared on the RYA-1127 identity key, because
    (ion, band, holding, treatment) is NOT unique over this feed."""
    def ident(q):
        return tuple(str(q.get(k) or "") for k in
                     ("element", "ion", "band", "instrument", "holding", "tier",
                      "selector", "route", "treatment"))
    was = {ident(p): p["A"] for p in published["products"]}
    assert len(was) == len(published["products"]), "identity key is not unique"
    moved = [f"{ident(p)} {was[ident(p)]} -> {p['A']}"
             for p in feed["products"] if p["A"] != was[ident(p)]]
    assert not moved, moved


def _sha(p):
    return (p.get("provenance") or {}).get("sha256")


def _cross_tier_groups(feed):
    """Artifacts published at more than one tier, keyed on the artifact HASH.

    ⚠️ INTERSECT THE LEDGERS, NEVER ASK "CONTAINS" (RYA-1194's lesson). The pairing is
    established by `provenance.sha256`, which is a statement about bytes; tier, line_set and
    grade are all labels that differ BECAUSE of the tier and would pair nothing.
    """
    by = defaultdict(list)
    for p in feed["products"]:
        if _sha(p):
            by[_sha(p)].append(p)
    return {s: g for s, g in by.items() if len({q["tier"] for q in g}) > 1}


# ── scope A: the floor is derived from the runs, and disagreement is refused ───
def test_the_floor_is_read_from_the_runs_and_every_run_declares_it():
    """The POSITIVE control. A guard that only ever refuses proves nothing (RYA-1178)."""
    floor = xi_min_paired()
    assert floor == 3
    for path in XI_BAND_RUNS:
        doc = json.loads(path.read_text())
        assert doc.get("min_paired") == floor, f"{path.name} disagrees with the derived floor"


def _runs(tmp_path, floors):
    out = []
    for i, mp in enumerate(floors):
        q = tmp_path / f"run{i}.json"
        body = {"ticket": f"RYA-TEST-{i}", "pools": []}
        if mp is not None:
            body["min_paired"] = mp
        q.write_text(json.dumps(body))
        out.append(q)
    return tuple(out)


def test_runs_that_disagree_on_the_floor_are_refused_not_reconciled(tmp_path, monkeypatch):
    """🔴 A per-run floor is not a floor. Picking the min would publish the laxest run's
    standard everywhere; picking the max would silently re-grade pools a run had passed."""
    monkeypatch.setattr("rya1178_emit_fe_schema.XI_BAND_RUNS", _runs(tmp_path, [3, 4]))
    with pytest.raises(SystemExit, match="disagree on `min_paired`"):
        xi_min_paired()


def test_a_run_that_declares_no_floor_is_refused_and_no_default_is_invented(tmp_path,
                                                                            monkeypatch):
    """🔴 AN OMITTED ARGUMENT IS NOT A NEUTRAL ONE (RYA-1196). Defaulting a missing floor to
    3 would let a run that never applied one pass as though it had."""
    monkeypatch.setattr("rya1178_emit_fe_schema.XI_BAND_RUNS", _runs(tmp_path, [3, None]))
    with pytest.raises(SystemExit, match="declares no `min_paired`"):
        xi_min_paired()


# ── scope B: the rule, on the live feed, at every tier ────────────────────────
def test_no_sub_floor_pool_carries_a_derivative_at_any_tier(feed):
    """THE TICKET'S RULE. `n_paired <= n_lines` always, so `n_lines < floor` proves the
    product's own pool cannot clear it -- without needing a run on a pool no run answered."""
    floor = xi_min_paired()
    sub = [p for p in feed["products"]
           if p.get("n_lines") is not None and p["n_lines"] < floor]
    #: 🔴 THE CONTROL THAT KEEPS THIS FROM GOING VACUOUS (RYA-1080's lesson: a guard whose
    #: two sides converge proves nothing). If the feed ever holds no sub-floor pool, this
    #: test is asserting over an empty set and must say so rather than pass quietly.
    assert sub, "no sub-floor pool in the feed — this guard would be vacuous"
    for p in sub:
        assert p["xi_state"] == "UNMEASURED", f"{p['tier']}/{p['holding']}/{p['treatment']}"
        assert p.get("sigma_xi") is None
        assert "dA_dxi_dex_per_kms" not in p, (
            f"a sub-floor pool must not publish a readable derivative: "
            f"{p['tier']}/{p['holding']}/{p['treatment']}")
        assert f"below min_paired={floor}" in p["xi_note"]


def test_the_floor_is_tier_uniform_and_not_a_fe_ii_patch(feed):
    """🔴 THE TICKET'S CRITICAL: the rule applied only to Fe II rather than uniformly.
    Every tier that HOLDS a sub-floor pool must hold it, and the states must not depend on
    which tier label the same pool was published under."""
    floor = xi_min_paired()
    sub = [p for p in feed["products"]
           if p.get("n_lines") is not None and p["n_lines"] < floor]
    tiers = {p["tier"] for p in sub}
    assert len(tiers) > 1, f"only one tier holds a sub-floor pool ({tiers}) — cannot show uniformity"
    assert {p["xi_state"] for p in sub} == {"UNMEASURED"}
    #: and it is not an Fe-II rule: the Fe I NIR ENGINE-A pool is held on the same grounds.
    assert {p["ion"] for p in sub} == {"I", "II"}, {p["ion"] for p in sub}


def test_no_product_publishes_a_borrowed_derivative(feed):
    """ALIASED is retired by MEASUREMENT, not by renaming. Nothing may carry the campaign's
    'measured on a pool of N lines; this product has M' note, at any state."""
    for p in feed["products"]:
        assert p.get("xi_state") != "ALIASED", f"{p['tier']}/{p['holding']}/{p['treatment']}"
        assert "was measured on a pool of" not in (p.get("xi_note") or ""), (
            f"a borrowed-pool note survived on {p['tier']}/{p['holding']}/{p['treatment']}")


def test_one_artifact_published_at_two_tiers_gets_one_xi_verdict(feed):
    """🔴 THE DEFECT ITSELF. `provenance.sha256` EQUAL means one CSV, one measurement --
    `line_set_resolved` differs only because RYA-1127 derives it from `tier` at read time.
    Two xi verdicts on one artifact is the inconsistency this ticket exists to remove."""
    groups = _cross_tier_groups(feed)
    assert len(groups) >= 18, f"only {len(groups)} cross-tier artifacts — the guard has lost its subject"
    for sha, g in groups.items():
        states = {q["xi_state"] for q in g}
        assert len(states) == 1, (
            f"{sha[:12]} {g[0]['ion']} {g[0]['band']} {g[0]['holding']} {g[0]['treatment']}: "
            f"one artifact, {len(states)} xi verdicts " +
            ", ".join(f"{q['tier']}={q['xi_state']}" for q in g))


def test_a_measured_product_traces_to_a_run_on_its_own_pool(feed):
    """THE TICKET'S SMOKE TEST: n_paired matches its n_lines, not a sibling's."""
    for p in feed["products"]:
        note = p.get("xi_note") or ""
        if p.get("xi_state") == "MEASURED" and "OWN POOL" in note:
            assert f"all {p['n_lines']} of its {p['n_lines']} lines paired" in note, note
            assert "`provenance.sha256` is EQUAL" in note, note


# ── scope C: the cross-tier route's gate is IDENTITY, and it mutation-tests ───
def _deep_1dlte(feed):
    return next(p for p in feed["products"]
                if p["ion"] == "II" and p["band"] == "VIS" and p["tier"] == "DEEPGRADED"
                and p["treatment"] == "1D-LTE")


def test_the_cross_tier_route_serves_the_real_case(feed, idx):
    """The POSITIVE control for scope C, so the three refusals below are not vacuous."""
    e = xi_same_artifact_entry(_deep_1dlte(feed), idx, feed)
    assert e is not None and e["xi_state"] == "MEASURED"
    assert e["_sibling_tier"] == "REFERENCE"
    assert e["n_paired"] == _deep_1dlte(feed)["n_lines"]


@pytest.mark.parametrize("field,value,why", [
    ("A", 9.99, "a different abundance is a different measurement"),
    ("n_lines", 7, "a different pool size is a different pool"),
    ("n_excluded", 4, "a different exclusion set is a different pool"),
])
def test_the_cross_tier_route_declines_when_the_two_rows_are_not_one_product(
        feed, idx, field, value, why):
    """🔴 MATCHING NUMBERS IS A PROXY; THE HASH IS THE CLAIM. RYA-1112 has a Fe II cell
    carrying the Fe I anchor's own digits on a different scale, so resemblance has to be
    refused even when the hash agrees."""
    f = copy.deepcopy(feed)
    prod = _deep_1dlte(f)
    sib = next(q for q in f["products"]
               if _sha(q) == _sha(prod) and q["tier"] != prod["tier"])
    sib[field] = value
    assert xi_same_artifact_entry(prod, idx, f) is None, why


def test_the_cross_tier_route_declines_on_a_different_artifact(feed, idx):
    """The hash is the gate. Break it and the route must fall through to the campaign."""
    f = copy.deepcopy(feed)
    prod = _deep_1dlte(f)
    for q in f["products"]:
        if q is not prod and _sha(q) == _sha(prod):
            q["provenance"]["sha256"] = "0" * 64
    assert xi_same_artifact_entry(prod, idx, f) is None


def test_the_cross_tier_route_declines_a_partially_paired_pool(feed, idx, monkeypatch):
    """🔴 STRICTER THAN THE SAME-TIER PATH, DELIBERATELY. A run that paired only some of
    the pool is the run's own business within its tier; reached ACROSS the tier key it is
    no longer provably this product's derivative, so it is declined."""
    prod = _deep_1dlte(feed)
    band = dict(idx["band"])
    k = next(key for key, e in band.items()
             if key[0] == "II" and key[4] == "VIS" and key[2] == "REFERENCE"
             and key[3] == "1D-LTE" and key[1] == prod["holding"])
    band[k] = {**band[k], "n_paired": band[k]["n_paired"] - 1}
    assert xi_same_artifact_entry(prod, {**idx, "band": band}, feed) is None


# ── scope D: ordering — an exemption is not a hold ────────────────────────────
def test_a_sub_floor_full_3d_product_stays_not_applicable(idx):
    """🔴 THE ORDER IS THE PHYSICS. A term that does not apply cannot be 'unmeasured'.
    Were the floor tested first, a two-line full-3D product would be relabelled from an
    exemption into an owed measurement."""
    prod = {"ion": "I", "band": "VIS", "holding": "solar_iag", "tier": "REFERENCE",
            "treatment": "ENGINE-A-3DNLTE", "route": "SYNTH", "n_lines": 2}
    out = xi_terms(prod, idx)
    assert out["xi_state"] == "NOT_APPLICABLE"
    assert out["sigma_xi"] == 0.0


def test_the_hoisted_disposition_changed_no_live_verdict(feed):
    """RYA-1224 moved `xi_disposition` above the band read. No band-keyed run covers a
    NOT_APPLICABLE treatment, so the move must be inert on the live feed -- measured here
    rather than asserted in a comment."""
    na = {p["treatment"] for p in feed["products"] if p["xi_state"] == "NOT_APPLICABLE"}
    covered = {k[3] for k in xi_index(load_sources()[3], feed)["band"]}
    assert not (na & covered), f"a band-keyed run covers a NOT_APPLICABLE treatment: {na & covered}"


# ── scope E: the index contract ───────────────────────────────────────────────
def test_a_partial_index_is_refused_by_name(feed, idx):
    """🔴 A DEFAULT HERE WOULD BE THE VACUOUS-GUARD FAILURE MODE. A caller that hand-rolls
    `{"band": ...}` must not quietly skip the floor and the identity route."""
    with pytest.raises(KeyError, match="min_paired"):
        xi_terms(_deep_1dlte(feed), {"band": idx["band"]})


def test_the_hold_is_a_function_of_the_pool_and_the_floor_alone():
    """The unit, so the rule can be read without the feed."""
    floor = 3
    assert xi_min_paired_hold({"n_lines": 3}, floor) is None
    assert xi_min_paired_hold({"n_lines": 99}, floor) is None
    assert xi_min_paired_hold({"n_lines": None}, floor) is None
    for n in (0, 1, 2):
        out = xi_min_paired_hold({"n_lines": n}, floor)
        assert out["xi_state"] == "UNMEASURED" and out["sigma_xi"] is None
        assert "dA_dxi" not in out
