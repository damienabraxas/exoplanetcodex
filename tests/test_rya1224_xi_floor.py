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


# ── scope F: the legacy exemption keys on science, and ONLY the clock is ignored ──
#
# 🔴 RYA-1224 FOLD-IN. `assert_publication_feed` exempted an existing row via whole-dict
# equality, but `enrich` restamps `generated_at`/`code_commit` on every product every run.
# All 160 Fe rows therefore lost the exemption to a ride-along and the feed became
# unwritable by its own emitter. The comparison now drops exactly those two fields.
#
# ⚠️ THE DANGEROUS FAILURE IS THE OPPOSITE ONE: an exclusion list that grows until the
# guard cannot see a real move. Every test below the positive control is a REFUSAL test.
from pipeline.uncertainty_contract import (LEGACY_RIDE_ALONG_STAMPS,  # noqa: E402
                                           UncertaintyError,
                                           assert_publication_feed)


def _legacy_row(published):
    """A real published row, so these tests bind to the true product shape."""
    return copy.deepcopy(next(p for p in published["products"]
                              if p["ion"] == "II" and p["band"] == "VIS"
                              and p["tier"] == "DEEPGRADED" and p["treatment"] == "1D-LTE"))


def _feeds(published, **changes):
    was = _legacy_row(published)
    now = copy.deepcopy(was)
    now.update(changes)
    return {"products": [now]}, {"products": [was]}


def test_the_ride_along_list_is_exactly_the_two_clock_fields():
    """🔴 PINNED BY VALUE. The list is the whole safety argument: these two are readings
    about the emit, everything else is a statement about the science. A third name here
    would blind the guard, so it is asserted rather than described."""
    assert LEGACY_RIDE_ALONG_STAMPS == ("generated_at", "code_commit")


def test_a_clock_tick_alone_is_not_a_product_change(published):
    """THE POSITIVE CONTROL, and the actual bug: without this the whole feed is refused.
    A guard that only ever refuses proves nothing (RYA-1178)."""
    doc, prev = _feeds(published, generated_at="2030-01-01T00:00:00Z",
                       code_commit="0" * 40)
    assert_publication_feed(doc, previous=prev)      # must not raise


def test_both_clock_fields_are_ignored_individually(published):
    """Each on its own, so a passing pair cannot hide one field still tripping it."""
    for field, value in (("generated_at", "2030-01-01T00:00:00Z"),
                         ("code_commit", "0" * 40)):
        doc, prev = _feeds(published, **{field: value})
        assert_publication_feed(doc, previous=prev)  # must not raise


#: 🔴 RYAN'S GUARD-STILL-WORKS PROOF, WIDENED FROM A TO EVERY FIELD A CHANGE CAN RIDE IN.
#: 0.001 dex is deliberately below anything that reads as significant: the guard must catch
#: the smallest real move, not just an obvious one.
#:
#: ⚠️ THE xi-LAYER FIELDS ARE HERE TOO, AND THEY ARE NOT EXEMPT -- they are judged by the
#: scope-G route instead of waved through, and each value below is one that route must also
#: refuse. `sigma_reported`/`sigma_xi` break the recomputed quadrature; `xi_state` claims a
#: verdict the row's own evidence contradicts. The ONE xi-layer change the route does accept
#: on a fully evidenced row is prose-only (`xi_note`), which has its own scope-G test --
#: keeping it in this list would assert the opposite of the ruling.
@pytest.mark.parametrize("field,value,why", [
    ("A", 7.910, "a 0.001 dex nudge to a published abundance"),
    ("n_lines", 7, "a pool that lost a line"),
    ("n_excluded", 1, "a changed exclusion set"),
    ("sigma_reported", 0.190173, "a 1e-6 nudge to the published bar"),
    ("sigma_stat", 0.1640, "a changed random term"),
    ("sigma_syst", 0.0650, "a changed systematic term"),
    ("sigma_xi", 0.071345, "a changed xi term that breaks the quadrature"),
    ("xi_state", "ALIASED", "a verdict the row's own sigma_xi contradicts"),
    ("grade", "Reference Grade", "a relabelled tier"),
    ("tier", "REFERENCE", "a re-identified product"),
    ("A", None, "a withdrawn abundance"),
])
def test_a_real_change_still_trips_the_guard(published, field, value, why):
    doc, prev = _feeds(published, **{field: value})
    with pytest.raises(UncertaintyError, match="refuses incomplete live products"):
        assert_publication_feed(doc, previous=prev)


def test_a_real_change_trips_it_even_alongside_a_clock_tick(published):
    """🔴 THE COMBINATION IS THE REALISTIC CASE -- a re-emit ALWAYS ticks the clock, so a
    guard that only catches a lone change would never fire in production."""
    doc, prev = _feeds(published, A=7.910, generated_at="2030-01-01T00:00:00Z",
                       code_commit="0" * 40)
    with pytest.raises(UncertaintyError, match="refuses incomplete live products"):
        assert_publication_feed(doc, previous=prev)


def test_a_brand_new_product_is_never_exempt(published):
    """The exemption is RETENTION of an existing row, not a blanket pass."""
    doc, prev = _feeds(published)
    doc["products"][0]["holding"] = "solar_some_new_holding"
    with pytest.raises(UncertaintyError, match="refuses incomplete live products"):
        assert_publication_feed(doc, previous=prev)


def test_a_dropped_field_is_a_change_not_a_match(published):
    """Excluding a field from the comparison is not the same as tolerating its removal."""
    doc, prev = _feeds(published)
    del doc["products"][0]["sigma_xi"]
    with pytest.raises(UncertaintyError, match="refuses incomplete live products"):
        assert_publication_feed(doc, previous=prev)


def test_a_product_claiming_its_own_budget_is_always_checked(published):
    """`"uncertainty" not in product` is half the exemption: a row that DOES claim a
    canonical budget must be validated even if it is otherwise byte-identical to legacy."""
    doc, prev = _feeds(published, uncertainty={"schema": "not-the-real-one"})
    with pytest.raises(UncertaintyError, match="refuses incomplete live products"):
        assert_publication_feed(doc, previous=prev)


# ── scope G: the xi-layer correction route (Ryan's ruling, 2026-09-18) ────────
#
# 🔴 WHY THE ROUTE EXISTS. RYA-1224 corrected `xi_state`/`sigma_xi` on real products -- and
# those are precisely the fields the legacy exemption must refuse a change to without
# evidence. So the gate was preserving a KNOWN DEFECT: the demonstrably wrong number (a
# derivative borrowed from a pool the product does not have) stayed live because it was
# already in the file, while the honest one could not be published. The route lets an
# xi-CONFINED correction through on named evidence, and nothing else.
#
# ⚠️ EVERY TEST BELOW THE TWO POSITIVE CONTROLS IS A REFUSAL. An exemption route is only as
# good as what it turns away.
from pipeline.uncertainty_contract import (XI_IMMUTABLE_FIELDS,  # noqa: E402
                                           XI_LAYER_FIELDS,
                                           xi_layer_correction_problems)


def _corrected_hold(row):
    """`row` re-published as an honest HOLD: the RYA-1224 ALIASED -> UNMEASURED move."""
    out = copy.deepcopy(row)
    out["xi_state"] = "UNMEASURED"
    out["sigma_xi"] = None
    out.pop("dA_dxi_dex_per_kms", None)
    out.pop("xi_source", None)
    out["xi_note"] = "below min_paired=3, no derivative is published for it"
    out["sigma_syst_components"] = {"published_syst": out.get("sigma_syst"), "sigma_xi": None}
    out["sigma_syst_complete"] = round(out["sigma_syst"], 6)
    out["sigma_reported"] = round(
        (out["sigma_stat"] ** 2 + out["sigma_syst_complete"] ** 2) ** 0.5, 6)
    return out


def _corrected_measured(row, sigma_xi=0.05, dadxi=-0.2):
    """`row` re-published ASSERTING a derivative, naming a real artifact."""
    out = copy.deepcopy(row)
    out["xi_state"] = "MEASURED"
    out["sigma_xi"] = sigma_xi
    out["dA_dxi_dex_per_kms"] = dadxi
    out["xi_source"] = "data/results/rya1213/reference_xi_dadxi.json"
    out["xi_note"] = "measured on this product's own pool"
    out["sigma_syst_components"] = {"published_syst": out.get("sigma_syst"),
                                    "sigma_xi": sigma_xi}
    out["sigma_syst_complete"] = round((out["sigma_syst"] ** 2 + sigma_xi ** 2) ** 0.5, 6)
    out["sigma_reported"] = round(
        (out["sigma_stat"] ** 2 + out["sigma_syst_complete"] ** 2) ** 0.5, 6)
    return out


def test_a_withdrawal_needs_no_artifact_because_nothing_is_claimed(published):
    """POSITIVE CONTROL 1 -- the ALIASED -> UNMEASURED direction. Evidence is required to
    ASSERT a number, not to stop asserting one; demanding an artifact to withdraw would
    keep the borrowed value live, which is the defect."""
    was = _legacy_row(published)
    assert xi_layer_correction_problems(_corrected_hold(was), was) == []


def test_an_assertion_qualifies_when_it_names_a_resolvable_artifact(published):
    """POSITIVE CONTROL 2 -- the ALIASED -> MEASURED direction."""
    was = _legacy_row(published)
    assert xi_layer_correction_problems(_corrected_measured(was), was) == []


def test_the_route_has_no_opinion_on_a_change_outside_the_xi_layer(published):
    """🔴 `None`, NOT `[]`. The route must DECLINE to judge a non-xi change so the full
    contract still decides it -- returning `[]` there would be a blanket pass."""
    was = _legacy_row(published)
    for field, value in (("A", 7.910), ("n_lines", 7), ("n_excluded", 1),
                         ("grade", "Reference Grade"), ("sigma_stat", 0.9),
                         ("provenance", {"sha256": "0" * 64})):
        now = _corrected_hold(was)
        now[field] = value
        assert xi_layer_correction_problems(now, was) is None, field


def test_an_identical_row_is_not_an_xi_correction(published):
    """No change at all is the plain legacy exemption's business, not this route's."""
    was = _legacy_row(published)
    assert xi_layer_correction_problems(copy.deepcopy(was), was) is None


@pytest.mark.parametrize("break_it,expect", [
    (lambda r: r.pop("xi_source"), "must name its artifact"),
    (lambda r: r.update(xi_source="   "), "must name its artifact"),
    (lambda r: r.update(xi_source="data/results/rya1213/no_such_file.json"), "does not resolve"),
    (lambda r: r.pop("dA_dxi_dex_per_kms"), "no readable derivative"),
    (lambda r: r.update(xi_state="ALIASED"), "xi_state is 'ALIASED'"),
    (lambda r: r.update(sigma_reported=0.5), "not quadrature(sigma_stat"),
    (lambda r: r.update(sigma_syst_complete=0.5), "not quadrature(sigma_syst, sigma_xi)"),
])
def test_an_unevidenced_assertion_is_refused(published, break_it, expect):
    was = _legacy_row(published)
    now = _corrected_measured(was)
    break_it(now)
    problems = xi_layer_correction_problems(now, was)
    assert problems, "the route accepted an unevidenced assertion"
    assert any(expect in p for p in problems), problems


def test_a_hold_that_keeps_the_borrowed_number_is_refused(published):
    """🔴 THE RYA-1224 DEFECT ITSELF, as a refusal: the honest label with the dishonest
    number still attached. `update_xi_budget` really did leave this behind."""
    was = _legacy_row(published)
    now = _corrected_hold(was)
    now["dA_dxi_dex_per_kms"] = -0.245          # the borrowed float, left readable
    problems = xi_layer_correction_problems(now, was)
    assert any("readable dA_dxi_dex_per_kms" in p for p in problems), problems


def test_a_hold_may_not_claim_to_be_measured(published):
    was = _legacy_row(published)
    now = _corrected_hold(was)
    now["xi_state"] = "MEASURED"
    problems = xi_layer_correction_problems(now, was)
    assert any("not a hold" in p for p in problems), problems


def test_the_immutable_fields_are_checked_BY_NAME_not_only_by_the_subset(published,
                                                                         monkeypatch):
    """🔴 THE SECOND LINE OF DEFENCE, MUTATION-TESTED. `A` is outside `XI_LAYER_FIELDS`, so
    the subset test already stops it and the by-name check is unreachable today -- which is
    exactly why it must be proven live. Widen the layer to include `A` (the plausible
    future mistake) and the named RYA-161 check must still refuse."""
    was = _legacy_row(published)
    now = _corrected_hold(was)
    now["A"] = round(was["A"] + 0.001, 6)
    monkeypatch.setattr("pipeline.uncertainty_contract.XI_LAYER_FIELDS",
                        XI_LAYER_FIELDS + XI_IMMUTABLE_FIELDS)
    problems = xi_layer_correction_problems(now, was)
    assert problems, "widening the layer blinded the guard to an A move"
    assert any("RYA-161" in p and "A moved" in p for p in problems), problems


def test_the_xi_layer_list_never_contains_a_firewall_field():
    """A static guard on the two lists, so a widening cannot be silent."""
    assert not set(XI_LAYER_FIELDS) & set(XI_IMMUTABLE_FIELDS)
    assert XI_IMMUTABLE_FIELDS == ("A", "n_lines", "n_excluded")


def test_end_to_end_an_xi_correction_writes_and_an_A_nudge_does_not(published):
    """The route through the real gate, both directions."""
    was = _legacy_row(published)
    assert_publication_feed({"products": [_corrected_hold(was)]},
                            previous={"products": [was]})        # must not raise
    bad = _corrected_hold(was)
    bad["A"] = round(was["A"] + 0.001, 6)
    with pytest.raises(UncertaintyError, match="refuses incomplete live products"):
        assert_publication_feed({"products": [bad]}, previous={"products": [was]})


def test_a_refused_xi_correction_says_it_was_the_xi_route_that_refused(published):
    """Diagnostics: a row that tried this route and failed must not be reported as though
    it had simply never had a budget."""
    was = _legacy_row(published)
    now = _corrected_measured(was)
    now["xi_source"] = "data/results/rya1213/no_such_file.json"
    with pytest.raises(UncertaintyError, match="xi-layer correction refused"):
        assert_publication_feed({"products": [now]}, previous={"products": [was]})


def test_an_ambiguous_predecessor_disqualifies_the_route(published):
    """🔴 NO UNIQUE PREDECESSOR, NO ROUTE. Two legacy rows on one identity key means the
    correction cannot be attributed, and picking one is the neighbour-matching error
    RYA-1206 exists to refuse."""
    was = _legacy_row(published)
    twin = copy.deepcopy(was)
    twin["sigma_stat"] = 0.999          # same identity key, different content
    with pytest.raises(UncertaintyError, match="refuses incomplete live products"):
        assert_publication_feed({"products": [_corrected_hold(was)]},
                                previous={"products": [was, twin]})


def test_a_prose_only_note_change_is_permitted_on_an_evidenced_row(published):
    """🔴 DOCUMENTED, NOT ACCIDENTAL. This is what let the five already-held rows recover the
    run attribution the floor had swallowed. It is the ONLY xi-layer change the route accepts
    with nothing else moving -- and it is still gated: the row must keep a resolvable
    `xi_source`, a readable derivative, and sigma arithmetic that recomputes. A note is prose;
    every number in the row is checked."""
    was = _legacy_row(published)
    now = copy.deepcopy(was)
    now["xi_note"] = was["xi_note"] + " (attribution restored)"
    assert xi_layer_correction_problems(now, was) == []
    assert_publication_feed({"products": [now]}, previous={"products": [was]})


def test_a_note_change_does_not_carry_a_broken_number_with_it(published):
    """The companion refusal: prose is permitted, prose PLUS a silent number is not."""
    was = _legacy_row(published)
    now = copy.deepcopy(was)
    now["xi_note"] = "reworded"
    now["sigma_reported"] = round(was["sigma_reported"] + 0.000001, 6)
    problems = xi_layer_correction_problems(now, was)
    assert problems and any("sigma_reported" in p for p in problems), problems
