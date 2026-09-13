"""The live tracker must SEE every committed band product (RYA-990).

RYA-935 built the tracker to derive its status from products on disk precisely so it
could not go stale the way a hand-typed dashboard does. But deriving from disk only
helps if the reader recognises what is on disk: RYA-984 began tagging artifacts with a
SELECTOR (`_DEEPGRADED`, `_FROMEW[-GRADED|-UNGRADED]`) and the tracker's filename
pattern had no place for one, so it silently dropped every product carrying a tag.

Two merged VIS Fe I legs were invisible for that reason -- RYA-984's Kitt Peak deep run
and RYA-991's HARPS deep run -- while the tracker displayed the 55-line shallow run as
the only VIS synth product. Nothing failed; the page simply under-reported reality,
which is the exact failure the derive-don't-type design existed to prevent.

So these tests pin the INVARIANT (no committed band product is unseen), not the two
filenames that happened to break it.
"""
from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BAND_PRODUCTS = ROOT / "data" / "results" / "band_products"


def _tracker():
    """Import the generator without running it (it reads the registries in main())."""
    spec = importlib.util.spec_from_file_location(
        "rya935_live_status", ROOT / "scripts" / "rya935_live_status.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# The real registry keys (scripts/measure_band_ew._INSTRUMENT_HOLDINGS). Held literally
# rather than imported: importing that module resolves the Kitt Peak atlas and exits when
# it is absent, and this is a test of the PARSER, which takes these as plain inputs.
def _registries() -> tuple[set[str], set[str]]:
    """The real keys, READ FROM THE SOURCE rather than hand-copied.

    🔴 RYA-1009 — the hand-copied version DRIFTED and hid real products. It listed
    three instruments and four holdings; `_INSTRUMENT_HOLDINGS` had FOUR and EIGHT,
    and the missing one was `iag_fts_solar_atlas` / `solar_iag`. Every committed IAG
    band product was therefore unparseable to this test — i.e. invisible to the page —
    and the test reported that as the PRODUCT's fault. Same shape as
    `gf_rung.LAB_GRADED_SPECIES`: a second copy of a registry that silently went stale.

    Still not an import: importing `measure_band_ew` resolves the Kitt Peak atlas and
    exits when it is absent, and this is a test of the PARSER. Reading the literals out
    of the source keeps that property while making drift impossible.
    """
    src = (ROOT / "scripts" / "measure_band_ew.py").read_text()
    blk = src[src.index("_INSTRUMENT_HOLDINGS: dict"):]
    blk = blk[:blk.index("\n}")]
    return (set(re.findall(r'^    "([a-z_0-9]+)": \(', blk, re.M)),
            set(re.findall(r'HoldingSpec\("([a-z_0-9]+)"', blk)))


INSTRUMENTS, HOLDINGS = _registries()


# ── the parser keeps the selector rather than dropping the product ────────────────────
@pytest.mark.parametrize("name,selector", [
    ("FeI_4200_6910_kpno_solar_atlas_SYNTH_products.csv", "default"),
    ("FeI_4200_6910_kpno_solar_atlas_SYNTH_DEEPGRADED_products.csv", "DEEPGRADED"),
    ("FeI_4200_6910_harps_solar_harps_molecfit_corrected_SYNTH_DEEPGRADED_products.csv",
     "DEEPGRADED"),
    ("FeI_3800_6910_kpno_solar_atlas_PROFILEFIT_FROMEW_products.csv", "FROMEW"),
    ("FeI_3800_6910_kpno_solar_atlas_PROFILEFIT_FROMEW-GRADED_products.csv",
     "FROMEW-GRADED"),
    ("FeI_3800_6910_kpno_solar_atlas_PROFILEFIT_FROMEW-UNGRADED_products.csv",
     "FROMEW-UNGRADED"),
])
def test_a_selector_tagged_product_is_parsed_and_keeps_its_tag(name, selector):
    meta = _tracker().parse_stem(name, INSTRUMENTS, HOLDINGS)
    assert meta is not None, f"tracker cannot see {name} -- it would be silently dropped"
    assert meta["selector"] == selector


def test_the_handler_is_not_swallowed_by_the_selector():
    """`SYNTH` must stay the handler when a tag follows it, not be absorbed into it."""
    meta = _tracker().parse_stem(
        "FeI_4200_6910_kpno_solar_atlas_SYNTH_DEEPGRADED_products.csv",
        INSTRUMENTS, HOLDINGS)
    assert meta["handler"] == "SYNTH"
    assert meta["instrument"] == "kpno_solar_atlas"
    # The Kitt Peak stems name the INSTRUMENT and no holding, so the parser must say the
    # holding is absent rather than invent one -- an instrument can serve a corrected and
    # an uncorrected holding, and guessing is the collapse RYA-933/934 prevented.
    assert meta["holding"] is None

    harps = _tracker().parse_stem(
        "FeI_4200_6910_harps_solar_harps_molecfit_corrected_SYNTH_DEEPGRADED_products.csv",
        INSTRUMENTS, HOLDINGS)
    assert harps["handler"] == "SYNTH"
    assert harps["instrument"] == "harps"
    assert harps["holding"] == "solar_harps_molecfit_corrected"
    assert harps["holding_source"] == "filename"


# ── the invariant: nothing committed goes unseen ──────────────────────────────────────
def test_every_committed_band_product_is_visible_to_the_tracker():
    """🔴 A tracker that cannot see a merged product under-reports without failing.

    Scoped to `data/results/band_products/`, which is where the band-product convention
    (RYA-933/934) applies. Other result trees hold different artifact families whose
    names were never claimed by this pattern.
    """
    if not BAND_PRODUCTS.is_dir():
        pytest.skip("no band_products/ directory in this checkout")
    tracker = _tracker()
    unseen = [p.name for p in sorted(BAND_PRODUCTS.glob("*_products.csv"))
              # A name carrying a HANDLER token is claiming the band-product convention;
              # pre-RYA-933/934 artifacts (no band, no handler) never did and are out of
              # scope for this pattern.
              if ("_PROFILEFIT" in p.name or "_SYNTH" in p.name)
              and tracker.parse_stem(p.name, INSTRUMENTS, HOLDINGS) is None]
    assert unseen == [], (
        "committed band products the tracker cannot parse, so they never reach the "
        f"page: {unseen}")


def test_two_selectors_on_one_holding_stay_two_products():
    """RYA-946 firewall: differing line sets must not collapse into one cell."""
    tracker = _tracker()
    shallow = tracker.parse_stem(
        "FeI_4200_6910_kpno_solar_atlas_SYNTH_products.csv", INSTRUMENTS, HOLDINGS)
    deep = tracker.parse_stem(
        "FeI_4200_6910_kpno_solar_atlas_SYNTH_DEEPGRADED_products.csv",
        INSTRUMENTS, HOLDINGS)
    assert shallow["selector"] != deep["selector"]
    # Everything else about them is identical, so the selector is the ONLY thing keeping
    # them apart -- which is why it has to be carried through to the emitted row.
    assert {k: v for k, v in shallow.items() if k != "selector"} == \
           {k: v for k, v in deep.items() if k != "selector"}


# ── the same invariant, one layer in: the FEED path (RYA-1097) ────────────────────────
# 🔴 EVERY TEST ABOVE GUARDS A CODE PATH THE LIVE PAGE NO LONGER USES. RYA-1097 made the
# element feed the source and left `parse_stem` behind `--from-results`, so the invariant
# "nothing merged goes unseen" was still pinned -- against the scanner. It held. All 70
# live products in `data/products/solar/Fe.json` v1.94 are read and either rendered or
# explicitly withheld with a reason.
#
# What went unseen instead was the FIELDS. `_feed_row` projected the record through a
# hand-written allowlist written when the feed was smaller, and by v1.94 the feed
# published `grade`, the folded-xi block, a per-product `telluric_state` and
# `sigma_syst_complete` on all 70 -- none of which the allowlist named, so none of which
# reached the page. The page drew its systematic wireframe from the PUBLISHED `sigma_syst`
# while the feed sat there carrying `sigma_syst_complete`, larger on 51 of 70.
#
# So the invariant is pinned where it can actually fail now: a field the feed publishes on
# every live product must survive the projection.

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"

#: Dropped from the emitted row ON PURPOSE, each with the reason it is not a loss.
#: Named individually rather than pattern-matched, and every exemption is ASSERTED below
#: rather than merely declared -- an exemption nobody checks is a hole with a comment.
_FLATTENED = {
    "provenance": "flattened into measured_at / ingested_at / source",
}


def _live_records() -> list[dict]:
    import json
    if not FEED.exists():
        pytest.skip("no solar Fe feed in this checkout")
    return json.loads(FEED.read_text()).get("products") or []


def _universal_fields(records: list[dict]) -> set[str]:
    """Fields the feed publishes on EVERY live product.

    Not "any field on any product": `stat_basis` sits on 12 of 70 and `line_set` on 4,
    and a projection that omits an optional field is a different (smaller) question than
    one that omits a field the publisher puts on everything.
    """
    return set.intersection(*(set(r) for r in records)) if records else set()


def test_the_feed_projection_drops_no_universally_published_field():
    records = _live_records()
    assert records, "the feed publishes no live products -- the invariant is vacuous"
    tracker = _tracker()
    universal = _universal_fields(records)
    # A guard whose two sides converge proves nothing: if the feed ever publishes a bare
    # handful of fields, this test passes by describing nothing.
    assert len(universal) >= 20, (
        f"only {len(universal)} fields are universal; this invariant has gone vacuous")

    row = tracker._feed_row(records[0], star="solar", feed=FEED, pool="products",
                            instruments=set(), holdings=set())
    dropped = sorted(f for f in universal if f not in row and f not in _FLATTENED)
    assert dropped == [], (
        "fields the feed publishes on every live product that never reach the page: "
        f"{dropped}")


def test_every_named_flattening_exemption_still_applies():
    """An exemption must keep earning itself. `provenance` is exempt because it is
    UNPACKED, not discarded -- so assert the unpacked fields are actually there."""
    records = _live_records()
    tracker = _tracker()
    rec = next(r for r in records if (r.get("provenance") or {}).get("copied_to"))
    row = tracker._feed_row(rec, star="solar", feed=FEED, pool="products",
                            instruments=set(), holdings=set())
    assert set(_FLATTENED) == {"provenance"}, "a new exemption was added without a check"
    assert row["source"] == rec["provenance"]["copied_to"]
    assert row["measured_at"] == rec["provenance"].get("artifact_mtime")
    assert row["ingested_at"] == rec["provenance"].get("ingested_at")


def test_the_dropped_field_check_has_teeth():
    """The control. A projection that DOES drop a universal field must be caught --
    otherwise the test above is a green light with no mechanism (RYA-1080)."""
    records = _live_records()
    universal = _universal_fields(records)
    assert "grade" in universal
    crippled = {k: v for k, v in records[0].items() if k != "grade"}
    dropped = sorted(f for f in universal
                     if f not in crippled and f not in _FLATTENED)
    assert dropped == ["grade"], "the check cannot see a dropped field"


@pytest.mark.parametrize("field", ["grade", "xi_state", "xi_value_kms", "delta_xi_kms",
                                   "sigma_xi", "telluric_state", "sigma_syst_complete",
                                   "sigma_reported"])
def test_the_specific_regressions_reach_the_row(field):
    """The four facts the page was asked for and could not show: the published grade, the
    folded-xi state, the per-holding telluric state, and the COMPLETE systematic."""
    records = _live_records()
    tracker = _tracker()
    for rec in records:
        row = tracker._feed_row(rec, star="solar", feed=FEED, pool="products",
                                instruments=set(), holdings=set())
        # PRESENCE, not non-nullness. `sigma_xi` is legitimately null on the 11
        # UNMEASURED products -- xi was not measured, so it carries no sigma -- and
        # demanding a value there would be asking the projection to invent one. What
        # the projection owes is that the record's answer, null included, arrives
        # unchanged.
        assert field in row, (
            f"{field} is published on this product and absent from the emitted row")
        assert row[field] == rec[field], f"{field} was altered in transit"


def test_stat_basis_is_resolved_through_the_eligibility_module():
    """🔴 Read raw, `stat_basis` is non-null on 12 of 70 and the page renders 58 blanks.
    `pe.stat_basis_of` resolves all 70 through its route fallback. A blank basis reads as
    'unknown', which is absence-as-conclusion, and a second hand-rolled copy of a rule
    that already has a module is how a ratified rule quietly stops being enforced."""
    from pipeline import product_eligibility as pe
    records = _live_records()
    tracker = _tracker()
    raw = sum(1 for r in records if r.get("stat_basis"))
    rows = [tracker._feed_row(r, star="solar", feed=FEED, pool="products",
                              instruments=set(), holdings=set()) for r in records]
    resolved = sum(1 for r in rows if r["stat_basis"])
    assert resolved == len(records), (
        f"{len(records) - resolved} products render with a blank sigma_stat basis")
    assert resolved > raw, (
        "the module resolves no more than the raw field -- this test has gone vacuous")
    for rec, row in zip(records, rows):
        assert row["stat_basis"] == pe.stat_basis_of(rec)
        assert row["stat_basis_declared"] == rec.get("stat_basis")


def test_the_telluric_basis_quotes_the_products_own_stated_state():
    """RYA-1194's shape: the per-HOLDING registry was answering a question the per-PRODUCT
    record answers. `telluric_state_of` returned the constant 'named holding' over the top
    of five distinct stated states."""
    records = _live_records()
    tracker = _tracker()
    for rec in records:
        row = tracker._feed_row(rec, star="solar", feed=FEED, pool="products",
                                instruments=set(), holdings=set())
        state = tracker.telluric_state_of(row, row.get("committed"))
        assert state["telluric_basis"] == rec["telluric_state"]
        assert state["telluric_basis"] != "named holding"


def test_the_complete_systematic_actually_differs_from_the_published_one():
    """The reason the previous point matters, measured rather than asserted. If these two
    ever agree everywhere, the finding is spent and this test says so instead of passing
    quietly on a distinction that stopped existing."""
    records = _live_records()
    differ = [r for r in records
              if abs(float(r["sigma_syst_complete"]) - float(r["sigma_syst"])) > 5e-4]
    assert differ, ("sigma_syst_complete no longer differs from sigma_syst anywhere -- "
                    "re-check whether the page still needs to distinguish them")
    worst = max(differ, key=lambda r: float(r["sigma_syst_complete"]) - float(r["sigma_syst"]))
    assert float(worst["sigma_syst_complete"]) > float(worst["sigma_syst"]), (
        "the complete budget is SMALLER than the published one -- that inverts the finding")


# ── the display gate must not draw a SUBSET of the published grades ───────────────────
# 🔴 THIS GATE HAS NARROWED SILENTLY TWICE. It drew `selector.startswith("GRADED")`,
# which withheld every DEEPGRADED product and rendered the entire near-UV band blank --
# all 8 of its live products are deep, because the near-UV lab pool is deep almost to a
# line. Corrected to the two measured families, it still held back the reference-set
# reproduction. Both times the page showed a subset and said nothing a reader could use.
#
# Ryan's ruling (2026-08-24, carried verbatim in the feed's own quarantine record): "the
# product grid is GRADED / DEEPGRADED / CONSISTENT and nothing else". Three grades ship
# as Reference / Codex / Deep, and every one of them is a product. What belongs in the
# appendix is `tier=ALL`, which the eligibility gate already withdraws upstream.

def test_every_published_grade_is_drawn():
    import json
    out = ROOT / "data" / "results" / "rya935" / "live_status.json"
    if not out.exists():
        pytest.skip("tracker has not been generated in this checkout")
    status = json.loads(out.read_text())
    drawn = {str(r.get("grade")) for r in status["products"]}
    live = {str(r.get("grade")) for r in _live_records()}
    assert live - drawn == set(), (
        f"grades published in the feed that the page never draws: {sorted(live - drawn)}")
    # Not vacuous: there really are several distinct grades to lose.
    assert len(live) >= 3, f"only {len(live)} grades in the feed; this check has no teeth"


def test_no_live_product_is_withheld_by_the_display_policy():
    """The display policy's remaining job is the telluric gate, not a tier filter. If a
    live product IS withheld, the reason must be the holding's telluric state -- never a
    selector, which is what silently removed two thirds of the products."""
    import json
    out = ROOT / "data" / "results" / "rya935" / "live_status.json"
    if not out.exists():
        pytest.skip("tracker has not been generated in this checkout")
    status = json.loads(out.read_text())
    held = [r for r in status["products_withheld"] if r.get("pool") == "products"]
    by_selector = [r for r in held if "selector" in str(r.get("not_displayed_because"))]
    assert by_selector == [], (
        "live products withheld on a SELECTOR: "
        f"{[(r.get('band'), r.get('selector')) for r in by_selector]}")


def test_every_band_with_a_live_product_reaches_the_plot():
    """🔴 The failure this whole section exists for, stated as the invariant. near-UV had
    8 live products and drew none of them, and the page rendered the band as absent."""
    import json, collections
    out = ROOT / "data" / "results" / "rya935" / "live_status.json"
    if not out.exists():
        pytest.skip("tracker has not been generated in this checkout")
    status = json.loads(out.read_text())
    live = collections.Counter(r["band"] for r in _live_records())
    drawn = collections.Counter(r["band"] for r in status["products"])
    missing = {b: n for b, n in live.items() if not drawn.get(b)}
    assert missing == {}, (
        f"bands with live products and nothing on the plot: {missing}")
