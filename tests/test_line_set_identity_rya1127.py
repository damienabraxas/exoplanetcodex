"""
RYA-1127 — `line_set` is part of the product identity, and the migration did not
re-identify anything.

🔴 THE FAILURE MODE HERE IS NOT A WRONG NUMBER. No abundance moves in this ticket. Adding a
field to `KEY_FIELDS` re-identifies every product at once, and the way that goes wrong is a
wrong SPLIT or a wrong MERGE -- records that were one identity becoming two, or two
collapsing into one, with nobody deciding. So the tests below are about the PARTITION of
the feed, not about values.

The motivating case is RYA-1106: the Amarsi 3D-NLTE method measured on Asplund's own AGSS21
line set produced four products whose key was identical to the our-graded Amarsi products'
-- same everything except the pool of lines, which the key did not carry. Worse, the two
had been kept apart only by a MISLABEL (`route=EW-3D`, the stranded ProfileFitHandler's
route that RYA-1104 refuted), so correcting the label is what exposed the latent defect.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import product_eligibility as pe        # noqa: E402
from pipeline import reference_lineset as rls         # noqa: E402
from pipeline import model_registry as mr             # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
AUDIT = ROOT / "scripts" / "rya1127_key_migration_audit.py"


@pytest.fixture(scope="module")
def doc():
    return json.loads(FEED.read_text())


# ── the key itself ────────────────────────────────────────────────────────────────

def test_line_set_is_in_the_identity_key():
    assert "line_set" in pe.KEY_FIELDS


def test_the_two_KEY_FIELDS_definitions_still_agree():
    """`publish_product` keeps its own literal copy; they must not drift."""
    sys.path.insert(0, str(ROOT / "scripts"))
    import publish_product as pp
    assert tuple(pp.KEY_FIELDS) == tuple(pe.KEY_FIELDS)


def test_two_products_differing_ONLY_in_line_set_get_distinct_keys(doc):
    """🔴 THE RYA-1106 CASE, as a test. This is the whole ticket in four lines."""
    #: 🔴 RYA-1185 -- "OUR OWN" MEANS THE ONES THAT DO NOT STORE THE AXIS, NOT EVERY
    #: ENGINE-A-3DNLTE PRODUCT. This selected on treatment alone, which was unambiguous
    #: while the four RYA-1106 Asplund replications were unpublished. RYA-1185 published
    #: them, and they are ENGINE-A-3DNLTE too -- so the clone below became
    #: `dict(asplund_row, line_set="asplund")`, i.e. the row compared against ITSELF, and
    #: the test reported a collision that does not exist (the live feed has 70 products
    #: and 70 distinct keys). A replication is identified by STORING `line_set`; our own
    #: products derive it from `tier` and never store it (RYA-1111).
    amarsi = [p for p in doc["products"]
              if p.get("treatment") == "ENGINE-A-3DNLTE" and not p.get("line_set")]
    assert amarsi, "the our-graded Amarsi products should be live"
    for ours in amarsi:
        replication = dict(ours, line_set="asplund")
        assert pe.key_of(ours) != pe.key_of(replication), (
            f"{ours['holding']}: the Asplund replication still collides with the "
            f"our-graded leg")
        assert pe.key_of(ours).endswith("|our-graded")
        assert pe.key_of(replication).endswith("|asplund")


def test_the_collision_was_REAL_under_the_old_key(doc):
    """⚠️ CONTROL — the test above is only meaningful if the pair genuinely collided.

    A test that two keys differ proves nothing if they differed before as well. This
    recomputes the PREVIOUS key and asserts it was the same string for both, so the fix is
    shown to have fixed something.
    """
    old_fields = tuple(f for f in pe.KEY_FIELDS if f != "line_set")
    def old_key(p):
        return "|".join(str(p.get(k) or "") for k in old_fields)
    # same scoping as above (RYA-1185): our own products are the ones that store nothing.
    amarsi = [p for p in doc["products"]
              if p.get("treatment") == "ENGINE-A-3DNLTE" and not p.get("line_set")]
    for ours in amarsi:
        assert old_key(ours) == old_key(dict(ours, line_set="asplund"))


# ── resolution is loud, never defaulted (RYA-869) ─────────────────────────────────

def test_line_set_resolves_for_EVERY_record_in_EVERY_pool(doc):
    """Putting the field in the key means every pool must answer, not just `products[]`.

    Nine quarantined records could not answer when this ticket started -- six tier
    UNGRADED and three tier ALL -- which is why the vocabulary gained `our-ungraded` and
    `our-all` rather than `key_of` gaining a default.
    """
    for pool in ("products", "superseded", "archive", "quarantine"):
        for p in doc.get(pool) or []:
            assert rls.line_set_for_product(p) in mr.LINE_SETS, (pool, p.get("treatment"))


def test_an_unresolvable_tier_RAISES_rather_than_defaulting():
    """POSITIVE CONTROL. A resolver that never refuses is a default wearing a disguise."""
    with pytest.raises(rls.ReferenceLineSetError):
        rls.line_set_for_product({"tier": "NOT-A-TIER"})
    with pytest.raises(rls.ReferenceLineSetError):
        rls.line_set_for_product({})


def test_CONSISTENT_still_refuses_after_the_vocabulary_widened():
    """⚠️ The widening must be narrow. RYA-1105 retires the Consistent tier, so it must
    NOT have acquired a name when UNGRADED and ALL did."""
    assert "consistent" not in mr.LINE_SETS
    with pytest.raises(rls.ReferenceLineSetError):
        rls.line_set_for_product({"tier": "CONSISTENT"})


def test_key_of_RESOLVES_line_set_rather_than_reading_the_field(doc):
    """🔴 The bug this avoids: our own products do not STORE `line_set`.

    A `.get("line_set")` would have returned "" for all of them -- a key column that is
    blank everywhere, which looks like it works and still lets the RYA-1106 pair collide.
    """
    ours = next(p for p in doc["products"] if p.get("tier") == "GRADED")
    assert "line_set" not in ours, "this test is stale: the field is now stored"
    assert pe.key_of(ours).endswith("|our-graded")


# ── the migration ─────────────────────────────────────────────────────────────────

def test_the_migration_audit_reports_no_unintended_identity_change():
    """🔴 THE GATE. Runs the audit in --check mode over every committed feed."""
    r = subprocess.run([sys.executable, str(AUDIT), "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-2000:]


def test_the_migration_split_nothing_that_is_already_published(doc):
    """The partition of the LIVE feed is unchanged: same number of distinct identities.

    The RYA-1106 products are not in the feed yet -- they are what this unblocks -- so
    today the schema widens and nothing already published is re-identified. If a future
    edit makes a live product split, that is a decision and this test makes it visible.
    """
    old_fields = tuple(f for f in pe.KEY_FIELDS if f != "line_set")
    old = {"|".join(str(p.get(k) or "") for k in old_fields) for p in doc["products"]}
    new = {pe.key_of(p) for p in doc["products"]}
    assert len(old) == len(new), (
        f"{len(old)} identities became {len(new)} -- a live product was re-identified")


def test_the_published_plot_grid_keys_carry_the_new_axis(doc):
    """The feed publishes `key_fields` and every cell's `product_key`; both must have
    migrated, or the site joins cells to products on a key the products no longer have."""
    grid = doc["plot_grid"]
    assert list(grid["key_fields"]) == list(pe.KEY_FIELDS)
    live = {pe.key_of(p) for p in doc["products"]}
    for section in grid["sections"]:
        for cell in section["cells"]:
            if cell.get("product_key"):
                assert cell["product_key"] in live, cell["product_key"]


def test_the_four_RYA1106_products_PUBLISH_without_colliding(doc):
    """🔴 SPEC ITEM 5, demonstrated end to end on a COPY of the feed.

    Showing the keys differ is necessary but not sufficient: publishing also has to clear
    the duplicate-cell check and the eligibility gate, and the rendered grid has to keep
    resolving. This stages the four replication products beside the live ones IN MEMORY --
    nothing is written -- and asserts the feed is still well-formed.

    ⚠️ The staged records borrow the live records' shape and carry the RYA-1106 measured
    values, because a demo built from a hand-typed record would prove the schema works on
    a record that does not exist.
    """
    import pandas as pd
    from pipeline import plot_grid

    res = ROOT / "data" / "results" / "rya1106"
    hold_by_dir = {"kpno_kurucz2005": "solar_kpno_kurucz2005_corrected",
                   "kpno_molecfit": "solar_kpno_molecfit_corrected",
                   "harps_molecfit": "solar_harps_molecfit_corrected",
                   "iag": "solar_iag"}
    live = {p["holding"]: p for p in doc["products"]
            if p.get("treatment") == "ENGINE-A-3DNLTE"}

    staged = []
    for d, holding in hold_by_dir.items():
        csv = res / d / "asplund_lines_products.csv"
        if not csv.exists():
            pytest.skip(f"RYA-1106 artifact not present: {csv.relative_to(ROOT)}")
        r = pd.read_csv(csv).iloc[0]
        staged.append(dict(live[holding],
                           line_set="asplund",
                           A=round(float(r["value"]), 4),
                           sigma_stat=float(r["stat_dex"]),
                           sigma_syst=float(r["syst_dex"]),
                           n_lines=int(r["n_lines"]),
                           stat_basis=str(r["stat_basis"])))

    merged = dict(doc, products=list(doc["products"]) + staged)

    # 1. no identity is duplicated once line_set separates them
    assert pe.duplicate_live_cells(merged) == {}
    # 2. and the count of identities grew by exactly the four we added
    assert len({pe.key_of(p) for p in merged["products"]}) == \
           len({pe.key_of(p) for p in doc["products"]}) + 4
    # 3. every staged product clears the eligibility gate
    bad = {k: [f"{r.code}: {r.detail}" for r in v]
           for k, v in pe.evaluate_feed(merged).items() if v}
    assert not bad, bad
    # 4. and the grid still builds, with the our-graded cells still resolving
    grid = plot_grid.build(merged["products"])
    keys = {pe.key_of(p) for p in merged["products"]}
    for section in grid["sections"]:
        for cell in section["cells"]:
            if cell.get("product_key"):
                assert cell["product_key"] in keys


def test_the_audits_split_detection_is_NOT_vacuous():
    """⚠️ CONTROL — the audit reports zero splits on the real feed, so nothing has ever
    exercised the code that finds one. A detector that has never fired is not evidence.

    Staged: the four Asplund products beside the live ones, which is the RYA-1106 case.
    The audit must find exactly four splits and attribute each to line_set.
    """
    import importlib.util, json as _json, tempfile
    spec = importlib.util.spec_from_file_location("aud", AUDIT)
    aud = importlib.util.module_from_spec(spec); spec.loader.exec_module(aud)

    feed = _json.loads(FEED.read_text())
    amarsi = [p for p in feed["products"]
              if p.get("treatment") == "ENGINE-A-3DNLTE" and not p.get("line_set")]
    staged = dict(feed, products=list(feed["products"])
                  + [dict(p, line_set="asplund") for p in amarsi])
    d = Path(tempfile.mkdtemp()) / "solar"; d.mkdir(parents=True)
    (d / "Fe.json").write_text(_json.dumps(staged, indent=2))

    r = aud.audit_feed(d / "Fe.json")
    assert len(r["splits"]) == len(amarsi) == 4, r["splits"]
    for s in r["splits"]:
        assert s["line_sets"] == ["asplund", "our-graded"]
    assert r["merges"] == []
    assert r["unresolved"] == []


def test_the_refinement_invariant_is_what_the_audit_actually_asserts():
    """🔴 The check with teeth. While `line_set` is APPENDED, a merge is structurally
    impossible and every split is explained by line_set for free -- so those two are not
    evidence. The invariant that makes them true is asserted instead.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("aud", AUDIT)
    aud = importlib.util.module_from_spec(spec); spec.loader.exec_module(aud)
    r = aud.audit_feed(FEED)
    assert r["refinement_holds"], r["refinement_broken"]
    assert tuple(pe.KEY_FIELDS[:len(aud.OLD_KEY_FIELDS)]) == tuple(aud.OLD_KEY_FIELDS), (
        "line_set is no longer a pure suffix on KEY_FIELDS -- the audit's merge and "
        "split reasoning no longer follows, and both need to become real checks")
