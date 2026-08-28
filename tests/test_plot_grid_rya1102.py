"""RYA-1102 — every section shows the SAME rows, and a missing product says N/A.

🔴 THE DEFECT. `solar-report.js::plotRows` emitted one row per product from a build-time
artifact. Three consequences, all visible on the live page:

  * a section's row count was whatever the data contained -- Kitt Peak rendered 19;
  * a model with no product for that instrument simply VANISHED, indistinguishable from
    a model that does not exist;
  * the artifact was static, so re-running a product could not change the plot. Its Kitt
    Peak rows were pre-continuum-fix values repeated three times each.

Ryan: "each instrument section on the error plot should have 8 lines total, no more no
less. If the data product is not available for that instrument at that band we simply put
N/A... when we redo the products, they should simply update live."
"""
import ast
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data" / "products" / "solar" / "Fe.json"

from pipeline import plot_grid as pg                    # noqa: E402
from pipeline import product_eligibility as pe          # noqa: E402
from pipeline import treatment_axes as tx               # noqa: E402


@pytest.fixture(scope="module")
def feed():
    return json.loads(FEED.read_text())


@pytest.fixture(scope="module")
def grid(feed):
    assert "plot_grid" in feed, "the feed publishes no plot_grid — the site has nothing to render"
    return feed["plot_grid"]


# ── the axis ────────────────────────────────────────────────────────────────

def test_the_axis_is_DERIVED_not_typed():
    """Names come from `display_for`, so the axis cannot drift from the labels products
    actually carry. A typed axis is the RYA-1100 defect one level up."""
    for row in tx.plot_row_axis():
        assert row["name"] == tx.display_for(row["treatment"], route_token=row["route"])


def test_the_axis_is_the_ratified_ladder(grid):
    """Ryan's order: 1D-LTE, then 1D-NLTE, then <3D>, then full 3D; EW before Synth."""
    assert grid["axis"] == [
        "EW · 1D-LTE",
        "Synth · 1D-LTE",
        "Synth · 1D-LTE · Gerber",
        "EW · 1D-NLTE · Bergemann",
        "Synth · 1D-NLTE · Bergemann",
        "Synth · 1D-NLTE · Gerber",
        "Synth · <3D>-LTE · Gerber · stagger",
        "Synth · <3D>-NLTE · Gerber · stagger",
        "EW · 3D-NLTE · Amarsi",
    ]


def test_the_codex_deck_is_DECLARED_but_not_rendered(grid):
    """The 10th model is the deck we solve ourselves. It is on the axis so the ladder is
    complete and turning it on is one flag -- but a row for a model that does not exist
    yet would read as a measurement we failed to make."""
    pending = [r for r in tx.plot_row_axis(include_pending=True) if r["pending"]]
    assert [r["treatment"] for r in pending] == ["synth-mean3D-NLTE-gerber-codex"]
    assert len(grid["axis"]) == 9
    assert "codex" not in " ".join(grid["axis"])


# ── the grid: same rows everywhere, N/A where absent ────────────────────────

def test_EVERY_section_has_exactly_the_axis_rows_in_order(grid):
    """The headline requirement. No more, no less, same order, every section."""
    for s in grid["sections"]:
        got = [c["row"] for c in s["cells"]]
        assert got == grid["axis"], (
            f"{s['instrument']} {s['band']} {s['holding']} rendered {len(got)} rows: {got}")


def test_an_absent_product_is_an_N_A_CELL_not_a_missing_row(grid):
    """The distinction the old plot could not make."""
    empty = [c for s in grid["sections"] for c in s["cells"] if c["product_key"] is None]
    assert empty, "no N/A cells at all — the grid is not exercising the empty case"
    for c in empty:
        assert c["row"] in grid["axis"]


def test_every_product_key_RESOLVES_to_a_live_product(feed, grid):
    """A cell pointing at nothing renders blank and looks like an N/A, which would hide a
    join bug as a coverage gap."""
    keys = {pe.key_of(p) for p in feed["products"]}
    for s in grid["sections"]:
        for c in s["cells"]:
            if c["product_key"]:
                assert c["product_key"] in keys, f"dangling key {c['product_key']}"


def test_the_grid_stores_a_KEY_not_a_COPY_of_the_product(grid):
    """Duplicating the numbers into the grid would create a second place they can be
    wrong, and they would drift on the first re-publish (RYA-353)."""
    for s in grid["sections"]:
        for c in s["cells"]:
            for banned in ("A", "sigma_stat", "sigma_syst", "n_lines"):
                assert banned not in c, f"the grid copied {banned} instead of referencing it"


# ── the display rule ────────────────────────────────────────────────────────

def test_GRADED_only_except_where_only_DEEPGRADED_exists(feed, grid):
    """Ryan: DEEPGRADED is a secondary product documented in its own section, not
    showcased -- unless a band has nothing else, where "no graded product here" and "no
    product at all" are different facts."""
    by_key = {pe.key_of(p): p for p in feed["products"]}
    for s in grid["sections"]:
        for c in s["cells"]:
            if not c["product_key"]:
                continue
            p = by_key[c["product_key"]]
            assert p["tier"] == "GRADED" or s["only_deepgraded"], (
                f"{p['tier']} product rendered in a section that has GRADED products")


def test_the_near_UV_exception_actually_FIRES(grid):
    """A rule with no live case is untested. near-UV is DEEPGRADED-only, so it must be
    rendering -- if this ever goes empty the exception has silently stopped applying."""
    nuv = [s for s in grid["sections"] if s["band"] == "near-UV"]
    assert nuv, "no near-UV section at all"
    assert all(s["only_deepgraded"] for s in nuv)
    assert any(c["product_key"] for s in nuv for c in s["cells"])


def test_EW_and_SYNTH_are_SEPARATE_ROWS_never_collapsed(grid):
    """A profile fit and a flux fit measure different line pools. Collapsing them would
    make one of two real measurements disappear."""
    assert "EW · 1D-LTE" in grid["axis"] and "Synth · 1D-LTE" in grid["axis"]
    assert "EW · 1D-NLTE · Bergemann" in grid["axis"]
    assert "Synth · 1D-NLTE · Bergemann" in grid["axis"]


# ── slot ties, and the promise that nothing vanishes ────────────────────────

def test_a_deprecated_alias_NEVER_wins_a_slot_from_its_canonical_label(feed, grid):
    """🔴 MEASURED DEFECT. Under RYA-1100's derived names the two ENGINE-B rows render
    "EW · 1D-LTE" and were BEATING the genuine product -- harps showed 7.431/n=5 instead
    of 7.498/n=6, kpno 7.479/n=9 instead of 7.445/n=13."""
    by_key = {pe.key_of(p): p for p in feed["products"]}
    for s in grid["sections"]:
        for c in s["cells"]:
            if c["product_key"]:
                assert by_key[c["product_key"]]["treatment"] not in tx.DEPRECATED_ALIASES


def test_the_two_known_alias_rows_are_the_ONLY_slot_losers(feed, grid):
    """Pinned exactly. A new slot loser is a real product silently leaving the plot."""
    placed = {c["product_key"] for s in grid["sections"] for c in s["cells"] if c["product_key"]}
    only = pg._sections_with_only_deepgraded(feed["products"])
    lost = [p for p in feed["products"]
            if pg.is_displayable(p, only) and pe.key_of(p) not in placed]
    assert {p["treatment"] for p in lost} == {"ENGINE-B"}, \
        f"unexpected slot losers: {[(p['treatment'], p['display']) for p in lost]}"
    assert len(lost) == 2


def test_EVERY_live_product_is_ACCOUNTED_FOR(feed, grid):
    """🔴 RYA-711 applied to a renderer: if it is not shown, it must still be counted.
    Every live product is placed, or lost a slot tie, or was excluded by the graded-only
    rule -- and the three must sum to the whole pool with nothing left over."""
    P = feed["products"]
    placed = {c["product_key"] for s in grid["sections"] for c in s["cells"] if c["product_key"]}
    only = pg._sections_with_only_deepgraded(P)
    shown = [p for p in P if pg.is_displayable(p, only)]
    lost = [p for p in shown if pe.key_of(p) not in placed]
    excluded = [p for p in P if not pg.is_displayable(p, only)]
    assert len(placed) + len(lost) + len(excluded) == len(P)


def test_an_OFF_AXIS_product_type_is_reported_not_dropped(grid):
    """A product whose display name is not on the axis must appear in `off_axis`, never
    just disappear. Today the axis covers everything displayable, so this asserts the
    channel exists and is empty rather than that it is unused."""
    for s in grid["sections"]:
        assert "off_axis" in s
    assert not [n for s in grid["sections"] for n in s["off_axis"]]


# ── the source ──────────────────────────────────────────────────────────────

def test_the_grid_is_built_from_the_FEED_not_a_static_report():
    """The 'updates live' requirement. `build` takes the product list; nothing in this
    module reads a generated report."""
    src = (ROOT / "pipeline" / "plot_grid.py").read_text()
    tree = ast.parse(src)
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    assert "SOLAR_REPORT" not in names and "feEvidence" not in src


# ── the grid cannot go stale, because it is rebuilt on every write ──────────

def test_every_feed_write_goes_through_write_feed():
    """🔴 THE STALENESS FIX. The plot's old source was a build artifact nobody
    regenerated. If a write site bypasses `write_feed`, the grid it leaves behind is the
    same defect wearing a new filename, so the AST is checked rather than trusting that
    three call sites stay in step."""
    src = (ROOT / "scripts" / "publish_product.py").read_text()
    tree = ast.parse(src)
    bad = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)):
            continue
        if node.func.attr != "write_text":
            continue
        fn = next((f for f in ast.walk(tree)
                   if isinstance(f, ast.FunctionDef) and node in ast.walk(f)), None)
        if fn is None or fn.name != "write_feed":
            bad.append(getattr(fn, "name", "<module>"))
    assert not bad, f"feed written outside write_feed() in: {bad}"


def test_write_feed_REBUILDS_the_grid_rather_than_carrying_it_forward(tmp_path):
    """A behavioural check, not a structural one: change the products and the grid must
    follow without anyone asking it to."""
    import json as _json
    from scripts.publish_product import write_feed
    doc = _json.loads(FEED.read_text())
    before = sum(1 for s in doc["plot_grid"]["sections"] for c in s["cells"] if c["product_key"])
    doc["products"] = [p for p in doc["products"] if p["band"] != "NIR"]
    out = tmp_path / "Fe.json"
    write_feed(out, doc)
    after = _json.loads(out.read_text())["plot_grid"]
    filled = sum(1 for s in after["sections"] for c in s["cells"] if c["product_key"])
    assert filled < before, "the grid did not follow the products — it was carried forward"


def test_the_feed_is_written_ASCII_ESCAPED_like_the_rest_of_the_repo():
    """⚠️ Every display name contains '·'. Writing it literally reformats every line that
    carries one — a whole-file diff that hides the real change and that the next publish
    reverts. Measured against the committed file, not against intent."""
    raw = FEED.read_bytes()
    assert b"\\u00b7" in raw, "the feed was written with ensure_ascii=False"
    assert "·".encode() not in raw, "a literal '·' reached the committed feed"
