"""RYA-1100 — a display name is DERIVED, and a label that names nothing must not survive.

🔴 THE DEFECT. `scripts/publish_product.py` carried a hand-typed `_DISPLAY` dict under a
comment that claimed the name was "derived, never stored by hand". It published three
real things to the live element page:

  1. `1D-LTE` and `ENGINE-B` both mapped to "Synth · 1D-LTE" -- read as a naming
     collision, actually one product under two labels (identical axes in LEGACY).
  2. `ENGINE-A` mapped to "EW · 1D-NLTE · Bergemann" unconditionally, so all 22 live
     ENGINE-A products on route=SYNTH displayed as EW measurements. RYA-1002 had already
     un-pinned that route in `_ROUTE_BY_LABEL`; the private map re-asserted it.
  3. The three axis-native treatments had NO entry, so they rendered as raw tokens
     ("synth-mean3D-NLTE-gerber-stagger") on the page.

All three are the same failure: a name typed beside the physics instead of read from it.
"""
import ast
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
PUBLISH = ROOT / "scripts" / "publish_product.py"

from pipeline import treatment_axes as tx                      # noqa: E402
from scripts.publish_product import display_name               # noqa: E402


@pytest.fixture(scope="module")
def live():
    return json.loads(FEED.read_text())


# ── the premise: ENGINE-B is not an engine ───────────────────────────────────

def test_ENGINE_B_declares_the_same_axes_as_1D_LTE():
    """The whole basis for retiring the label. If this ever fails, ENGINE-B has become a
    real treatment and RYA-1100's withdrawal must be REVERSED, not re-run."""
    assert tx.LEGACY["ENGINE-B"] == tx.LEGACY["1D-LTE"]


def test_the_alias_registry_names_1D_LTE_as_the_canonical_label():
    assert tx.DEPRECATED_ALIASES["ENGINE-B"] == ("1D-LTE", "synth")


def test_the_alias_still_RESOLVES_because_history_must_stay_readable():
    """Dual-label forever (Ryan). Rewriting the stored `treatment` column would make every
    historical product incomparable to its own past -- the RYA-874 lesson. The alias is
    deprecated for NEW products, never unreadable for old ones."""
    assert tx.display_for("ENGINE-B", route_token="SYNTH") == "Synth · 1D-LTE"


# ── the fix: the name is derived, not typed ──────────────────────────────────

def test_publish_product_has_NO_hand_typed_display_map():
    """An AST test, not a substring test: this file's own docstring names `_DISPLAY`, and
    a `grep` for the token goes red on the prose describing the bug it fixed."""
    tree = ast.parse(PUBLISH.read_text())
    assigned = {t.id for n in ast.walk(tree) if isinstance(n, ast.Assign)
                for t in n.targets if isinstance(t, ast.Name)}
    assert "_DISPLAY" not in assigned, (
        "publish_product re-grew a hand-typed display map; derive via "
        "treatment_axes.display_for instead (RYA-906/RYA-1100)")


def test_display_name_delegates_to_the_axis_registry():
    tree = ast.parse(PUBLISH.read_text())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "display_name")
    called = {ast.unparse(c.func) for c in ast.walk(fn) if isinstance(c, ast.Call)}
    assert any("display_for" in c for c in called), \
        f"display_name must call treatment_axes.display_for; it calls {called}"


@pytest.mark.parametrize("treatment,route,expected", [
    # 🔴 the 22-product defect: ENGINE-A on the synthesis route is NOT an EW measurement.
    ("ENGINE-A",                         "SYNTH",      "Synth · 1D-NLTE · Bergemann"),
    ("ENGINE-A",                         "PROFILEFIT", "EW · 1D-NLTE · Bergemann"),
    ("1D-LTE",                           "SYNTH",      "Synth · 1D-LTE"),
    ("1D-LTE",                           "PROFILEFIT", "EW · 1D-LTE"),
    ("ENGINE-B",                         "SYNTH",      "Synth · 1D-LTE"),
    ("ENGINE-A-3DNLTE",                  "EW-3D",      "EW · 3D-NLTE · Amarsi"),
    ("ENGINE-B-NLTE",                    "SYNTH",      "Synth · 1D-NLTE · Gerber"),
    # the axis-native trio had NO map entry and rendered as bare tokens on the page.
    ("synth-1D-LTE-gerber",              "SYNTH",      "Synth · 1D-LTE · Gerber"),
    ("synth-mean3D-LTE-gerber-stagger",  "SYNTH",      "Synth · <3D>-LTE · Gerber · stagger"),
    ("synth-mean3D-NLTE-gerber-stagger", "SYNTH",      "Synth · <3D>-NLTE · Gerber · stagger"),
])
def test_the_published_display_name_matches_the_derived_axes(treatment, route, expected):
    assert display_name(treatment, gf="kurucz", route=route) == expected


def test_no_published_product_renders_as_a_bare_treatment_token(live):
    """A label is not a name. If the derived name equals the raw token, the axes did not
    produce a name and the product is shipping under a string nobody owns."""
    for p in live["products"]:
        got = display_name(p["treatment"], gf=p.get("gf"), route=p["route"])
        assert got != p["treatment"], f"{p['treatment']} still renders as its own token"


def test_the_dominant_ERROR_TERM_is_never_read_as_a_gf_POOL():
    """⚠️ THE TRAP. `dominant` carries the string "gf scale (cited lab)" -- an error
    budget, not a pedigree -- and the old call site fell back to it when `gf` was empty,
    which would stamp `· lab-gf` on a Kurucz product. Measured as never having fired (all
    89 committed rows carry gf='kurucz'), so this pins the fix, not a past defect."""
    assert display_name("1D-LTE", gf="gf scale (cited lab)", route="SYNTH") \
        == "Synth · 1D-LTE"
    assert display_name("1D-LTE", gf="lab", route="SYNTH") == "Synth · 1D-LTE · lab-gf"


# ── the route-token witness ──────────────────────────────────────────────────

def test_route_token_is_a_WITNESS_and_an_unknown_one_does_not_guess():
    assert tx.resolve_route("ENGINE-A", route_token="SYNTH") == ("synth", "route_token")
    # ENGINE-A pins no route (RYA-1002), so with a token it does not recognise the answer
    # is "unknown" -- stated, never defaulted.
    assert tx.resolve_route("ENGINE-A", route_token="WAT") == (None, "unknown")


def test_the_handler_still_outranks_the_route_token():
    """RYA-869 put `handler` on the product to be the authoritative witness. A stored
    route field must not be able to overrule the class that did the work."""
    assert tx.resolve_route("1D-LTE", handler="SynthesisHandler",
                            route_token="PROFILEFIT") == ("synth", "handler")


# ── the feed: no duplicate cell survives, and the open case is NOT hidden ────

def test_no_live_cell_carries_both_ENGINE_B_and_an_identical_1D_LTE(live):
    """The user-visible symptom: two rows, same name, same numbers, one page."""
    from scripts.rya1100_retire_engine_b_duplicates import classify
    dupes, _, _ = classify(live)
    assert not dupes, (
        f"{len(dupes)} ENGINE-B products still duplicate a live 1D-LTE twin; "
        f"run scripts/rya1100_retire_engine_b_duplicates.py --apply")


def test_the_two_route_CONTRADICTIONS_are_still_live_and_still_flagged(live):
    """🔴 NOT SILENTLY WITHDRAWN. These are a different line selection, not a duplicate
    (n=5 vs 6, n=9 vs 13, both from the same CSV and sha256), and their route contradicts
    the label. One of each pair is wrong; which is a science call Ryan has not made.
    Withdrawing them to tidy the page would be picking the answer by convenience (RYA-161)."""
    from scripts.rya1100_retire_engine_b_duplicates import classify
    _, contra, _ = classify(live)
    assert len(contra) == 2, f"expected the 2 known contradictions, found {len(contra)}"
    for p, twin in contra:
        assert p["route"] == "PROFILEFIT"
        assert p["route"] != tx._ROUTE_BY_LABEL["ENGINE-B"].upper()
        assert (p["A"], p["n_lines"]) != (twin["A"], twin["n_lines"])


def test_every_withdrawal_is_REVERSIBLE_and_says_why(live):
    """RYA-711: quarantined, not culled. A withdrawal that does not name what replaced it
    cannot be undone by anyone but its author."""
    for rec in live.get("superseded", []):
        if rec.get("superseded_reason_code") != "DUPLICATE_RETIRED_LABEL":
            continue
        assert rec["treatment"] == "ENGINE-B"
        assert rec["superseded_by"]["treatment"] == "1D-LTE"
        assert rec["A"] is not None and rec["n_lines"] is not None, "record not kept in full"
        assert "RYA-711" in rec["superseded_reason"]


# ── RYA-1100 now OWNS the `display` field (handed over from RYA-1080's pin) ───

def test_every_published_display_is_DERIVABLE_from_its_own_axes(live):
    """🔴 THE CHECK THAT REPLACED RYA-1080's BASELINE PIN, and is stricter than it.

    RYA-1080 pinned `display` to a baseline commit, which made the field impossible to
    correct in code -- it blocked fixing 22 ENGINE-A products published as "EW" while
    running on route=SYNTH. Ownership moved here in exchange for an assertion a pin
    cannot make: the label must EQUAL what the axis registry derives for that row.

    A pin is satisfied by never touching the field. This cannot be satisfied by a hand
    edit at all -- the only way to change a display name is to change the physics it is
    derived from, or the deriver, and the deriver is tested above against the ratified
    table. Every product, every pool.
    """
    from scripts.publish_product import display_name
    for pool in ("products", "superseded", "quarantine", "archive"):
        for p in live.get(pool) or []:
            # RYA-1106: the product's own pool, not a hardcoded one -- see the note in
            # test_feed_repo_reconciliation_rya1080.
            want = display_name(p["treatment"], gf=p.get("gf"), route=p.get("route"))
            assert p["display"] == want, (
                f"{pool}: {p['treatment']} on route={p.get('route')} publishes "
                f"{p['display']!r} but its axes derive {want!r}")


def test_no_two_live_products_share_a_display_name_AND_a_cell(live):
    """The symptom the user actually reported: two rows, same name, same instrument.

    Scoped to a CELL rather than to the whole feed, because sharing a name across
    different bands or holdings is correct -- `Synth · 1D-LTE` on VIS and on NIR are the
    same model measured twice, and the page distinguishes them by their own columns.
    """
    import collections
    seen = collections.defaultdict(list)
    for p in live["products"]:
        seen[(p["ion"], p["band"], p["instrument"], p["holding"],
              p["tier"], p["selector"], p["display"])].append(p)
    dupes = {k: v for k, v in seen.items() if len(v) > 1}

    #: 🔴 THE TWO THAT REMAIN, AND WHY THEY ARE NOT FIXED HERE. Both are the ENGINE-B
    #: rows whose published route is wrong (see the module docstring). Correcting the
    #: route DOES NOT REMOVE THE COLLISION -- measured: moved to SYNTH they land beside
    #: the 1D-LTE SYNTH product of the SAME cell (harps 7.431/n=5 vs 7.512/n=67; kpno
    #: 7.479/n=9 vs 7.442/n=67), same model, same name, different line POOLS, and both
    #: carry selector="GRADED".
    #:
    #: So this is not a labelling defect at all: two real measurements over different
    #: pools have no identity field that distinguishes them. Inventing a selector to make
    #: the page tidy would be manufacturing provenance. Recorded as debt, tied to the
    #: ticket, and asserted EXACTLY -- this test fails if a new collision appears AND if
    #: these two are resolved, so neither can pass unnoticed.
    KNOWN_OPEN = {("I", "VIS", "harps", "solar_harps_molecfit_corrected",
                   "GRADED", "GRADED", "EW · 1D-LTE"),
                  ("I", "VIS", "kpno_solar_atlas", "solar_kpno_kurucz2005_corrected",
                   "GRADED", "GRADED", "EW · 1D-LTE")}
    new_dupes = {k: v for k, v in dupes.items() if k not in KNOWN_OPEN}
    assert not new_dupes, f"{len(new_dupes)} NEW cells show two rows with one name: " + \
        "; ".join(f"{k[6]!r} x{len(v)} in {k[2]} {k[1]} {k[0]}" for k, v in new_dupes.items())
    assert set(dupes) == KNOWN_OPEN, (
        f"the known-open set moved: expected {len(KNOWN_OPEN)} documented collisions, "
        f"found {sorted(dupes)}. If these were fixed, delete KNOWN_OPEN (RYA-1100).")


# ── the publisher reads the ROW's route, not the filename's ──────────────────

def test_the_route_comes_from_the_ROW_not_the_file_wide_flag():
    """🔴 A products CSV IS NOT SINGLE-ROUTE. The PROFILEFIT-named HARPS and KPNO GRADED
    files each carry two ProfileFitHandler rows AND one SynthesisHandler row. Publishing
    `--route PROFILEFIT` over all three stamped a synthesis product with an EW route --
    and `route` is in KEY_FIELDS, so that is corrupted IDENTITY, not a cosmetic label.

    Read against the committed artifact, not a hand-built frame."""
    import pandas as pd
    from scripts.publish_product import normalise
    src = ROOT / ("data/results/band_products/FeI_4200_6910_harps_"
                  "solar_harps_molecfit_corrected_PROFILEFIT_GRADED_products.csv")
    if not src.exists():
        pytest.skip(f"{src.name} not committed here")
    rows = {r["treatment"]: r for r in normalise(
        pd.read_csv(src), holding="h", tier="GRADED", route="PROFILEFIT", selector="GRADED")}
    assert rows["1D-LTE"]["route"] == "PROFILEFIT"
    assert rows["ENGINE-A"]["route"] == "PROFILEFIT"
    assert rows["ENGINE-B"]["route"] == "SYNTH", \
        "the SynthesisHandler row still inherits the file-wide PROFILEFIT flag"
    assert rows["ENGINE-B"]["display"] == "Synth · 1D-LTE"


def test_a_row_that_declares_NO_route_falls_back_to_the_flag():
    """The other half. A row predating the `handler` column must not be made to look like
    it said something -- the flag is the stated fallback, never a guess."""
    import pandas as pd
    from scripts.publish_product import normalise
    df = pd.DataFrame([{"element": "Fe", "ion": "I", "band": "VIS", "instrument": "x",
                        "treatment": "1D-LTE", "A": 7.5, "n_lines": 10, "n_excluded": 0,
                        "stat_dex": 0.01, "syst_dex": 0.02, "gf": "kurucz"}])
    got = normalise(df, holding="h", tier="GRADED", route="PROFILEFIT", selector="GRADED")
    assert got[0]["route"] == "PROFILEFIT"
