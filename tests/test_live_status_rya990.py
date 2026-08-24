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
