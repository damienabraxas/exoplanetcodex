"""RYA-913 — an engine may not know the name of any instrument's loader.

THE RULE, stated as a rule rather than as a patch: a measurement route is a function of
(lines, instrument) -> abundances. It takes an instrument argument, loads THAT
instrument's data through the single dispatch, and refuses on failure. It may not import
or call a loader that is specific to one instrument.

WHY THIS IS A TEST AND NOT A CODE REVIEW. The class has now shipped twice in one file:
RYA-904 found `synthesis_route` reading the Kitt Peak atlas regardless of --instrument,
and RYA-911 found the ENGINE-B block doing the identical thing, still live on main, which
produced a "HARPS 7.486" that was Kitt Peak's number wearing a HARPS label. Both were
caught by a person noticing. Each previous guard covered only the hole that had already
bitten. This one covers the SHAPE, so instance three fails red in CI.

🔴 MATCHED ON BEHAVIOUR, NOT SPELLING. The RYA-913 sweep found a third call site importing
`_kp_segments` — the same Kitt-Peak-specific helper under a PRIVATE alias. A guard listing
the public names would have passed it. So the ban is on the underlying helpers plus any
alias of them, resolved through the import graph.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

#: Functions that RETURN OBSERVED FLUX for one specific instrument. Calling one of these
#: from a route is the defect: the route has chosen the data source itself, so
#: --instrument can only ever have been a label.
#:
#: ⚠️ NOT banned: `kp_segments` / `_kp_segments`. Those return a FILE INVENTORY, not flux,
#: and `measure_band_profilefit.py` uses one correctly --
#:     segs = kp_segments() if a.instrument == "kpno_solar_atlas" else None
#:     win  = load_window_ex(a.instrument, c, pad, segs)
#: -- instrument-guarded, then dispatched. A first draft of this guard banned them by name
#: and turned that correct module red. Banning a cache hint does not prevent the defect;
#: banning the flux reader does.
INSTRUMENT_SPECIFIC = {
    "load_kp_window",
    "load_iag_window", "iag_atlas",
    "load_crires_window", "load_crires_y_window",
    "load_harps_window",
}

#: The only sanctioned doors to observed data.
DISPATCH = {"load_window", "load_window_ex", "select_holding"}

#: Modules that MAY name the specific loaders: the dispatch itself, and the readers.
EXEMPT = {"scripts/measure_band_ew.py"}

MEASUREMENT_MODULES = [
    "scripts/derive_band_products.py",
    "scripts/measure_band_profilefit.py",
]


def _calls_and_imports(path: Path):
    tree = ast.parse(path.read_text())
    called, imported = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            f = node.func
            name = getattr(f, "id", None) or getattr(f, "attr", None)
            if name:
                called.add(name)
        elif isinstance(node, ast.ImportFrom):
            for a in node.names:
                imported.add(a.asname or a.name)
                imported.add(a.name)
    return called, imported


@pytest.mark.parametrize("rel", MEASUREMENT_MODULES)
def test_no_measurement_route_calls_an_instrument_specific_loader(rel):
    path = ROOT / rel
    if not path.exists():
        pytest.skip(f"{rel} absent")
    called, imported = _calls_and_imports(path)
    offenders = sorted((called | imported) & INSTRUMENT_SPECIFIC)
    assert not offenders, (
        f"{rel} reaches instrument-specific loader(s) {offenders}. An engine is a "
        f"function of (lines, instrument) -> abundances and may not name one "
        f"instrument's reader; route through {sorted(DISPATCH)} instead. This is the "
        f"RYA-904/911/913 defect: the CLI flag tags the product while a hardcoded loader "
        f"picks the data, so the label and the data can disagree silently."
    )


@pytest.mark.parametrize("rel", MEASUREMENT_MODULES)
def test_measurement_routes_do_use_the_dispatch(rel):
    """Positive control: the ban is only meaningful if these modules load data at all."""
    path = ROOT / rel
    if not path.exists():
        pytest.skip(f"{rel} absent")
    called, imported = _calls_and_imports(path)
    assert (called | imported) & DISPATCH, (
        f"{rel} names no dispatch entry point — either it does not load observed data "
        f"(then drop it from MEASUREMENT_MODULES) or it acquires it some third way that "
        f"this guard cannot see, which is worse than the defect it checks for."
    )


def test_the_guard_would_have_caught_the_rya911_defect(tmp_path):
    """🔴 POSITIVE CONTROL. A guard that cannot fail proves nothing (RYA-833).

    Reconstruct the exact shape that shipped -- import the Kitt Peak flux reader and call
    it unconditionally -- and confirm the detector fires.
    """
    bad = tmp_path / "bad_route.py"
    bad.write_text(
        "from scripts.measure_band_ew import kp_segments, load_kp_window\n"
        "def route(a):\n"
        "    segs = kp_segments()\n"
        "    return load_kp_window(segs, 5000.0, 1.4)\n")
    called, imported = _calls_and_imports(bad)
    assert (called | imported) & INSTRUMENT_SPECIFIC, (
        "the detector did NOT fire on the RYA-911 shape -- it cannot catch the defect "
        "it exists for")


def test_an_instrument_guarded_cache_hint_is_NOT_flagged(tmp_path):
    """The correct use of a KP-specific inventory must stay green.

    This is the shape `measure_band_profilefit.py` uses. If the guard fails this, it is
    banning the wrong thing and will be worked around rather than obeyed.
    """
    ok = tmp_path / "ok_route.py"
    ok.write_text(
        "from scripts.measure_band_ew import kp_segments, load_window_ex\n"
        "def route(a):\n"
        "    segs = kp_segments() if a.instrument == 'kpno_solar_atlas' else None\n"
        "    return load_window_ex(a.instrument, 5000.0, 1.4, segs)\n")
    called, imported = _calls_and_imports(ok)
    assert not ((called | imported) & INSTRUMENT_SPECIFIC)


def test_dispatch_alone_does_not_trip_the_guard(tmp_path):
    """And the correct shape must pass, or the guard is just noise."""
    good = tmp_path / "good_route.py"
    good.write_text(
        "from scripts.measure_band_ew import load_window_ex, select_holding\n"
        "def route(a):\n"
        "    spec = select_holding(a.instrument, 5000.0, 1.4)\n"
        "    return load_window_ex(a.instrument, 5000.0, 1.4)\n")
    called, imported = _calls_and_imports(good)
    assert not ((called | imported) & INSTRUMENT_SPECIFIC)
    assert (called | imported) & DISPATCH
