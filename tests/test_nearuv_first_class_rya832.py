"""
tests/test_nearuv_first_class_rya832.py — RYA-832
=================================================
The near-UV Fe product is derived by SYNTHESIS, because equivalent-width measurement is
undefined in that band — `pipeline.band_policy` forbids profile-fit and
interval-integration in 3000-3800 A, and RYA-759 falsified profile fitting there in
practice (901 candidates, 0 measurable). Wiring it into the EW-keyed driver therefore
needed a second route, and this file pins the three things that route can get wrong
without anyone noticing:

  * it takes the wrong branch — the near-UV silently falls into the EW path, which then
    fails on a missing EW file that is missing BY DESIGN;
  * it emits the wrong SCHEMA — the cell reaches the matrix carrying `n=40` and a NaN
    abundance, which reads as a measurement that failed rather than a wiring mismatch.
    That happened on the first run of this route;
  * it drifts from RYA-759's configuration, so the wired value stops being the published
    one. The route defends against this by CALLING 759's functions rather than copying
    them, and the constants below are asserted to still match what 759 used.

The synthesis itself is not exercised here — it needs iSpec and Turbospectrum, which live
on Sirius. What is exercised is the routing and the output contract, which is where the
wiring bugs live.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.band_policy import resolve as resolve_band  # noqa: E402


# ── the routing decision ──────────────────────────────────────────────────────

def test_the_near_uv_band_forbids_the_ew_route_and_permits_synthesis():
    """The routing predicate the driver branches on. If the policy ever stopped
    forbidding profile-fit here, the driver would silently go back to the EW path."""
    pol = resolve_band(0.5 * (3000.0 + 3780.0))
    assert "profile-fit" in pol.forbidden_methods
    assert "interval-integration" in pol.forbidden_methods
    assert "synthesis" in pol.permitted_methods


def test_the_optical_and_ir_bands_still_take_the_ew_route():
    """The branch must not swallow the bands it was never meant to touch — this is the
    positive control for the routing predicate."""
    for mid in (0.5 * (3800.0 + 6910.0), 0.5 * (6910.0 + 9199.0)):
        pol = resolve_band(mid)
        assert "profile-fit" not in getattr(pol, "forbidden_methods", ())


# ── the output contract ───────────────────────────────────────────────────────

def test_the_products_schema_matches_what_the_matrix_reads():
    """THE BUG THIS FILE EXISTS FOR. `products_frame` emits `value`; the matrix reads
    `A`. The first run of this route wrote the wrong one and the near-UV cell landed in
    the matrix with n=40 and a NaN abundance — present, and empty."""
    import derive_band_products as dbp
    src = Path(dbp.__file__).read_text()
    route = src[src.index("def synthesis_route"):src.index("def asdict_line")]
    for col in ("A=", "stat_dex=", "syst_dex=", "n_lines=", "n_excluded=",
                "band=", "treatment="):
        assert col in route, f"the synthesis route does not emit {col!r}"
    assert "products_frame([product]).to_csv" not in route, (
        "the route writes products_frame's schema again — the matrix cannot read it")


def test_the_route_carries_the_759_pseudo_continuum_systematic():
    """RYA-759's 0.100 dex does NOT average down and is NOT in the line scatter, so it
    has to be in the budget rather than in prose."""
    import derive_band_products as dbp
    assert dbp.NEARUV_PSEUDO_CONTINUUM_DEX == 0.100


def test_the_route_still_uses_759s_own_configuration():
    """Drift guard on the numbers the ticket forbids moving. These are 759's defaults;
    if any changes, the wired value stops being the published 7.487 and the change must
    be deliberate and stated."""
    import derive_band_products as dbp
    assert dbp.NEARUV_HALF_WIDTH_A == 0.40
    assert dbp.NEARUV_MIN_SEP_A == 4.0
    assert dbp.NEARUV_N_LINES == 40


def test_the_route_reuses_759s_functions_rather_than_copying_them():
    """A hand-adapted copy keeps its source's identity in places nobody looks (RYA-701:
    one Ba->Al copy produced 13 defects) and would let the value drift silently. The
    route must IMPORT the validated pieces."""
    import derive_band_products as dbp
    src = Path(dbp.__file__).read_text()
    route = src[src.index("def synthesis_route"):src.index("def asdict_line")]
    assert "from rya759_nearuv_fe_product import" in route
    assert "select_lines" in route and "fit_one" in route
    assert "from pipeline.nearuv_synth import" in route


def test_the_route_loud_fails_rather_than_emitting_an_empty_product():
    """A missing line list must stop the run. Emitting a 0-line product would reach the
    matrix as a real cell reporting nothing, which is indistinguishable from a physics
    result. (Verified live: the first Sirius run refused exactly this way.)"""
    import derive_band_products as dbp
    src = Path(dbp.__file__).read_text()
    route = src[src.index("def synthesis_route"):src.index("def asdict_line")]
    assert "raise SystemExit" in route
    assert "Refusing to emit an" in route


def test_the_matrix_report_no_longer_claims_the_near_uv_is_underivable():
    """The report used to print 'NOT derivable by this route ... RYA-759, which is not
    merged'. Every clause is now false. A report describing a state the repo left behind
    is the defect class this project keeps paying for."""
    txt = (ROOT / "scripts" / "rya783_fe_matrix_report.py").read_text()
    # Check what the report PRINTS, not what the file contains. The removal comment
    # legitimately quotes the old wording so a reader knows what changed and why, and a
    # cruder assertion over the whole file flagged that comment — which would have
    # pressured me to delete the explanation to make a test pass.
    printed = "\n".join(l for l in txt.splitlines() if l.strip().startswith("print("))
    assert "which is not merged" not in printed
    assert "NOT derivable by this route" not in printed
    assert "is DERIVED and first-class" in printed


def test_the_near_uv_cell_is_never_coadded_with_another_band():
    """RYA-712: the atomic unit is per-band. `band_products` has no combine() and must
    not grow one."""
    import pipeline.band_products as bp
    assert not hasattr(bp, "combine")
    src = Path(bp.__file__).read_text()
    assert "def combine(" not in src
