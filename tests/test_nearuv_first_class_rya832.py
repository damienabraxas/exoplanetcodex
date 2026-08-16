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


def test_the_route_does_not_re_add_the_pseudo_continuum_term():
    """🔴 THE TEST THAT USED TO LIVE HERE ASSERTED `NEARUV_PSEUDO_CONTINUUM_DEX == 0.100`
    AND SO PINNED A BUG (RYA-845).

    Asserting a magic constant locks the value. It does not check that the value is
    correct, and — the failure here — it does not check that the value is applied ONCE.
    The near-UV systematic was counted twice for as long as that test was green:
    `error_budget.build()` already adds the term for any band whose policy calls for a
    pseudo-continuum, and this route added it again in quadrature, publishing 0.2211
    where the budget's own answer is 0.1972.

    The replacement asserts the RELATIONSHIP instead: the budget owns the term, and the
    route must not touch `syst` after asking for it. That is a property a double-add
    cannot satisfy.
    """
    import derive_band_products as dbp
    src = Path(dbp.__file__).read_text()
    route = src[src.index("def synthesis_route"):src.index("def asdict_line")]
    assert "np.hypot(syst" not in route, (
        "the synthesis route is modifying `syst` after build_budget returned it — that "
        "is how the pseudo-continuum term came to be counted twice (RYA-845)")
    assert "PSEUDO_CONTINUUM_DEX" not in route, (
        "the route re-declares the pseudo-continuum systematic; it belongs to "
        "pipeline/error_budget.py alone, so that no caller can add what the budget "
        "already holds")


def test_the_budget_carries_the_pseudo_continuum_term_exactly_once():
    """The other half: the term must still be THERE. Removing the double-add must not
    quietly remove the term, which would be the opposite error and just as invisible."""
    from pipeline.error_budget import build

    b = build("Fe", 3390.0, 40, scatter_dex=0.413, gf_graded=False,
              harness_residual_dex=0.0, handler="SynthesisHandler")
    hits = [t for t in b.terms if "pseudo" in t.name.lower()]
    assert len(hits) == 1, f"expected exactly one pseudo-continuum term, found {len(hits)}"


def test_the_pseudo_continuum_term_is_near_uv_only():
    """RYA-841/845: the branch fires on the band POLICY, so it must not reach any band
    whose continuum is actually observed. If it ever did, every VIS and IR cell would
    silently inflate."""
    from pipeline.error_budget import build

    for wave_A in (5500.0, 8000.0, 15000.0):
        b = build("Fe", wave_A, 40, scatter_dex=0.10, gf_graded=False,
                  harness_residual_dex=0.0, handler="ProfileFitHandler")
        assert not [t for t in b.terms if "pseudo" in t.name.lower()], (
            f"a pseudo-continuum term appeared at {wave_A} A — that band observes its "
            f"continuum and must not carry this systematic")


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


def test_the_published_nearuv_cells_reconstruct_with_the_term_applied_once():
    """REGENERATE-AND-DIFF, as a test: the matrix's own near-UV systematics must equal
    what `error_budget.build()` produces — no more, no less.

    This is the check the constant-assertion could never be. `syst` does not depend on
    the line scatter at all (the scatter term averages down and lands in `stat`), so the
    budget's systematic is fully determined by the band and whether the pool is on graded
    gf. That makes the published number exactly reconstructible, and a re-added
    pseudo-continuum term shows up as a mismatch rather than as a green test.
    """
    import numpy as np
    import pandas as pd
    from pipeline.error_budget import build

    matrix = ROOT / "data" / "results" / "rya783" / "fe_product_matrix.csv"
    if not matrix.exists():
        pytest.skip("Fe product matrix absent")
    df = pd.read_csv(matrix)
    near = df[df.band == "near-UV"]
    assert len(near) >= 2, "expected both near-UV cells (1D-LTE and 1D-LTE-LABGF)"

    for _, row in near.iterrows():
        # The LABGF pool is the one on primary laboratory gf; the Kurucz pool is not.
        graded = row.treatment.endswith("LABGF")
        b = build("Fe", 3390.0, int(row.n_lines), scatter_dex=0.4,
                  gf_graded=graded, harness_residual_dex=0.0,
                  handler="SynthesisHandler")
        _, syst = b.total()

        assert round(float(syst), 4) == pytest.approx(float(row.syst_dex), abs=5e-4), (
            f"{row.treatment}: published syst {row.syst_dex} does not reconstruct from "
            f"the budget ({syst:.4f}). If the difference is a factor of the "
            f"pseudo-continuum term, something is adding it twice again (RYA-845)")

        # and state the failure mode explicitly, so the test names the bug it guards
        doubled = float(np.hypot(syst, 0.10))
        assert not np.isclose(float(row.syst_dex), doubled, atol=5e-4), (
            f"{row.treatment}: published syst equals the DOUBLE-ADDED value {doubled:.4f} "
            f"— the RYA-845 defect has returned")
