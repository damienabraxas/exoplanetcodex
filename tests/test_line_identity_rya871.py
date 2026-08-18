"""
tests/test_line_identity_rya871.py — RYA-871
============================================
A WAVELENGTH DOES NOT IDENTIFY A LINE. The EW-route per-line artifact carried none of the
excitation potential its own source row had, so `gf_rung.resolve_lines` matched a measured
line back to the loaded line list on wavelength alone and 16 of 152 VIS Fe I lines did not
resolve at all — 14 with no row inside the 0.005 A window and 2 with two rows inside it.

WHY BOTH HALVES ARE THE SAME MISSING KEY
-----------------------------------------
`measure_band_profilefit` draws its candidates from `data/audit/line_accounting/
per_line.csv`, whose rows are FEATURES rather than lines: `line_accounting_rya709.
features()` groups line-list rows within 0.05 A and reports the group MEAN wavelength.
So a blended feature sits BETWEEN its components by construction, and the 0.006-0.02 A
offsets are that geometry, not measurement error. Reaching the component means widening
the window; widening the window without a second key means CHOOSING between the rows it
now contains, which is what RYA-855 refused. The EP is the second key.

WHAT THESE PIN
--------------
Relationships and refusals, not the counts — the counts are a property of today's pools
and are measured in `data/results/rya871/`, not asserted here.

1. the two tolerances travel WITH the key: a line carrying no EP is never widened;
2. an EP that disagrees REFUSES a match the wavelength alone would have accepted —
   the guard is not decoration;
3. a genuine (wavelength AND EP) degeneracy stays UNRESOLVED, never `iloc[0]`;
4. the parallel arrays cannot slip — a length mismatch raises rather than keying a line
   on its neighbour's EP;
5. the emitters RAISE on a source row with no EP instead of emitting a blank identity;
6. `ep_eV` reaches the artifact — a field that exists on the object and not in the CSV is
   the RYA-843 defect (`red_chi2` lived its whole life that way).
"""
from __future__ import annotations

import sys
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import gf_rung                                       # noqa: E402
from pipeline.band_products import LineMeasurement                 # noqa: E402


#: Two REAL Fe I transitions close in wavelength and far apart in EP — the shape RYA-855
#: met at 3125.651 / 3125.683 A (0.990 and 2.404 eV). 0.015 A apart, so BOTH fall inside
#: the EP-keyed window and NEITHER falls inside the wavelength-only one at a midpoint.
_PAIR = pd.DataFrame({
    "species": ["Fe 1", "Fe 1"],
    "wavelength_air_A": [5000.000, 5000.015],
    "ep_eV": [0.990, 2.404],
    "log_gf": [-1.100, -2.200],
})


def _resolve(waves, eps=None, ll=_PAIR):
    with mock.patch.object(gf_rung, "linelist_frame", lambda _: ll):
        return gf_rung.resolve_lines("Fe", "I", waves, None, measured_ep_eV=eps)


# ── 1. the tolerance travels with the key ────────────────────────────────────────────

def test_a_line_with_no_ep_is_not_widened():
    """THE REGRESSION GUARD. A wider window with no second key buys a CHOICE, not an
    identification — measured: at 0.020 A with no EP, 7 of the 136 VIS Fe I lines that
    already resolved change WHICH ROW they resolve to. So the wide tolerance must be
    unreachable without an EP."""
    r = _resolve([5000.0075])                       # midway between the pair
    assert not r.resolved.iloc[0]
    assert f"{gf_rung.LINELIST_MATCH_TOL_A}" in r.unresolved_why.iloc[0]
    assert "wavelength alone" in r.unresolved_why.iloc[0]


def test_the_same_line_resolves_once_it_carries_its_ep():
    """The ticket, in one assertion: the identity was recoverable all along."""
    r = _resolve([5000.0075], eps=[2.404])
    assert bool(r.resolved.iloc[0])
    assert r.ep_eV.iloc[0] == pytest.approx(2.404)
    assert r.log_gf.iloc[0] == pytest.approx(-2.200)   # the RIGHT row, not the nearer one


def test_the_ep_key_picks_the_further_row_when_that_is_the_right_one():
    """A wavelength-only rule would take the nearest row. The EP must OVERRULE proximity,
    or it is a tie-breaker rather than an identification."""
    near, far = 5000.000, 5000.015
    w = near + 0.002                                  # much closer to `near`
    assert abs(w - near) < abs(w - far)
    r = _resolve([w], eps=[2.404])                    # ...but the EP says `far`
    assert bool(r.resolved.iloc[0])
    assert r.log_gf.iloc[0] == pytest.approx(-2.200)


# ── 2. the guard refuses, it does not only permit ────────────────────────────────────

def test_a_disagreeing_ep_refuses_a_match_wavelength_alone_would_accept():
    """VERIFY THE TEST DISCRIMINATES (RYA-805/845). If a wrong EP still matched, the key
    would be decoration and every assertion above would pass on a no-op."""
    exact = 5000.000
    assert bool(_resolve([exact]).resolved.iloc[0])            # wavelength alone: yes
    r = _resolve([exact], eps=[7.777])                          # nothing sits at 7.777 eV
    assert not r.resolved.iloc[0]
    assert "absent" in r.unresolved_why.iloc[0]
    assert "wavelength+EP" in r.unresolved_why.iloc[0]


# ── 3. a real degeneracy is still refused ────────────────────────────────────────────

def test_two_rows_identical_in_wavelength_and_ep_stay_unresolved():
    """HFS components share a lower level and sit ~0.001 A apart. Nothing can separate
    them, so the honest answer is UNGRADEABLE — never `iloc[0]`, which is how RYA-853
    manufactured 12-dex defects."""
    hfs = pd.DataFrame({"species": ["Fe 1", "Fe 1"],
                        "wavelength_air_A": [5000.000, 5000.004],
                        "ep_eV": [1.500, 1.500], "log_gf": [-1.0, -1.4]})
    r = _resolve([5000.002], eps=[1.500], ll=hfs)
    assert not r.resolved.iloc[0]
    assert r.unresolved_why.iloc[0].startswith("ambiguous")
    assert np.isnan(r.log_gf.iloc[0]), "an unresolved line must carry NO gf"


# ── 4. the parallel arrays cannot slip ───────────────────────────────────────────────

def test_a_length_mismatch_raises_rather_than_keying_on_a_neighbour():
    with pytest.raises(ValueError, match="per-line parallel arrays"):
        _resolve([5000.0, 5000.015], eps=[1.5])


def test_a_generator_of_wavelengths_still_works():
    """`for_lines` builds these from comprehensions today, but sizing the EP list by
    consuming the wavelengths would empty the loop for any caller that passes an
    iterator — a silent zero-line pool, which `decide` reports as rung 1."""
    r = _resolve((w for w in [5000.000]), eps=[0.990])
    assert len(r) == 1 and bool(r.resolved.iloc[0])


# ── 5. no silent fallback at the emitters ────────────────────────────────────────────

def test_the_emitter_raises_on_a_source_row_with_no_ep():
    """A blank identity column is worse than no column: a consumer cannot tell "this route
    carries no EP" from "the EP is missing for this line" once it exists but is empty."""
    from pipeline.band_products import carried_ep
    for bad in ({"wave_air_A": 5000.0, "ep_eV": np.nan}, {"wave_air_A": 5000.0}):
        with pytest.raises(ValueError, match="no ep_eV"):
            carried_ep(pd.Series(bad), wavelength_A=5000.0, element="Fe", ion="I")
    ok = pd.Series({"wave_air_A": 5000.0, "ep_eV": 2.404})
    assert carried_ep(ok, wavelength_A=5000.0, element="Fe", ion="I") == 2.404


@pytest.mark.parametrize("mod", ["measure_band_profilefit", "measure_band_ew"])
def test_both_emitters_call_the_shared_carry_and_define_no_copy(mod):
    """ONE HOME. The first cut of this fix wrote `carried_ep` into both drivers, which is
    the RYA-845/855/869 shape three tickets deep: a rule at two call sites drifts between
    them and each copy is internally consistent while the pair is wrong. Checked in the
    source because importing either driver needs the Kitt Peak atlas."""
    src = (ROOT / "scripts" / f"{mod}.py").read_text()
    assert "def carried_ep" not in src, f"{mod} defines its own copy of the carry"
    assert "carried_ep(" in src, f"{mod} does not carry the EP at all"
    assert "carried_ep" in src[:src.index("def ")], f"{mod} must import it"


def test_the_accounting_table_actually_has_the_column_the_emitter_reads():
    """POSITIVE CONTROL. The guard above proves the emitter refuses a missing EP; this
    proves the EP is THERE, so the refusal is not the normal path (RYA-833: an absence
    needs a positive control)."""
    acc = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"
    d = pd.read_csv(acc, usecols=["element", "ion", "wave_air_A", "ep_eV"])
    assert "ep_eV" in d.columns
    fe = d[(d.element == "Fe") & (d.ion == "I")]
    assert len(fe) > 100, "no Fe I rows — the control cannot discriminate"
    assert fe.ep_eV.notna().all(), "the emitter would raise on real rows"


# ── 6. the field reaches the artifact ────────────────────────────────────────────────

def test_ep_reaches_the_per_line_csv():
    """RYA-843: `red_chi2` was returned by the fitter for its whole life and never reached
    the per-line CSV, because there was no field for it and no check that would notice."""
    import derive_band_products as D
    lm = LineMeasurement(element="Fe", ion="I", wavelength_air_A=5000.0,
                         instrument="kpno_solar_atlas", ew_mA=40.0, ew_method="t",
                         abundance=7.5, ep_eV=2.404)
    assert D.asdict_line(lm)["ep_eV"] == pytest.approx(2.404)
    assert "ep_eV" not in D.ASDICT_LINE_OMITTED


def test_a_route_that_does_not_carry_an_ep_still_works_and_stays_narrow():
    """Default None means "this route does not carry one", and must not become 0 eV —
    which would key every such line on a level none of them has."""
    assert LineMeasurement(element="Fe", ion="I", wavelength_air_A=5000.0,
                           instrument="i", ew_mA=1.0, ew_method="t").ep_eV is None
    ms = [LineMeasurement(element="Fe", ion="I", wavelength_air_A=5000.000,
                          instrument="i", ew_mA=1.0, ew_method="t", abundance=7.5,
                          ew_inversion=False)]
    with mock.patch.object(gf_rung, "linelist_frame", lambda _: _PAIR):
        r = gf_rung.for_lines("Fe", "I", ms, linelist=None)
    assert r.rung == 1 or r.n_lines == 1        # it resolved narrowly, or not at all
    assert r.n_unresolved in (0, 1)


def test_for_lines_hands_the_resolver_the_measurements_own_ep():
    """The wiring, asserted rather than assumed: a product whose lines carry an EP must
    reach the EP-keyed rule, or every emitter change above is inert."""
    seen = {}

    def _spy(element, ion, wavelengths, linelist, measured_ep_eV=None):
        seen["eps"] = list(measured_ep_eV or [])
        return pd.DataFrame([{"wavelength_air_A": 5000.0, "ep_eV": 0.99,
                              "log_gf": -1.1, "resolved": True, "unresolved_why": ""}])

    ms = [LineMeasurement(element="Fe", ion="I", wavelength_air_A=5000.0, instrument="i",
                          ew_mA=1.0, ew_method="t", abundance=7.5, ep_eV=2.404,
                          ew_inversion=False)]
    with mock.patch.object(gf_rung, "resolve_lines", _spy):
        gf_rung.for_lines("Fe", "I", ms, linelist=None)
    assert seen["eps"] == [2.404]


# ── the constants are a plateau, not a knob ──────────────────────────────────────────

def test_the_wide_tolerance_is_only_reachable_with_an_ep():
    """The wide constant must be READ once, inside the branch that has an EP.

    Checked on the AST rather than the text: a docstring naming the constant is not a
    use of it, and the first cut of this guard counted both and failed on its own prose.
    """
    import ast
    tree = ast.parse((ROOT / "pipeline" / "gf_rung.py").read_text())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "resolve_lines")
    reads = [n for n in ast.walk(fn)
             if isinstance(n, ast.Name) and n.id == "LINELIST_MATCH_TOL_EPKEY_A"]
    assert len(reads) == 1, f"read {len(reads)} times; it belongs in one branch"
    # ...and that one read is the false arm of a conditional testing for a missing EP,
    # so it is UNREACHABLE without one.
    ifexps = [n for n in ast.walk(fn) if isinstance(n, ast.IfExp)
              and isinstance(n.orelse, ast.Name)
              and n.orelse.id == "LINELIST_MATCH_TOL_EPKEY_A"]
    assert len(ifexps) == 1, "the wide tolerance is not guarded by an EP test"
    assert isinstance(ifexps[0].body, ast.Name)
    assert ifexps[0].body.id == "LINELIST_MATCH_TOL_A"


def test_the_ep_window_is_far_narrower_than_the_separation_it_must_resolve():
    """The EP tolerance is a ROUNDING tolerance (the accounting table stores 4 dp), not a
    physical one. If it ever grew to the scale of a real level separation it would stop
    discriminating — the 3125.65 pair is 1.414 eV apart."""
    assert gf_rung.EP_MATCH_TOL_EV <= 0.01
    assert (2.404 - 0.990) / gf_rung.EP_MATCH_TOL_EV > 100
