"""
tests/test_gf_rung_rya855.py — RYA-855
======================================
`derive_band_products.py` passed `gf_graded=False` to `error_budget.build()`
unconditionally at BOTH of its call sites, so a three-rung gf ladder that has existed
since RYA-850 could only ever be reached at rung 1 and every regenerated band product
charged the ungraded Kurucz 0.17 regardless of what its lines were.

WHAT THESE TESTS PIN, AND WHAT THEY DELIBERATELY DO NOT
--------------------------------------------------------
They do NOT assert any published number. The fix does not move a value and the bars it
can move belong to pools whose composition is data — asserting `syst == 0.1972` would
pin whatever the pool happens to be today, which is exactly how the RYA-832 test pinned
RYA-845's double-count (see `test_nearuv_first_class_rya832`). What is pinned is the
RELATIONSHIP: which rung a given pool is entitled to, that the routes ASK rather than
assert, and that the three refusals hold.

Every "must not" here is a failure this ladder can produce silently:

* a MIXED pool claiming the graded rung on the strength of a subset, which would
  attribute a laboratory pedigree to lines that do not have one;
* a non-Fe-I species graded through the Fe I lab table on wavelength alone, which would
  MANUFACTURE a graded pool out of a coincidence;
* the cited sigma clamped to the 0.041 bound, which would turn a measurement back into
  the assumption it supersedes.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import gf_grades, gf_rung                          # noqa: E402
from pipeline.error_budget import (GRADED_GF_SYSTEMATIC_DEX,     # noqa: E402
                                   UNGRADED_GF_SYSTEMATIC_DEX, build)

LAB_CSV = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"


# ── fixtures: a line list that is the pool, so the "value the pool used" is exact ──

def _linelist(rows: pd.DataFrame, species: str = "Fe 1"):
    """An iSpec-shaped structured array from (wavelength, EP, log gf).

    Shaped like the real thing rather than mocked: `linelist_frame` reads `wave_A`,
    `element`, `lower_state_eV` and `loggf` off the loaded list, and a dict standing in
    for it would stop testing that read.
    """
    a = np.zeros(len(rows), dtype=[("element", "U8"), ("wave_A", "f8"),
                                   ("loggf", "f8"), ("lower_state_eV", "f8")])
    a["element"] = species
    a["wave_A"] = rows.wavelength_air_A.to_numpy(dtype=float)
    a["loggf"] = rows.loggf.to_numpy(dtype=float)
    a["lower_state_eV"] = rows.elo_eV.to_numpy(dtype=float)
    return a


@pytest.fixture(scope="module")
def lab() -> pd.DataFrame:
    if not LAB_CSV.exists():
        pytest.skip(f"primary-lab Fe I gf table absent at {LAB_CSV}")
    return pd.read_csv(LAB_CSV)


@pytest.fixture(scope="module")
def lab_pool(lab) -> pd.DataFrame:
    """A dozen real Fe I lab lines in the VIS. Real, so the grader's two-key match,
    its value-agreement check and the paper citations are all genuinely exercised."""
    d = lab[(lab.wavelength_air_A > 3800) & (lab.wavelength_air_A < 6910)]
    if len(d) < 12:
        pytest.skip("too few VIS lab lines to build a pool")
    return d.head(12).reset_index(drop=True)


def _rung(pool: pd.DataFrame, element="Fe", ion="I", species="Fe 1"):
    ll = _linelist(pool, species=species)
    lines = gf_rung.resolve_lines(element, ion, pool.wavelength_air_A, ll)
    return gf_rung.decide(element, ion, lines)


# ── 1. the defect: all three rungs must be REACHABLE ─────────────────────────────

def test_all_three_rungs_are_reachable(lab_pool):
    """THE BUG THIS FILE EXISTS FOR. The ladder had three rungs and the deriver could
    reach exactly one of them. Reaching one is indistinguishable, from the outside, from
    a ladder that only has one — which is why this asserts all three, not just the fix."""
    seen = set()

    seen.add(_rung(lab_pool).rung)                                    # all lab + sigmas

    mixed = lab_pool.copy()
    mixed.loc[0, "loggf"] = float(mixed.loc[0, "loggf"]) + 0.5        # one line off-scale
    seen.add(_rung(mixed).rung)

    # Rung 2 is the graded pool whose cited sigmas do not cover it. Constructed by
    # blanking the sigmas rather than by inventing lines, so the pool is still real.
    thin = lab_pool.copy()
    r = _rung(thin)
    assert r.rung == 3
    seen.add(_no_sigma_rung(thin))

    assert seen == {1, 2, 3}, (
        f"only rungs {sorted(seen)} are reachable; error_budget has carried three since "
        f"RYA-850 and the deriver could reach one of them for its whole life")


def _no_sigma_rung(pool: pd.DataFrame) -> int:
    """Rung 2: every line graded, too few citable sigmas to describe the pool.

    Driven through the real `decide` with the grader's sigma monkeypatched to NaN, so
    the coverage refusal is exercised on the actual code path rather than restated here.
    """
    real = gf_grades.grade_line

    def blank(w, ep, gf):
        v = real(w, ep, gf)
        return gf_grades.GradeVerdict(v.gf_grade, v.gf_grade_source, float("nan"),
                                      v.gf_reference_tag, v.gf_ref_loggf,
                                      v.gf_delta_dex, v.note)
    gf_grades.grade_line = blank
    try:
        return _rung(pool).rung
    finally:
        gf_grades.grade_line = real


# ── 2. a mixed pool must NOT inherit the best rung it contains ───────────────────

def test_one_ungraded_line_drops_the_whole_pool_to_ungraded(lab_pool):
    """The rule the ticket states in as many words. A single line whose gf is not the
    laboratory value makes the pool's bar a claim about lines it does not describe."""
    graded = _rung(lab_pool)
    assert graded.gf_graded and graded.n_graded == len(lab_pool)

    mixed = lab_pool.copy()
    mixed.loc[0, "loggf"] = float(mixed.loc[0, "loggf"]) + 0.5
    r = _rung(mixed)
    assert r.rung == 1 and r.gf_graded is False
    assert r.cited_sigma_dex is None
    # ...and it must say so with the counts, not merely return False.
    assert f"{len(lab_pool) - 1} of {len(lab_pool)}" in r.reason
    assert "MIXED POOL" in r.reason


def test_the_mixed_rule_is_all_and_not_any(lab_pool):
    """VERIFY THE TEST DISCRIMINATES. A pool with ONE graded line and eleven ungraded
    ones is the case an `any()` rule would wave through, and it is the worst one: the
    bar would be a laboratory sigma describing eleven Kurucz lines."""
    mostly_bad = lab_pool.copy()
    mostly_bad.loc[1:, "loggf"] = mostly_bad.loc[1:, "loggf"] + 0.5
    r = _rung(mostly_bad)
    assert r.n_graded == 1 and r.n_lines == len(lab_pool)
    assert r.rung == 1, "one graded line out of twelve bought the pool a graded bar"


# ── 3. the species gate: a graded pool must never be MANUFACTURED ────────────────

def test_a_non_fe1_species_is_never_refereed_against_the_fe1_lab_table(lab_pool):
    """`gf_grades` is Fe I — its lab table is Fe I and `canonical_fe1()` filters to
    'Fe I' — and neither checks the species of the line handed to it. Al lines sitting
    at Fe I laboratory wavelengths would therefore match on wavelength and EP alone.

    This is the manufactured-absence rule run the other way: a manufactured PRESENCE.
    The pool below is byte-identical to the one that grades GF-LAB above; only the
    species label differs, so a pass here cannot come from the lines being different.
    """
    assert _rung(lab_pool).rung == 3, "control: as Fe I this same pool IS graded"
    r = _rung(lab_pool, element="Al", ion="I", species="Al 1")
    assert r.rung == 1 and r.gf_graded is False
    assert "no primary-laboratory gf table exists for Al I" in r.reason
    assert r.n_graded == 0


def test_fe_ii_has_no_lab_pool_and_says_so(lab_pool):
    """RYA-850/852: no primary-lab-gf pool exists for Fe II in any band. That must come
    out of the decider as a STATED absence, not as a silent False."""
    r = _rung(lab_pool, element="Fe", ion="II", species="Fe 2")
    assert r.rung == 1 and "Fe II" in r.reason


# ── 4. the cited sigma is a MEASUREMENT and is not clamped ───────────────────────

def test_the_cited_sigma_is_not_clamped_to_the_generic_bound(lab, lab_pool):
    """RYA-850 amended: for the Fe I pools the cited sigma is LARGER than the 0.041
    bound, so the bound was optimistic — and a grade-A pool would legitimately fall
    BELOW it. Both directions are checked, because a clamp in either direction turns the
    measurement back into the assumption it supersedes."""
    tight = lab[lab.e_loggf_dex <= 0.02]
    wide = lab[lab.e_loggf_dex >= 0.05]
    if len(tight) < 5 or len(wide) < 5:
        pytest.skip("lab table lacks both a tight and a wide sub-pool")

    r_tight = _rung(tight.head(8).reset_index(drop=True))
    r_wide = _rung(wide.head(8).reset_index(drop=True))
    assert r_tight.rung == r_wide.rung == 3
    assert r_tight.cited_sigma_dex < GRADED_GF_SYSTEMATIC_DEX < r_wide.cited_sigma_dex, (
        f"cited sigma {r_tight.cited_sigma_dex:.4f} / {r_wide.cited_sigma_dex:.4f} did "
        f"not straddle the {GRADED_GF_SYSTEMATIC_DEX} bound — it is being clamped")


def test_the_cited_source_names_the_papers(lab_pool):
    """`error_budget.cited_gf_term` refuses an unsourced sigma; this refuses a source
    that is not a citation. Built from the lines actually in the pool, so it cannot name
    a paper the pool did not use."""
    r = _rung(lab_pool)
    assert r.rung == 3
    assert "DOI" in r.cited_source and "assumed" not in r.cited_source.lower()
    for tok in r.cited_source.split(";"):
        assert any(ch.isdigit() for ch in tok), f"citation fragment {tok!r} names no year"


# ── 5. a line we cannot price is UNGRADEABLE, not ungraded-by-default ────────────

def test_a_line_absent_from_the_loaded_linelist_forces_rung_one(lab_pool):
    """A grade describes the log gf the pool ACTUALLY used (RYA-799). If the line is not
    in the list the run loaded, there is no such value — so the pool cannot be graded,
    and the count must be visible rather than absorbed."""
    ll = _linelist(lab_pool.iloc[1:], species="Fe 1")     # drop one line from the LIST
    lines = gf_rung.resolve_lines("Fe", "I", lab_pool.wavelength_air_A, ll)
    r = gf_rung.decide("Fe", "I", lines)
    assert r.n_unresolved == 1
    assert r.rung == 1
    assert "UNRESOLVED" in r.grade_counts


def test_an_ambiguous_wavelength_is_unresolved_rather_than_iloc_zero(lab_pool):
    """RYA-853: a wavelength window is not a unique line ID, and `iloc[0]` there is how
    that ticket manufactured 12-dex 'defects'. Two same-species rows inside the
    tolerance must come back unresolved."""
    dup = pd.concat([lab_pool, lab_pool.iloc[[0]]], ignore_index=True)
    ll = _linelist(dup)
    lines = gf_rung.resolve_lines("Fe", "I", lab_pool.wavelength_air_A, ll)
    assert int((~lines.resolved).sum()) == 1
    assert "same-species rows" in lines[~lines.resolved].iloc[0].unresolved_why


def test_the_lab_table_itself_is_unambiguous_under_the_graders_match_window(lab):
    """POSITIVE CONTROL on `gf_grades._nearest`, which resolves a multi-row match with
    `argmin` — a silent choice. It is only safe while no two lab rows sit inside BOTH
    the 0.02 A and the 0.02 eV window of each other. That is true today (0 of 465); this
    test exists so that extending the table can never quietly make it untrue."""
    w = lab.wavelength_air_A.to_numpy(dtype=float)
    e = lab.elo_eV.to_numpy(dtype=float)
    twins = [i for i in range(len(lab))
             if ((np.abs(w - w[i]) <= gf_grades.WAVE_TOL_A)
                 & (np.abs(e - e[i]) <= gf_grades.EP_TOL_EV)).sum() > 1]
    assert not twins, (
        f"{len(twins)} lab row(s) now have a twin inside the grader's match window, so "
        f"`_nearest`'s argmin has become a CHOICE — decide it explicitly")


# ── 6. the decision reaches the budget, and only through one door ────────────────

@pytest.mark.parametrize("kwargs,expect", [
    ({"gf_graded": False}, "gf scale (UNGRADED)"),
    ({"gf_graded": True}, "gf scale (NIST-graded)"),
    ({"gf_graded": True, "cited_gf_sigma_dex": 0.0522,
      "cited_gf_source": "Ruffoni et al. 2014, MNRAS 441, 3127"}, "gf scale (cited lab)"),
])
def test_budget_kwargs_select_the_intended_term(kwargs, expect):
    b = build("Fe", 5000.0, 40, scatter_dex=0.10, harness_residual_dex=0.0129,
              handler="ProfileFitHandler", **kwargs)
    names = [t.name for t in b.terms]
    assert expect in names
    assert len([n for n in names if n.startswith("gf scale")]) == 1, (
        "more than one gf term reached the budget — the rungs describe the SAME "
        "quantity, so carrying two double-counts the oscillator strengths")


def test_budget_kwargs_round_trip_through_build(lab_pool):
    """The decider's output must be directly consumable. Passing the graded flag without
    its cited sigma is a caller error `error_budget.build` raises on, and keeping the
    pair inside one mapping is what makes that unreachable."""
    r = _rung(lab_pool)
    b = build("Fe", 5000.0, r.n_lines, scatter_dex=0.10, harness_residual_dex=0.0129,
              handler="ProfileFitHandler", **r.budget_kwargs())
    gf = [t for t in b.terms if t.name.startswith("gf scale")]
    assert len(gf) == 1 and gf[0].name == r.term_name


# ── 7. the wiring itself: both call sites ASK, neither ASSERTS ───────────────────

def _code(path: Path) -> str:
    """Source with `#` comment lines removed — these checks are about what the module
    DOES. The comments explaining the defect necessarily quote it, and a check that
    cannot tell the fix from its own explanation is not checking anything."""
    return "\n".join(l for l in path.read_text().splitlines()
                      if not l.strip().startswith("#"))


def test_neither_route_hardcodes_the_rung():
    """THE DEFECT, AS SOURCE. `gf_graded=False` appeared literally at both call sites;
    a fix that left either one behind would leave the ladder unreachable on that route
    and the two routes disagreeing — the RYA-847 shape, where one of two routes bypassed
    the decider for its whole life."""
    src = _code(ROOT / "scripts" / "derive_band_products.py")
    assert "gf_graded=False" not in src and "gf_graded=True" not in src, (
        "derive_band_products still states a gf rung of its own")
    assert src.count("rung.budget_kwargs()") == 2, (
        "expected exactly two budget call sites to consume the decision")


def test_both_routes_call_the_same_decider():
    """Do NOT re-implement the rung decision separately in the two call sites — the
    ticket's own instruction, and the defect shape RYA-845/847 keep finding. Checked as
    'both routes reach the module', which a second in-file implementation cannot fake."""
    src = _code(ROOT / "scripts" / "derive_band_products.py")
    synth = src[src.index("def synthesis_route"):src.index("def asdict_line")]
    ew = src[src.index("def main()"):]
    assert "gf_rung.for_lines(" in synth
    assert "gf_rung.for_product(" in ew
    for route, name in ((synth, "synthesis"), (ew, "EW")):
        assert "GRADED_GF_SYSTEMATIC" not in route and "UNGRADED_GF_SYSTEMATIC" not in route, (
            f"the {name} route names a gf systematic constant directly, which is a "
            f"second declaration of the ladder")


def test_the_coverage_threshold_is_declared_once():
    """RYA-845: a number with two homes drifts. `rya850_graded_products` applies the same
    refusal and must import it rather than restate it."""
    src = _code(ROOT / "scripts" / "rya850_graded_products.py")
    assert "from pipeline.gf_rung import CITED_COVERAGE_MIN" in src
    assert "CITED_COVERAGE_MIN = 0.90" not in src


# ── 8. aggregate membership: an excluded line must not price the bar ─────────────

def test_only_lines_in_the_aggregate_decide_the_rung(lab_pool):
    """A quarantined line is not part of the measurement, so its pedigree is not part of
    the measurement's uncertainty. The inverse would let a line the product REFUSED drag
    an otherwise-graded pool to rung 1 — a bar set by a line that contributed nothing."""
    from pipeline.band_products import LineMeasurement

    ms = []
    for i, r in lab_pool.iterrows():
        ms.append(LineMeasurement(element="Fe", ion="I",
                                  wavelength_air_A=float(r.wavelength_air_A),
                                  instrument="kpno_solar_atlas", ew_mA=float("nan"),
                                  ew_method="synthesis flux-fit", abundance=7.5,
                                  treatment="1D-LTE", ew_inversion=False))
    ll = _linelist(lab_pool)
    assert gf_rung.for_lines("Fe", "I", ms, linelist=ll).rung == 3

    # Now quarantine one line AND take it off the list it was priced from: excluded, so
    # it must not be consulted at all.
    ms[0].in_aggregate = False
    ms[0].excluded_reason = "quarantined for this test"
    ll2 = _linelist(lab_pool.iloc[1:])
    r = gf_rung.for_lines("Fe", "I", ms, linelist=ll2)
    assert r.n_lines == len(lab_pool) - 1 and r.n_unresolved == 0 and r.rung == 3
