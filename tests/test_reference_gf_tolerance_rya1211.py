"""RYA-1211 — the grading window is a property of the SOURCE, not a constant.

`pipeline.reference_lineset` has said so in its docstring since RYA-1109 and every join it
owns takes a DERIVED per-set `match_tol_A`. `pipeline.gf_grades` is the join that never
asked: it graded AGSS21's 0.1 A-printed wavelengths through a 0.02 A window, 15 of the 21
Reference Grade Fe I lines missed their own `canonical_gf` row by 0.024-0.050 A, and
because `grade_line` has no "not found" verdict the misses fell through to the blanket
Kurucz systematic and were REPORTED as `systematic:K07 x15` — a source none of them uses.

These tests pin the fix and, just as much, the three things that must stay true for the
wider window to be legitimate rather than the RYA-1109 move of widening until a count
improves: a zero displaced null, zero ambiguity, and a global constant that does NOT move.
"""
from __future__ import annotations

import ast
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline import gf_grades as gg
from pipeline import gf_rung
from pipeline import reference_lineset as rls

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
PER_LINE = ROOT / "data/results/rya1106/kpno_kurucz2005/asplund_lines_per_line.csv"
DISPLACE_A = 0.35


@pytest.fixture(scope="module")
def pool() -> pd.DataFrame:
    d = pd.read_csv(PER_LINE)
    s = d[d["a_3dnlte"].notna()].copy()
    assert len(s) == 21, f"the Reference Grade pool is 21 lines, got {len(s)}"
    return s


def _ties(s, wtol, shift=0.0):
    """(canonical ties, lab ties, ambiguous) at this window."""
    can, lab = gg.canonical_species("Fe I"), gg.lab_lines("Fe I")
    n_can = n_lab = n_amb = 0
    for _, r in s.iterrows():
        w, ep = float(r.wavelength_air_A) + shift, float(r.elo_eV)
        for d, wc, ec, is_lab in ((can, "wavelength_air_A", "excitation_potential_eV", False),
                                  (lab, "wavelength_air_A", "elo_eV", True)):
            m = d[(np.abs(d[wc] - w) <= wtol) & (np.abs(d[ec] - ep) <= gg.EP_TOL_EV)]
            if len(m) > 1:
                n_amb += 1
            elif len(m) == 1:
                n_lab += is_lab
                n_can += not is_lab
    return n_can, n_lab, n_amb


# ── the generalisation is additive ────────────────────────────────────────────

def test_the_default_window_is_unchanged():
    """Every pre-RYA-1211 caller must grade EXACTLY as before. A fix that silently
    re-graded the whole Fe matrix would move published bars nobody asked about."""
    assert gg.WAVE_TOL_A == 0.02


def test_omitting_the_argument_grades_identically_to_passing_the_default(pool):
    for _, r in pool.iterrows():
        a = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV), float(r.loggf_asplund))
        b = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV), float(r.loggf_asplund),
                          wave_tol_A=gg.WAVE_TOL_A)
        assert a.gf_grade == b.gf_grade and a.gf_reference_tag == b.gf_reference_tag


# ── the defect itself ─────────────────────────────────────────────────────────

def test_the_agss21_pool_ties_to_canonical_gf_at_its_own_derived_window(pool):
    """🔴 THE DEFECT. At 0.02 A only 6 of 21 tie; at AGSS21's own derived window all 21 do.

    This is the assertion that fails on the pre-fix source — `grade_line` had no
    `wave_tol_A` at all, so the derived window could not be asked for.
    """
    tol = rls.grading_tol_A("asplund")
    assert _ties(pool, gg.WAVE_TOL_A)[0] == 6, "the pre-fix window ties 6 of 21"
    assert _ties(pool, tol)[0] == 21, "every line must find its own canonical_gf row"


def test_not_one_reference_grade_line_actually_carries_a_kurucz_gf(pool):
    """The premise the ticket was written on, tested: `systematic:K07 x15` described a
    LOOKUP MISS, not a Kurucz source. Graded at AGSS21's own window, the tags are
    NIST-C+, NIST ASD, primary laboratory, VALD3 and BK+BWL — and K07 appears nowhere."""
    tol = rls.grading_tol_A("asplund")
    tags = {gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV),
                          float(r.loggf_asplund), wave_tol_A=tol).gf_reference_tag
            for _, r in pool.iterrows()}
    assert "" not in tags, "a line with no row found is the defect, not a grade"
    assert not [t for t in tags if gg.K07_TAG in str(t)], f"expected no Kurucz tag, got {tags}"


def test_the_ep_axis_is_what_discriminates(pool):
    """The window is widened on WAVELENGTH only because wavelength is what AGSS21 rounded.
    EP is printed to 5 dp and agrees to better than 0.001 eV on every line — which is why
    a 2.5x wider wavelength window admits no second candidate."""
    can = gg.canonical_species("Fe I")
    tol = rls.grading_tol_A("asplund")
    for _, r in pool.iterrows():
        m = can[(np.abs(can.wavelength_air_A - float(r.wavelength_air_A)) <= tol)
                & (np.abs(can.excitation_potential_eV - float(r.elo_eV)) <= gg.EP_TOL_EV)]
        assert len(m) == 1
        assert abs(float(m.iloc[0].excitation_potential_eV) - float(r.elo_eV)) < 1e-3


# ── the three controls ────────────────────────────────────────────────────────

def test_the_displaced_null_finds_nothing(pool):
    """🔴 THE CONTROL THAT MAKES THE WIDER WINDOW LEGITIMATE. Displace the pool by 0.35 A
    and the same join must find NOTHING at any tolerance out to 0.15 A. If these were
    proximity accidents the null would collect them too."""
    for tol in (0.02, 0.05, 0.06, 0.10, 0.15):
        can, lab, _ = _ties(pool, tol, shift=DISPLACE_A)
        assert (can, lab) == (0, 0), f"displaced null found {can}/{lab} at {tol} A"


def test_the_census_plateaus_and_never_goes_ambiguous(pool):
    """A count read off a single window is an artifact waiting to happen. The census
    saturates at 21 canonical / 6 laboratory and stays there, with zero ambiguity
    throughout — so nothing here is a match picked by proximity (RYA-1034)."""
    seen = [_ties(pool, t) for t in (0.06, 0.08, 0.10, 0.15)]
    assert {s[:2] for s in seen} == {(21, 6)}, f"no plateau: {seen}"
    assert {s[2] for s in seen} == {0}, "an ambiguous match must never be counted as one"


def test_a_global_widening_is_refused_and_the_reason_is_measured():
    """Why this is a per-source argument and not a new constant: at 0.06 A the canonical
    Fe I table's own self-ambiguous rows go 26 -> 210, so every pool stated at full
    precision would start REFUSING lines that grade cleanly today."""
    can = gg.canonical_species("Fe I")
    w = np.sort(can["wavelength_air_A"].to_numpy(float))
    order = np.argsort(can["wavelength_air_A"].to_numpy(float))
    e = can["excitation_potential_eV"].to_numpy(float)[order]

    def self_ambiguous(tol):
        lo = np.searchsorted(w, w - tol, "left")
        hi = np.searchsorted(w, w + tol, "right")
        return sum(np.count_nonzero(np.abs(e[lo[i]:hi[i]] - e[i]) <= gg.EP_TOL_EV) > 1
                   for i in range(len(w)))

    narrow, wide = self_ambiguous(0.02), self_ambiguous(0.06)
    assert narrow < wide, "if widening cost nothing, a constant would have been fine"
    assert wide > 5 * narrow, f"expected a large jump, got {narrow} -> {wide}"


# ── it moves no published number ──────────────────────────────────────────────

def test_the_correction_changes_no_published_value(pool):
    """The re-grade fixes the STATED REASON, not the bar. 5 of 21 lines are GF-LAB either
    way, so the mixed-pool rule binds either way: rung 1, `gf_graded=False`, 0.17 dex.
    16 of the 21 have no primary-laboratory Fe I measurement in the repo AT ANY tolerance
    (the lab census plateaus at 6/21), so that is a real absence, not a lookup failure."""
    lines = pd.DataFrame({"wavelength_air_A": pool.wavelength_air_A.to_numpy(float),
                          "ep_eV": pool.elo_eV.to_numpy(float),
                          "log_gf": pool.loggf_asplund.to_numpy(float)})
    before = gf_rung.decide("Fe", "I", lines)
    after = gf_rung.decide("Fe", "I", lines, wave_tol_A=rls.grading_tol_A("asplund"))
    assert before.rung == after.rung == 1
    assert before.budget_kwargs() == after.budget_kwargs() == {"gf_graded": False}
    assert after.n_graded == 5 and after.n_lines == 21
    assert _ties(pool, 0.15)[1] == 6, "no more laboratory gf exists to be found"


def test_the_cited_sigma_fraction_still_cannot_price_the_pool(pool):
    """The upgrade is real but does not reach the bar: 17 of 21 lines carry a citable
    per-line sigma, and `CITED_COVERAGE_MIN` is 90%. An RMS over 81% of a pool does not
    describe the pool, so the graded route stays refused — for a stated reason."""
    tol = rls.grading_tol_A("asplund")
    v = [gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV), float(r.loggf_asplund),
                       wave_tol_A=tol) for _, r in pool.iterrows()]
    cited = [x for x in v if x.has_cited_sigma]
    assert len(cited) == 17
    assert len(cited) / len(v) < gf_rung.CITED_COVERAGE_MIN


# ── the tolerance keeps ONE home ──────────────────────────────────────────────

def test_the_agss21_caller_derives_its_window_and_does_not_write_one():
    """🔴 A LITERAL HERE WOULD BE THE DEFECT AGAIN, one layer up. The tolerance is a
    property of the source, so the caller must ASK the registry. Asserted on the AST, so
    a later hand-tuned float cannot pass by being numerically equal today."""
    src = ROOT / "scripts" / "rya1106_asplund_replication.py"
    tree = ast.parse(src.read_text())
    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call)
             and getattr(n.func, "id", getattr(n.func, "attr", None)) == "budget_from_pool"]
    assert calls, "budget_from_pool call not found — did the entry point move?"
    for c in calls:
        kw = {k.arg: k.value for k in c.keywords}
        assert "wave_tol_A" in kw, "the AGSS21 budget must state its own grading window"
        node = kw["wave_tol_A"]
        # RYA-1212 hoisted the window into a local so the empirical-sigma call can share
        # it, so follow one level of assignment. A NAME is only acceptable if the thing
        # bound to it is itself the derived call — a literal bound to `tol` must still
        # fail, which is the whole point of the check.
        if isinstance(node, ast.Name):
            bound = [t for st in ast.walk(tree) if isinstance(st, ast.Assign)
                     for tgt in st.targets
                     if isinstance(tgt, ast.Name) and tgt.id == node.id
                     for t in [st.value]]
            assert bound, f"wave_tol_A is {node.id!r} but nothing assigns it"
            node = bound[-1]
        assert isinstance(node, ast.Call), \
            f"wave_tol_A must be DERIVED, not written: got {ast.dump(node)[:80]}"
        assert getattr(node.func, "attr", None) == "grading_tol_A", \
            "derive it from reference_lineset.grading_tol_A — the tolerance has one home"


def test_the_derived_window_matches_agss21s_printed_resolution(pool):
    """0.05 A is not a taste: AGSS21 prints lambda in nm to 2 dp = 0.1 A, and a value
    printed to a step S lies within +/-S/2 of truth. The registry states that basis."""
    spec = rls.SETS["asplund"]
    assert spec.match_tol_A == 0.05
    assert "printed" in spec.tol_basis.lower()
    assert max(len(str(w).split(".")[-1]) for w in pool.wavelength_air_A) == 1
    # closed at the edge: 5247.0 sits EXACTLY 0.05 A from the GES row at 5247.05
    assert rls.grading_tol_A("asplund") > spec.match_tol_A
    assert _ties(pool.head(1), spec.match_tol_A)[0] == 0, "5247.0 is the closed-edge case"
    assert _ties(pool.head(1), rls.grading_tol_A("asplund"))[0] == 1
