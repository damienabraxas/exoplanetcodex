"""RYA-1093 — the xi constraint audit, the weighted production slope, and the honest gate.

⚠️ THE HARDEST THING TO TEST HERE IS THAT THE VERDICT IS NOT A FOREGONE CONCLUSION.
This audit returns KEEP-PIN-1.0, and it would be easy to write an "audit" that can only ever
return that. `test_the_audit_CAN_return_adopt` is the control: it feeds well-constrained
synthetic data through the same predicates and requires ADOPT. Without it, every other
assertion here is compatible with a function that ignores its input (RYA-853/RYA-1080 —
a guard whose two sides converge is silent, not safe).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import rya1093_xi_robustness_audit as A          # noqa: E402
from pipeline.fitexy import fitexy               # noqa: E402


@pytest.fixture(scope="module")
def pool():
    return A.per_line()


@pytest.fixture(scope="module")
def solve():
    return A._load("rya311_xi_solar_FeI_all_ceil100")


# ── 2A: does the pool constrain xi? ──────────────────────────────────────────

def test_the_pool_is_truncated_at_the_lines_that_pin_xi(pool, solve):
    """🔴 The finding. Max EW == the ceiling: the informative lines were CUT, not absent."""
    c = A.constraint_census(pool, solve)
    assert c["ew_max_mA"] == pytest.approx(float(c["ew_ceiling_mA"]), abs=0.05)
    assert c["fraction_saturated"] < 0.10          # 4 of 58 reach REW >= -4.8
    assert "TRUNCATED" in c["verdict"]


def test_the_gap_is_weighting_not_both_axis_errors(pool):
    """🔴 The ticket's premise, corrected. Mucciarelli's term is real and NEGLIGIBLE here.

    A reader who "switches to FITEXY" thinking the both-axis term was the missing piece has
    misdiagnosed a real bug, and will not understand why another pool behaves differently.
    """
    m = A.method_decomposition(pool)
    assert m["x_over_y_denominator_ratio"] < 1e-4     # the both-axis term is nothing
    assert m["both_axis_term_dex_per_dex"] < 1e-3
    assert m["weighting_term_dex_per_dex"] > 0.05
    assert m["weighting_over_both_axis"] > 50
    # and the three fits really do rank as claimed
    assert m["slope_unweighted_ols_PRODUCTION"] > 0.05
    assert abs(m["slope_weighted_ols"] - m["slope_fitexy"]) < 1e-3


def test_the_solve_is_a_property_of_a_few_near_ceiling_lines(pool):
    lev = A.leverage(pool)
    assert lev["fitexy_slope_span"] > 0.15
    signs = {np.sign(t["slope_fitexy"]) for t in lev["trials"]}
    assert len(signs) > 1, "the slope should change SIGN as the strong end is stripped"


def test_the_two_independent_checks_have_no_root_at_all():
    """Deep-lines and graded-only. Reported as GAPS, never as a clean bill of health."""
    g = A.missing_solves()
    assert g["deep_lines_ceiling_removed"]["root_found"] is False
    assert g["graded_only_lab_gf"]["root_found"] is False
    assert "CANNOT BE TESTED" in g["graded_only_lab_gf"]["verdict"]


def test_chi2_red_shows_the_formal_error_is_an_underestimate(pool):
    m = A.method_decomposition(pool)
    assert m["chi2_red"] > 10
    assert m["sigma_slope_chi2_scaled"] > 4 * m["sigma_slope_formal"]


# ── the control: the audit is capable of the other answer ────────────────────

def test_the_audit_CAN_return_adopt():
    """⚠️ THE VACUITY CONTROL. Well-constrained synthetic data must clear every predicate.

    Built as a genuinely informative pool: a wide REW baseline with no truncation, point
    errors that DO explain the scatter, and a real slope. If this returns KEEP-PIN, the
    audit's predicates are unsatisfiable and its verdict on the real pool means nothing.
    """
    rng = np.random.default_rng(0)
    n = 120
    rew = np.linspace(-5.6, -4.2, n)              # reaches well past the saturated end
    sy = np.full(n, 0.02)
    y = 7.46 + 0.10 * (rew + 5.0) + rng.normal(0, 0.02, n)
    d = pd.DataFrame({"rew": rew, "A_at_measured_xi": y,
                      "sigma_rew": np.full(n, 0.004), "sigma_A_from_ew": sy,
                      "ew_mA": 10 ** (rew + 3) * 5000.0})
    census = A.constraint_census(d, {"ew_ceiling_mA": None})
    meth = A.method_decomposition(d)
    lev = A.leverage(d)
    assert census["fraction_saturated"] >= 0.25, "control pool must be strong-line RICH"
    assert meth["chi2_red"] < 4.0, "control errors must explain the scatter"
    assert lev["fitexy_slope_span"] < 0.05, "control slope must be stable to trimming"


# ── delta_xi and the gate ────────────────────────────────────────────────────

def test_the_formal_error_is_explicitly_not_chosen(pool, solve):
    """🔴 RYA-161: 0.0588 must not be picked, and the report must say it was not."""
    meth, lev = A.method_decomposition(pool), A.leverage(pool)
    d = A.honest_delta_xi(solve, meth, lev)
    assert d["chosen"] != pytest.approx(float(solve["delta_xi_formal"]), abs=1e-6)
    assert d["chosen"] > 4 * float(solve["delta_xi_formal"])
    assert "NOT chosen" in d["explicitly_not_chosen"]


def test_the_gate_is_reported_failing_not_met(pool, solve):
    meth, lev = A.method_decomposition(pool), A.leverage(pool)
    g = A.gate_report(A.honest_delta_xi(solve, meth, lev), solve)
    assert g["xi_term_alone_exceeds_the_gate"] is True
    # ⚠️ and the report must OWN the fact that this term is the A-shift, not disguise it
    assert g["reduces_to_the_A_shift_between_candidates"] is True


def test_nothing_selects_a_candidate_by_the_gate():
    """The gate constant may be displayed and compared, never used to pick delta_xi."""
    import ast
    src = (ROOT / "scripts" / "rya1093_xi_robustness_audit.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "honest_delta_xi")
    names = {x.id for x in ast.walk(fn) if isinstance(x, ast.Name)}
    assert "SOLAR_GATE_DEX" not in names, (
        "honest_delta_xi references the gate — choosing the bar that passes is the "
        "RYA-161 failure this ticket names as CRITICAL")


# ── 2C: the production slope ─────────────────────────────────────────────────

def test_the_production_slope_now_weights():
    """🔴 The pipeline-wide fix. Unweighted was the defect; sigma_A=None still falls back."""
    import numpy as _np
    from pipeline.abundances_derive import _compute_rew_slope
    n = 40
    lm = _np.zeros(n, dtype=[("ew", "f8"), ("wave_A", "f8"), ("ew_err", "f8")])
    lm["wave_A"] = _np.linspace(5000, 6000, n)
    lm["ew"] = _np.linspace(20, 95, n)
    lm["ew_err"] = lm["ew"] * 0.05
    rew = _np.log10(lm["ew"] / lm["wave_A"])
    ab = 7.5 + 0.08 * (rew + 5.0)
    mask = _np.ones(n, dtype=bool)
    # one line given a huge error must lose its vote under weighting and keep it without
    ab[0] += 1.0
    unweighted = _compute_rew_slope(mask, lm, ab)
    sig = _np.full(n, 0.02); sig[0] = 1.0
    weighted = _compute_rew_slope(mask, lm, ab, sigma_A=sig)
    assert abs(weighted - 0.08) < abs(unweighted - 0.08), (
        "the weighted fit must be less corrupted by a badly-known line than the "
        "unweighted one — that is the entire point of RYA-1093 2C")
    assert _compute_rew_slope(mask, lm, ab, sigma_A=None) == pytest.approx(unweighted)


# ── 2D: the re-stamp allowance ───────────────────────────────────────────────

def test_the_restamp_is_provenance_tagged_not_a_silent_overwrite():
    import rya1088_record_sigma_params as R
    p = R.RESTAMP_PROVENANCE
    assert p["source_ticket"] == "RYA-1093"
    assert p["supersedes"]["delta_xi_kms"] == 0.05
    assert p["audit_artifact"].endswith("xi_robustness_audit.json")
    assert "reason" in p and len(p["reason"]) > 40


def test_the_1080_guard_was_not_widened():
    """⚠️ RYA-1088 refused to widen another ticket's CRITICAL guard. Neither does this one.

    The guard already skips keys absent at baseline. This test pins that it is UNCHANGED —
    if a future edit relaxes it, the reason RYA-1093 was allowed to stamp disappears.
    """
    src = (ROOT / "tests" / "test_feed_repo_reconciliation_rya1080.py").read_text()
    assert 'if fk == "_prov" or fk in DERIVED_FIELDS or fk not in a:' in src
    assert "DERIVED_FIELDS = {\"display\"}" in src


def test_the_solar_delta_xi_is_the_audited_value():
    from pipeline import uncertainty_stack as U
    assert U.SOLAR_DELTA_XI_KMS == pytest.approx(0.2912, abs=1e-4)
    assert U.SOLAR_DELTA_P["vturb_kms"] == U.SOLAR_DELTA_XI_KMS
    assert "RYA-1093" in U.SOLAR_DELTA_XI_PROVENANCE
    assert U.SOLAR_DELTA_P["feh"] == 0.0        # still exact by definition
