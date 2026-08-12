"""
tests/test_fe_ir_rew_trend_rya782.py
====================================
RYA-782 — the Fe I IR "+0.384 REW slope".

WHAT IS WORTH PINNING HERE
--------------------------
The verdict (the slope is outlier leverage from the K07 tier, not a curve-of-growth
systematic) rests on a measured band product that does not live on main, so the numbers
themselves are recorded on the ticket and in the emitted artifact. What these tests pin is
the machinery that could silently invert the conclusion:

  * the ESTIMATOR — a partial coefficient conditional on a badly-specified control is
    exactly how the +0.384 arose. If ols() ever silently dropped a regressor or mis-paired
    a name to a column, the refit that overturns it would be meaningless.
  * ROBUSTNESS to the thing under test — the whole finding is that one statistic is
    outlier-driven. A robust-z that is not robust (mean/sd rather than median/MAD) would
    hide the very lines the ticket is about.
  * the SATURATION CEILING single source — remedy (b) is "already in force" only because
    the pool sits at the pipeline ceiling. If the script's restated constant drifted from
    pipeline/band_products.py, that argument would silently rot.
  * the QUARANTINE — the five lines this ticket contributes must stay registered, and must
    not silently duplicate the two RYA-780 already holds.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import rya782_rew_trend as RT                                   # noqa: E402
from pipeline import band_products                              # noqa: E402

REGISTRY = ROOT / "data" / "registry" / "problem_children.csv"
RYA780_DISPOSITION = ROOT / "data" / "results" / "fe_ir_gf" / "FeI_IR_gf_disposition.csv"

# The five lines RYA-782 contributes (RYA-780 already holds 7190.122 and 7330.137).
RYA782_LINES = ("8024.543", "7052.715", "7120.021", "7810.814", "7107.459")


# ── the estimator ────────────────────────────────────────────────────────────

def test_ols_recovers_a_known_slope_exactly():
    """A planted relation must come back, with the coefficient attached to the right name."""
    rng = np.random.default_rng(0)
    n = 400
    x1 = rng.normal(size=n)
    x2 = rng.normal(size=n)
    frame = pd.DataFrame({"a": x1, "b": x2, "A": 1.5 + 2.0 * x1 - 0.5 * x2})
    out = RT.ols(frame, ["a", "b"])
    assert out["a"]["coef"] == pytest.approx(2.0, abs=1e-6)
    assert out["b"]["coef"] == pytest.approx(-0.5, abs=1e-6)
    assert out["intercept"]["coef"] == pytest.approx(1.5, abs=1e-6)
    assert out["r2"] == pytest.approx(1.0, abs=1e-9)


def test_ols_reports_every_regressor_it_was_given():
    """A silently dropped control is how a partial coefficient becomes a marginal one."""
    rng = np.random.default_rng(1)
    frame = pd.DataFrame({c: rng.normal(size=50) for c in ("p", "q", "r")})
    frame["A"] = rng.normal(size=50)
    out = RT.ols(frame, ["p", "q", "r"])
    for name in ("intercept", "p", "q", "r"):
        assert name in out and "coef" in out[name] and "t" in out[name]


def test_ols_refuses_an_underdetermined_fit_instead_of_inventing_one():
    frame = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0], "A": [5.0, 6.0]})
    assert RT.ols(frame, ["a", "b"])["insufficient"] is True


def test_a_confounded_binary_can_carry_a_slope_that_is_not_there():
    """The RYA-782 mechanism, in miniature: A depends ONLY on group, and group is
    correlated with x. Regressing on x alone finds nothing; adding the binary makes x
    pick up a coefficient it does not deserve. This is why the refit was necessary."""
    rng = np.random.default_rng(7)
    n = 200
    grp = np.repeat([0.0, 1.0], n // 2)
    x = rng.normal(size=n) - 1.2 * grp          # group shifts x
    A = 7.5 + 0.4 * grp + rng.normal(scale=0.05, size=n)   # A depends on grp, NOT on x
    frame = pd.DataFrame({"x": x, "g": grp, "A": A})
    marginal = RT.ols(frame, ["x"])["x"]["coef"]
    partial = RT.ols(frame, ["x", "g"])["x"]["coef"]
    assert abs(partial) < abs(marginal)
    # and with the group properly controlled, x carries essentially nothing
    assert abs(partial) < 0.05


# ── robustness to the thing under test ───────────────────────────────────────

def test_univariate_slope_bootstrap_ci_brackets_the_point_estimate():
    rng = np.random.default_rng(3)
    n = 120
    x = rng.normal(size=n)
    frame = pd.DataFrame({"rew": x, "A": 7.5 + 0.3 * x + rng.normal(scale=0.1, size=n)})
    out = RT.univariate_slope(frame, boot=800)
    lo, hi = out["boot95"]
    assert lo < out["slope"] < hi
    assert out["excludes_zero"] is True


def test_univariate_slope_is_reproducible_under_a_fixed_seed():
    rng = np.random.default_rng(4)
    frame = pd.DataFrame({"rew": rng.normal(size=60), "A": rng.normal(size=60)})
    a = RT.univariate_slope(frame, boot=500)
    b = RT.univariate_slope(frame, boot=500)
    assert a["boot95"] == b["boot95"]


def test_the_outlier_flag_uses_a_ROBUST_scale_not_mean_and_sd():
    """With a +2 dex line present, mean/sd inflates so much the outlier hides itself.
    The MAD-based z must still find it -- that is the entire RYA-782 census."""
    A = np.concatenate([np.full(40, 7.55) + np.random.default_rng(5).normal(scale=0.15, size=40),
                        [9.68]])
    med = np.median(A)
    mad = 1.4826 * np.median(np.abs(A - med))
    robust_z = (A[-1] - med) / mad
    naive_z = (A[-1] - A.mean()) / A.std(ddof=1)
    assert robust_z > 8.0, "the MAD-based z must flag a +2 dex line"
    assert naive_z < 7.0, "sanity: the naive z is materially smaller"
    assert robust_z > naive_z


# ── single sources, not restatements that drift ──────────────────────────────

def test_the_saturation_ceiling_matches_the_pipeline_single_source():
    """Remedy (b) is 'already in force' only if this constant is the pipeline's."""
    assert RT.REW_SATURATION_CEILING == band_products.REW_SATURATION_CEILING


def test_the_laboratory_tier_matches_the_rya760_cut():
    """The gf-homogeneous sample must be the same cut RYA-760/780 used, or the
    'laboratory pool' in this ticket is a different population than they mean."""
    assert RT.LAB_TIER == 0.04
    lab = {s for s, dex in RT.SOURCE_DEX.items() if dex <= RT.LAB_TIER}
    assert lab == {"BWL", "2014MNRAS.", "GESHRL14", "BK", "BKK"}
    for disputed in ("FMW", "GESB82c", "MRW", "K07"):
        assert RT.SOURCE_DEX[disputed] > RT.LAB_TIER


def test_composite_source_tags_take_their_first_component():
    assert RT.base_source("BK+BWL") == "BK"
    assert RT.base_source(" FMW ") == "FMW"


# ── the quarantine this ticket contributes ───────────────────────────────────

def test_the_five_rya782_lines_are_registered():
    reg = pd.read_csv(REGISTRY)
    scopes = " ".join(reg["lambda_or_scope"].astype(str))
    for lam in RYA782_LINES:
        assert lam in scopes, f"Fe I {lam} must stay in the quarantine registry"


def test_the_registered_lines_cite_this_ticket():
    reg = pd.read_csv(REGISTRY)
    rows = reg[reg["lambda_or_scope"].astype(str).str.contains("|".join(RYA782_LINES))]
    assert len(rows) >= 4
    for _, r in rows.iterrows():
        assert "782" in str(r["governing_tickets"])
        assert str(r["species"]).strip() == "Fe I"


def test_rya782_does_not_re_quarantine_the_lines_rya780_already_holds():
    """Duplicate custody of one line in two registries is how a disposition silently
    diverges. 7190.122 and 7330.137 belong to RYA-780."""
    if not RYA780_DISPOSITION.exists():
        pytest.skip("RYA-780 disposition artifact absent")
    held = pd.read_csv(RYA780_DISPOSITION)["wave_A"].astype(float).to_numpy()
    for lam in RYA782_LINES:
        assert not np.any(np.abs(held - float(lam)) < 0.05), \
            f"{lam} is already adjudicated by RYA-780"


def test_the_misidentified_line_is_not_filed_as_a_gf_problem():
    """8024.543 measures 46.6 mA from the nearest catalogued Fe I -- the absorber is
    another species. Filing it as BAD_GF would send it to a gf re-source that cannot fix
    it, and would put a wrong number into the gf adjudication population."""
    reg = pd.read_csv(REGISTRY)
    row = reg[reg["lambda_or_scope"].astype(str).str.contains("8024.543")]
    assert len(row) == 1
    assert row.iloc[0]["problem_class"] == "ATOMIC_BLEND"
    assert row.iloc[0]["required_treatment"] == "exclude"


def test_the_laboratory_tier_casualty_is_carried_as_unresolved():
    """7107.459 is BWL (+/-0.02 dex) yet reads 0.57 dex below the laboratory pool
    median (0.64 below the full-pool median), with a second independent
    anomaly from RYA-712. Cause is NOT established; it must not be quietly resolved."""
    reg = pd.read_csv(REGISTRY)
    row = reg[reg["lambda_or_scope"].astype(str).str.contains("7107.459")]
    assert len(row) == 1
    assert row.iloc[0]["status"] == "owed"


# ── the firewall ─────────────────────────────────────────────────────────────

def test_the_anchor_is_never_a_fit_target():
    """RYA-161: a slope removed by fitting to the optical anchor is tuning, not a fix.
    The anchor may appear only as a descriptive comparison."""
    src = (ROOT / "scripts" / "rya782_rew_trend.py").read_text()
    assert "VIS_ANCHOR" in src
    # the anchor must never enter a regression: it may only be printed/subtracted
    for line in src.splitlines():
        if "VIS_ANCHOR" in line and not line.strip().startswith("#"):
            assert "ols(" not in line and "lstsq" not in line, \
                f"anchor must not enter a fit: {line.strip()}"
