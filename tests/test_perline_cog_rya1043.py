"""Saturation is a PER-LINE property, and stage 2 gates on it — RYA-1043.

The scalar `rew_max` is retired. It could not be derived: an ensemble curve of growth
for solar Fe I has NO resolvable knee, and that is physics rather than method. The
clean-and-untruncated 1805-line intersection IS a genuine COG (corr +0.362, medians
monotone -5.22 -> -4.77) and its local dREW/dX still carries three sign changes, because
damping, blending and core depth vary line to line and smear the transition. A
single-line COG is sharp; an ensemble COG is broad BY CONSTRUCTION.

So stage 2 asks "has THIS line's own response flattened", using dREW/dA measured per line
by `scripts/rya1043_perline_cog.py`.

The property that matters most here is the third test: a line with NO measured derivative
must be UNKNOWN, never a pass. The per-line COG is expensive and incomplete, and treating
"not yet computed" as "linear" would admit exactly the saturated lines the stage exists
to catch (RYA-833).
"""
from __future__ import annotations

import numpy as np
import pytest

import pipeline.gf_empirical as G


def _th(**over):
    t = G.Thresholds(rew_min=-6.5, min_d_rew_dA=0.5, min_depth=0.03, max_red_chi2=5.0,
                     admission_tol_dex=0.20, min_anchor_lines=5,
                     consistency_bound_dex=0.30, fallback_sigma_dex=0.50)
    for k, v in over.items():
        setattr(t, k, v)
    return t


def test_the_scalar_rew_max_is_gone():
    """🔴 Keeping the field would invite someone to set it to Gray's -4.9 again — a
    borrowed textbook constant, which is the class RYA-968 3.1 exists to reject."""
    assert not hasattr(G.Thresholds(), "rew_max"), (
        "rew_max is retired: saturation cannot be expressed as one number for all lines")


def test_a_linear_line_passes_and_a_saturated_one_does_not():
    """Measured over 617 VIS lines, dREW/dA runs 1.044 down to -0.129 and falls
    monotonically with line strength (corr -0.807)."""
    assert G.stage2_saturation(-5.2, _th(), d_rew_dA=0.92)[0] == G.YES
    v, rc = G.stage2_saturation(-4.4, _th(), d_rew_dA=0.37)
    assert (v, rc) == (G.NO, G.RC_SATURATION)


def test_a_line_with_no_measured_slope_is_UNKNOWN_never_a_pass():
    """🔴 The rule that keeps an unfinished computation from reading as a passed test."""
    assert G.stage2_saturation(-5.2, _th(), d_rew_dA=None)[0] == G.UNKNOWN
    assert G.stage2_saturation(-5.2, _th(), d_rew_dA=float("nan"))[0] == G.UNKNOWN
    assert G.stage2_saturation(-5.2, _th())[0] == G.UNKNOWN      # not supplied at all


def test_the_weak_end_floor_still_bites():
    """`rew_min` survives, because "too weak to measure" is a DIFFERENT question from
    "saturated" and a curve of growth does not answer it."""
    v, rc = G.stage2_saturation(-7.0, _th(rew_min=-6.5), d_rew_dA=0.99)
    assert (v, rc) == (G.NO, G.RC_SATURATION), (
        "a line below the noise floor must fail even though it is perfectly linear")


def test_the_flattening_threshold_is_declared_or_nothing():
    """min_d_rew_dA joins the raise-until-set contract: the number is Ryan's (RYA-1041
    Step 2), derived from the measured distribution, never a default."""
    with pytest.raises(G.ThresholdNotDeclared):
        G.Thresholds().require("min_d_rew_dA")


def test_the_measured_distribution_is_on_disk_and_spans_the_transition():
    """Guard against the gate being wired to data that does not exist or cannot
    discriminate. The artifact must reach the linear limit AND the saturated floor,
    otherwise no threshold placed in it means anything."""
    import pandas as pd
    from pathlib import Path
    p = (Path(__file__).resolve().parents[1] / "data" / "audit" / "rya1041"
         / "rya1041_perline_cog_4200_6910.csv")
    if not p.exists():
        pytest.skip("per-line COG artifact not present in this checkout")
    d = pd.read_csv(p)
    d = d[d.status == "ok"]
    assert len(d) > 300, f"only {len(d)} lines measured — too few to set a threshold from"
    lo, hi = d.d_rew_dA.min(), d.d_rew_dA.max()
    assert hi > 0.9, f"nothing reaches the linear limit (max {hi:.3f}) — the blend "
    assert lo < 0.2, f"nothing reaches saturation (min {lo:.3f})"
    r = np.corrcoef(d.rew, d.d_rew_dA)[0, 1]
    assert r < -0.5, (
        f"corr(REW, dREW/dA) = {r:+.3f}; the derivative must FALL with line strength or "
        "it is not measuring saturation. It was -0.333 while blend-contaminated and "
        "-0.807 after the blend baseline was subtracted.")
