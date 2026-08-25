"""Stage 4 is ROUTE-AWARE, and a synth line never meets the red_chi2 cap — RYA-1041.

`stage4_fit` gated EVERY line on `red_chi2 <= max_red_chi2`. On the synthesis route that
statistic is disqualified, and it is disqualified by measurement rather than by argument:

  * the GRADED laboratory anchor -- lines whose gf is KNOWN-good -- has synth `red_chi2`
    median 67 and max 3404 (RYA-1041 D1). Any conventional cap deletes the reference.
  * RYA-981 measured corr(red_chi2, |A - meanA|) = 0.025 on this route: it does not
    predict abundance error at all.
  * the two routes differ by five and a half orders of magnitude on the SAME lines --
    EW median 2.2e-4 against synth 67 -- so one number cannot serve both.

The synth branch gates on RYA-992's `frac_rise_weaker` at its PER-ARM cut, already
derived and validated on the 163-line fixture (PRs #352/#363).

The regression that matters is the third test: the graded anchor is known-good BY
CONSTRUCTION, so a gate that rejects it is wrong no matter how principled it looks.
"""
from __future__ import annotations

import numpy as np
import pytest

from pipeline.gf_empirical import NO, RC_BAD_FIT, RC_OK, RouteNotDeclared, Thresholds, YES, stage4_fit
from pipeline.synth_gof import SYNTH_GOF_REF_CUT, reference_arm, synth_gof_cut


def _th(cap=5.0):
    t = Thresholds()
    object.__setattr__(t, "max_red_chi2", cap) if getattr(type(t), "__setattr__", None) is None else None
    t.max_red_chi2 = cap
    return t


def test_route_is_required_not_defaulted():
    """🔴 A missing route must REFUSE. Defaulting to the EW branch would put synth lines
    against a cap that deletes the anchor -- a silent, catastrophic mis-grade."""
    with pytest.raises(RouteNotDeclared):
        stage4_fit(70.0, None, None, 5000.0, _th())


def test_a_synth_line_never_touches_the_red_chi2_cap():
    """red_chi2 = 3404 (the anchor's own worst) must NOT fail on the synth route."""
    arm = reference_arm()
    good = synth_gof_cut(arm) * 2.0
    v, rc = stage4_fit(3404.0, None, None, 5000.0, _th(cap=5.0),
                       route="synth", instrument=arm, frac_rise_weaker=good)
    assert (v, rc) == (YES, RC_OK), (
        "a synth line was judged by red_chi2 — that is the RYA-1041 regression; the cap "
        "is EW-route-only.")


def test_the_graded_anchor_survives_stage4_on_synth():
    """The anchor is known-good by construction. A gate that deletes it is wrong."""
    arm = reference_arm()
    cut = synth_gof_cut(arm)
    # the anchor's measured synth red_chi2 spread (RYA-1041 D1): median 67, max 3404
    anchor_red_chi2 = [2.5, 8.5, 67.1, 118.0, 435.0, 3404.0]
    kept = sum(stage4_fit(r, None, None, 5000.0, _th(cap=5.0), route="synth",
                          instrument=arm, frac_rise_weaker=cut * 1.5)[0] == YES
               for r in anchor_red_chi2)
    assert kept == len(anchor_red_chi2), (
        f"{len(anchor_red_chi2) - kept} of {len(anchor_red_chi2)} anchor-like lines were "
        "rejected on the synth route")


def test_the_ew_branch_is_unchanged():
    """The EW route keeps the cap: it is the one place red_chi2 still means something."""
    assert stage4_fit(0.0002, None, None, 5000.0, _th(cap=5.0), route="ew")[0] == YES
    assert stage4_fit(9.0, None, None, 5000.0, _th(cap=5.0), route="ew") == (NO, RC_BAD_FIT)


def test_the_synth_branch_still_rejects_a_genuinely_unconstrained_fit():
    """Guard the guard: the route split must not become a pass-through. A fit whose chi2
    barely rises as A is driven weaker did not determine an abundance (RYA-992)."""
    arm = reference_arm()
    v, rc = stage4_fit(1.0, None, None, 5000.0, _th(cap=5.0), route="synth",
                       instrument=arm, frac_rise_weaker=synth_gof_cut(arm) * 0.1)
    assert (v, rc) == (NO, RC_BAD_FIT)


def test_the_synth_branch_refuses_without_an_arm():
    """The frac_rise cut is PER-ARM (RYA-992). Defaulting it would apply Kitt Peak's cut
    to another arm silently -- the RYA-913/922 no-module-instrument-constant rule."""
    with pytest.raises(RouteNotDeclared):
        stage4_fit(1.0, None, None, 5000.0, _th(), route="synth", frac_rise_weaker=2.0)


def test_admission_tol_is_set_from_a_measurement_not_a_literal():
    """RYA-1041 Step 3: every wired threshold traces to a measured input. 0.0119 is the
    anchor's own differential scatter, and it must never be tightened (RYA-161)."""
    t = Thresholds()
    assert t.admission_tol_dex == 0.0119
    src = __import__("inspect").getsource(Thresholds)
    assert "RYA-1041" in src and "differential" in src.lower(), (
        "the value carries no inline citation to the diagnostic it came from")
    assert t.admission_tol_dex < 0.1290, (
        "0.1290 is the SOLO scatter — the wrong input (RYA-898)")
