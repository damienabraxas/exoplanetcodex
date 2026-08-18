"""RYA-847 — the constraint metrics, and the plumbing that must carry them.

The metric tests assert RELATIONSHIPS and closed forms, never magic constants (RYA-845).

The plumbing tests are the ones that would have caught the original defect. `red_chi2`
was returned by `_fit_synth_flux` for its entire life and never reached the per-line CSV,
because `LineMeasurement` had no field for it and `asdict_line` is a hand-written
projection that silently drops whatever it does not mention. Nothing failed. A quantity
with nowhere to live is a quantity nobody can gate on, and that is how two unconstrained
NIR lines reached a published aggregate (RYA-843).
"""
from __future__ import annotations

import math
from dataclasses import fields

import numpy as np
import pytest

from pipeline.fit_constraint import (CURVATURE_PROBE_STEP_DEX as STEP,
                                     ConstraintMetrics, curvature_sigma,
                                     measure_constraint)


def _parabola(a0, curv, floor=0.0):
    return lambda a: floor + curv * (a - a0) ** 2


# ── sigma_A: the bracket-free metric ──────────────────────────────────────────

@pytest.mark.parametrize("red_chi2", [1.0, 4.0, 25.0, 66.553, 1226.1])
def test_sigma_A_scales_as_sqrt_red_chi2(red_chi2):
    """THE RELATIONSHIP. The same curvature against a chi2 N times too large must give a
    sigma sqrt(N) times larger — rescaling to red_chi2 == 1 is the entire correction."""
    f = _parabola(7.5, 1000.0)
    kw = dict(a_best=7.5, chi2_min=0.0, a_lo=4.5, a_hi=12.5, edge=False)
    unit = curvature_sigma(f, red_chi2=1.0, **kw)
    got = curvature_sigma(f, red_chi2=red_chi2, **kw)
    assert got == pytest.approx(unit * math.sqrt(red_chi2), rel=1e-12)


def test_sigma_A_is_independent_of_the_bracket():
    """🔴 THE PROPERTY THAT MAKES IT UNIVERSAL, and the one `frac_rise` lacks.

    847 item 3 wants a criterion that transfers across bands. Handlers bracket
    differently, so any metric that moves when only the bracket moves cannot be compared
    between them. sigma_A must not.
    """
    f = _parabola(7.5, 800.0, floor=250.0)
    wide = measure_constraint(f, a_best=7.5, chi2_min=250.0, red_chi2=9.0,
                              a_lo=1.0, a_hi=20.0)
    narrow = measure_constraint(f, a_best=7.5, chi2_min=250.0, red_chi2=9.0,
                                a_lo=7.0, a_hi=8.0)
    assert wide.sigma_A == pytest.approx(narrow.sigma_A, rel=1e-12)
    # and the contrast: frac_rise moves by orders of magnitude on the same fit
    assert not (wide.frac_rise_weaker ==
                pytest.approx(narrow.frac_rise_weaker, rel=0.1))


def test_sigma_A_is_independent_of_an_additive_baseline():
    """Delta-chi2 is a DIFFERENCE, so the depressed-baseline offset that makes an
    ABSOLUTE red_chi2 untransferable cancels exactly. Measured on RYA-843's own numbers:
    a good NIR line sits at red_chi2 = 72 and another good one at 2.6."""
    curv, a0 = 500.0, 7.5
    bare = curvature_sigma(_parabola(a0, curv), a_best=a0, chi2_min=0.0, red_chi2=4.0,
                           a_lo=4.5, a_hi=12.5, edge=False)
    OFF = 1.0e6
    shifted = curvature_sigma(_parabola(a0, curv, floor=OFF), a_best=a0, chi2_min=OFF,
                              red_chi2=4.0, a_lo=4.5, a_hi=12.5, edge=False)
    assert shifted == pytest.approx(bare, rel=1e-12)


def test_a_railed_fit_is_not_measurable_never_a_number():
    """RYA-848: the old arithmetic returned 0.000 here — the tightest possible bar for
    the least constrained possible fit."""
    a_hi = 12.496
    got = curvature_sigma(lambda a: 1000.0 - 10.0 * a, a_best=a_hi,
                          chi2_min=1000.0 - 10.0 * a_hi, red_chi2=50.0,
                          a_lo=4.496, a_hi=a_hi, edge=True)
    assert math.isnan(got)


def test_a_flat_objective_is_not_measurable():
    """No curvature to invert. NaN says so; 1.000 would read as a measured 1 dex."""
    got = curvature_sigma(lambda a: 4242.0, a_best=7.5, chi2_min=4242.0, red_chi2=10.0,
                          a_lo=4.5, a_hi=12.5, edge=False)
    assert math.isnan(got)


# ── the full metric bundle ────────────────────────────────────────────────────

def test_measure_constraint_matches_closed_form_and_costs_four_evals():
    f = _parabola(7.5, 1000.0, floor=500.0)
    m = measure_constraint(f, a_best=7.5, chi2_min=500.0, red_chi2=25.0,
                           a_lo=4.5, a_hi=12.5)
    assert m.sigma_A == pytest.approx(STEP / math.sqrt(1000.0 * STEP ** 2 / 25.0), rel=1e-12)
    assert m.frac_rise_lo == pytest.approx(1000.0 * 9.0 / 500.0, rel=1e-12)
    assert m.frac_rise_hi == pytest.approx(1000.0 * 25.0 / 500.0, rel=1e-12)
    assert m.frac_rise_weaker == pytest.approx(m.frac_rise_lo, rel=1e-12)
    assert m.n_probe_evals == 4, "two curvature probes and two bracket ends"


def test_no_threshold_is_applied_anywhere_in_this_module():
    """🔴 THE FIREWALL (847 item 1 as refined, RYA-161). This module MEASURES. The metric
    shape and the cut are chosen from a sweep across every synthesis band, and a
    threshold appearing here would be a cut chosen from whatever product came first."""
    import inspect

    from pipeline import fit_constraint
    src = inspect.getsource(fit_constraint)
    for banned in ("MAX_SIGMA", "SIGMA_GATE", "CONSTRAINT_GATE", "THRESHOLD"):
        assert banned not in src, f"{banned} — a cut does not belong in the measurer"
    assert not any(f.name in ("passed", "accepted", "verdict")
                   for f in fields(ConstraintMetrics)), "no verdict field"


# ── the plumbing that dropped red_chi2 for its whole life ─────────────────────

def test_asdict_line_emits_every_field_or_declares_why():
    """🔴 THE DEFECT THAT LET THIS HAPPEN, now caught.

    `asdict_line` is a hand-written projection — a SECOND declaration of the per-line
    schema. Adding a field to `LineMeasurement` and forgetting it here produces an
    artifact that is silently missing the column, which is exactly how `red_chi2` never
    reached a CSV. Any new field must be emitted or explicitly declared as omitted.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from derive_band_products import ASDICT_LINE_OMITTED, asdict_line

    from pipeline.band_products import LineMeasurement

    emitted = set(asdict_line(LineMeasurement(
        element="Fe", ion="I", wavelength_air_A=5000.0, instrument="x",
        ew_mA=float("nan"), ew_method="m")).keys())
    declared = {f.name for f in fields(LineMeasurement)}
    missing = declared - emitted - ASDICT_LINE_OMITTED
    assert not missing, (
        f"LineMeasurement fields reach no artifact and are not declared omitted: "
        f"{sorted(missing)}. Add them to asdict_line or to ASDICT_LINE_OMITTED "
        f"with a reason.")


def test_the_constraint_metrics_are_carried_on_the_line():
    """They must exist as fields, and default to None rather than NaN.

    None means "not applicable" — an EW inversion has no chi2 surface, so the question
    does not apply to it. NaN in a CSV column reads like a measurement that failed, which
    would put every EW-route line under permanent suspicion.
    """
    from pipeline.band_products import LineMeasurement
    lm = LineMeasurement(element="Fe", ion="I", wavelength_air_A=5000.0,
                         instrument="x", ew_mA=10.0, ew_method="EW")
    for name in ("sigma_A", "frac_rise_weaker", "edge_distance_dex", "red_chi2"):
        assert hasattr(lm, name)
        assert getattr(lm, name) is None


def test_the_759_harness_forwards_every_metric_the_fitter_returns():
    """The harness sits between `_fit_synth_flux` and the band product, and anything it
    does not forward is invisible downstream — `red_chi2`'s original fate."""
    import inspect
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import rya759_nearuv_fe_product as h

    src = inspect.getsource(h.fit_one)
    for k in ("sigma_A", "frac_rise_weaker", "edge_distance_dex", "red_chi2"):
        assert k in src, f"fit_one drops {k}"


def test_cno_uses_the_shared_definition_not_its_own():
    """RYA-847 hoisted `curvature_sigma`; three copies of an accept/reject rule is how
    the three paths drifted into disagreeing."""
    from pipeline import cno_synthesis
    assert cno_synthesis.curvature_sigma is curvature_sigma
    assert cno_synthesis.CURVATURE_PROBE_STEP_DEX == STEP
    import inspect
    src = inspect.getsource(cno_synthesis)
    assert "def curvature_sigma" not in src, "cno_synthesis redefines it again"


def test_the_handler_has_one_line_constructor_that_owns_the_ew_optout():
    """RYA-847 — structural, not remembered.

    Both handler exits used to build a `LineMeasurement` independently and each had to
    remember `ew_inversion=False`; RYA-770's guard checks that by COUNTING literal
    occurrences, which is a proxy that a refactor can satisfy or break without the
    property changing. `_mk_line` now owns the flag, so no call site can forget it, and
    the metrics reach rejected and accepted lines alike.
    """
    import inspect
    import re

    from pipeline.measure import synthesis as syn
    src = inspect.getsource(syn)
    mk = src[src.index("def _mk_line("):]
    mk = mk[:mk.index("\n        def quarantine(")]
    assert "ew_inversion=False" in mk, "_mk_line must own the EW opt-out"
    for metric in ("sigma_A", "frac_rise_weaker", "edge_distance_dex", "red_chi2"):
        assert metric in mk, f"_mk_line must carry {metric}"
    # no exit may hand ew_inversion in itself any more -- that would be a second owner
    assert not re.search(r"_mk_line\([^)]*ew_inversion", src, re.S)


def test_every_post_fit_rejection_states_the_value_the_fit_returned():
    """RYA-847 item 6 — the reason string is the ONLY place the value survives.

    `quarantine` builds its line through `_mk_line`, which sets no abundance, so EVERY
    line this handler rejects lands in the artifact with `abundance` blank. That is fine
    for a fit that never converged and wrong for one that did: the reader of a
    NON-MINIMUM exclusion has to be able to check that a plausible-looking number was cut
    on the shape of chi2 rather than on being implausible — 3617.318 fitted at a
    perfectly solar 7.477, and no plausibility check could ever have found it.

    CHI2-GATE and FIT-EDGE-PINNED already stated their value; the constraint gate did not,
    so the VIS artifact recorded an exclusion with no number anywhere on the row. This
    asserts the PROPERTY over the region rather than one reason's wording: every
    rejection that happens after the fit has produced `a_best` must report it.
    """
    import inspect

    from pipeline.measure import synthesis as syn
    src = inspect.getsource(syn)
    start = src.index("# FIT-QUALITY GATE (RYA-713)")
    end = src.index("synth_line_accepted", start)
    region = src[start:end]

    calls = [i for i in range(len(region)) if region.startswith("return quarantine(", i)]
    assert len(calls) >= 3, "expected the post-fit rejections to still live here"
    for i in calls:
        stmt = region[i:region.index("\n", region.index(")", i))]
        assert "a_best" in stmt, (
            "a rejection after a converged fit must state the value it returned — "
            "`quarantine` blanks `abundance`, so this string is the only record: " + stmt)


# ── the appendix must name what was removed (RYA-844 / 847 items 5 and 7) ──────

def test_the_appendix_names_excluded_lines_with_their_reason():
    """🔴 RYA-844 puts the burden of proof on the appendix, and it was not carrying it.

    The appendix reported band verdicts and never said which lines had been removed to
    reach them, so a reader could see A(Fe) and its bars with no indication that a line
    had been dropped. An exclusion and a tuning are INDISTINGUISHABLE in the output —
    only the stated reason separates them — which is exactly why the reason has to be
    printed where the number is.
    """
    from pipeline.validate_element import excluded_lines_section
    text = "\n".join(excluded_lines_section("Fe"))
    assert "11119.795" in text, "the RYA-847 item 7 exclusion is not named"
    assert "TELLURIC_ADJACENT" in text
    # the physical evidence, not the consequence
    assert "0.872" in text and "11083.4" in text
    # ⚠️ THIS GUARD WAS ORIGINALLY A BLANKET BAN ON THE WORD "scatter", and it fired on
    # a reason that mentions dispersion as evidence AGAINST tuning: 3617.318's removal
    # makes the near-UV scatter WORSE, which is precisely what shows the cut is not
    # self-serving. A word ban cannot tell a justification from a counter-example, so it
    # bans the JUSTIFYING PHRASES instead. Narrower, and it still fails the case it
    # exists for.
    low = text.lower()
    for tuning in ("inflated the scatter", "reduce the scatter", "reduces the scatter",
                   "improve the scatter", "improves the scatter", "tighten the scatter",
                   "to lower the dispersion", "inflated the dispersion"):
        assert tuning not in low, (
            f"an exclusion reason justifies itself by dispersion ({tuning!r}) — that is "
            f"tuning, not a physical cause (RYA-161/844)")


def test_the_appendix_separates_excluded_from_kept_and_flagged():
    """RYA-807: the discriminator is `status`, not `required_treatment`. A row that is
    `owed` is KEPT and FLAGGED — removing it on an undiagnosed cause would be tuning."""
    from pipeline.validate_element import excluded_lines_section
    text = "\n".join(excluded_lines_section("Fe"))
    assert "Excluded from the aggregate" in text
    assert "Kept and flagged" in text
    assert text.index("Excluded from the aggregate") < text.index("Kept and flagged")


def test_a_missing_registry_makes_no_claim_rather_than_claiming_none():
    """An absent input must not render as an empty exclusion set — that would turn a
    missing file into the confident statement "nothing was excluded" (RYA-833)."""
    import pipeline.validate_element as ve
    from pathlib import Path
    orig = ve._PROBLEM_CHILDREN
    try:
        ve._PROBLEM_CHILDREN = Path("/nonexistent/problem_children.csv")
        text = "\n".join(ve.excluded_lines_section("Fe"))
        assert "NO claim" in text and "missing input" in text
    finally:
        ve._PROBLEM_CHILDREN = orig


def test_an_exclusion_without_a_reason_is_called_out():
    """A registry row with an empty `notes` must render loudly, not silently."""
    import pipeline.validate_element as ve
    import inspect
    src = inspect.getsource(ve.excluded_lines_section)
    assert "NO REASON RECORDED" in src


# ── the decider (pipeline/constraint_gate.py) ─────────────────────────────────

def _metrics(**kw):
    from pipeline.fit_constraint import ConstraintMetrics
    base = dict(sigma_A=0.05, frac_rise_lo=10.0, frac_rise_hi=20.0,
                frac_rise_weaker=10.0, edge_distance_dex=3.0, red_chi2=50.0,
                n_probe_evals=4)
    base.update(kw)
    return ConstraintMetrics(**base)


def test_no_metric_threshold_is_applied_and_the_gate_says_so():
    """🔴 UPDATED WHEN THE DECISION CHANGED, not because it was inconvenient.

    This test used to assert the gate "gates nothing" and announces itself as NOT
    RATIFIED. That was the honest description while the sweep was pending. The sweep has
    since SETTLED it: no metric threshold is defensible, so none is applied — but the
    non-minimum correctness check IS. Both halves must be stated, or the product cannot
    tell a deliberate absence from a missing feature (RYA-843's actual defect).
    """
    from pipeline import constraint_gate as g
    # a line with an appalling sigma_A but a real minimum still passes: no metric cut
    assert g.verdict(_metrics(sigma_A=999.0, frac_rise_weaker=5.0)).ok is True
    d = g.describe()
    assert "NON-MINIMUM check is applied" in d
    assert "NO METRIC THRESHOLD IS APPLIED" in d


def test_a_cut_is_never_defaulted_from_the_chi2_gate():
    """SYNTH_CHI2_GATE exists and is tempting and does NOT transfer. Inheriting it would
    refuse a good NIR line at red_chi2 = 72 and all 40 near-UV lines."""
    from config import constants
    assert getattr(constants, "SYNTH_CONSTRAINT", "missing") is None
    assert constants.SYNTH_CHI2_GATE == 10.0, "unchanged; a different question"


def test_a_ratified_cut_excludes_and_explains(monkeypatch):
    from config import constants
    from pipeline import constraint_gate as g
    monkeypatch.setattr(constants, "SYNTH_CONSTRAINT", ("sigma_A", 0.5), raising=False)
    ok = g.verdict(_metrics(sigma_A=0.05))
    bad = g.verdict(_metrics(sigma_A=0.90))
    assert ok.ok and ok.ratified and ok.reason == ""
    assert not bad.ok and bad.ratified
    assert "UNCONSTRAINED FIT" in bad.reason and "0.9" in bad.reason
    assert "RYA-711" in bad.reason, "a quarantined line is retained, never dropped"
    assert "NOT RATIFIED" not in g.describe()


def test_not_measurable_is_the_STRONGEST_failure_not_a_pass(monkeypatch):
    """🔴 NaN sigma_A means the curvature could not be inverted — railed, or flat. Letting
    it through as 'no data' would admit exactly the fits this gate exists to catch."""
    from config import constants
    from pipeline import constraint_gate as g
    monkeypatch.setattr(constants, "SYNTH_CONSTRAINT", ("sigma_A", 0.5), raising=False)
    for v in (float("nan"), None):
        r = g.verdict(_metrics(sigma_A=v) if v is not None else {"sigma_A": None})
        assert not r.ok, f"{v!r} must fail the gate"
        assert "not measurable" in r.reason


def test_a_lower_is_worse_metric_is_applied_in_the_right_direction(monkeypatch):
    """frac_rise fails LOW (flat chi2), sigma_A fails HIGH. Getting this backwards would
    exclude precisely the well-constrained lines."""
    from config import constants
    from pipeline import constraint_gate as g
    monkeypatch.setattr(constants, "SYNTH_CONSTRAINT",
                        ("frac_rise_weaker", 1.0), raising=False)
    assert g.verdict(_metrics(frac_rise_weaker=10.0)).ok
    assert not g.verdict(_metrics(frac_rise_weaker=0.02)).ok


def test_an_unknown_metric_fails_loudly_rather_than_passing_everything(monkeypatch):
    from config import constants
    from pipeline import constraint_gate as g
    monkeypatch.setattr(constants, "SYNTH_CONSTRAINT", ("wibble", 1.0), raising=False)
    with pytest.raises(g.ConstraintGateError):
        g.verdict(_metrics())


def test_measure_and_decide_stay_in_separate_modules():
    """The firewall. A threshold in the measurer is how the cut gets chosen from whatever
    product is in front of whoever is editing (RYA-161)."""
    import inspect
    from pipeline import constraint_gate, fit_constraint
    assert "SYNTH_CONSTRAINT" not in inspect.getsource(fit_constraint)
    assert "SYNTH_CONSTRAINT" in inspect.getsource(constraint_gate)


def test_both_serialisers_of_a_line_carry_the_metrics():
    """🔴 THE SAME OBJECT HAS TWO SERIALISERS AND ONLY ONE IS SAFE.

    `Product.to_frame()` uses `dataclasses.asdict`, so it cannot drop a field. The
    band-product route's `asdict_line` is HAND-WRITTEN and drops whatever it does not
    mention — that is the one that lost `red_chi2` for the entire life of the fitter.

    The Engine-B cells (eight of the ten synthesis products) serialise through
    `to_frame`, the near-UV/NIR route through `asdict_line`. If either drops the metrics
    the sweep produces a band that looks like it measured nothing, so both are pinned.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from derive_band_products import asdict_line

    from pipeline.band_products import LineMeasurement, build_product

    lm = LineMeasurement(element="Fe", ion="I", wavelength_air_A=5000.0,
                         instrument="kpno_solar_atlas", ew_mA=40.0,
                         ew_method="SYNTHESIS", abundance=7.5, ew_inversion=False,
                         sigma_A=0.0731, frac_rise_weaker=12.5,
                         edge_distance_dex=3.02, red_chi2=48.6)
    wanted = {"sigma_A": 0.0731, "frac_rise_weaker": 12.5,
              "edge_distance_dex": 3.02, "red_chi2": 48.6}

    row = asdict_line(lm)
    for k, v in wanted.items():
        assert row.get(k) == v, f"asdict_line dropped {k}"

    frame = build_product("Fe", "I", "kpno_solar_atlas", "VIS", "ENGINE-B",
                          [lm]).to_frame()
    for k, v in wanted.items():
        assert k in frame.columns, f"Product.to_frame dropped {k}"
        assert float(frame.iloc[0][k]) == v


# ── the settled gate: a correctness check, and NO threshold (RYA-847 sweep) ────

def test_the_non_minimum_check_fires_with_no_cut_ratified():
    """🔴 THE ADOPTED CRITERION. It is a CORRECTNESS check, so it must not depend on a
    ratified threshold — with SYNTH_CONSTRAINT = None the gate still refuses a fit whose
    chi2 never rose away from its reported minimum."""
    from config import constants
    from pipeline import constraint_gate as g
    assert getattr(constants, "SYNTH_CONSTRAINT", "missing") is None
    v = g.verdict(_metrics(frac_rise_weaker=-0.117))
    assert not v.ok
    assert "NON-MINIMUM" in v.reason and "-0.117" in v.reason
    assert "RYA-711" in v.reason, "quarantined, never dropped"


def test_the_sign_test_uses_le_zero_not_lt_zero():
    """A perfectly FLAT objective rises by exactly zero and is just as unconstrained as
    one that dips. The sweep found no line tied at zero, so the two agree on today's
    data — but 'chi2 did not rise' is the semantics that stays correct when a tie
    appears, and a `< 0` implementation would silently pass a flat fit."""
    from pipeline import constraint_gate as g
    assert not g.verdict(_metrics(frac_rise_weaker=0.0)).ok
    assert not g.verdict(_metrics(frac_rise_weaker=-1e-9)).ok
    assert g.verdict(_metrics(frac_rise_weaker=1e-9)).ok


def test_a_real_minimum_passes_however_bad_its_value_or_its_chi2():
    """🔴 THE FIREWALL. The gate must NOT reject on implausibility or on model adequacy.

    Fe I 11593.588 fits at A = 10.559 with red_chi2 = 1226 and a GENUINE minimum. It is a
    non-measurement too — but the sweep proved no defensible red_chi2 cut exists
    (per-band spread x489), so auto-gating it would mean inventing the threshold the
    sweep refuted. It goes to problem_children under RYA-844 instead, by a human, with a
    stated reason.
    """
    from pipeline import constraint_gate as g
    absurd = _metrics(frac_rise_weaker=0.08, sigma_A=0.267, red_chi2=1226.0)
    assert g.verdict(absurd).ok, (
        "the gate rejected a real minimum on its value or its chi2 — that is the "
        "threshold the RYA-847 sweep refuted, re-entering through the back door")


def test_no_threshold_is_ratified_and_the_reason_is_recorded():
    """`SYNTH_CONSTRAINT = None` is the SETTLED answer, not a pending one. If someone
    sets it later it must be earned the same way, so the evidence lives in the code."""
    import inspect
    from pipeline import constraint_gate as g
    src = inspect.getsource(g)
    for evidence in ("x39.5", "x489.4", "581", "continuous"):
        assert evidence in src, f"the sweep evidence for None is not recorded ({evidence})"
    d = g.describe()
    assert "NON-MINIMUM check is applied" in d
    assert "NO METRIC THRESHOLD IS APPLIED" in d
    assert "measured conclusion" in d


def test_the_excluded_lines_are_reproducible_from_the_registry_alone():
    """RYA-844's skeptic test: both RYA-847 exclusions must be recoverable from their
    stated reasons, and neither reason may appeal to the dispersion it produced."""
    from pipeline.validate_element import excluded_lines_section
    t = "\n".join(excluded_lines_section("Fe"))
    assert "3617.318" in t and "NON_MINIMUM" in t
    assert "11119.795" in t and "TELLURIC_ADJACENT" in t
    # 3617.318's reason may mention that the scatter got WORSE -- that is evidence
    # against tuning -- but must not justify the cut by tidiness.
    assert "inflated the scatter" not in t.lower()
