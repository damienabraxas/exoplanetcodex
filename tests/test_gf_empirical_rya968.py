"""
RYA-968 — the empirical gf-grade framework.

The firewall is the point of this file. RYA-161 says grade on a line's own PRECISION, never on
whether it agrees with an expected answer, and "we were careful" is not a mechanism. Three of
these tests are the mechanism:

  F1  the module cannot import a reference abundance (checked by reading the source);
  F2  stages 1-5 are INVARIANT under a constant offset (checked by sweeping one);
  F4  the Cr canary is unmoved.

F3 is the absence of defaults: every threshold raises until declared, and a test asserts that
rather than trusting the dataclass.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

from pipeline import gf_empirical as G  # noqa: E402
from pipeline.error_budget import empirical_gf_term, zero_point_term  # noqa: E402


def _th(**over):
    """Fully-declared thresholds for tests. Values are TEST fixtures, not project defaults."""
    t = G.Thresholds(rew_min=-6.0, rew_max=-4.8, min_depth=0.03, max_red_chi2=5.0,
                     admission_tol_dex=0.20, min_anchor_lines=5,
                     consistency_bound_dex=0.30, fallback_sigma_dex=0.50)
    for k, v in over.items():
        setattr(t, k, v)
    return t


def _row(**over):
    # RYA-1041 step zero: the per-line record now carries its ROUTE. These fixtures were
    # written against the EW semantics (stage 4's red_chi2 cap), so `ew` is what they have
    # always meant -- the field makes that explicit rather than changing it. A synth line
    # is graded on RYA-992's per-arm frac_rise instead and never meets the cap.
    r = dict(wavelength_air_A=5500.0, observed_depth=0.30, rew=-5.1, red_chi2=1.0,
             ew_mA=45.0, problem_class="", excluded_reason="", gf_tier="KURUCZ",
             abundance=7.50, gf_sigma_dex=None, cross_measurements=None,
             inferred_sigma_dex=None, route="ew", instrument=None,
             frac_rise_weaker=None)
    r.update(over)
    return r


# ── F3: declared, or nothing ─────────────────────────────────────────────────────────
#: Thresholds still awaiting a DERIVATION. `admission_tol_dex` has left this list: RYA-1041
#: D2 measured the setting input (the anchor's own differential scatter) and it is now set
#: WITH its citation. The rule has not been relaxed -- a value still may not appear without
#: a measurement behind it, which `test_a_set_threshold_carries_its_measurement` enforces
#: for anything that leaves this list.
@pytest.mark.parametrize("name", [
    "rew_min", "rew_max", "min_depth", "max_red_chi2",
    "min_anchor_lines", "consistency_bound_dex", "fallback_sigma_dex"])
def test_every_undelivered_threshold_refuses_until_declared(name):
    with pytest.raises(G.ThresholdNotDeclared, match="has not been declared"):
        G.Thresholds().require(name)


def test_a_set_threshold_carries_its_measurement():
    """🔴 THE RULE THAT REPLACES raise-until-set FOR A DELIVERED THRESHOLD.

    Removing a field from the refusing list is only safe if something else stops an
    unsourced constant appearing. A number with no provenance is exactly the blanket 0.17
    this framework exists to retire, wearing a smaller value.

    Every threshold that is SET must name, in the source beside it, the ticket and the
    measurement it came from -- so `git blame` is not the only record of why it is that
    number.
    """
    import inspect
    src = inspect.getsource(G.Thresholds)
    t = G.Thresholds()
    for name in ("rew_min", "rew_max", "min_depth", "max_red_chi2", "admission_tol_dex",
                 "min_anchor_lines", "consistency_bound_dex", "fallback_sigma_dex"):
        if getattr(t, name) is None:
            continue                                    # still undelivered, covered above
        # walk BACK from the field's own line, collecting the `#:` block above it
        lines = src.split("\n")
        i = next(k for k, l in enumerate(lines) if l.strip().startswith(f"{name}:"))
        block = []
        for l in reversed(lines[:i]):
            st = l.strip()
            if st.startswith("#:") or st.startswith("#"):
                block.append(st)
            elif not st:
                continue
            else:
                break
        text = " ".join(block)
        assert "RYA-" in text, f"{name} is set but names no ticket"
        assert any(w in text.lower() for w in ("measured", "diagnostic", "derived")), (
            f"{name} is set but does not say it came from a measurement")


def test_admission_tol_is_the_differential_not_the_solo_scatter():
    """RYA-898's distinction, pinned as a number. Solo scatter on the same pool is
    0.1290 dex -- an order of magnitude wider -- and sizing admission from it would admit
    lines on a spread the differential never pays."""
    t = G.Thresholds()
    assert t.admission_tol_dex == 0.0119, "the RYA-1041 D2 same-instrument differential"
    assert t.admission_tol_dex < 0.05, (
        "0.05 is GBS's borrowed tolerance, which admits only 28.3% of our own anchor")


def test_the_refusal_explains_why_a_borrowed_constant_is_not_a_default():
    """§3.1 measured this: the GBS window admits 100% of our pool."""
    with pytest.raises(G.ThresholdNotDeclared) as e:
        G.Thresholds().require("rew_min")
    assert "100%" in str(e.value) and "no default" in str(e.value)


# ── F1: the reference abundance is not reachable ─────────────────────────────────────
def test_the_module_cannot_reach_a_reference_abundance():
    src = (ROOT / "pipeline" / "gf_empirical.py").read_text(encoding="utf-8")
    code = "\n".join(l for l in src.splitlines()
                     if not l.lstrip().startswith("#"))
    banned = ("gold_reference", "solar_gold", "STAR_PARAMS", "asplund",
              "reference_abundance", "phase_c_verdict", "get_gold")
    for b in banned:
        assert b.lower() not in code.lower(), (
            f"gf_empirical references {b!r}. A grading function that can see the expected "
            f"answer can grade toward it — RYA-161 F1.")
    assert "import" in code


def test_the_physical_stages_are_not_even_given_an_abundance():
    """F1 by signature: stages 1-5 take no abundance argument at all."""
    import inspect
    for fn in (G.stage1_depth, G.stage2_saturation, G.stage3_purity, G.stage4_fit, G.stage5_gf):
        params = set(inspect.signature(fn).parameters)
        assert not (params & {"abundance", "anchor", "reference"}), (
            f"{fn.__name__} accepts an abundance; it must not be able to see one")


# ── F2: offset invariance — the checkable form of "precision, not accuracy" ──────────
@pytest.mark.parametrize("delta", [-0.5, -0.2, -0.05, 0.05, 0.2, 0.5])
def test_stages_1_to_5_are_invariant_under_a_constant_offset(delta):
    th = _th()
    rows = [_row(wavelength_air_A=5000 + 10 * i, abundance=7.4 + 0.05 * (i % 7),
                 rew=-5.6 + 0.05 * (i % 9), observed_depth=0.05 + 0.02 * (i % 8),
                 red_chi2=0.5 + 0.4 * (i % 6),
                 gf_tier=["LAB", "NIST-C+", "KURUCZ", "VALD3"][i % 4]) for i in range(40)]
    base = [G.grade_line(r, th) for r in rows]
    moved = [G.grade_line({**r, "abundance": r["abundance"] + delta}, th) for r in rows]
    for a, b in zip(base, moved):
        assert a.stages == b.stages, "a physical stage moved when every abundance shifted"
        assert a.tier == b.tier
        assert a.reason_code == b.reason_code
        assert a.sigma_dex == b.sigma_dex


def test_stage6_is_deliberately_NOT_offset_invariant_and_says_so():
    """It is the one stage that compares abundances — that is what it is for."""
    th = _th()
    anchor = G.build_anchor([7.50, 7.52, 7.47, 7.55, 7.49, 7.51], th)
    v1, d1 = G.stage6_anchor_consistency(7.50, anchor, th)
    v2, d2 = G.stage6_anchor_consistency(7.50 + 1.0, anchor, th)
    assert v1 == G.YES and v2 == G.NO
    assert d1 != d2, "the distance must move; it is the audit for the non-invariant stage"


def test_stage6_returns_its_distance_so_the_one_uncheckable_stage_is_auditable():
    th = _th()
    anchor = G.build_anchor([7.50, 7.52, 7.47, 7.55, 7.49], th)
    pool = G.grade_pool([_row(abundance=7.9)], th, anchor)
    assert pool.verdicts[0].anchor_distance_dex is not None
    assert any(r["anchor_distance_dex"] is not None for r in G.audit_table(pool))


# ── the anchor and the zero-point ────────────────────────────────────────────────────
def test_a_thin_anchor_is_refused_rather_than_used():
    th = _th(min_anchor_lines=12)
    with pytest.raises(ValueError, match="below the declared minimum"):
        G.build_anchor([7.5, 7.4, 7.6], th)


def test_the_refusal_forbids_lowering_the_minimum_to_pass():
    th = _th(min_anchor_lines=12)
    with pytest.raises(ValueError) as e:
        G.build_anchor([7.5] * 3, th)
    assert "do NOT lower the minimum" in str(e.value)


def test_the_zero_point_does_not_shrink_when_ungraded_lines_are_admitted():
    """🔴 The central constraint. Admission buys statistical precision, never scale."""
    th = _th()
    lab = [7.50, 7.52, 7.47, 7.55, 7.49, 7.51, 7.45]
    anchor = G.build_anchor(lab, th)
    thin = G.grade_pool([_row(abundance=7.50, wavelength_air_A=5000 + i) for i in range(3)],
                        th, anchor)
    fat = G.grade_pool([_row(abundance=7.50, wavelength_air_A=5000 + i) for i in range(300)],
                       th, anchor)
    assert thin.zero_point_dex == fat.zero_point_dex, (
        "the zero-point moved when more ungraded lines were admitted — they carry no "
        "independent information about the scale")
    assert fat.n_in_value > thin.n_in_value


def test_the_zero_point_matches_the_measured_vis_cell():
    """RYA-959 VIS: 7 laboratory lines, sd 0.157 -> 0.059 dex."""
    th = _th()
    rng = np.random.default_rng(968)
    a = 7.498 + 0.157 * (rng.standard_normal(7) / np.std(rng.standard_normal(7), ddof=1))
    anchor = G.build_anchor(list(a), th)
    assert anchor.n == 7
    assert anchor.zero_point_dex == pytest.approx(anchor.sd / np.sqrt(7))


# ── tiers, and the boundary that keeps quarantine from becoming delete ───────────────
def test_a_physical_failure_is_invalid_with_a_named_cause():
    th = _th()
    v = G.grade_line(_row(observed_depth=0.001), th)
    assert v.tier == G.TIER_INVALID
    assert v.reason_code == G.RC_CONTINUUM
    assert not v.in_value


def test_scattered_is_never_a_reason_for_invalid():
    """The boundary: a wide line is documented, not called an artifact."""
    th = _th()
    anchor = G.build_anchor([7.50, 7.52, 7.47, 7.55, 7.49], th)
    v = G.grade_line(_row(abundance=9.0), th, anchor)
    assert v.tier == G.TIER_SCATTERED
    assert v.tier != G.TIER_INVALID
    assert not v.in_value, "a scattered line must not enter the value"


def test_a_lab_line_is_golden_and_keeps_its_pedigree():
    th = _th()
    v = G.grade_line(_row(gf_tier="LAB", gf_sigma_dex=0.03), th)
    assert v.tier == G.TIER_GOLDEN and v.sigma_source == "cited" and v.in_value


def test_every_reason_code_is_from_the_existing_registry_vocabulary():
    """RYA-809's taxonomy, not a parallel one."""
    import csv
    reg = ROOT / "data" / "registry" / "problem_children.csv"
    with open(reg, encoding="utf-8") as fh:
        known = {r["problem_class"].strip() for r in csv.DictReader(fh)}
    known |= {"BAD_FIT", ""}
    ours = {getattr(G, n) for n in dir(G) if n.startswith("RC_")}
    assert ours <= known, f"reason codes not in the registry vocabulary: {ours - known}"


# ── the pool: a mixed pool no longer collapses to its worst member ───────────────────
def test_an_admitted_line_with_no_measured_sigma_carries_the_fallback_and_still_counts():
    """The fallback is ignorance, not a measurement of width — it must not exclude."""
    th = _th()
    anchor = G.build_anchor([7.50, 7.52, 7.47, 7.55, 7.49], th)
    v = G.grade_line(_row(gf_tier="KURUCZ", abundance=7.50), th, anchor)
    assert v.tier == G.TIER_CONSISTENT and v.in_value
    assert v.sigma_source == "fallback" and v.sigma_dex == th.fallback_sigma_dex
    assert "no measured width" in v.note


def test_a_measured_sigma_above_the_bound_is_scattered():
    th = _th()
    anchor = G.build_anchor([7.50, 7.52, 7.47, 7.55, 7.49], th)
    v = G.grade_line(_row(gf_tier="KURUCZ", abundance=7.50,
                          cross_measurements=[7.1, 7.5, 7.9]), th, anchor)
    assert v.sigma_source.startswith("self-reported")
    assert v.tier == G.TIER_SCATTERED and v.reason_code == G.RC_OUTLIER


def test_a_mixed_pool_is_a_quadrature_not_a_collapse():
    """🔴 The reason RYA-855 moved 0 of 36 bars."""
    th = _th()
    anchor = G.build_anchor([7.50, 7.52, 7.47, 7.55, 7.49, 7.51, 7.45], th)
    rows = ([_row(gf_tier="LAB", gf_sigma_dex=0.03, abundance=7.50,
                  wavelength_air_A=5000 + i) for i in range(7)]
            + [_row(gf_tier="KURUCZ", abundance=7.50, wavelength_air_A=6000 + i)
               for i in range(50)])
    pool = G.grade_pool(rows, th, anchor)
    assert pool.tier_counts[G.TIER_GOLDEN] == 7
    assert pool.gf_sigma_dex < th.fallback_sigma_dex, (
        "the laboratory lines surrendered their pedigree to their neighbours — this is exactly "
        "the all-or-nothing behaviour this framework replaces")
    assert pool.gf_sigma_dex > 0.03


def test_sigma_precedence_is_cited_then_self_then_inferred_then_fallback():
    th = _th()
    assert G.resolve_sigma(cited_sigma_dex=0.03, cross_measurement_abundances=[1, 2, 3],
                           inferred_sigma_dex=0.4, th=th)[1] == "cited"
    s, src = G.resolve_sigma(cross_measurement_abundances=[7.5, 7.6, 7.4],
                             inferred_sigma_dex=0.4, th=th)
    assert src.startswith("self-reported")
    assert G.resolve_sigma(inferred_sigma_dex=0.4, th=th)[1] == "inferred"
    assert G.resolve_sigma(th=th) == (0.50, "fallback")


def test_self_reported_needs_three_independent_measurements():
    th = _th()
    assert G.resolve_sigma(cross_measurement_abundances=[7.5, 7.6], th=th)[1] == "fallback"


# ── budget wiring ────────────────────────────────────────────────────────────────────
def test_the_budget_terms_refuse_an_unsourced_or_impossible_value():
    with pytest.raises(ValueError):
        empirical_gf_term(0.0, n_lines=10, provenance="x")
    with pytest.raises(ValueError, match="never 'assumed'"):
        empirical_gf_term(0.2, n_lines=10, provenance="")
    with pytest.raises(ValueError):
        zero_point_term(-1.0, n_anchor=7)


def test_neither_budget_term_averages_down():
    a = empirical_gf_term(0.2, n_lines=10, provenance="cited x7, fallback x3")
    z = zero_point_term(0.059, n_anchor=7)
    assert a.averages_down is False and z.averages_down is False
    assert "RYA-968" in a.source and "cannot improve it" in z.source


# ── F4: the Cr canary ────────────────────────────────────────────────────────────────
def test_the_cr_canary_is_not_something_this_framework_can_see():
    """+0.402 is a standing residual. Grading must not be able to move it, because grading
    never sees an expected answer — so the canary's value cannot appear in this module."""
    # Behaviour, not a grep: the docstring names the canary (that is documentation), so the
    # real check is that shifting an abundance BY the canary changes nothing the framework
    # decides — because nothing it decides can see an abundance except stage 6.
    th = _th()
    anchor = G.build_anchor([7.50, 7.52, 7.47, 7.55, 7.49], th)
    base = G.grade_pool([_row(abundance=7.5)], th, anchor)
    shifted = G.grade_pool([_row(abundance=7.5 + 0.402)], th, anchor)
    assert base.verdicts[0].stages["depth"] == shifted.verdicts[0].stages["depth"]
    assert base.verdicts[0].sigma_dex == shifted.verdicts[0].sigma_dex
