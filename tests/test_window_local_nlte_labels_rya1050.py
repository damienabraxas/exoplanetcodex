"""RYA-1050 — the NLTE-label guard must ask the WINDOW question, not the element question.

🔴 THE DEFECT. `assert_linelist_supports_nlte` accepts optional wavelength bounds, and all
three call sites omitted them. Without bounds it asks "is this element labelled ANYWHERE
in the list", and the guard's own docstring says that is the form that misses:

    bsyn applies departures per LINE and falls back to departure = 1 for any line whose
    levels are unidentified ... Fe answers the global question yes 15706 times while
    every line in the Fe II 6910-9199 window is unlabelled.

So a fully unlabelled BAND passes, provided the element is labelled somewhere else.

⚠️ AND PARTIAL COVERAGE NEVER RAISES AT ALL -- IT DILUTES. The fallback is per line, so a
half-labelled band emits a well-formed NLTE product carrying roughly half the effect, with
nothing anywhere recording that. That is not something an assert can catch, which is why
the labelled-line COUNT now travels in the product's provenance rather than living only
inside a check.

⚠️ NO COVERAGE FLOOR IS IMPOSED. A "refuse below X% labelled" rule would be a chosen cut,
and RYA-161 says a threshold has to be derived or not applied. Coverage is STATED and the
reader decides.
"""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DRIVER = (ROOT / "scripts" / "derive_band_products.py").read_text()
GUARD = (ROOT / "pipeline" / "gerber_nlte.py").read_text()


def test_every_call_site_passes_the_band_bounds():
    """The whole fix. An element-global call is the version that cannot see the defect."""
    n_calls = DRIVER.count("assert_linelist_supports_nlte(")
    n_bounded = DRIVER.count("wave_lo_A=float(a.lo)")
    assert n_calls >= 3, f"expected at least 3 call sites, found {n_calls}"
    assert n_bounded == n_calls, (
        f"{n_calls} call sites but only {n_bounded} pass the window -- an unbounded call "
        f"asks the question the guard's docstring says is the one that misses")


def test_the_guard_actually_supports_a_window():
    """Pinning the parameter, because the fix is worthless if the signature drops it."""
    assert "wave_lo_A" in GUARD and "wave_hi_A" in GUARD
    # the banner now names BOTH narrowings, since the stage turned out to matter more
    assert "NOT ELEMENT-GLOBAL" in GUARD
    assert "WINDOW-LOCAL AND STAGE-LOCAL" in GUARD


def test_the_guard_narrows_to_the_window_before_counting():
    """It must filter to the band and THEN look for labels, not the reverse.

    The narrowing now lives in `_select_species_rows`, which the assert calls before it
    counts anything -- so the ordering is checked across the two, not inside one body.
    """
    seg = GUARD[GUARD.index("def assert_linelist_supports_nlte"):]
    i_select = seg.index("_select_species_rows(linelist, Z, ion")
    i_count = seg.index('flagged["nlte_label_low"]')
    assert i_select < i_count, "rows must be narrowed before labels are counted"
    # and the selector itself narrows on BOTH axes before returning
    sel = GUARD[GUARD.index("def _select_species_rows("):]
    sel = sel[:sel.index("\n\ndef ")]
    assert sel.index("species == Z") < sel.index("return linelist[sel]")
    assert sel.index("wave_lo_A is not None") < sel.index("return linelist[sel]")


def test_coverage_is_STATED_in_provenance_not_just_asserted():
    """🔴 Partial coverage dilutes silently and no assert can catch it. The count has to
    reach the product, or a half-labelled band looks exactly like a fully labelled one."""
    # `pool_label_coverage(` contains `label_coverage(`, so count the window calls by
    # their distinguishing prefix rather than by the shared substring.
    n_window = (DRIVER.count("_gn.label_coverage(")
                + DRIVER.count("gnlte.label_coverage("))
    assert n_window == 3, (
        f"every NLTE leg must measure its own window coverage (found {n_window})")
    assert DRIVER.count("_b_source += _cov_txt") + DRIVER.count("eb_prov += _cov_txt") \
        + DRIVER.count("+ _cov_txt)") == 3, (
        "and every one must put it in the PRODUCT's provenance, not just print it")


def test_coverage_is_reported_as_a_FRACTION_not_a_bare_count():
    """🔴 9,048 looks reassuring until it is set beside 18,366 — which is the real number
    for Fe in our own 4200-6910 band, i.e. 49% coverage. A bare count hides that."""
    GUARD_SRC = (ROOT / "pipeline" / "gerber_nlte.py").read_text()
    assert "def label_coverage(" in GUARD_SRC
    assert "of {_cov_tot}" in DRIVER, "the TOTAL must appear beside the labelled count"
    assert "total not measurable" in DRIVER, (
        "0/0 must be distinguishable from zero coverage")


def test_no_coverage_FLOOR_is_imposed():
    """A 'refuse below X%' rule would be a chosen cut. RYA-161: derive it or omit it."""
    for bad in ("coverage < 0.5", "coverage_floor", "MIN_LABEL_FRACTION"):
        assert bad not in DRIVER, f"{bad!r} looks like a picked threshold"


def test_the_LTE_comparand_is_never_asked_for_NLTE_labels():
    """The check belongs only to the NLTE member. Demanding labels of a run that
    deliberately applies none would refuse the comparand for lacking what it does not use."""
    i = DRIVER.index("_gn.assert_linelist_supports_nlte(")
    # search back to the nearest enclosing guard; the RYA-1050 comment block sits between
    # the guard and the call, so a short window would miss it and test nothing.
    before = DRIVER[max(0, i - 2000):i]
    assert "if _nlte:" in before, "the synthesis_route call must sit under `if _nlte:`"
    # and nothing re-opens the branch between the guard and the call
    assert "else:" not in before[before.rindex("if _nlte:"):]


# ── the ionisation stage, which is where the real damage was ────────────────

def test_the_selection_is_ION_AWARE_not_just_element_aware():
    """🔴 THE BIGGER HALF OF THIS DEFECT. `turbospectrum_species` is the ONE field in the
    list that collapses the stages -- every Fe row carries 26.0 -- and both callers were
    built on exactly that field. Measured on GESv6 420-920 nm:

        4200-6910 A   Fe 1   8194 of 8243   99.4%
        4200-6910 A   Fe 2    854 of 8870    9.6%
                      ALL Fe 9048 of 18366  49.3%   <- what a Z-only filter sees

    Fe I is essentially COMPLETE, Fe II essentially ABSENT, and the 49% is an average of
    two populations with nothing to do with each other. An Fe II NLTE product would pass a
    Z-only check ON Fe I's LABELS and synthesise ~90% of its own lines in LTE under an
    NLTE label.
    """
    assert "def _select_species_rows(" in GUARD
    assert 'if ion is not None and "element" in names:' in GUARD, (
        "the stage must come from a field that CARRIES it -- `element` is 'Fe 1'/'Fe 2'")
    assert "parse_ion" in GUARD, "use the existing helper, do not add a sixth roman map"


def test_both_the_assert_and_the_coverage_report_share_one_selection():
    """If they filtered separately they could disagree, and the number asserted on would
    not be the number reported."""
    # three consumers now: the assert, the window report, and the pool report -- all
    # narrowing through the same function so their numbers cannot disagree.
    assert GUARD.count("_select_species_rows(linelist, Z, ion") == 3


def test_every_call_site_passes_the_ion():
    n_calls = DRIVER.count("assert_linelist_supports_nlte(") + DRIVER.count("label_coverage(")
    assert DRIVER.count("ion=a.ion)") >= n_calls, (
        "every guard and coverage call must name the ionisation stage it means")


def test_ion_None_still_means_element_global_on_purpose():
    """Not every caller wants a stage. The old behaviour must remain REACHABLE, just not
    the default thing a band run silently gets."""
    assert "`ion=None` keeps the old element-global behaviour deliberately" in GUARD


# ── the POOL, which is what actually gets synthesised ───────────────────────

def test_the_POOL_is_checked_not_only_the_window():
    """🔴 THE THIRD AND SHARPEST NARROWING. A run does not synthesise "the window", it
    synthesises a POOL, and the two can disagree completely. Measured:

        Fe II 4200-6910 window        : 854 of 8870 labelled (9.6%)  -> window check PASSES
        RYA-877's measured Fe II pool :   0 of   11 labelled (0.0%)  -> every line is LTE

    So an Fe II NLTE run passes a window check on 854 lines it will never touch, then
    synthesises all eleven it does touch in LTE under an NLTE label.
    """
    assert "def pool_label_coverage(" in GUARD
    n_pool = (DRIVER.count("_gn.pool_label_coverage(")
              + DRIVER.count("gnlte.pool_label_coverage("))
    assert n_pool == 3, f"every NLTE leg must check its pool (found {n_pool})"


def test_a_pool_with_ZERO_labelled_lines_REFUSES():
    """Dilution is reported; total absence is refused. A product where NOT ONE fitted line
    carries a label is not a weakened NLTE result, it is an LTE result wearing an NLTE
    name."""
    assert "_pool_tot and _pool_lab == 0" in DRIVER
    assert DRIVER.count("_pool_tot and _pool_lab == 0") == 3
    assert "LTE under" in DRIVER and "an NLTE label" in DRIVER


def test_the_refusal_explains_why_the_window_check_passed():
    """Otherwise the next reader sees a window count of 854 and a refusal and concludes
    the guard is broken."""
    i = DRIVER.index("_pool_tot and _pool_lab == 0")
    seg = DRIVER[i:i + 1400]
    assert "they are not the lines being fitted" in seg


def test_the_window_number_is_still_reported_because_it_answers_BLENDS():
    """The pool check does not replace the window check -- unlabelled neighbours inside
    the fitting window contribute to the profile in LTE whether or not they are targets."""
    assert "label_coverage(" in DRIVER and "pool_label_coverage(" in DRIVER
    assert "Both numbers are real and they answer different" in GUARD


# ── the FOURTH call site, found 2026-08-30 ──────────────────────────────────
#
# 🔴 The three sites above are all in the DRIVER. There is a fourth, inside the fit
# itself, and it could not be fixed the same way: `_fit_synth_flux` had no `ion`
# parameter, so the caller's stage died at its signature. `SynthesisHandler.measure_line`
# had the stage all along and dropped it there. That path is reached with ion="II" from
# OUTSIDE the driver (scripts/rya783_nlte_engagement_probe.py), where none of the three
# driver guards run -- so on that path the weakest form of the check was the only form.

import ast
import inspect

import numpy as np

FIT_SRC = (ROOT / "pipeline" / "abundances_derive.py").read_text()
HANDLER_SRC = (ROOT / "pipeline" / "measure" / "synthesis.py").read_text()


def _calls_named(src, func_name):
    """Every Call node whose callee ends in `func_name` — AST, not grep.

    Read as CODE: a docstring or comment naming the call must not satisfy a test that
    the call is MADE, and a `grep` for 'ion=' cannot tell an argument from prose.
    """
    out = []
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Call):
            f = node.func
            name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
            if name == func_name:
                out.append(node)
    return out


def test_the_fit_can_be_told_which_stage_it_is_fitting():
    """The fix is worthless if the signature drops it -- this is where it died before."""
    from pipeline.abundances_derive import _fit_synth_flux
    sig = inspect.signature(_fit_synth_flux)
    assert "ion" in sig.parameters, (
        "_fit_synth_flux takes no `ion`, so the guard inside it can only ask the "
        "element-global question no matter what the caller knows")
    assert sig.parameters["ion"].default is None, (
        "`ion=None` must stay the default so every existing caller is unchanged")


def test_the_fourth_call_site_passes_the_ion():
    calls = _calls_named(FIT_SRC, "assert_linelist_supports_nlte")
    assert calls, "the guard call inside the fit has gone missing"
    for c in calls:
        kw = {k.arg for k in c.keywords}
        assert "ion" in kw, (
            "the guard call inside `_fit_synth_flux` does not name the stage. "
            "`turbospectrum_species` collapses Fe I and Fe II onto 26.0, so a Z-only "
            "check clears an Fe II synthesis on Fe I's labels -- an LTE spectrum under "
            "an NLTE label, which is the outcome this guard exists to refuse")


def test_the_handler_forwards_the_stage_it_already_has():
    """`measure_line(..., ion=...)` knows the answer; the fit one layer down did not."""
    calls = [c for c in _calls_named(HANDLER_SRC, "_fit")]
    assert calls, "SynthesisHandler no longer delegates to the validated fit"
    for c in calls:
        kw = {k.arg: k.value for k in c.keywords}
        assert "ion" in kw, "the handler drops `ion` on the way into the fit"
        assert isinstance(kw["ion"], ast.Name) and kw["ion"].id == "ion", (
            "forward the `ion` the handler was CALLED with, do not re-derive it")


def _two_stage_linelist():
    """Fe I labelled, Fe II unlabelled, both in one window — the measured GESv6 shape.

    Real numbers this mimics: 4200-6910 A, Fe 1 8194/8243 labelled (99.4%), Fe 2
    854/8870 (9.6%). A Z-only filter averages them to 49.3% and sees nothing wrong.
    """
    n = 8
    ll = np.zeros(n, dtype=[("turbospectrum_species", "f8"), ("element", "U8"),
                            ("wave_A", "f8"), ("nlte", "U2"),
                            ("nlte_label_low", "U16"), ("nlte_label_up", "U16")])
    ll["turbospectrum_species"] = 26.0            # collapsed: identical for BOTH stages
    ll["wave_A"] = np.linspace(4300.0, 6800.0, n)
    ll["nlte"] = "T"
    ll["element"][:4] = "Fe 1"
    ll["nlte_label_low"][:4] = "a5D"
    ll["nlte_label_up"][:4] = "z5F"
    ll["element"][4:] = "Fe 2"
    ll["nlte_label_low"][4:] = "none"             # every Fe II line unidentified
    ll["nlte_label_up"][4:] = "none"
    return ll


def test_a_stage_blind_guard_CLEARS_an_FeII_run_on_FeI_labels():
    """🔴 THE CONTROL. This is what the fourth site permitted, stated as a passing test.

    Without a stage the guard counts Fe I's four labels and returns happily for a run
    that will synthesise four Fe II lines, every one of them in LTE. The assertion here
    is that the OLD behaviour does not raise -- if this ever starts raising, the
    element-global escape hatch has changed meaning and the tests below are no longer
    measuring what they claim.
    """
    from pipeline.gerber_nlte import assert_linelist_supports_nlte
    ll = _two_stage_linelist()
    n = assert_linelist_supports_nlte(ll, 26, "Fe", 4200.0, 6910.0, ion=None)
    assert n == 4, f"the stage-blind count is Fe I's labels, got {n}"


def test_with_the_stage_the_FeII_run_is_REFUSED():
    """🔴 RYA-1135 CHANGED *WHY* THIS REFUSES, NOT WHETHER. The assertion used to pin the
    message "NONE carry NLTE level labels", i.e. the LABEL-COUNT refusal. Fe II is now
    refused one step earlier, on the ATOM's reach: `atom.fe607a` has no Fe II bound-bound
    transitions at all, so label coverage cannot decide the case — a pool with labelled
    Fe II lines would have sailed past the label count and emitted LTE under an NLTE label.

    Pinning the old string would now pass only by accident of ordering, so the assertion
    moves to what actually matters and keeps its teeth: it is refused, AND the reason names
    the atom rather than the list."""
    from pipeline.gerber_nlte import assert_linelist_supports_nlte, GerberDeckError
    ll = _two_stage_linelist()
    with pytest.raises(GerberDeckError, match="NO LINE LIST CAN FIX THIS") as e:
        assert_linelist_supports_nlte(ll, 26, "Fe", 4200.0, 6910.0, ion=2)
    assert "ZERO Fe II bound-bound transitions" in str(e.value)


def test_the_FeII_refusal_does_not_depend_on_the_LIST_at_all():
    """⚠️ CONTROL for the above. If the refusal still keyed on label coverage, handing it a
    list where every Fe II line IS labelled would let it through — which is precisely the
    hole RYA-1135 found (854 of 8,870 Fe II lines in the VIS window carry labels). Passing
    `None` for the linelist is the sharpest form: nothing can be counted, and it still
    refuses, so the verdict cannot be coming from the list."""
    from pipeline.gerber_nlte import assert_linelist_supports_nlte, GerberDeckError
    with pytest.raises(GerberDeckError, match="NO LINE LIST CAN FIX THIS"):
        assert_linelist_supports_nlte(None, 26, "Fe", 4200.0, 6910.0, ion=2)


def test_and_the_FeI_run_on_the_same_list_still_passes():
    """The refusal must be about the STAGE, not about the list being hard to read."""
    from pipeline.gerber_nlte import assert_linelist_supports_nlte
    ll = _two_stage_linelist()
    assert assert_linelist_supports_nlte(ll, 26, "Fe", 4200.0, 6910.0, ion=1) == 4


def test_the_direct_caller_that_runs_no_guard_of_its_own_forwards_it_too():
    """`scripts/rya759_nearuv_fe_product.py` calls `_fit_synth_flux` DIRECTLY and runs no
    NLTE guard itself, so the guard inside the fit is the only one on that path. It has
    no ion axis today (ctx comes from `build_solar_context('Fe', ...)`), which is why
    this is latent rather than live -- but the wiring must exist before an axis is added.
    """
    src = (ROOT / "scripts" / "rya759_nearuv_fe_product.py").read_text()
    calls = _calls_named(src, "_fit_synth_flux")
    assert calls, "rya759 no longer calls the fit"
    for c in calls:
        assert "ion" in {k.arg for k in c.keywords}, (
            "rya759 forwards `nlte_deck` but not the stage, and runs no guard of its own")
