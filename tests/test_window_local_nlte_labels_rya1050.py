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
    assert DRIVER.count("label_coverage(") == 3, (
        "every NLTE leg must measure its own window coverage")
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
    assert GUARD.count("_select_species_rows(linelist, Z, ion") == 2


def test_every_call_site_passes_the_ion():
    n_calls = DRIVER.count("assert_linelist_supports_nlte(") + DRIVER.count("label_coverage(")
    assert DRIVER.count("ion=a.ion)") >= n_calls, (
        "every guard and coverage call must name the ionisation stage it means")


def test_ion_None_still_means_element_global_on_purpose():
    """Not every caller wants a stage. The old behaviour must remain REACHABLE, just not
    the default thing a band run silently gets."""
    assert "`ion=None` keeps the old element-global behaviour deliberately" in GUARD
