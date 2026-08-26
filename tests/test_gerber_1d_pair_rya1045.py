"""RYA-1045 — the 1D rung's paired comparand, and the atmosphere its deck demands.

🔴 TWO DEFECTS, ONE SHAPE: the RYA-1040 paired-comparand discipline was built for the
⟨3D⟩ rung and never carried to the 1D one.

1. `synthesis_route`'s Engine-B leg passed `ctx["atmosphere"]` to the Gerber 1D deck. That
   deck is computed on MARCS.GES (RYA-798) and pairs index-for-index with it and nothing
   else; the route's own atmosphere is ATLAS9.Castelli. `assert_depth_match` refused the
   run — 56 layers against 72 — so `--engine-b-deck gerber-nlte` had never been reachable
   on this route at all. Caught rather than applied, which matters: iSpec overwrites the
   departure tau with the atmosphere's, so a mismatch would have applied every departure
   at the wrong depth SILENTLY.

2. There was no LTE product on that deck's atmosphere, so the only available difference
   (ENGINE-B-NLTE minus the route's 1D-LTE) reported the ATLAS9→MARCS.GES ATMOSPHERE
   change as non-LTE physics — the RYA-542 confound, on the rung nobody checked. Measured
   on solar Fe I: that mixed delta is +0.050, and it is NOT the 1D NLTE effect.
"""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DRIVER = (ROOT / "scripts" / "derive_band_products.py").read_text()

from pipeline import band_products as bp                       # noqa: E402
from pipeline import treatment_axes as taxes                    # noqa: E402


def test_the_1d_deck_gets_its_own_atmosphere_not_the_routes():
    """🔴 DEFECT 1. The deck decides the atmosphere; the route does not get a vote."""
    assert 'model_grid="MARCS.GES")' in DRIVER
    assert "_layers1d" in DRIVER
    assert "gnlte.assert_depth_match(_dep, _layers1d)" in DRIVER, (
        "the depth match must be asserted against the deck's OWN atmosphere")


def test_the_depth_match_is_asserted_before_any_fit():
    """The guard is what caught this. It must stay ahead of the synthesis."""
    i_assert = DRIVER.index("gnlte.assert_depth_match(_dep, _layers1d)")
    i_fit = DRIVER.index("eb_lines = _fit_lines(eb_treatment")
    assert i_assert < i_fit


def test_the_comparand_exists_and_is_paired_to_the_nlte_member():
    """🔴 DEFECT 2. An NLTE product with no same-atmosphere LTE product cannot yield an
    NLTE effect — only an atmosphere-plus-NLTE mixture."""
    tok = taxes.GERBER1D_LTE_MARCS.token
    assert taxes.comparand_for("ENGINE-B-NLTE") == tok


def test_the_comparand_is_on_the_same_atmosphere_as_the_nlte_member():
    """The whole point: ONE atmosphere, departures the only difference (RYA-542)."""
    nlte = taxes.axes_for("ENGINE-B-NLTE", handler="SynthesisHandler")
    lte = taxes.GERBER1D_LTE_MARCS
    assert nlte.atmos == lte.atmos == "marcs-ges"
    assert nlte.model == lte.model == "gerber"
    assert nlte.is_nlte is True and lte.is_nlte is False


def test_the_comparand_does_not_claim_the_stagger_deck():
    """⚠️ `deck` records ⟨3D⟩ SOLVE PROVENANCE. The 1D MARCS deck is not a ⟨3D⟩ deck, and
    labelling it `stagger` would claim it came from the STAGGERmean3D solve."""
    assert taxes.GERBER1D_LTE_MARCS.deck == "none"
    assert "stagger" not in taxes.GERBER1D_LTE_MARCS.token


def test_the_treatment_is_registered_or_the_product_dies_after_the_synthesis():
    """🔴 THE RYA-798 FAILURE, HIT AGAIN WHILE FIXING THIS. A treatment the tuple has
    never heard of raises in `build_product` — AFTER every line has been fitted. Cost ~40
    minutes of Sirius time. The token is IMPORTED here, never retyped."""
    assert taxes.GERBER1D_LTE_MARCS.token in bp.TREATMENTS
    # Stronger than importing one name: TREATMENTS now ENUMERATES the axis-native
    # registry, so a member registered there cannot be unknown here. Naming members
    # individually is what let this one slip through in the first place.
    assert "_AXIS_NATIVE_TREATMENTS = tuple(treatment_axes.AXIS_NATIVE)" in (
        ROOT / "pipeline" / "band_products.py").read_text()
    assert set(taxes.AXIS_NATIVE) <= set(bp.TREATMENTS)


def test_the_driver_uses_the_derived_token_never_a_typed_label():
    """RYA-906/RYA-1040: new products are DEFINED as axes, not coined as labels."""
    assert "taxes.GERBER1D_LTE_MARCS.token" in DRIVER
    assert "ENGINE-B-LTE-MARCS" not in DRIVER, (
        "the coined label revives the ENGINE-* vocabulary RYA-906 retired")


def test_gerber_1d_lte_is_an_lte_run_so_it_withholds_departures():
    """Same deck, same atmosphere, same lines — departures OFF is the only difference."""
    assert '"gerber-1d-lte"' in DRIVER
    i = DRIVER.index('elif _nlte or a.engine_b_deck == "gerber-1d-lte":')
    seg = DRIVER[i:i + 3000]
    assert "if _nlte:" in seg and 'nlte_deck="gerber"' in seg, (
        "the deck must be handed over only on the NLTE member")


def test_gerber_1d_lte_is_not_counted_as_nlte():
    """`_nlte` gates the departure wiring; a comparand caught by it would not be one."""
    assert '_nlte = a.engine_b_deck in ("gerber-nlte", "gerber-mean3d")' in DRIVER
