"""RYA-1040: the ⟨3D⟩ treatment pair, defined as AXES rather than coined as a label.

Ryan ratified the naming on 2026-08-25:

    "DROP the label ENGINE-B-MEAN3D — it revives the ENGINE-A/B vocabulary RYA-906
     retired. The science-facing label follows route.scale.model.deck-provenance …
     PRODUCE ITS PAIRED COMPARAND, same deck, departures off — so the NLTE effect is
     (⟨3D⟩-NLTE minus ⟨3D⟩-LTE) on ONE atmosphere."

So there are two new products, they are defined by their axes, both names are derived, and
the pair is mandatory. These tests hold each of those to be true — and one of them holds a
defect this change EXPOSED in the module it extends.
"""
import pytest

from pipeline import band_products as BP
from pipeline import treatment_axes as T


# ── the defect adding an LTE scale exposed ──────────────────────────────────

def test_is_nlte_is_not_outgrown_by_an_LTE_variant():
    """🔴 THE MODULE'S OWN INVARIANT, FAILING IN THE DIRECTION NOBODY WATCHED.

    `is_nlte` was `scale != "1D-LTE"` and its docstring claimed it "cannot be outgrown by
    a new variant". That is a spelling list of length one wearing a field test's clothes:
    it assumed the only LTE scale would ever be the 1D one. RYA-1040's REQUIRED comparand
    is `<3D>-LTE`, and under the old implementation it answered True — an LTE product
    counted as NLTE by the property that exists to answer "everything NLTE"."""
    assert T.MEAN3D_LTE_STAGGER.is_nlte is False
    assert T.MEAN3D_NLTE_STAGGER.is_nlte is True
    # every declared scale is classified, and the two vocabularies do not overlap
    assert set(T.LTE_SCALES) | set(T.NLTE_SCALES) == set(T.SCALES)
    assert not (set(T.LTE_SCALES) & set(T.NLTE_SCALES))


def test_an_unclassified_scale_RAISES_rather_than_defaulting_to_NLTE():
    """A silent default is how RYA-869 published four wrong systematics. The old code
    defaulted anything that was not `1D-LTE` to NLTE, which is the wrong direction to be
    wrong in: it over-claims."""
    rogue = T.Axes(route="synth", scale="4D-SIDEWAYS", model="none",
                   atmos="atlas9", gf="kurucz")
    with pytest.raises(ValueError, match="in neither LTE_SCALES nor NLTE_SCALES"):
        _ = rogue.is_nlte


# ── the naming Ryan ratified ────────────────────────────────────────────────

def test_no_ENGINE_letter_survives_in_the_new_names():
    """🔴 THE POINT OF THE RATIFICATION. `ENGINE-B-MEAN3D` would have been the seventh
    product named after a letter, in a vocabulary RYA-906 retired for causing three real
    defects in one cycle."""
    for ax in (T.MEAN3D_NLTE_STAGGER, T.MEAN3D_LTE_STAGGER):
        assert "ENGINE" not in ax.token.upper()
        assert "ENGINE" not in ax.display.upper()


def test_the_label_carries_all_four_ratified_axes():
    """route · scale · model · deck-provenance. The deck axis is what makes `stagger`
    (the vendor's deck) distinguishable from `codex` (ours, later) in the NAME rather than
    only in a provenance file nobody reads at a glance."""
    d = T.MEAN3D_NLTE_STAGGER.display
    assert "Synth" in d and "<3D>-NLTE" in d and "Gerber" in d and "stagger" in d
    assert T.MEAN3D_LTE_STAGGER.display.count("<3D>-LTE") == 1


def test_mean3D_is_not_3D_and_the_vocabulary_keeps_them_apart():
    """🔴 `<3D>` AND `3D` ARE DIFFERENT PHYSICS. `3D-NLTE` is Amarsi+2022 solving through
    the inhomogeneous cube; `<3D>-NLTE` averages the cube on constant-τ500 surfaces FIRST,
    which removes the very term that distinguishes 3D from 1D. Reporting a mean-3D number
    under `3D-NLTE` claims the expensive result while having run the cheap one."""
    assert "<3D>-NLTE" in T.SCALES and "3D-NLTE" in T.SCALES
    assert T.MEAN3D_NLTE_STAGGER.scale != "3D-NLTE"
    # the existing full-3D product is untouched and still says 3D
    assert T.axes_for("ENGINE-A-3DNLTE").scale == "3D-NLTE"


def test_the_token_is_filesystem_safe_and_the_label_is_not_expected_to_be():
    """🔴 WHY THERE ARE TWO DERIVED NAMES. `derive_band_products` writes
    `f"{stem}_{prod.treatment}_lines.csv"`, so the token lands in a FILENAME — and `<`/`>`
    are shell redirection characters. The science-facing label keeps Ryan's `<3D>` form;
    the key is a slug of the SAME axes."""
    for ax in (T.MEAN3D_NLTE_STAGGER, T.MEAN3D_LTE_STAGGER):
        assert not (set(ax.token) & set('<>|&;:*?"\'$ /\\')), ax.token
        assert "mean3D" in ax.token
    assert "<3D>" in T.MEAN3D_NLTE_STAGGER.display


def test_the_token_and_the_label_cannot_drift_apart():
    """Both are DERIVED from one Axes object, so a change to one that is not a change to
    the other is impossible by construction — this asserts the relationship, not the two
    strings (RYA-845)."""
    for ax in T.AXIS_NATIVE.values():
        assert ax.scale.replace("<3D>", "mean3D") in ax.token
        assert ax.scale in ax.display
        # `token` deliberately omits a "none" deck (there is no provenance to record),
        # so the relationship is conditional -- RYA-1045's 1D comparand is the first
        # axis-native product with no ⟨3D⟩ deck behind it.
        if ax.deck != "none":
            assert ax.deck in ax.token and ax.deck in ax.display
        assert (ax.model in ax.token) and (T.DISPLAY_MODEL[ax.model] in ax.display)


# ── the pair is mandatory, and the physics says why ─────────────────────────

def test_the_NLTE_product_declares_its_LTE_comparand():
    """🔴 THE NLTE EFFECT IS A DIFFERENCE OF TWO PRODUCTS ON ONE ATMOSPHERE. Differencing
    ⟨3D⟩-NLTE against 1D-LTE instead folds the 1D→mean-3D ATMOSPHERE shift into the
    reported non-LTE effect — the confound RYA-542 disentangled for Ti, where MARCS at
    +0.203 reproduced the deck's +0.221 and the difference WAS the atmosphere."""
    assert T.comparand_for(T.MEAN3D_NLTE_STAGGER.token) == T.MEAN3D_LTE_STAGGER.token


def test_the_pair_differs_in_SCALE_ALONE():
    """That is what makes the difference attributable to non-LTE physics. Same deck, same
    atmosphere, same line list, same gf pool, same route — departures off."""
    a, b = T.MEAN3D_NLTE_STAGGER, T.MEAN3D_LTE_STAGGER
    differing = {f for f in ("route", "model", "atmos", "gf", "deck")
                 if getattr(a, f) != getattr(b, f)}
    assert differing == set(), f"the pair must differ in scale alone, but also in {differing}"
    assert a.scale != b.scale


def test_the_LTE_member_keeps_the_deck_identity():
    """Calling it `model="none"` would make it indistinguishable from an unrelated LTE run
    and destroy the pairing — its identity is 'the LTE limit OF THAT DECK'S SETUP'."""
    assert T.MEAN3D_LTE_STAGGER.model == "gerber"
    assert T.MEAN3D_LTE_STAGGER.deck == "stagger"
    assert T.MEAN3D_LTE_STAGGER.atmos == "stagger-mean3d"


def test_the_atmosphere_is_the_3D_one_and_that_breaks_the_legacy_inference():
    """⚠️ The legacy table infers `gerber -> marcs-ges`. That correlation is now FALSE for
    these two, which is exactly the day the module's docstring predicted — and the reason
    Ryan ratified storing `atmos` independently rather than deriving it."""
    assert T.MEAN3D_NLTE_STAGGER.atmos == "stagger-mean3d"
    assert T._LEGACY_ATMOS["gerber"] == "marcs-ges"      # still right for the 1D decks
    assert T.axes_for("ENGINE-B-NLTE").atmos == "marcs-ges"


# ── the registries agree, by construction ───────────────────────────────────

def test_axes_for_resolves_an_axis_native_token_without_touching_LEGACY():
    """A new product must not depend on the legacy label table at all."""
    for token, ax in T.AXIS_NATIVE.items():
        assert token not in T.LEGACY
        assert T.axes_for(token) is ax
        assert T.axes_for(token).route_basis == "axis-native"


def test_band_products_imports_the_tokens_rather_than_retyping_them():
    """RYA-798 emitted a treatment `TREATMENTS` had never heard of, and the product died
    at `build_product` AFTER the synthesis had run. Importing makes that disagreement
    impossible rather than merely tested-for."""
    for ax in T.AXIS_NATIVE.values():
        assert ax.token in BP.TREATMENTS
    assert set(BP._AXIS_NATIVE_TREATMENTS) == set(T.AXIS_NATIVE), (
        "TREATMENTS must ENUMERATE the axis-native registry, not name its members: "
        "naming them is how RYA-1045's comparand was registered in one place and "
        "unknown in the other, and build_product raised after the synthesis had run")
    assert set(BP._MEAN3D_TREATMENTS) <= set(T.AXIS_NATIVE)


def test_the_handler_table_stays_exhaustive_over_the_grown_vocabulary():
    """The RYA-869 guarantee, re-checked with two more members in the tuple."""
    from pipeline.harness_residual import _BANKED_PROFILEFIT_HANDLER as H
    assert set(H) == set(BP.TREATMENTS)
    for ax in T.AXIS_NATIVE.values():
        assert H[ax.token] == "SynthesisHandler"


def test_every_legacy_product_is_unchanged():
    """The whole vocabulary change must be a pure ADDITION — a moved legacy label makes
    every historical product incomparable to its own past (the RYA-874 lesson)."""
    assert T.axes_for("1D-LTE").display == "1D-LTE"
    assert T.axes_for("ENGINE-B-NLTE").display == "Synth · 1D-NLTE · Gerber"
    # RYA-1106 -- the ONE legacy display this ticket deliberately moves, and it moves
    # on measured grounds: the route was a stranded handler's (RYA-1104) and the gf
    # axis said kurucz over a 67/67 primary-laboratory pool. Everything else above is
    # still asserted unchanged, which is what makes this line a statement rather than
    # a loosened test.
    assert T.axes_for("ENGINE-A-3DNLTE").display == "3D-NLTE · Amarsi · lab-gf"
    for label in T.LEGACY:
        assert T.axes_for(label).deck == "none", label
