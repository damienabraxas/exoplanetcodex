"""RYA-798 — the Gerber NLTE wiring, pinned where it is silent-wrong rather than loud.

Every failure mode this ticket had to defend against returns an ANSWER rather than an
error: iSpec synthesises in LTE when an element is unlabelled, pairs departures with the
wrong depths when the layer counts disagree, and accepts an ATLAS9 atmosphere whose tau
column is all zeros. So these tests are about the guards, not the happy path.
"""
from __future__ import annotations

import numpy as np
import pytest

from pipeline import gerber_nlte as gn

STAGED = gn.Path(gn.GT).exists()
sirius = pytest.mark.skipif(not STAGED, reason="Gerber decks are Sirius-only")


def test_lte_path_is_untouched_by_default():
    """RYA-770 stabilised the LTE flux fit at -0.026 dex against the banked answer. The
    NLTE parameter defaults to None and the LTE branch must remain the same call."""
    import inspect

    from pipeline.abundances_derive import _fit_synth_flux, _synth_flux_at_abund

    assert inspect.signature(_synth_flux_at_abund).parameters["nlte_departures"].default is None
    assert inspect.signature(_fit_synth_flux).parameters["nlte_deck"].default is None
    src = inspect.getsource(_synth_flux_at_abund)
    # the LTE call must not pass the NLTE kwarg at all
    lte_branch = src.split("if nlte_departures is not None:")[-1].split("return ispec")[-1]
    assert "nlte_departure_coefficients" not in lte_branch


def test_unregistered_element_raises_rather_than_running_lte():
    """⚠️ RYA-1005 re-pointed this AT THE INVARIANT. It used Al as its example of an
    unregistered deck; Al then got one, and the test failed for the best possible
    reason. The rule is "a species with no deck RAISES rather than silently running
    LTE" — so the example is chosen from the registry at run time and asserted absent,
    which cannot expire the same way."""
    unregistered = next(e for e in ("Ti", "Cr", "Ba", "Zn")
                        if e not in gn.DECKS)
    with pytest.raises(gn.GerberDeckError, match="no TS-native Gerber deck registered"):
        gn.for_node(unregistered, 5772.0, 4.44, 0.0)


@sirius
def test_non_solar_node_refuses_instead_of_returning_solar_departures():
    """The departure path passes ONE MARCS model eight times — a degenerate interpolation
    box. Asking for Procyon must fail, not silently hand back the Sun's departures."""
    with pytest.raises(gn.GerberDeckError, match="solar-node-only"):
        gn.for_node("Fe", 6554.0, 4.0, 0.0)


@sirius
def test_departures_do_not_depend_on_the_abundance():
    """The finding that shapes the whole design: the deck has NO abundance axis, so the
    coefficients are cached per node and only the STAMP follows each trial abundance.
    If this ever becomes false, the chi2 loop must re-interpolate per evaluation."""
    p = gn.for_node("Fe", 5772.0, 4.44, 0.0)
    t_lo = gn.as_ispec_tuple(p, 7.20)
    t_hi = gn.as_ispec_tuple(p, 7.80)
    assert np.array_equal(t_lo[0], t_hi[0]), "departures must not move with the stamp"
    assert t_lo[5] == 7.20 and t_hi[5] == 7.80, "the stamp must follow the trial value"


@sirius
def test_departures_are_physical_and_oriented_for_ispec():
    p = gn.for_node("Fe", 5772.0, 4.44, 0.0)
    assert p["nk"] == 607, "atom.fe607a carries 607 levels"
    assert p["departures"].shape == (p["ndep"], p["nk"]), "file layout is (ndep, nk)"
    # b -> 1 in the deep layers is the format check RYA-402 uses on these grids
    assert p["departures"][-1].mean() == pytest.approx(1.0, abs=0.02)
    # ...and iSpec wants the transpose, which it checks and raises on
    assert gn.as_ispec_tuple(p, 7.46)[0].shape == (p["nk"], p["ndep"])


@sirius
def test_depth_mismatch_raises_because_ispec_overwrites_the_tau():
    """iSpec replaces the departure tau with `atmosphere_layers[:, 7]`. If the layer counts
    differ, depth i of the departures is paired with a different physical depth — and it
    still returns a spectrum. ATLAS9 (72 layers) vs this deck (56) is exactly that case."""
    p = gn.for_node("Fe", 5772.0, 4.44, 0.0)
    with pytest.raises(gn.GerberDeckError, match="depth mismatch"):
        gn.assert_depth_match(p, np.zeros((72, 10)))
    gn.assert_depth_match(p, np.zeros((56, 11)))       # MARCS.GES matches


def test_unlabelled_linelist_raises_rather_than_synthesising_lte():
    """THE RYA-764 TRAP. iSpec appends the element to `nlte_ignored` and proceeds in LTE
    without raising, and never reports `nlte_available` back to the caller — so the check
    has to happen before the call or not at all."""
    ll = np.array([(26.0, "T", "none", "none"), (26.0, "F", "none", "none")],
                  dtype=[("turbospectrum_species", "f8"), ("nlte", "U1"),
                         ("nlte_label_low", "U8"), ("nlte_label_up", "U8")])
    with pytest.raises(gn.GerberDeckError, match="NONE carry NLTE level labels"):
        gn.assert_linelist_supports_nlte(ll, 26, "Fe")

    ok = np.array([(26.0, "T", "a5D0", "z7D1")],
                  dtype=[("turbospectrum_species", "f8"), ("nlte", "U1"),
                         ("nlte_label_low", "U8"), ("nlte_label_up", "U8")])
    assert gn.assert_linelist_supports_nlte(ok, 26, "Fe") == 1


def test_absent_element_in_window_also_raises():
    ll = np.array([(20.0, "T", "a5D0", "z7D1")],
                  dtype=[("turbospectrum_species", "f8"), ("nlte", "U1"),
                         ("nlte_label_low", "U8"), ("nlte_label_up", "U8")])
    with pytest.raises(gn.GerberDeckError, match="no Fe lines"):
        gn.assert_linelist_supports_nlte(ll, 26, "Fe")


@sirius
def test_deck_abundance_comes_from_provenance_never_a_reference_value():
    """Fe's grid was computed at 7.50 — which IS the value in the atom header — and bsyn
    STOPs on a mismatch. It must be read from the committed provenance record.

    🔴 CORRECTED RYA-1035. This docstring said "7.46, NOT the 7.50 in the atom header",
    and had it exactly backwards: the atom header is the grid's own declaration
    (`atom.fe607a` line 2, `7.50  55.85`, md5-matched to the staged copy), while the 7.46
    was OUR OWN INPUT coming back — `abu_ref` is read from stdin
    (interpol_modeles_nlte.f:206), written verbatim into the departure file (:761), loaded
    as `abundance_nlte` and printed by bsyn (:988), and `gerber_nlte` fed that stdin from
    `deck_abundance()` itself. A closed loop with no external referee, and 7.46 survived
    because it is the Asplund solar A(Fe) — the number that looks right on arrival.

    The test's PREMISE stands and is the reason it exists: never a reference value, always
    the record. RYA-1035 additionally has `deck_abundance` cross-examine that record
    against the deck's own aux, so the loop cannot re-close."""
    assert gn.deck_abundance("Fe") == pytest.approx(7.50)
    with pytest.raises(gn.GerberDeckError, match="no provenance record"):
        gn.deck_abundance("Xx")


def test_engine_b_nlte_is_a_separate_product_never_a_correction():
    """RYA-712: engines are separate products and are never combined."""
    src = (gn.ROOT / "scripts" / "derive_band_products.py").read_text()
    assert "ENGINE-B-NLTE" in src
    assert "MARCS.GES" in src, "the deck is MARCS; ATLAS9's index 7 is not a tau column"
