"""RYA-821: a <3D> deck is read DIRECTLY; no vendor binary can consume one.

These tests do NOT touch the deck binary -- it is Sirius-only and several GB. They pin
the CONTRACT and the REFUSALS, which is what a wrong answer here would come through.
The numbers quoted in the comments were measured against the real deck on Sirius.
"""
import numpy as np
import pytest

from pipeline import gerber_nlte as G


# ── the deck is registered, and carries its own atmosphere ───────────────────

def test_the_mean3d_deck_is_registered_separately_from_the_marcs_one():
    """Two decks for the SAME element differing only in the atmosphere their departures
    were solved on. A registry keyed by element alone cannot hold both."""
    assert "Al" in G.DECKS and "Al@mean3D" in G.DECKS
    assert G.DECKS["Al@mean3D"]["grid"] != G.DECKS["Al"]["grid"]


def test_each_deck_carries_the_atmosphere_it_was_solved_on():
    """Applying <3D> departures on a 1D MARCS structure produces a perfectly well-formed
    file, which is exactly why this is a guard and not a comment."""
    assert G.deck_atmosphere("Al") == G.MARCS_SOLAR
    assert G.deck_atmosphere("Al@mean3D") != G.MARCS_SOLAR
    assert "stagger" in G.deck_atmosphere("Al@mean3D").lower()


def test_the_mean3d_deck_is_flagged_for_the_direct_reader():
    """The MARCS decks must keep the interpolator path byte-for-byte."""
    assert G.DECKS["Al@mean3D"].get("read_via") == "direct"
    assert G.DECKS["Al"].get("read_via") is None
    assert G.DECKS["Fe"].get("read_via") is None


# ── the name parser: two conventions for one model ───────────────────────────

def test_both_naming_conventions_resolve_to_the_same_node():
    """🔴 THE TRAP THAT BROKE THE VENDOR BINARY. The record carries STAGGER's own name
    while the aux carries a MARCS-style alias -- which is why the aux file is called
    `_marcs_names`. A STRING comparison rejects the CORRECT record, and that is exactly
    what the first version of this reader did."""
    assert G._node_from_model_name("t5777g44m00") == (5777.0, 4.4, -0.0)
    assert G._node_from_model_name(
        "p5777_g+4.4_m0.0_t02_st_z+0.00_a+0.00") == (5777.0, 4.4, -0.0)


def test_the_metallicity_letter_is_a_SIGN_not_a_digit():
    """⚠️ In the STAGGER form `m10` is [Fe/H] = -1.0, NOT +1.0. Reading it unsigned is how
    a metallicity sign gets lost -- a defect this project has already met once, in the
    STAGGER GRID_STATUS table."""
    assert G._node_from_model_name("t5777g44m10")[2] == pytest.approx(-1.0)
    assert G._node_from_model_name("t4500g25m20")[2] == pytest.approx(-2.0)
    assert G._node_from_model_name("t5777g44m00")[2] == pytest.approx(0.0)


def test_an_unparseable_name_is_not_silently_accepted():
    assert G._node_from_model_name("garbage") is None


# ── the tolerances are DERIVED, and they matter ──────────────────────────────

def test_the_node_tolerance_admits_our_solar_star_but_not_a_neighbour():
    """MEASURED node spacing: the deck's Teff axis is 4000/4500/5000/5500/**5777**/6000/
    6500/7000 -- an otherwise 500 K grid with STAGGER's solar member dropped in, so
    5777's nearest neighbour is 223 K away. Our solar star is 5772/4.438 (IAU) against
    the deck's 5777/4.44 (STAGGER): a 5 K CONVENTION difference for the same physical
    star. The tolerance must clear the first and never reach the second."""
    from config.constants import STAR_PARAMS
    solar = STAR_PARAMS["solar"]
    assert abs(float(solar["teff"]) - 5777.0) <= G._NODE_TOL_TEFF
    assert abs(float(solar["logg"]) - 4.44) <= G._NODE_TOL_LOGG
    assert G._NODE_TOL_TEFF < 223.0 / 2, "tolerance could reach a neighbouring node"


# ── the refusals, which are the safety property ──────────────────────────────

def test_reading_an_unregistered_deck_refuses():
    with pytest.raises(G.GerberDeckError):
        G.read_deck_node("Xx@mean3D", 5777.0, 4.44, 0.0, 6.43)


def test_the_direct_route_is_taken_only_for_flagged_decks():
    """Pinned on the source: a silent fall-through to the interpolator for a <3D> deck
    would return all-zero departures (measured) that still look like a valid file."""
    import inspect
    src = inspect.getsource(G.for_node)
    assert 'read_via' in src and 'read_deck_node' in src


def test_the_reader_returns_the_same_contract_as_the_vendor_parser():
    """Both paths must hand downstream the same keys, or the caller has to know which
    route ran -- and the whole point is that it does not."""
    import inspect
    src = inspect.getsource(G.read_deck_node)
    for key in ("abundance", "ndep", "nk", "tau", "departures", "corners"):
        assert f"{key}=" in src, f"{key} missing from the direct reader's return"


def test_the_deck_tau_is_converted_to_LOG():
    """🔴 THE DECK STORES LINEAR TAU (measured 1e-5..1e5); `read_departure_file` returns
    LOG tau. Nothing downstream reads these values today -- iSpec overwrites the departure
    tau with the atmosphere's own and `assert_depth_match` checks only the COUNT -- which
    is precisely why a 5-orders-of-magnitude mismatch would have sat here undetected until
    the first caller that did read them."""
    import inspect
    assert "np.log10(tau)" in inspect.getsource(G.read_deck_node)


# ── the matrix must not be able to claim a wiring the code lacks ─────────────

def test_the_matrix_derives_the_wiring_claim_from_the_code():
    """🔴 THE DEFECT THAT STARTED THIS. `_MEAN3D_WIRED` was a hard-coded
    `{"Al": "Al@mean3D"}` and it LIED: main rendered "WIRED via gerber_nlte deck
    'Al@mean3D'" for days while that deck existed only on an unmerged branch. A status
    surface whose job is to say what the pipeline can do must READ the pipeline."""
    from pipeline.model_availability_matrix import _MEAN3D_WIRED, _mean3d_wired
    assert _MEAN3D_WIRED == _mean3d_wired()
    assert set(_MEAN3D_WIRED.values()) <= set(G.DECKS), (
        "the matrix names a deck that gerber_nlte does not register")


# ── the abundance axis: whether departures may be hoisted is a DECK property ─

def test_hoisting_out_of_the_chi2_loop_is_decided_PER_DECK():
    """🔴 THE BUG THIS CLOSES. `abundances_derive` computed the departures ONCE, outside
    the fit, on a comment that is TRUE OF Fe AND FALSE OF Al. Fe's deck holds one A(X);
    Al's resolves 31, and RYA-1005 MEASURED its departures genuinely differing between
    them (adjacent 0.1-dex nodes differ by up to 10.3 in b).

    Worse than wrong: `for_node` REFUSES an axis deck given no abundance, so Engine-B NLTE
    could not run Al AT ALL -- on either the MARCS or the <3D> deck -- while the matrix
    reported the deck as staged and wired.
    """
    # Read as TEXT, not imported: `abundances_derive` uses `X | None` annotations that
    # some interpreters evaluate at import, so importing it makes this pass or fail on the
    # interpreter rather than on the code -- an environment check wearing a policy check's
    # label, the same reason the RYA-1026 tests parse `measure_band_ew` with `ast`.
    from pathlib import Path as _P
    src = (_P(__file__).resolve().parent.parent
           / "pipeline" / "abundances_derive.py").read_text()
    assert "has_abundance_axis" in src, (
        "the call site must ask the DECK whether departures depend on abundance")
    assert "abundance=float(a_x)" in src, (
        "an axis deck must be evaluated at each TRIAL abundance inside the fit")


def test_interpolation_degenerates_to_an_exact_read_on_a_node():
    """If it did not, a run that happened to land on a node would disagree with a direct
    read of the same node -- two answers to one question."""
    import inspect
    src = inspect.getsource(G.departures_at_abundance)
    assert "read_deck_node(element, teff, logg, feh, exact[0])" in src


def test_off_axis_abundance_is_refused_not_extrapolated():
    """Extrapolating departures off the grid is not a correction, it is an invention."""
    import inspect
    src = inspect.getsource(G.departures_at_abundance)
    assert "outside this deck's abundance axis" in src


def test_bracketing_nodes_must_agree_in_shape_before_blending():
    """Interpolating across differently-shaped nodes would pair unrelated depths and
    levels while still producing a well-formed block."""
    import inspect
    assert "different shapes" in inspect.getsource(G.departures_at_abundance)
