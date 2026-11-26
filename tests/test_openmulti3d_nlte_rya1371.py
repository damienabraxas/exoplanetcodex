"""RYA-1371: the OpenMULTI3D Engine-B deck provider and its treatment pair.

No engine run is needed: the deck is replaced by a synthetic one where the physics is not the
point, so these run anywhere CI does.
"""
import numpy as np
import pytest

from pipeline import openmulti3d_nlte as om
from pipeline import treatment_axes as ta

KEEPER_HASHES = {  # Sirius Keeper snapshot, SHA256SUMS 2026-11-18
    "atom.o41f": "5cc8a84d7ebc6ba23b22d882e4cf0e7a38c4a707bc441e323350236bfed21af5",
    "atom.o41j": "9dc06ba050bfa0db2af144bdddd5d796b18add30aff075bf9030be0eefaf749d",
    "atom.ni538qm": "dcbebdf513f9a39643604ee2610b99fbfd0a942fdbe051f64d21f8714d3a0872",
}


def test_decks_pin_the_keeper_snapshot_hashes():
    for key, spec in om.DECKS.items():
        assert spec["atom_sha256"] == KEEPER_HASHES[spec["atom"]], key


def test_the_same_api_surface_as_the_gerber_provider():
    """abundances_derive swaps the module, so every name it calls must exist on both."""
    from pipeline import gerber_nlte
    for name in ("has_abundance_axis", "abundance_axis", "for_node", "as_ispec_tuple",
                 "assert_linelist_supports_nlte", "assert_depth_match"):
        assert callable(getattr(om, name)) and callable(getattr(gerber_nlte, name)), name


def test_axes_are_uniform_and_cover_the_solar_value():
    solar_abundance = {"O": 8.69, "Ni": 6.20}
    for key, spec in om.DECKS.items():
        ax = np.array(spec["abundance_axis"])
        assert np.allclose(np.diff(ax), ax[1] - ax[0]), key
        assert ax[0] < solar_abundance[om.element_of(key)] < ax[-1], key


def test_variant_keys_resolve_to_their_element():
    assert om.element_of("O") == "O" and om.element_of("O@o41j") == "O"


@pytest.fixture
def synthetic_deck(monkeypatch):
    ax = np.array(om.DECKS["O"]["abundance_axis"])
    dep = np.stack([np.full((56, 41), 1.0 + 0.1 * i) for i in range(len(ax))])
    monkeypatch.setitem(om._CACHE, "O", dict(axis=ax, dep=dep, lgTau5=np.linspace(-5, 2, 56), deck_dir="synthetic"))
    monkeypatch.setattr(om, "deck_spec", lambda key: dict(om.DECKS[key], atom_path="/x/atom.o41f", atmosphere="/x/m.mod",
                                                           abundance_axis=list(om.DECKS[key]["abundance_axis"])))
    return ax


def test_departures_are_linear_between_axis_nodes(synthetic_deck):
    ax = synthetic_deck
    p = om.for_node("O", 5772, 4.438, 0.0, abundance=float(ax[2] + 0.025))
    assert p["departures"].shape == (56, 41) and p["ndep"] == 56 and p["nk"] == 41
    assert np.allclose(p["departures"], 1.0 + 0.1 * 2.25)
    assert p["Z"] == 8 and p["deck_key"] == "O" and p["provider"] == "openmulti3d"


def test_refuses_another_star(synthetic_deck):
    with pytest.raises(om.OpenMULTI3DDeckError, match="solar MARCS node only"):
        om.for_node("O", 4800, 4.6, 0.0, abundance=8.7)


def test_refuses_a_missing_or_off_axis_abundance(synthetic_deck):
    with pytest.raises(om.OpenMULTI3DDeckError, match="abundance is required"):
        om.for_node("O", 5772, 4.438, 0.0)
    with pytest.raises(om.OpenMULTI3DDeckError, match="outside"):
        om.for_node("O", 5772, 4.438, 0.0, abundance=9.5)


def test_unregistered_key_is_refused():
    with pytest.raises(om.OpenMULTI3DDeckError, match="no OpenMULTI3D deck"):
        om.deck_spec("Xx")


def test_ispec_tuple_is_level_major(synthetic_deck):
    p = om.for_node("O", 5772, 4.438, 0.0, abundance=8.7)
    t = om.as_ispec_tuple(p, 8.7)
    assert t[0].shape == (41, 56) and t[4] == 8 and t[5] == 8.7


def test_treatment_pair_is_axis_native_and_paired():
    nl, lt = ta.OM1D_NLTE_MARCS, ta.OM1D_LTE_MARCS
    assert nl.token in ta.AXIS_NATIVE and lt.token in ta.AXIS_NATIVE
    assert ta.comparand_for(nl.token) == lt.token
    assert nl.is_nlte and not lt.is_nlte
    assert nl.atmos == lt.atmos == "marcs-ges" and nl.model == lt.model == "openmulti3d"
    assert "ENGINE" not in nl.token and "OpenMULTI3D" in nl.display


def test_the_fit_names_the_provider_and_refuses_unknown_ones():
    import inspect
    from pipeline import abundances_derive
    src = inspect.getsource(abundances_derive._fit_synth_flux)
    assert "'openmulti3d'" in src and "unknown NLTE deck" in src


def test_the_departure_tuple_declares_abundance_dependence(synthetic_deck):
    """synth_cache must step aside for it without consulting the Gerber registry: O has no
    Gerber deck (the lookup raised: 'no deck registered for O'), and an OpenMULTI3D Fe deck
    would be mis-cached as node-fixed because Gerber's Fe deck is."""
    t = om.as_ispec_tuple(om.for_node("O", 5772, 4.438, 0.0, abundance=8.7), 8.7)
    assert getattr(t, "abundance_dependent", False) is True
    a, b, *_ = t                                  # iSpec unpacks it as a plain tuple
    assert a.shape == (41, 56)


def test_synth_cache_steps_aside_for_an_abundance_dependent_tuple(synthetic_deck, monkeypatch):
    import inspect
    from pipeline import synth_cache
    src = inspect.getsource(synth_cache)
    assert 'getattr(t, "abundance_dependent", False)' in src


def test_an_x_label_is_no_label():
    """GESv6 marks O I 8446 nlte=T with labels "x"; bsyn finds no level "x" and runs it in
    LTE. Every guard counted it labelled, so it entered an NLTE aggregate (RYA-1371)."""
    from pipeline import gerber_nlte as gn
    assert list(gn.is_label(["3s5S2*", "x", "none", "", " x "])) == [True, False, False, False, False]
    rec = np.array([("O 1", 8.0, 7771.944, "T", "3s5S2*", "3p5P3"),
                    ("O 1", 8.0, 8446.359, "T", "x", "x")],
                   dtype=[("element", "U8"), ("turbospectrum_species", float),
                          ("wave_A", float), ("nlte", "U2"),
                          ("nlte_label_low", "U20"), ("nlte_label_up", "U20")])
    m = gn.pool_label_mask(rec, 8, [7771.944, 8446.359], ion=None)
    assert list(m) == [True, False]


def test_the_orch_dispatches_the_pair_and_holds_what_it_cannot_run():
    from pipeline import run_matrix as rm
    from pipeline.run_descriptor import RunDescriptor
    assert rm.DECK_EMITS["openmulti3d-nlte"] == (ta.OM1D_NLTE_MARCS.token,)
    assert rm.DECK_EMITS["openmulti3d-1d-lte"] == (ta.OM1D_LTE_MARCS.token,)
    assert rm.DECK_PAIRS["openmulti3d-nlte"] == "openmulti3d-1d-lte"
    rm.assert_registry_covered()
    cell = dict(ion="I", instrument="iag_fts_solar_atlas", holding="solar_iag", lo_A=6910.0, hi_A=9199.0)
    assert rm._deck_hold(RunDescriptor(element="O", engine_deck="openmulti3d-nlte", **cell)) == ""
    assert "no OpenMULTI3D deck for Zn" in rm._deck_hold(
        RunDescriptor(element="Zn", engine_deck="openmulti3d-nlte", **cell))
    assert "solar MARCS node only" in rm._deck_hold(
        RunDescriptor(element="O", engine_deck="openmulti3d-1d-lte", star="procyon", **cell))
