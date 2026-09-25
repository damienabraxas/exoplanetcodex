"""RYA-1220 -- molecular identity is declared, and the 1D->3D term is measured.

Two defects this pins:

1. `molecule = 'NH' if diag.key == 'NH_AX' else '12C14N'` stamped OH_H and CO_K as
   CN A-X (0-0), a band with no transitions in either window. The label is hashed into
   indicator_id, the key the uncertainty contract joins on.
2. Every molecular route carried `3d_nlte_model: HOLD_MOLECULAR` while Amarsi 2021
   Table 2 already prices the 1D->3D shift per molecule.
"""
from __future__ import annotations

import pytest

from pipeline import molecular_identity as mi
from pipeline import molecular_model_form as mf


# --- identity ------------------------------------------------------------------------

def test_unknown_diagnostic_raises_instead_of_defaulting_to_cn():
    """The fallback is the bug. Absence must be loud."""
    with pytest.raises(mi.MolecularIdentityError):
        mi.identity_for("OH_K_not_declared")


def test_oh_and_co_are_not_cn():
    for key, molecule, element in [("OH_H", "OH", "O"), ("CO_K", "12C16O", "C")]:
        ident = mi.identity_for(key)
        assert ident.molecule == molecule
        assert ident.element == element
        assert ident.system == "X-X", "these are vibration-rotation, not electronic"
        assert ident.molecule != "12C14N"


def test_no_route_invents_a_vibrational_band():
    """Where the assignment is not established, say so rather than writing (0-0)."""
    for key in ("OH_H", "CO_K", "CN_red"):
        assert mi.identity_for(key).band == mi.UNSPECIFIED


def test_the_cn_ax_00_band_is_only_claimed_where_it_exists():
    for key in ("CN_AX_J", "CN_AX_IR"):
        ident = mi.identity_for(key)
        assert (ident.molecule, ident.system, ident.band) == ("12C14N", "A-X", "(0-0)")


def test_pool_shape_is_preserved():
    pool = mi.pool_for("CN_AX_J", [[11640.31, 11643.31]])
    assert set(pool) == {"molecule", "system", "band", "diagnostic", "windows_air_A"}


# --- model form ----------------------------------------------------------------------

def test_term_is_measured_on_the_reference_band():
    got = mf.model_form_term("CN", 11640, 12076)
    assert got["state"] == "MEASURED"
    assert got["sigma_dex"] == pytest.approx(0.061, abs=1e-3)
    assert got["sigma_dex"] > 0, "a model-form term is never zero"


def test_offband_routes_hold_with_a_cited_reason():
    """RYA-1226: 'unpriceable' with no reason re-buries measured work."""
    for molecule, lo, hi in [("CN", 6125, 6200), ("NH", 3358, 3373)]:
        got = mf.model_form_term(molecule, lo, hi)
        assert got["state"] == "HOLD"
        assert got["sigma_dex"] is None
        assert "outside the reference used-line range" in got["reason"]
        assert str(int(lo)) in got["reason"]


def test_hold_raises_rather_than_returning_a_zero_term():
    with pytest.raises(mf.ModelFormError):
        mf.as_term("NH", 3358, 3373)


def test_the_term_is_molecule_specific():
    """A single uniform molecular 3D offset would be wrong by up to 0.18 dex."""
    ch = mf.model_form_term("CH", 4254, 37955)["published_3D_minus_MARCS_dex"]
    co = mf.model_form_term("CO", 22954, 63292)["published_3D_minus_MARCS_dex"]
    assert ch > 0 > co
    assert abs(ch - co) > 0.15


def test_an_unknown_molecule_does_not_borrow_a_neighbours_value():
    with pytest.raises(mf.ModelFormError):
        mf.model_form_term("SiH", 4000, 5000)
