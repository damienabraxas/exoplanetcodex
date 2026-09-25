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
    """Where the assignment is not established, say so rather than writing (0-0).

    OH_H and CO_K were UNSPECIFIED here until the diagnostic registry's own reference
    strings supplied their bands -- sourced, so they are now declared. CN_red still has
    no published assignment and must stay honest about it.
    """
    assert mi.identity_for("CN_red").band == mi.UNSPECIFIED
    for key in ("OH_H", "CO_K"):
        band = mi.identity_for(key).band
        assert band != mi.UNSPECIFIED and band != "(0-0)"


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


# --- red_chi2 must never regain the status of a gate -----------------------------------

def test_red_chi2_is_never_used_as_a_gate():
    """RYA-847/RYA-1152: this guard was tried and refuted on this data.

    red_chi2 > 10 flags 1262 of 1366 in-aggregate lines because the HARPS ERR column is
    all NaN, so the absolute chi2 is meaningless. Any field that reads as a verdict
    (fit_constrained, fit_ok) derived from it would re-introduce the refuted guard.
    """
    import json
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[1] / "data/output/rya1220"
    for name in ("molecular_route_model_form.json", "cno_hold_readjudication.json"):
        path = root / name
        if not path.exists():
            continue
        doc = json.loads(path.read_text())
        blob = json.dumps(doc)
        assert "fit_constrained" not in blob, (
            f"{name} derives a verdict from red_chi2; that guard is refuted")
        for route in doc.get("routes", []):
            assert "red_chi2_max" not in route, (
                "red_chi2 must be carried only under an UNCALIBRATED name")


def test_no_frac_rise_threshold_is_applied():
    """RYA-847 swept 9 cells / 581 lines and refuted every candidate cut."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1]
           / "scripts/rya1220_readjudicate_holds.py").read_text()
    assert "ranked_by_frac_rise_weaker" in src
    assert "no_threshold_applied" in src


# --- completeness against the live diagnostic registry --------------------------------

def test_every_molecular_diagnostic_in_the_registry_is_declared():
    """The table must cover what the pipeline can actually run, not what I happened to see.

    The first version declared 6 identities -- the ones that appeared in the RYA-1220 route
    outputs -- and the Sirius run refused on CH_Gband. That was the no-fallback guard doing
    its job, but the table should have been derived from the registry in the first place.

    cno_synthesis needs ispec, so this can only run on a synthesis host; it skips elsewhere
    rather than passing vacuously.
    """
    try:
        from pipeline.cno_synthesis import REGION_DIAGNOSTICS
    except Exception:
        pytest.skip("cno_synthesis needs ispec; this check belongs on a synthesis host")

    required = {d.key for diags in REGION_DIAGNOSTICS.values() for d in diags
                if getattr(d, "kind", None) == "molecular_band"}
    missing = sorted(required - set(mi.IDENTITIES))
    assert not missing, f"molecular diagnostics with no declared identity: {missing}"


def test_the_nine_known_diagnostics_stay_declared():
    """Runs everywhere -- pins the set even where ispec is absent."""
    assert set(mi.IDENTITIES) == {
        "C2_Swan", "CH_Gband", "CN_AX_IR", "CN_AX_J", "CN_red",
        "CO_K", "NH_AX", "OH_AX", "OH_H"}


def test_vibrational_bands_are_taken_from_the_registry_not_invented():
    """OH_H and CO_K were UNSPECIFIED until the registry's own reference strings were read."""
    assert mi.identity_for("OH_H").band == "(2-0)/(3-1)/(4-2)"
    assert mi.identity_for("CO_K").band == "(2-0)/(3-1)"
    # CN_red genuinely has no published assignment, so it stays honest
    assert mi.identity_for("CN_red").band == mi.UNSPECIFIED


def test_the_two_near_uv_electronic_routes_are_distinguished_from_their_ir_namesakes():
    """OH_AX (3064 A, electronic) is not OH_H (1.5 um, vibration-rotation)."""
    assert mi.identity_for("OH_AX").system == "A-X"
    assert mi.identity_for("OH_H").system == "X-X"
    assert mi.identity_for("OH_AX").molecule == mi.identity_for("OH_H").molecule
