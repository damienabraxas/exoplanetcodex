"""RYA-1135 — Fe II <3D>-LTE, and the phantom-departure gate that makes it safe.

Fe II NLTE is STRUCTURALLY UNAVAILABLE: `atom.fe607a` declares 12,635 bound-bound
transitions and not one involves an Fe II level (RYA-1055), so bsyn applies departure = 1
to every Fe II line whatever the line list says. An Fe II product labelled NLTE would be
LTE wearing an NLTE name — the class found three times already (RYA-783's two Fe II
ENGINE-B-NLTE cells, RYA-1113's near-UV n=7).

🔴 THE TICKET ASKS FOR A model5-vs-model6 COMPARISON AND THAT COMPARISON CANNOT BE BUILT.
`derive_band_products --engine-b-deck gerber-mean3d` on the Fe II VIS pool REFUSES before
synthesising, because not one of its 9 pooled lines carries an NLTE label. So there is no
model 6 artifact to difference against model 5 — and that refusal is the correct outcome,
not a failure. The tests below therefore assert the REFUSAL and the reachability of the
hole it left open, rather than differencing two products one of which cannot exist.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline import gerber_nlte as gn

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "data" / "audit" / "rya1135_fe2_mean3d_lte"


def test_the_atom_cannot_serve_Fe_II_and_the_registry_says_so():
    capable, limit = gn.nlte_ion_capability("Fe", "II")
    assert capable is False
    assert "ZERO Fe II bound-bound transitions" in limit
    assert "no line list can enable it" in limit


def test_Fe_I_is_still_capable_so_the_gate_is_not_a_blanket():
    """⚠️ CONTROL. A gate that refused every ion would pass the test above and be useless."""
    capable, _ = gn.nlte_ion_capability("Fe", "I")
    assert capable is True


def test_the_gate_refuses_Fe_II_BEFORE_counting_labels():
    """🔴 THE FIX. The capability verdict used to be consulted only INSIDE the `n == 0`
    branch, i.e. only when no line was labelled — so a pool with even one labelled Fe II
    line skipped it entirely and emitted LTE under an NLTE label. It now gates the
    function. `linelist=None` proves the point: the refusal happens before anything is
    read, so no label count can reach it."""
    with pytest.raises(gn.GerberDeckError, match="NO LINE LIST CAN FIX THIS"):
        gn.assert_linelist_supports_nlte(None, 26, "Fe", 4200, 6910, ion="II")


def test_the_gate_does_NOT_short_circuit_a_capable_ion():
    """⚠️ CONTROL, the other direction. Fe I must still reach the label machinery — if it
    were refused here too, the gate would be hiding the RYA-764 check rather than adding
    to it. It gets past the capability step and fails on the absent linelist instead."""
    with pytest.raises(gn.GerberDeckError, match="cannot read NLTE labels"):
        gn.assert_linelist_supports_nlte(None, 26, "Fe", 4200, 6910, ion="I")


def test_the_hole_this_closed_was_REACHABLE_not_theoretical():
    """The refusal only fired for RYA-1135's pool because none of its 9 lines is labelled.
    The window holds 854 labelled Fe II lines, so a differently-selected pool would have
    passed every pre-existing guard. That count is recorded, not recalled."""
    v = json.loads((AUDIT / "verdict.json").read_text())
    h = v["🔴 phantom_hole_closed"]
    assert h["labelled_Fe_II_lines_in_window"] == 854
    assert h["pooled_lines_labelled"] == 0
    assert h["was_reachable"] is True


def test_the_Fe_II_mean3D_product_is_labelled_LTE_and_never_NLTE():
    """RYA-1050's rule. The <3D>-mean atmosphere is ion-agnostic, so Fe II synthesised on
    it in LTE is a genuine <3D>-LTE number — but it is LTE and must say so."""
    v = json.loads((AUDIT / "verdict.json").read_text())
    p = v["product"]
    assert p["scale"] == "<3D>-LTE"
    assert "NLTE" not in p["treatment"].replace("mean3D-LTE", "")
    assert p["treatment"] == "synth-mean3D-LTE-gerber-stagger"
    assert "departures WITHHELD" in p["provenance"]


def test_the_product_carries_the_RYA_1099_over_correction_caveat():
    """Frontier/experimental, never adopted (RYA-1123 Rule 6). The <3D>-mean route inflates
    A(Fe); this number inherits that and must be reported with it."""
    v = json.loads((AUDIT / "verdict.json").read_text())
    c = v["product"]["⚠️ rya1099_over_correction_caveat"]
    assert "+0.08" in c and "0.11" in c
    assert v["product"]["status"] == "EXPERIMENTAL_FRONTIER_NOT_ADOPTED"


def test_model6_was_refused_and_the_refusal_is_recorded_verbatim():
    v = json.loads((AUDIT / "verdict.json").read_text())
    m6 = v["model6_mean3D_NLTE"]
    assert m6["emitted"] is False
    assert "NOT ONE of the 9 pooled lines carries an NLTE label" in m6["refusal"]
