"""RYA-731 — an element with NO registered Engine-A source (Ni) gets no Engine-A axis, not a crash.

The Ni IR diagnostic cells measured their 1D-LTE products and then died in
`engine_a_model("Ni")`, which was called eagerly as an argument to `_emit_departure` even
though no Engine-A product was going to be written. The guard is kept where it matters: a
registered element with no mappable axis still raises, and `_emit_departure` still refuses to
write an Engine-A product whose axis is None.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


def test_an_unregistered_element_has_no_engine_a_axis():
    import derive_band_products as dbp
    from config.constants import NLTE_CORRECTION_ELEMENTS

    assert "Ni" not in NLTE_CORRECTION_ELEMENTS
    assert dbp.engine_a_model("Ni") is None
    assert dbp.engine_a_source("Ni") == "no registered Engine-A source"


def test_registered_elements_keep_their_axes():
    import derive_band_products as dbp

    assert dbp.engine_a_model("Fe") == "bergemann"
    assert dbp.engine_a_model("Al") == "amarsi"


def test_a_registered_element_without_a_mappable_axis_still_raises(monkeypatch):
    import derive_band_products as dbp
    import config.constants as cc

    monkeypatch.setitem(cc.NLTE_CORRECTION_ELEMENTS, "Zz", {"grid": "somebody", "ref": "2099"})
    with pytest.raises(ValueError, match="no RYA-906 model-axis mapping for Zz"):
        dbp.engine_a_model("Zz")


def test_the_writer_refuses_an_engine_a_product_without_an_axis():
    src = (ROOT / "scripts" / "derive_band_products.py").read_text(encoding="utf-8")
    body = src[src.index("    def _emit_departure("):]
    body = body[:body.index("\n    def ", 10)] if "\n    def " in body[10:] else body
    assert 'if model is None and treatment.startswith("ENGINE-A")' in body
