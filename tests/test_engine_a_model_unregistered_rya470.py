"""RYA-470: an element with no registered Engine-A (NLTE) source gets model axis 'none' --
it must not raise after the 1D-LTE products are written (every solar V cell FAILED that way)."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _dbp():
    spec = importlib.util.spec_from_file_location("dbp", ROOT / "scripts" / "derive_band_products.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_unregistered_element_has_model_axis_none_and_matching_source():
    from config.constants import NLTE_CORRECTION_ELEMENTS
    from pipeline.treatment_axes import MODELS
    assert "V" not in NLTE_CORRECTION_ELEMENTS
    m = _dbp()
    assert m.engine_a_model("V") == "none" and "none" in MODELS
    assert m.engine_a_source("V") == "no registered Engine-A source"


def test_every_registered_element_still_maps_to_a_real_model():
    from config.constants import NLTE_CORRECTION_ELEMENTS
    m = _dbp()
    for el in NLTE_CORRECTION_ELEMENTS:
        try:
            model = m.engine_a_model(el)
        except ValueError:
            continue                      # a registered-but-unmapped source still refuses
        assert model != "none", el
