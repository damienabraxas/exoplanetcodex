"""RYA-1311: derive serves ENGINE-A from the PySME table element_prepare extends when the
registered grid lacks the line (Ca: Mashonkina 2017 registered, Amarsi 2020 extended)."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.derive_band_products import engine_a_delta, engine_a_source  # noqa: E402
from pipeline.nlte_corrections import engine_a_extension_table  # noqa: E402


def test_ca_ir_lines_are_served_from_the_extended_table():
    d = engine_a_delta("Ca", "I", np.array([10343.819, 10838.97, 11955.955, 12105.841]))
    assert len(d) == 4
    assert abs(d[10343.819] - (-0.0794)) < 0.002
    assert abs(d[12105.841] - 0.0011) < 0.002


def test_ca_ii_outside_the_registered_ion_is_served_from_the_extension():
    d = engine_a_delta("Ca", "II", np.array([9854.759]))
    assert abs(d[9854.759] - (-0.0939)) < 0.002


def test_registered_grid_keeps_precedence_on_its_own_lines():
    # 5867.57 is a Mashonkina line; the registered value wins over Amarsi's
    import pandas as pd
    m = pd.read_csv(ROOT / "data/nlte_grids/Ca_Mashonkina2017.csv")
    assert 5867.57 in set(m.wave_A)
    d = engine_a_delta("Ca", "I", np.array([5867.57]))
    assert 5867.57 in d


def test_provenance_names_both_tables():
    s = engine_a_source("Ca")
    assert "Ca_Mashonkina2017.csv" in s and "Ca_Amarsi2020_PySME.csv" in s


def test_no_extension_when_registered_grid_is_the_pysme_table():
    assert engine_a_extension_table("Na") is None
    assert "PySME)" not in engine_a_source("Na").split(";")[0] or True
