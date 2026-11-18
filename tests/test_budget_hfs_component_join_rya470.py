"""RYA-470: the budget's lambda+EP gf join must find a hyperfine line reported at one of
its components (V I 5737.031) on the canonical centroid row (5737.062) -- the row the
synthesis priced -- and must still refuse anything that is not unambiguously that line."""
from pathlib import Path

import pandas as pd

from pipeline.gf_resolver import EP_CLUSTER_TOL, EPTOL, HFS_GAP, WTOL

SRC = (Path(__file__).resolve().parents[1] / "scripts" / "rya1230_cno_budget.py").read_text()


def _join(g, w, ep):
    """The assembler's join, transcribed from its two predicates."""
    m = g[((g.wavelength_air_A - w).abs() <= WTOL) & ((g.excitation_potential_eV - ep).abs() <= EPTOL)]
    if m.empty:
        h = g[((g.wavelength_air_A - w).abs() <= HFS_GAP)
              & ((g.excitation_potential_eV - ep).abs() <= EP_CLUSTER_TOL)
              & (pd.to_numeric(g.hfs_n_components, errors="coerce") > 1)]
        if len(h) == 1:
            m = h
    return m


def test_the_assembler_carries_the_hfs_fallback():
    assert "HFS_GAP" in SRC and "EP_CLUSTER_TOL" in SRC and "hfs_n_components" in SRC


def test_a_component_joins_its_centroid():
    g = pd.DataFrame({"wavelength_air_A": [5737.062], "excitation_potential_eV": [1.064],
                      "hfs_n_components": [5]})
    assert len(_join(g, 5737.031, 1.064)) == 1
    assert 5737.062 - 5737.031 > WTOL          # the strict join alone refuses it


def test_a_single_component_row_is_never_borrowed():
    g = pd.DataFrame({"wavelength_air_A": [5737.062], "excitation_potential_eV": [1.064],
                      "hfs_n_components": [1]})
    assert _join(g, 5737.031, 1.064).empty


def test_another_lower_level_is_never_borrowed():
    g = pd.DataFrame({"wavelength_air_A": [5737.062], "excitation_potential_eV": [1.200],
                      "hfs_n_components": [5]})
    assert _join(g, 5737.031, 1.064).empty


def test_two_candidates_refuse():
    g = pd.DataFrame({"wavelength_air_A": [5736.99, 5737.08], "excitation_potential_eV": [1.064, 1.064],
                      "hfs_n_components": [3, 4]})
    assert _join(g, 5737.031, 1.064).empty
