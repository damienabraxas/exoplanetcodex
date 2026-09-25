"""RYA-1220 WP A -- the molecular literature matrix must keep two facts pinned.

1. The 1D->3D model-form term is MOLECULE-SPECIFIC. Applying one uniform molecular
   offset is wrong by up to 0.18 dex across the six species.
2. Two of our molecular N diagnostics are NOT on bands the reference analysis uses.
   If that ever silently becomes 'YES', someone has widened a window.
"""
from __future__ import annotations

import pathlib

import pandas as pd
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
REF = ROOT / "data/reference/molecular_cno_literature_rya1220/amarsi2021_molecular_reference_by_species.csv"
MATCH = ROOT / "data/reference/molecular_cno_literature_rya1220/our_diagnostics_vs_amarsi2021.csv"

pytestmark = pytest.mark.skipif(not REF.exists(), reason="literature matrix not generated")


def test_all_six_reference_molecules_are_present():
    d = pd.read_csv(REF)
    assert sorted(d.molecule) == ["12C16O", "C2", "CH", "CN", "NH", "OH"]
    assert d.n_lines.sum() == 408, "the published used-line count must stay 408"


def test_model_form_term_is_molecule_specific():
    d = pd.read_csv(REF).set_index("molecule")
    ch = d.loc["CH", "model_form_3D_minus_MARCS"]
    co = d.loc["12C16O", "model_form_3D_minus_MARCS"]
    assert ch > 0 and co < 0, "CH and CO must fall on opposite sides of zero"
    assert abs(ch - co) > 0.15, (
        "the spread across molecules is what forbids a single uniform molecular 3D offset")


def test_line_to_line_scatter_is_small_in_the_reference():
    """Our molecular routes run red_chi2 33-213; the published scatter is <=0.032 dex."""
    d = pd.read_csv(REF)
    assert d.line_to_line_sd_3D.max() <= 0.035


def test_cn_red_and_nh_uv_are_outside_the_reference_set():
    d = pd.read_csv(MATCH).set_index("our_key")
    assert d.loc["CN_red", "in_amarsi2021_used_set"] == "NO"
    assert d.loc["NH_AX", "in_amarsi2021_used_set"] == "NO"


def test_the_ir_cn_routes_are_inside_the_reference_set():
    d = pd.read_csv(MATCH).set_index("our_key")
    assert d.loc["CN_AX_J", "in_amarsi2021_used_set"] == "YES"
    assert d.loc["CN_AX_IR", "in_amarsi2021_used_set"] == "YES"
