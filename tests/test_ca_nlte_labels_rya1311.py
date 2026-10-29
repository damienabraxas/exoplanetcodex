"""RYA-1311: Ca NLTE labels and Ca's <3D> Gerber deck.

`atom.ca105b` stores Ca II energies above the Ca I ground and spells levels stage + orbital +
term (+J or 2J), term-averaged above its lowest levels. Before this, the labeller compared
raw energies and never generated that spelling: 0 Ca II lines labelled in any list although
the atom carries 224 Ca II bound-bound transitions, and the 2 near-UV Ca I REFERENCE cells
FAILED with no labelled line in their pool.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from pipeline import nlte_label_sidecar as nls
from pipeline.gerber_nlte import DECKS

ROOT = Path(__file__).resolve().parents[1]
SIDE = {k: ROOT / f"data/linelists/ispec_{k}/nlte_labels.tsv"
        for k in ("nearuv_3000_3780", "ir_9200_13000", "h_15007_17494", "k_19452_24846")}


def _side(k):
    return pd.read_csv(SIDE[k], sep="\t")


@pytest.mark.parametrize("des,j,ion,expect", [
    ("3p6.4s2 1S", 0.0, 1, "14S21S0"),
    ("3p6.4s.4p 3P*", 1.0, 1, "14P3P1"),
    ("3p6.3d 2D", 1.5, 2, "23D2D3"),      # half-integer J is written as 2J
    ("3p6.5p 2P*", 0.5, 2, "25P2P"),      # the term-averaged level, no J
])
def test_staged_candidates_spell_the_ca_atom(des, j, ion, expect):
    assert expect in nls._staged_candidates(des, j, ion)


def test_a_stage_is_matched_from_its_own_ground():
    lv = pd.DataFrame({"index": [1, 2, 3], "ion": [1, 2, 2], "term": ["a", "b", "c"],
                       "energy_eV": [0.0, 6.113, 7.806]})
    s = nls._stage_levels(lv, 2)
    assert s.energy_eV.round(3).tolist() == [0.0, 1.693]
    assert nls._stage_levels(lv, 1).energy_eV.tolist() == [0.0]          # no-op from 0


def test_ca_ii_is_labelled_in_every_vald_list():
    for k in SIDE:
        assert _side(k).element.eq("Ca 2").sum() > 0, k
    r = _side("ir_9200_13000").set_index("wave_A").loc[9854.759]
    assert (r.nlte_label_low, r.nlte_label_up) == ("25P2P", "26S2S")      # 5p 2P1/2 - 6s 2S


@pytest.mark.parametrize("wave", [3624.112, 3630.752, 3630.978, 3644.413, 3644.768, 3644.996])
def test_near_uv_lab_gf_ca_i_lines_are_labelled(wave):
    s = _side("nearuv_3000_3780")
    assert ((s.element == "Ca 1") & ((s.wave_A - wave).abs() < 0.01)).any()


def test_ca_mean3d_deck_is_registered_with_pinned_files():
    d = DECKS["Ca@mean3D"]
    assert d["atom"] == DECKS["Ca"]["atom"] == "atom.ca105b" and d["read_via"] == "direct"
    prov = json.loads((ROOT / "data/nlte_grids/gerber_ts/Ca_gerber2023.prov.json").read_text())
    assert prov["files"]["grid_3d_bin"]["name"] == d["grid"]
    assert prov["files"]["grid_3d_aux"]["name"] == d["aux"]
    assert len(prov["files"]["grid_3d_bin"]["md5"]) == 32


def test_ca105b_reach_is_measured_for_both_stages():
    r = json.loads((ROOT / "data/reference/nlte_atom_reach/atom.ca105b.json").read_text())
    assert {s["ion_stage"]: s["nlte_capable"] for s in r["stages"]} == {1: True, 2: True, 3: False}
