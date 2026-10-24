"""RYA-1250: NLTE level labels beside the VALD-built synthesis lists (Mg, Al, Fe K).

The sidecar is generated on Sirius (`scripts/rya1250_nlte_label_sidecars.py`); these tests
pin what it must say and that only the Gerber NLTE cells of the labelled element see it.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline import nlte_label_sidecar as nls

ROOT = Path(__file__).resolve().parents[1]
LISTS = {k: ROOT / f"data/linelists/ispec_{k}/atomic_lines.tsv"
         for k in ("ir_9200_13000", "h_15007_17494", "k_19452_24846")}


def _side(k):
    return pd.read_csv(nls.sidecar_path(LISTS[k]), sep="\t")


@pytest.mark.parametrize("k", LISTS)
def test_every_sidecar_row_is_a_row_of_its_list_and_both_endpoints_are_labelled(k):
    s = _side(k)
    d = pd.read_csv(LISTS[k], sep="\t", low_memory=False, usecols=nls.KEY)
    m = s.merge(d.drop_duplicates(), on=nls.KEY, how="left", indicator=True)
    assert (m._merge == "both").all()
    assert (s.nlte == "T").all() and (s.nlte_level_low > 0).all() and (s.nlte_level_up > 0).all()
    assert not s.duplicated(nls.KEY).any()


@pytest.mark.parametrize("wave,low,up", [
    (12083.346, "3d1D", "4f3F*"),     # PR17: 3d1D2 - 4f3F3, the weak intercombination line
    (12083.649, "3d1D", "4f1F*"),     # the strong 4f1F line, 0.3 A away
    (9255.778, "3d1D", "5f1F*"),
    (10811.084, "3d3D", "5f3F*"),
])
def test_mg_ir_lines_reach_the_right_term_of_the_term_resolved_atom(wave, low, up):
    r = _side("ir_9200_13000").set_index("wave_A").loc[wave]
    assert (r.nlte_label_low, r.nlte_label_up) == (low, up)


def test_finished_fe_ir_h_labels_are_not_extended_here():
    """RYA-1094 labelled Fe in the IR/H lists and published from them; the sidecar does not
    add Fe there (the extra 12 IR / 296 H lines are Ryan's call)."""
    for k in ("ir_9200_13000", "h_15007_17494"):
        assert not _side(k).element.str.startswith("Fe").any()
    assert _side("k_19452_24846").element.eq("Fe 1").sum() == 695


def test_mg_ii_is_not_labelled_against_a_ground_level_reservoir():
    for k in LISTS:
        assert not _side(k).element.eq("Mg 2").any()


@pytest.mark.parametrize("des,j,expect", [
    ("3s.7p 3P*", 1.0, "7p3P*"),
    ("3s.3p 3P*", 1.0, "3p3P1*"),
    ("3d7.(4F).4d f3D", 2.0, "f3D2"),
    ("3s2.4s 2S", 0.5, "3s2.4s2S"),
])
def test_label_candidates_spell_each_atoms_convention(des, j, expect):
    assert expect in nls.label_candidates(des, j)


def test_apply_merges_unlabelled_rows_keeps_labelled_ones_and_refuses_strangers(tmp_path):
    dt = [("element", "U8"), ("wave_A", "f8"), ("lower_state_eV", "f8"), ("lower_j", "f8"),
          ("upper_j", "f8"), ("nlte", "U1"), ("nlte_level_low", "i4"), ("nlte_level_up", "i4"),
          ("nlte_label_low", "U20"), ("nlte_label_up", "U20")]
    ll = np.array([("Mg 1", 15024.997, 5.1078, 1., 2., "F", 0, 0, "none", "none"),
                   ("Fe 1", 15025.000, 3.0, 1., 2., "T", 7, 9, "a5D1", "z5F2*")], dtype=dt)
    lst = tmp_path / "atomic_lines.tsv"
    lst.write_text("")
    rows = pd.DataFrame([["Mg 1", 15024.997, 5.1078, 1., 2., "T", 6, 9, "4s3S", "4p3P*", "x"],
                         ["Fe 1", 15025.000, 3.0, 1., 2., "T", 1, 1, "bad", "bad", "x"]],
                        columns=nls.KEY + nls.FIELDS + ["route"])
    rows.to_csv(nls.sidecar_path(lst), sep="\t", index=False)
    out, n = nls.apply_sidecar(ll.copy(), lst)
    assert n == 1 and out["nlte_label_up"][0] == "4p3P*" and out["nlte_label_low"][1] == "a5D1"
    rows.iloc[:1].assign(wave_A=15030.0).to_csv(nls.sidecar_path(lst), sep="\t", index=False)
    with pytest.raises(ValueError, match="matches no row"):
        nls.apply_sidecar(ll.copy(), lst)


def test_only_the_labelled_elements_gerber_nlte_cells_fingerprint_the_sidecar():
    from types import SimpleNamespace
    from pipeline import run_matrix as rm
    lst = LISTS["h_15007_17494"]
    assert nls.fingerprint_rows(lst, "Mg", 15000, 17500) is not None
    assert nls.fingerprint_rows(lst, "Fe", 15000, 17500) is None
    assert nls.fingerprint_rows(lst, "Si", 15000, 17500) is None
    resolved = SimpleNamespace(method="synthesis", steps=[])
    kinds = {}
    for deck in ("ts-lte", "gerber-1d-lte", "gerber-nlte"):
        d = rm.RunDescriptor(element="Mg", ion="I", instrument="crires_plus",
                             holding="solar_crires_plus_h_rya1094", lo_A=15007.11,
                             hi_A=17493.69, engine_deck=deck, method="synthesis", pool="reference")
        kinds[deck] = {r["kind"] for r in rm.input_fingerprints(d, resolved, manifest_path=None)}
    assert "nlte_labels" in kinds["gerber-nlte"]
    assert "nlte_labels" not in kinds["ts-lte"] | kinds["gerber-1d-lte"]
