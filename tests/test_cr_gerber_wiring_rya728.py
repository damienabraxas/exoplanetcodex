"""RYA-728 -- Cr's Gerber decks, and the ORCH gaps that kept them unwired since RYA-800/818.

The element is Cr; the rules are general: a level's term comes from the level (never the nearest
atom energy, RYA-818), routes that disagree leave a line unlabelled, a GES sidecar touches only
the cell's own element, and a staged-but-unregistered deck is named as a wiring gap.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import gerber_nlte as G                     # noqa: E402
from pipeline import nlte_label_sidecar as nls            # noqa: E402

TABLE = {("Cr 2", 7): (np.array([3.86441, 6.42099]), np.array(["a4F", "z4F*"]),
                       np.array(["3d4.(3F).4s a4F", "3d4.(5D).4p z4F*"])),
         ("Cr 1", 4): (np.array([3.0, 3.0005]), np.array(["a3F", "b3F"]),
                       np.array(["x a3F", "y b3F"]))}


def test_level_table_returns_the_levels_term():
    assert nls._level_term(TABLE, "Cr 2", 6.42099, 3.5) == "3d4.(5D).4p z4F*"


def test_level_table_refuses_a_level_two_terms_claim_and_never_takes_the_nearest():
    assert nls._level_term(TABLE, "Cr 1", 3.0, 2.0) == ""        # two terms within the key
    assert nls._level_term(TABLE, "Cr 2", 3.90, 3.5) == ""       # no level there: no nearest
    assert nls._level_term(TABLE, "Cr 2", 6.42099, 2.5) == ""    # wrong J


def _ges_like(tmp_path):
    d = tmp_path / "GESv6_atom_hfs_iso.420_920nm"
    d.mkdir(exist_ok=True)
    return d / "atomic_lines.tsv"


def _array():
    return np.array([("Cr 1", 5000.0, 1.0, 2.0, 3.0, "F", 7, 9, "bad", "bad"),
                     ("Fe 1", 5000.0, 1.0, 2.0, 3.0, "T", 1, 2, "keep", "keep")],
                    dtype=[("element", "U8"), ("wave_A", "f8"), ("lower_state_eV", "f8"),
                           ("lower_j", "f8"), ("upper_j", "f8"), ("nlte", "U1"),
                           ("nlte_level_low", "i4"), ("nlte_level_up", "i4"),
                           ("nlte_label_low", "U12"), ("nlte_label_up", "U12")])


def _write_ges_sidecar(monkeypatch, tmp_path):
    side = tmp_path / "side" / "nlte_labels.tsv"
    side.parent.mkdir()
    pd.DataFrame([{"element": "Cr 1", "wave_A": 5000.0, "lower_state_eV": 1.0, "lower_j": 2.0,
                   "upper_j": 3.0, "nlte": "T", "nlte_level_low": 11, "nlte_level_up": 62,
                   "nlte_label_low": "a3F", "nlte_label_up": "y3F*", "route": "term"}]
                 ).to_csv(side, sep="\t", index=False)
    monkeypatch.setattr(nls, "GES_SIDECAR", side)


def test_ges_sidecar_replaces_the_named_elements_labels(monkeypatch, tmp_path):
    _write_ges_sidecar(monkeypatch, tmp_path)
    ll, n = nls.apply_sidecar(_array(), _ges_like(tmp_path), element="Cr")
    assert n == 1
    assert (ll[0]["nlte_label_low"], ll[0]["nlte_level_up"]) == ("a3F", 62)
    assert ll[1]["nlte_label_low"] == "keep"


def test_ges_sidecar_leaves_every_other_cell_byte_identical(monkeypatch, tmp_path):
    """An Fe cell, or a caller naming no element, keeps its cache digest."""
    _write_ges_sidecar(monkeypatch, tmp_path)
    for element in ("Fe", None):
        ll, n = nls.apply_sidecar(_array(), _ges_like(tmp_path), element=element)
        assert n == 0 and ll.tobytes() == _array().tobytes()


def test_sidecar_path_keeps_vald_lists_beside_the_list_and_ges_in_the_repo(tmp_path):
    assert nls.sidecar_path(_ges_like(tmp_path)) == nls.GES_SIDECAR
    assert nls.sidecar_path(tmp_path / "ispec_x" / "atomic_lines.tsv").parent.name == "ispec_x"


def _levels():
    return pd.DataFrame({"index": [4, 1, 2, 3], "element": "Cr", "ion": 1,
                         "term": ["a7S", "a3F", "y3F*", "z9Z"], "g": [1.0, 1.0, 1.0, 1.0],
                         "energy_cm": [0.0, 1.0, 2.0, 3.0], "energy_eV": [0.0, 1.0, 3.0, 4.0],
                         "kind": "term"})


def _one_line_list(tmp_path):
    p = tmp_path / "list.tsv"
    pd.DataFrame([{"element": "Cr 1", "wave_A": 5000.0, "lower_state_eV": 1.0, "upper_state_eV": 3.0,
                   "lower_j": 2.0, "upper_j": 3.0}]).to_csv(p, sep="\t", index=False)
    raw = tmp_path / "raw.txt"
    raw.write_text("'Cr 1', 5000.000, 0.0, 1.000, 2.0, 3.000, 3.0, 0.0\n"
                   "'LS x a3F'\n'LS y y3F*'\n")
    return p, raw


def test_term_join_agreeing_with_the_vald_route_labels_the_line(tmp_path):
    p, raw = _one_line_list(tmp_path)
    rows, st = nls.label_list(p, "Cr", _levels(), [raw], term_join={(5000000, 100, 1): (1, 2)})
    assert st["Cr 1"]["labelled"] == 1 and st["Cr 1"]["conflict"] == 0
    assert (rows.iloc[0].nlte_level_low, rows.iloc[0].nlte_level_up) == (1, 2)


def test_routes_that_disagree_leave_the_line_unlabelled_and_counted(tmp_path):
    p, raw = _one_line_list(tmp_path)
    rows, st = nls.label_list(p, "Cr", _levels(), [raw], term_join={(5000000, 100, 1): (1, 3)})
    assert rows.empty and st["Cr 1"]["conflict"] == 1 and st["Cr 1"]["labelled"] == 0


def test_cr_decks_are_registered_on_the_term_resolved_atom():
    one, three = G.DECKS["Cr"], G.DECKS["Cr@mean3D"]
    assert one["atom"] == three["atom"] == "atom.cr374"
    assert three.get("read_via") == "direct" and "STAGGERmean3D" in three["grid"]
    assert three["atmos"].endswith("sun_avg3d_stagger.mod")


def test_atom_cr374_serves_both_ions():
    assert G.nlte_ion_capability("Cr", "I")[0] is True
    assert G.nlte_ion_capability("Cr@mean3D", "II")[0] is True


def test_a_staged_unregistered_deck_is_named(monkeypatch, tmp_path):
    (tmp_path / "auxData_Xx_MARCS_Jan-01-2025.dat").write_text("")
    (tmp_path / "auxData_Xx_STAGGERmean3D_Jan-01-2025_marcs_names.txt").write_text("")
    monkeypatch.setattr(G, "GT", str(tmp_path))
    assert G.staged_unregistered("Xx") == ["1D", "<3D>"]
    assert G.staged_unregistered("Yy") == []
    (tmp_path / "auxData_Cr_MARCS_Dec-16-2024.dat").write_text("")
    assert G.staged_unregistered("Cr") == []


# -- every element: the per-band literature sweep and the line inventory ----------------------
from types import SimpleNamespace                                  # noqa: E402

from pipeline import element_prepare as ep                         # noqa: E402
from pipeline import run_matrix as rm                              # noqa: E402

_BANDS = [("I", "VIS"), ("I", "H")]
_H_NEED = [{"species": "Cr I", "band": "H"}]


def _sweep(*entries):
    return {"literature_sweep": {"bands": list(entries)}}


def test_an_anchor_only_litscan_owes_every_band():
    st = ep.literature_sweep_status({"entries": [{"value": 5.65}]}, _BANDS, [])
    assert st["recorded"] == 0 and {o["band"] for o in st["owed"]} == {"VIS", "H"}


def test_an_entry_must_name_a_search_and_a_result():
    empty = {"ion": "I", "band": "VIS", "searched": [], "result": "none"}
    done = {"ion": "I", "band": "H", "searched": ["Smith+2021"], "result": "no lab gf", "diagnostic": "x"}
    st = ep.literature_sweep_status(_sweep(empty, done), _BANDS, _H_NEED)
    assert [o["band"] for o in st["owed"]] == ["VIS"]


def test_a_no_gf_band_is_owed_until_it_is_measured():
    searched = {"ion": "I", "band": "H", "searched": ["Smith+2021"], "result": "no lab gf"}
    st = ep.literature_sweep_status(_sweep(searched), [("I", "H")], _H_NEED)
    assert st["owed"][0]["why"].endswith("no measured diagnostic")
    st = ep.literature_sweep_status(_sweep({**searched, "diagnostic": "H: 5.6, n=4"}), [("I", "H")], _H_NEED)
    assert st["owed"] == []


def _report(owed):
    return dict(ready=False, scoped_readiness=True, nlte_missing=[], adopt=[], missing_lines=[],
                cull_candidates=[], holdings=[], blocking={"literature_sweep": len(owed)},
                literature={"litscan": True, "missing_bib_keys": []},
                literature_sweep={"recorded": 0, "owed": owed})


def _cell(ion, band):
    return SimpleNamespace(element="Cr", ion=ion, band=band, lo_A=1.0, hi_A=2.0, holding="h")


def test_an_owed_sweep_holds_only_its_own_ion_and_band():
    report = _report([{"ion": "I", "band": "H", "why": "no recorded search"}])
    assert "literature sweep owed for Cr I H" in rm.preparation_hold_for_cell(report, _cell("I", "H"))
    assert rm.preparation_hold_for_cell(report, _cell("I", "VIS")) == ""
    assert rm.preparation_hold_for_cell(report, _cell("II", "H")) == ""


def test_write_only_rows_leaves_every_other_line_byte_identical(monkeypatch, tmp_path):
    import scripts.rya1268_adopt_lab_gf as adopt
    p = tmp_path / "canon.csv"
    p.write_text("line_id,species,log_gf\na,Ti I,-1.2300000000000002\nb,Cr I,-1.0\nc,Ca I,1e-01\n")
    monkeypatch.setattr(adopt, "CANON", p)
    cg = pd.read_csv(p)
    cg.at[1, "log_gf"] = -2.0
    adopt.write_only_rows(cg, [1])
    assert p.read_text() == "line_id,species,log_gf\na,Ti I,-1.2300000000000002\nb,Cr I,-2.0\nc,Ca I,1e-01\n"


def _ir_list(tmp_path):
    p = tmp_path / "atomic_lines.tsv"
    p.write_text("element\twave_A\tloggf\nFe 1\t9500.000\t-1.0\nCr 1\t9571.813\t-2.3\n")
    return p


def test_a_scoped_addition_counts_only_in_its_own_elements_identity(tmp_path):
    p = _ir_list(tmp_path)
    (tmp_path / "element_scoped_additions.tsv").write_text("element\tion_element\twave_A\nCr\tCr 1\t9571.813\n")
    before = tmp_path / "before"
    before.mkdir()
    (before / "atomic_lines.tsv").write_text("element\twave_A\tloggf\nFe 1\t9500.000\t-1.0\n")
    fp = lambda path, el: rm._rows_fingerprint(path, "wave_A", 9000, 9900, sep="\t", for_element=el)  # noqa: E731
    assert fp(p, "Fe") == fp(before / "atomic_lines.tsv", "Fe")
    assert fp(p, "Cr") != fp(before / "atomic_lines.tsv", "Cr")


def test_a_row_outside_the_scoped_list_still_moves_every_element(tmp_path):
    p = _ir_list(tmp_path)
    fp = lambda el: rm._rows_fingerprint(p, "wave_A", 9000, 9900, sep="\t", for_element=el)  # noqa: E731
    first = fp("Fe")
    p.write_text(p.read_text() + "Cr 1\t9600.000\t-2.0\n")
    assert fp("Fe") != first


def test_appending_a_canonical_row_leaves_every_other_line_byte_identical(monkeypatch, tmp_path):
    from pipeline import published_line_add as pla
    p = tmp_path / "canon.csv"
    p.write_text("line_id,species,log_gf\na,Ti I,-1.2300000000000002\n")
    monkeypatch.setattr(pla, "CANONICAL", p)
    pla._append_canonical_rows([{"line_id": "b", "species": "Cr I", "log_gf": -2.21}], ["line_id", "species", "log_gf"])
    assert p.read_text() == "line_id,species,log_gf\na,Ti I,-1.2300000000000002\nb,Cr I,-2.21\n"


def test_a_transition_held_twice_never_looks_unique_after_one_row_is_adopted():
    import scripts.rya1268_adopt_lab_gf as adopt
    cg = pd.DataFrame({"species": ["Cr I", "Cr I"], "wavelength_air_A": [4351.762, 4351.7879],
                       "excitation_potential_eV": [1.03, 1.03], "log_gf": [-0.48, -0.1585]})
    open_rows = {"Cr I": cg.iloc[[1]]}
    assert adopt.match(cg, open_rows, "Cr I", 4351.7625, 1.03) is None


def test_the_wider_window_needs_the_gf_to_agree_and_an_open_row():
    import scripts.rya1268_adopt_lab_gf as adopt
    cg = pd.DataFrame({"species": ["Cr I"], "wavelength_air_A": [5224.925],
                       "excitation_potential_eV": [3.45], "log_gf": [0.07]})
    assert adopt.match_same_gf(cg, {"Cr I": cg}, "Cr I", 5224.97, 3.45, 0.07) == 0
    assert adopt.match_same_gf(cg, {"Cr I": cg}, "Cr I", 5224.97, 3.45, -1.0) is None
    assert adopt.match_same_gf(cg, {"Cr I": cg.iloc[[]]}, "Cr I", 5224.97, 3.45, 0.07) is None


def test_gurell_row_with_a_lower_level_reaching_no_known_level_takes_valds():
    import scripts.rya728_parse_gurell2010 as g
    vald = pd.DataFrame({"wavelength_air_A": [3617.273], "ep_eV": [4.9274], "ep_upper_eV": [8.3547], "j_lower": [2.5]})
    t2 = pd.DataFrame({"air": [3617.274], "A_1e7": [0.48], "unc_pct": [17.0]})
    t3 = pd.DataFrame({"vac": [3618.306], "sigma": [27637.24], "elow": [32603.40], "loggf": [-1.423]})
    out, counts = g.build(t2, t3, vald)
    assert counts["lower_level_from_vald"] == 1
    assert out.ep_eV.iloc[0] == 4.9274 and out.ep_upper_eV.iloc[0] == 8.3547
