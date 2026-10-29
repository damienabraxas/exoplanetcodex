"""RYA-1315 -- Ti ENGINE-A extension: levels from the grid's shipped label table, and one
unlabelled line held instead of aborting the batch."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import pysme_nlte as pn                     # noqa: E402

# Four rows of the Mallinson 2024 label_Ti.txt, verbatim fixed-width layout.
LABELS = """# level index, level species, level configuration, level term, level J, level energy (NIST)/eV
   1 Ti 1                        3d2.4s2        a3F   2.0          0.000000000
  44 Ti 1           3d2.(3F).4s.4p.(3P*)       z5D*   2.0          2.296720100
 120 Ti 1                 3d2.4s.(4F).5s        e5F   3.0          4.475600000
 600 Ti 2                     3d2.(3F).4s       a4F   1.5          0.000000000
"""


def test_label_table_gives_the_rya545_frozen_labels(monkeypatch, tmp_path):
    (tmp_path / "label_Ti.txt").write_text(LABELS)
    monkeypatch.setattr(pn, "_GRID_DIR", tmp_path)
    tl, tu, jl, ju, note = pn.labels_by_selection_rules("Ti", 2.297, 4.476, ion=1)
    want = pn.NLTE_LINES["Ti"][0]                           # 5689.460, frozen by RYA-545
    assert (tl, tu, jl, ju) == (want[6], want[7], want[3], want[5]) and note == ""


def test_label_table_keeps_each_ion_on_its_own_levels(monkeypatch, tmp_path):
    (tmp_path / "label_Ti.txt").write_text(LABELS)
    monkeypatch.setattr(pn, "_GRID_DIR", tmp_path)
    g = pn._label_levels("Ti")
    assert [pn.decode_grid_label(s) for s in g["spec"]] == ["Ti 1"] * 3 + ["Ti 2"]
    assert len(pn._ion_level_indices(g, "Ti", 2)) == 1


def test_an_unlabelled_line_is_held_not_fatal(monkeypatch):
    sys.path.insert(0, str(ROOT / "scripts"))
    import extend_nlte_grid as ext
    monkeypatch.setattr(ext, "ELEMENT", "Ti")
    monkeypatch.setattr(ext, "ION", 1)
    monkeypatch.setattr(ext, "UNLABELLED", {})
    monkeypatch.setattr(ext, "atomic", lambda wl: (-1.0, 1.0, 3.0, 0.0))

    def labels(el, ep, eup, ion=1):
        if eup == 3.0 and ep == 1.0 and labels.n == 0:
            labels.n += 1
            raise ValueError("Ti: no LS-allowed upper level")
        return "lo", "up", 1.0, 2.0, ""
    labels.n = 0
    monkeypatch.setattr(ext.pn, "labels_by_selection_rules", labels)
    out = ext.build_lines([5000.0, 5100.0])
    assert [l[0] for l in out] == [5100.0]
    assert "no LS-allowed" in ext.UNLABELLED[5000.0]


def test_a_near_uv_line_reads_the_near_uv_list(monkeypatch, tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    import extend_nlte_grid as ext
    nuv = tmp_path / "atomic_lines.tsv"
    nuv.write_text("element\twave_A\tlower_state_eV\twaals\trad\tstark\tupper_state_eV\tlower_j\tupper_j\n"
                   "Ti 2\t3008.311\t1.5718\t-7.85\t8.32\t-6.44\t5.69\t4.5\t5.5\n")
    monkeypatch.setattr(ext, "NEARUV", (3000.0, 3780.0, nuv))
    monkeypatch.setattr(ext, "ELEMENT", "Ti")
    monkeypatch.setattr(ext, "ION", 2)
    monkeypatch.setattr(ext, "_LL", {})
    assert list(ext._synth_rows(3008.311).wave_A) == [3008.311]
