"""RYA-1250: the Gerber interpolator's scratch directory is private to each process.

A host-wide `/tmp/rya798_gerber` let parallel derive processes delete or read each
other's `{node}_{element}_coef.dat` (the Mg IAG red-optical xi- leg lost its departures
under `--jobs 4`, and an abundance-axis deck could read a sibling's trial abundance).
"""
import inspect
import multiprocessing as mp

from pipeline import gerber_nlte as gn


def _child(q):
    q.put(gn.default_workdir())


def test_two_processes_get_two_workdirs():
    q = mp.get_context("spawn").Queue()
    p = mp.get_context("spawn").Process(target=_child, args=(q,))
    p.start()
    other = q.get(timeout=60)
    p.join(60)
    assert other != gn.default_workdir()


def test_for_node_has_no_shared_default_workdir():
    default = inspect.signature(gn.for_node).parameters["workdir"].default
    assert default is None
    assert "rya798_gerber" in gn.default_workdir() and str(__import__("os").getpid()) in gn.default_workdir()


def test_mg_ir_nist_grades_only_where_nist_is_the_rows_value():
    """RYA-1250 IR sweep: a NIST grade is written beside an unchanged gf, never instead of one."""
    import pandas as pd
    from pathlib import Path
    from pipeline.gf_grades import is_gf_graded
    root = Path(__file__).resolve().parents[1]
    c = pd.read_csv(root / "data/linelists/canonical_gf.csv", low_memory=False)
    s = pd.read_csv(root / "data/reference/nist/mg_ir_rya1250.csv")
    g = c[c.adjudication_status.eq("adopted_rya1250_nist_ir_grade")]
    assert len(g) == 49 and (g.gf_tier == "NIST-C+").all()
    assert set(g.species) <= {"Mg I", "Mg II"} and g.wavelength_air_A.between(9200, 24850).all()
    r = c[(c.species == "Mg I") & (c.wavelength_air_A.round(3) == 21458.863)].iloc[0]
    assert r.nist_grade == "B+" and r.log_gf == -1.319
    assert s[(s.species == "Mg I") & (s.wavelength_air_A.round(2) == 21458.86)].log_gf.iloc[0] == -1.319
    k = c[(c.species.isin(["Mg I", "Mg II"])) & c.wavelength_air_A.between(19450, 24850)]
    assert is_gf_graded(k).sum() == 23


def test_k_band_has_its_pr17_lab_line():
    """4d1D2-5f1F3 (PR17 24572.92 vac): the K band's one Reference Grade lab line."""
    import pandas as pd
    from pathlib import Path
    c = pd.read_csv(Path(__file__).resolve().parents[1] / "data/linelists/canonical_gf.csv",
                    low_memory=False).set_index("line_id")
    r = c.loc["gf_179143"]
    assert (r.gf_tier, r.lab_source_tag, r.log_gf, r.gf_sigma_dex) == ("LAB", "PR2017_exp", -0.391, 0.06)
