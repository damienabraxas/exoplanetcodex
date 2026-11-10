"""RYA-1350: the checklist reads our Asplund line-set offset with the settled flux-vs-disk
ruling (RYA-1311), so no run reopens it."""
from pipeline import coverage_report as cr


def _p(t, a, sel="SET-SI_AGSS21"):
    return {"selector": sel, "treatment": t, "A": a}


def test_an_offset_inside_the_flux_geometry_is_expected():
    cmp = cr.asplund_comparison([_p("ENGINE-A-3DNLTE", 7.57), _p("ENGINE-A-3DNLTE", 7.58),
                                 _p("1D-LTE", 7.52), _p("ENGINE-A", 7.49, "GRADED")], "Si")
    assert cmp["asplund"] == 7.51
    rows = {r["treatment"]: r for r in cmp["rows"]}
    assert set(rows) == {"ENGINE-A-3DNLTE", "1D-LTE"}            # line-set products only
    assert rows["ENGINE-A-3DNLTE"]["offset_dex"] == 0.065
    assert rows["ENGINE-A-3DNLTE"]["reading"].startswith("expected")
    assert "do not re-measure" in rows["1D-LTE"]["reading"]


def test_an_offset_beyond_the_geometry_points_at_the_fit():
    row = cr.asplund_comparison([_p("1D-LTE", 7.80)], "Si")["rows"][0]
    assert "look at the fit" in row["reading"]


def test_no_line_set_no_comparison():
    assert cr.asplund_comparison([_p("1D-LTE", 7.5, "GRADED")], "Si") is None
    assert cr.asplund_comparison([_p("1D-LTE", 7.5)], "Xx") is None


def test_the_checklist_shows_the_mu1_measurement_or_owes_it(tmp_path):
    import json
    asp = cr.asplund_comparison([_p("1D-LTE", 7.52)], "Si")
    cov = {"star": "solar", "element": "Si", "feed_version": "x", "counts": {}, "done": True,
           "checklist": [], "disagreements": [], "rows": [], "asplund": asp}
    assert "owed -- run_pipeline.py --element <El> --asplund-mu1" in cr.render(cov)
    (tmp_path / "asplund_mu1.json").write_text(json.dumps(
        {"mean_A_mu1": 7.50, "n_lines": 9, "mean_minus_A21_marcs_lte": 0.01, "n_paired": 9}))
    cov["asplund_mu1"] = cr._asplund_mu1(tmp_path)
    out = cr.render(cov)
    assert "[x] measured as Asplund does (mu=1, Liege): 7.500 (-0.010) on 9 lines" in out


def test_cno_mu1_pairs_only_the_indicators_a21_shares(tmp_path, monkeypatch):
    import json
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import asplund_mu1 as am
    (tmp_path / "cno_VIS").mkdir()
    (tmp_path / "cno_VIS" / "cno_mu1_per_band.json").write_text(json.dumps({"per_band": [
        {"key": "CH_Gband", "element": "C", "role": "primary", "A_X": 8.45},
        {"key": "CI_5052", "element": "C", "role": "cross_check", "A_X": 8.50},
        {"key": "C2_Swan", "element": "C", "role": "cross_check", "A_X": 8.40,
         "constrained": False},
        {"key": "CN_red", "element": "N", "role": "primary", "A_X": 7.90}]}))
    c = am.collect_cno("C", tmp_path)
    assert c["n_lines"] == 2 and c["n_paired"] == 2          # unconstrained fit is not a value
    a21 = {r["indicator"]: r["A21_marcs_lte"] for r in c["lines"]}
    assert c["mean_minus_A21_marcs_lte"] == round(((8.45 - a21["CH_Gband"])
                                                   + (8.50 - a21["CI_5052"])) / 2, 3)
    n = am.collect_cno("N", tmp_path)
    assert n["n_lines"] == 1 and n["n_paired"] == 0 and n["mean_minus_A21_marcs_lte"] is None
