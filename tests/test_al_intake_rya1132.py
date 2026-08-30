import json
import pandas as pd

from scripts.build_al_intake_rya1132 import OUT, build


def test_manifest_is_reproducible_and_never_derives_abundance(tmp_path):
    s = build(tmp_path)
    m = pd.read_csv(tmp_path / "al_line_manifest.csv")
    assert s["abundances_generated"] is False
    assert s["canonical_gf_mutated"] is False
    assert len(m) >= 500
    assert set(m.species) == {"Al I", "Al II"}
    assert not m.canonical_line_id.duplicated().any()
    assert {"FUV", "NUV", "near-UV", "VIS", "red-optical", "NIR", "H", "K"} <= set(m.band)
    for committed in OUT.iterdir():
        assert (tmp_path / committed.name).read_bytes() == committed.read_bytes()


def test_burheim_and_hfs_guards_are_explicit(tmp_path):
    build(tmp_path)
    m = pd.read_csv(tmp_path / "al_line_manifest.csv")
    vis = m[(m.wavelength_air - 6696.015).abs() < .03]
    wrong = m[(m.wavelength_air - 6696.185).abs() < .03]
    assert len(vis) == 1 and "BURHEIM" in vis.iloc[0].gf_source.upper()
    assert len(wrong) == 1 and "BURHEIM" not in wrong.iloc[0].gf_source.upper()
    blend = pd.read_csv(tmp_path / "conflict_ledger.csv")
    assert blend.rejection_problem_code.str.contains("BURHEIM_STRONG_COMPONENT").any()
    held = m[(m.wavelength_air - 11254.924).abs() < .08]
    assert len(held) == 1
    assert "BURHEIM" not in held.iloc[0].gf_source.upper()
    assert abs(held.iloc[0].loggf_adopted - 0.3538) < .001
    assert m[m.HFS_status == "COMPONENT_SUM_VERIFIED"].component_or_total.eq("TOTAL_TRANSITION_GF").all()


def test_gate_fails_closed(tmp_path):
    build(tmp_path)
    v=json.loads((tmp_path / "intake_verdict.json").read_text())
    assert v["UV"] == "CROSSMATCH_REVIEW"
    assert v["IR"] == "BLOCKED_PIPELINE_COVERAGE"
    assert v["measurement_unblocked"] is False


def test_web_followup_is_conservative(tmp_path):
    build(tmp_path)
    d = pd.read_csv(tmp_path / "web_source_followup.csv")
    assert {"Vujnovic et al. 2002", "Johnson, Smith & Parkinson 1986"} <= set(d.source)
    assert not d.disposition.str.contains("PROMOTE_NOW").any()
    assert d[(d.wavelength_A - 1670.78861).abs() < .001].iloc[0].disposition == "WAVELENGTH_ONLY_NOT_GF_LAB"
    assert d.evidence_class.eq("NEGATIVE_RESULT").any()
    bib = pd.read_csv(tmp_path / "source_bibliography.csv")
    actionable = bib[bib.source_id.isin(["Vujnovic2002", "Trabert1999", "Johnson1986"])]
    assert len(actionable) == 3
    assert actionable.article_url.str.startswith("https://").all()
    assert bib.loc[bib.source_id.eq("Vujnovic2002"), "download_url"].str.endswith(".pdf").all()
