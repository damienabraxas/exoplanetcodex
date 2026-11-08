"""RYA-1311: H/K-band Ca I gf from Civis et al. 2013 (A&A 554, A24, Table 5)."""
import csv
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.gf_grades import is_gf_graded  # noqa: E402


def _canon(species, w):
    rows = [r for r in csv.DictReader(open(ROOT / "data/linelists/canonical_gf.csv"))
            if r["species"] == species and abs(float(r["wavelength_air_A"]) - w) < 0.01]
    assert len(rows) == 1
    return rows[0]


def test_adopted_rows_are_codex_graded_but_not_lab():
    rows = pd.DataFrame([_canon("Ca I", 16150.7624), _canon("Ca I", 22651.1773)])
    assert is_gf_graded(rows).all()
    assert (rows.gf_tier == "C2013-QDT").all()
    assert float(rows.gf_sigma_dex.iloc[0]) == 0.11


def test_the_h_band_5p_5d_lines_take_the_weaker_civis_values():
    assert float(_canon("Ca I", 16150.7624)["log_gf"]) == -0.180
    assert float(_canon("Ca I", 16157.3637)["log_gf"]) == -0.161


def test_a_multiplet_our_rows_only_partly_carry_stays_on_kurucz():
    au = pd.read_csv(ROOT / "data/audit/rya1311_civis2013/adoption.csv")
    amb = set(au[au.action == "ambiguous"].canonical_A.round(2))
    assert {20494.04, 20495.48, 20497.45} <= amb
    assert _canon("Ca I", 20497.4507)["gf_tier"] == "KURUCZ"
    assert (au.action == "adopted").sum() == 33


def test_the_source_has_an_error_model_row():
    em = pd.read_csv(ROOT / "data/reference/gf_error_model.csv", comment="#")
    assert "C2013_qdt" in set(em.lab_source_tag)
