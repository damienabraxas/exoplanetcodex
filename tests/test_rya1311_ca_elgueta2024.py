"""RYA-1311: Y/J-band Ca gf from Elgueta et al. 2024 (MNRAS 532, 3694, Table 8)."""
import csv
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import rya1311_adopt_elgueta2024_ca as e24  # noqa: E402
from pipeline.gf_grades import is_gf_graded  # noqa: E402


def _canon(species, w):
    rows = [r for r in csv.DictReader(open(ROOT / "data/linelists/canonical_gf.csv"))
            if r["species"] == species and abs(float(r["wavelength_air_A"]) - w) < 0.01]
    assert len(rows) == 1
    return rows[0]


def test_adopted_rows_are_codex_graded_but_not_lab():
    rows = pd.DataFrame([_canon("Ca I", 10343.819), _canon("Ca II", 9854.759)])
    assert is_gf_graded(rows).all()
    assert (rows.gf_tier == "E2024-ASTRO").all()       # never the Reference (LAB) pool


def test_the_vald_calibration_is_preferred_and_mb99_is_the_fallback():
    assert float(_canon("Ca I", 10838.97)["log_gf"]) == -0.141          # VALD column
    assert float(_canon("Ca I", 12105.841)["log_gf"]) == -0.726         # MB99 only
    assert float(_canon("Ca I", 10838.97)["gf_sigma_dex"]) == 0.113


def test_lines_missing_from_our_list_are_reported_not_invented():
    au = pd.read_csv(ROOT / e24.AUDIT.relative_to(ROOT))
    assert set(au[au.action == "unmatched"].air_A) == {10516.14, 10791.45, 10846.79}
    assert (au.action == "adopted").sum() == 7


def test_the_source_has_an_error_model_row():
    em = pd.read_csv(ROOT / "data/reference/gf_error_model.csv", comment="#")
    assert "E2024_astro" in set(em.lab_source_tag)
