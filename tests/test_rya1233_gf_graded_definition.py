"""RYA-1233: Codex/Deep Grade lines are gf-graded = ANY published gf uncertainty."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.gf_grades import is_gf_graded  # noqa: E402


def _rows(**kw):
    base = dict(gf_tier="KURUCZ", gf_sigma_dex=float("nan"), nist_grade=None)
    return pd.DataFrame([{**base, **kw}])


def test_each_published_uncertainty_grades_a_line():
    assert is_gf_graded(_rows(gf_tier="LAB")).iloc[0]
    assert is_gf_graded(_rows(gf_tier="PR2024-CALC", gf_sigma_dex=0.05)).iloc[0]
    assert is_gf_graded(_rows(nist_grade="B")).iloc[0]
    assert is_gf_graded(_rows(gf_tier="NIST-C+")).iloc[0]


def test_no_published_uncertainty_is_not_graded():
    assert not is_gf_graded(_rows(gf_tier="VALD3")).iloc[0]
    assert not is_gf_graded(_rows(nist_grade="D")).iloc[0]
    assert not is_gf_graded(_rows()).iloc[0]
