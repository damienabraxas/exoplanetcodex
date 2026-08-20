from pathlib import Path

import pandas as pd
import pytest

from scripts.rya927_compare_website_ews import BASELINE, compare


def test_baseline_is_exactly_predeclared_o2_b_set():
    d = pd.read_csv(BASELINE, comment="#")
    assert len(d) == 7
    assert d.wavelength_air_A.between(6867.0, 6884.0, inclusive="both").all()
    assert (d.element == "Fe").sum() == 4


def test_comparison_reports_ew_change(tmp_path):
    d = pd.read_csv(BASELINE, comment="#")
    d["ew_mA"] *= 0.5
    p = tmp_path / "corrected.csv"
    d.to_csv(p, index=False)
    out = compare(p)
    assert (out.ratio_corrected_to_uncorrected == 0.5).all()
    assert (out.delta_ew_mA < 0).all()


def test_missing_predeclared_line_is_loud(tmp_path):
    d = pd.read_csv(BASELINE, comment="#").iloc[:-1]
    p = tmp_path / "corrected.csv"
    d.to_csv(p, index=False)
    with pytest.raises(ValueError, match="omitted predeclared"):
        compare(p)
