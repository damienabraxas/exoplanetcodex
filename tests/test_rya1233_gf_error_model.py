"""RYA-1233: gf errors are correlated the way each source's literature says, per element."""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.gf_error_model import covariance, load  # noqa: E402


def _term(sig, tags, src=None):
    cov, note, unrev = covariance(sig, src or ["s"] * len(sig), tags)
    w = np.full(len(sig), 1 / len(sig))
    return math.sqrt(w @ np.array(cov) @ w), note, unrev


def test_garz_averages_down_over_lines_as_amarsi_and_asplund_treat_it():
    """Si's gf term read 0.129 dex with Garz's 0.119 per-line error charged fully correlated;
    random per line + the O'Brian & Lawler lifetime scale (0.019) gives 0.048 on 7 lines."""
    t, note, unrev = _term([0.1186] * 7, ["AGSS21_Si_Amarsi2017"] * 7)
    assert 0.045 < t < 0.052 and not unrev
    assert "O'Brian" in load()["AGSS21_Si_Amarsi2017"]["basis"]


def test_an_unclassified_source_stays_conservative_and_is_named():
    t, note, unrev = _term([0.05] * 3, ["SomeLab2010"] * 3)
    assert abs(t - 0.05) < 1e-9
    assert unrev == ["SomeLab2010"] and "UNREVIEWED" in note


def test_one_scale_group_spans_experimental_and_calculated_pr2024():
    cov, _, _ = covariance([0.03, 0.05], ["a", "b"], ["PR2024_exp", "PR2024_calc"])
    assert cov[0][1] > 0


def test_every_registry_row_cites_its_basis():
    for tag, r in load().items():
        assert r["correlation"] in ("scale_only", "fully_correlated"), tag
        assert r["citation"].strip() and r["basis"].strip(), tag
