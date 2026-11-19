import csv
import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import rya1354_na_model as V61  # noqa: E402
import rya1354_onevalence_model as M  # noqa: E402
M.use("Na I")

OUT = ROOT / 'data/audit/rya1354/na_model_v62'


def test_sixj_sum_rule_and_d_line_ratio():
    # sum_{J,J'} (2J+1)(2J'+1) {l J 1/2; J' l' 1}^2 = 2S+1 = 2 for every allowed l -> l+1
    for l in range(3):
        tot = sum((j1 + 1) * (j2 + 1) * M.sixj(Fr(l), Fr(j1, 2), Fr(1, 2), Fr(j2, 2), Fr(l + 1), Fr(1)) ** 2
                  for j1 in (2 * l - 1, 2 * l + 1) if j1 > 0 for j2 in (2 * l + 1, 2 * l + 3) if abs(j1 - j2) <= 2)
        assert abs(tot - 2) < 1e-12
    d2 = 2 * 4 * M.sixj(Fr(0), Fr(1, 2), Fr(1, 2), Fr(3, 2), Fr(1), Fr(1)) ** 2
    d1 = 2 * 2 * M.sixj(Fr(0), Fr(1, 2), Fr(1, 2), Fr(1, 2), Fr(1), Fr(1)) ** 2
    assert abs(d2 / d1 - 2) < 1e-12


def test_solver_reproduces_hydrogen(monkeypatch):
    monkeypatch.setattr(M, 'potential', lambda r, l, p: -1 / r)
    w, z, k0 = M.solve(0, None, 3)
    assert k0 == 0
    assert np.allclose(w[:3], [-0.5, -0.125, -1 / 18], rtol=2e-4)
    w1, z1, _ = M.solve(1, None, 1)
    _, r, dx = M.grid()
    # <1s|r|2p> radial integral for hydrogen = 1.2903 a.u.
    assert abs(abs((z[:, 0] * z1[:, 0] * r * dx).sum()) - 1.29027) < 2e-3


def test_ff_parser_reads_both_gauges_of_the_d_lines():
    ff9 = V61.parse_ff('Na_11.9.U.lsj.db.gz')
    d2 = ff9[('3s', 1, '3p', 3)]
    assert abs(d2['gf_L'] - 1.28895) < 1e-5 and abs(d2['gf_V'] - 1.28710) < 1e-5
    assert len(ff9) == 385


def test_committed_score_is_consistent():
    s = json.load(open(OUT / 'summary.json'))
    assert s['status'] == 'CANDIDATE_HOLD_NOT_ADOPTED'
    # the model must beat FAC v61 on the experimental anchor and stay within 2% of it
    for v in s['experiment'].values():
        assert abs(v['dlogA']) < 0.009
    # classes must be ordered: class A agrees with the referee better than class C
    c = s['vs_ff9_by_class']
    assert c['A']['robust_sigma'] < c['C']['robust_sigma']
    assert s['sigma_coverage_vs_ff9']['frac_within_2sigma'] >= 0.85
    rows = list(csv.DictReader(open(OUT / 'solar_fill.csv')))
    assert all(r['model_loggf'] for r in rows if float(r['max_depth']) >= 0.01)
    # the Asplund Na lines: model within 0.015 dex of FF
    assert all(abs(v) < 0.015 for v in s['fill']['asplund_lines_vs_ff9'].values())


def test_mg2_committed_score_is_consistent():
    s = json.load(open(ROOT / 'data/audit/rya1354/mg2_model_v63/summary.json'))
    assert s['species'] == 'Mg II' and s['status'] == 'CANDIDATE_HOLD_NOT_ADOPTED'
    for v in s['experiment'].values():          # h & k within the class-A anchor's 3%
        assert abs(v['dlogA']) < 0.013
    assert s['sigma_coverage_vs_ff9']['frac_within_2sigma'] >= 0.9
    assert all(abs(v) < 0.02 for v in s['fill']['asplund_lines_vs_ff9'].values())
    # the near-cancelling 3s-4p line must not set the s-p family sigma (h & k stay better than class C)
    rows = {r['air_A']: r for r in csv.DictReader(open(ROOT / 'data/audit/rya1354/mg2_model_v63/solar_fill.csv'))}
    assert rows['2795.528']['model_class'] in 'AB'
