"""Discriminating controls for numerical-only replay comparisons."""
import copy
import json
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_pair_rejects_wrong_physics_case_and_relaxed_budget():
    from scripts.analyze_bride_thermal_refinement import verify_pair
    root = ROOT/'data/results/rya1281/thermal_slice_v1/runs'
    base = json.loads((root/'thermal_nominal/report.json').read_text())
    ref = json.loads((root/'refined_nominal/report.json').read_text())
    plan = json.loads((ROOT/'config/experiments/rya1281_mac_thermal_slice.json').read_text())
    verify_pair(base, ref, plan)
    for key, value in [('case', dict(ref['case'], teff=5822)), ('certificate', {}), ('mode', 'none')]:
        bad = copy.deepcopy(ref); bad[key] = value
        with pytest.raises(ValueError): verify_pair(base, bad, plan)
    bad = copy.deepcopy(ref); bad['coupled']['controls']['seconds_per_case'] += 1
    with pytest.raises(ValueError, match='budgets'): verify_pair(base, bad, plan)


def test_independent_balance_discriminates_populations_operator_and_nuclei():
    from scripts.analyze_bride_thermal_refinement import independent_balance
    # 0->1 = 2/s, 1->0 = 3/s, stationary n0:n1 = 3:2.
    off = np.array([[[0.], [3.]], [[2.], [0.]]])
    total = np.array([[[-2.], [3.]], [[2.], [-3.]]])
    pop = np.array([[3.], [2.]])
    result = independent_balance(pop, off, np.zeros_like(off), total, np.array([5.]), 1e-6)
    assert result == dict(residual=0., particle_conservation_error=0.)
    with pytest.raises(ValueError, match='balance'):
        independent_balance(pop[::-1], off, np.zeros_like(off), total, np.array([5.]), 1e-6)
    with pytest.raises(ValueError, match='conservation'):
        independent_balance(pop, off, np.zeros_like(off), total, np.array([6.]), 1e-6)
    with pytest.raises(AssertionError):
        independent_balance(pop, off, np.zeros_like(off), total*.9, np.array([5.]), 1e-6)
