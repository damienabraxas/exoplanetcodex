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
    from pipeline.bride_reference import fingerprint
    matched = copy.deepcopy(ref)
    matched['initialization'] = dict(mode='matched_physical_baseline', populations=base['coupled']['populations'],
        helper=fingerprint(ROOT/'pipeline/bride_thermal_initial.py'))
    matched['generator'] = fingerprint(ROOT/'scripts/run_bride_thermal_sweep.py')
    verify_pair(base, matched, plan)
    matched['initialization']['populations'] = {}
    with pytest.raises(ValueError, match='initialization'): verify_pair(base, matched, plan)


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
    with pytest.raises(ValueError, match='Disconnected'):
        independent_balance(pop, np.zeros_like(off), np.zeros_like(off), np.zeros_like(off), np.array([5.]), 1e-6)


def test_timeout_receipt_cannot_accept_a_report_or_different_physical_case(tmp_path):
    from scripts.analyze_bride_thermal_refinement import failed_attempt
    base = json.loads((ROOT/'data/results/rya1281/thermal_slice_v1/runs/thermal_T5572/report.json').read_text())
    plan = json.loads((ROOT/'config/experiments/rya1281_mac_thermal_slice.json').read_text())
    inputs = {k:base[k] for k in ('case', 'plan', 'sweep_adapter', 'certificate', 'atmosphere_sources',
        'local_binding', 'source_manifest', 'source_binding', 'lines', 'pressure_source', 'stark_source', 'generator')}
    controls = dict(plan['controls'], rays=5, qscale=2, flux_rays=7)
    inputs.update(schema='codex.bride_thermal_input/1', numeric_refinement=True,
        mode='cowley', population_mode='coupled', controls=controls)
    (tmp_path/'input_manifest.json').write_text(json.dumps(inputs))
    (tmp_path/'history.jsonl').write_text(json.dumps(base['coupled']['final_diagnostics'])+'\n')
    log = tmp_path/'attempt.log'
    log.write_text('TimeoutError: Fixed-population cap radiation budget exceeded\n')
    assert failed_attempt(tmp_path, base, plan, log, False)['status'] == 'BUDGET_EXHAUSTED_NOT_ACCEPTED'
    bad = copy.deepcopy(base); bad['case']['teff'] = 5822
    with pytest.raises(ValueError, match='physical/source'):
        failed_attempt(tmp_path, bad, plan, log, False)
    (tmp_path/'report.json').write_text('{}')
    with pytest.raises(ValueError, match='timeout'):
        failed_attempt(tmp_path, base, plan, log, False)
