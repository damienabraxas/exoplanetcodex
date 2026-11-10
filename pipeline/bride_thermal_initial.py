"""Verified same-physical-case starting populations for numerical refinement."""
import json
import numpy as np


def matched_start(directory, expected, target_nstar, target_total):
    """Load a converged baseline as an initial guess; never accept it as a solve."""
    from pipeline.bride_reference import fingerprint
    r = json.loads((directory/'report.json').read_text())
    if (r['schema'] != 'codex.bride_thermal_flux/1' or r['numeric_refinement']
            or r['mode'] != 'cowley' or r['population_mode'] != 'coupled'
            or r['scientific_admission'] != 'HOLD'):
        raise ValueError('Require coupled unrefined thermal baseline')
    for key in ('case', 'plan', 'sweep_adapter', 'certificate', 'atmosphere_sources',
                'local_binding', 'source_manifest', 'source_binding', 'lines',
                'pressure_source', 'stark_source'):
        if r[key] != expected[key]: raise ValueError('Warm-start physical/source mismatch: '+key)
    if r['source_populations'] != expected['initial_populations']:
        raise ValueError('Warm-start original population lineage mismatch')
    if r['coupled']['controls'] != expected['baseline_controls']:
        raise ValueError('Warm-start baseline controls mismatch')
    for name, key in (('populations.npz', 'populations'), ('rates.npz', 'rates'), ('history.jsonl', 'history')):
        if fingerprint(directory/name) != r['coupled'][key]: raise ValueError('Changed warm-start artifact: '+name)
    history = [json.loads(line) for line in (directory/'history.jsonl').read_text().splitlines()]
    if (len(history) != r['coupled']['iterations'] or len(history) < 3
            or history[-1] != r['coupled']['final_diagnostics']):
        raise ValueError('Wrong warm-start joint history')
    last = history[-1]; rtol = expected['baseline_controls']['rtol']
    if (max(last[k] for k in ('population_change', 'physical_rate_residual', 'J_change', 'cap_boundary_change')) >= rtol
            or last['particle_conservation_error'] > 1e-10
            or r['coupled']['fresh_physical_rate_residual'] >= rtol):
        raise ValueError('Unconverged warm-start baseline')
    with np.load(directory/'populations.npz') as saved:
        pop = saved['populations'].copy(); nstar = saved['nstar'].copy()
    np.testing.assert_array_equal(nstar, target_nstar)
    if (pop.shape != target_nstar.shape or not np.isfinite(pop).all() or np.any(pop <= 0)
            or not np.allclose(pop.sum(axis=0), target_total, rtol=1e-10, atol=0)):
        raise ValueError('Warm-start populations/nuclei mismatch')
    receipt = dict(mode='matched_physical_baseline', report=fingerprint(directory/'report.json'),
                   populations=fingerprint(directory/'populations.npz'), helper=fingerprint(__file__),
                   scope='Initial guess only; unchanged fresh joint convergence, physical rate and conservation gates required.')
    return pop, receipt
