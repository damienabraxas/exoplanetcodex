import copy
import json
from pathlib import Path
import numpy as np
import pytest


def test_matched_start_requires_case_fields_controls_and_exact_lte_state():
    from pipeline.bride_thermal_initial import matched_start
    d = Path(__file__).resolve().parents[1]/'data/results/rya1281/thermal_slice_v1/runs/thermal_T5572'
    r = json.loads((d/'report.json').read_text())
    expected = {k:r[k] for k in ('case', 'plan', 'sweep_adapter', 'certificate', 'atmosphere_sources',
        'local_binding', 'source_manifest', 'source_binding', 'lines', 'pressure_source', 'stark_source')}
    expected.update(initial_populations=r['source_populations'], baseline_controls=r['coupled']['controls'])
    with np.load(d/'populations.npz') as saved:
        nstar = saved['nstar'].copy(); pop = saved['populations'].copy()
    actual, receipt = matched_start(d, expected, nstar, pop.sum(axis=0))
    np.testing.assert_array_equal(actual, pop)
    assert receipt['mode'] == 'matched_physical_baseline'
    for key, value in (('case', dict(r['case'], teff=5822)), ('certificate', {}),
                       ('baseline_controls', dict(r['coupled']['controls'], rtol=.1))):
        bad = copy.deepcopy(expected); bad[key] = value
        with pytest.raises(ValueError): matched_start(d, bad, nstar, pop.sum(axis=0))
    badstar = nstar.copy(); badstar[0, 0] *= 1.01
    with pytest.raises(AssertionError): matched_start(d, expected, badstar, pop.sum(axis=0))
    with pytest.raises(ValueError, match='nuclei'):
        matched_start(d, expected, nstar, pop.sum(axis=0)*1.01)
