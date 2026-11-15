"""Strict opt-in retained-node native diagnostic grids; no default-grid change."""
import json
from pathlib import Path
import numpy as np


def load_native_frequency_grid(plan, *, resolver):
    from pipeline.reference_engines import ROOT, fingerprint
    from pipeline.native_reference_controls import validate_native_grid_plan
    qualified = ROOT / 'config/experiments/rya1371_openmulti3d_oxygen_shared16.json'
    control = validate_native_grid_plan(plan, json.loads(qualified.read_text()))
    source = Path(resolver(control['source_resource']))
    directory = Path(resolver(control['resource']))
    paths = {'source_report': source/'report.json', 'source_wave': source/'native/save/all_lam.bin',
             'producer': directory/'report.json', 'grid': directory/'frequency_grid_vac_A.txt'}
    for key, path in paths.items():
        if fingerprint(path) != control[key]:
            raise ValueError('Changed native grid identity: '+key)
    source_doc = json.loads(paths['source_report'].read_text())
    if not source_doc['computational_completed'] or source_doc['plan'] != fingerprint(qualified):
        raise ValueError('Native grid source is not qualified shared16 producer')
    entries = [r for r in source_doc['output_files'] if r['path'] == 'native/save/all_lam.bin']
    if len(entries) != 1 or entries[0]['receipt'] != control['source_wave']:
        raise ValueError('Native grid source output receipt mismatch')
    producer = json.loads(paths['producer'].read_text())
    if not producer['computational_completed'] or producer.get('all_original_nodes_exact') is not True or producer['grid'] != control['grid']:
        raise ValueError('Unqualified retained-node grid producer')
    for name in ['source_report', 'source_wave']:
        if producer['inputs'].get(str(paths[name])) != control[name]:
            raise ValueError('Native grid producer belongs to another source')
    original = np.fromfile(paths['source_wave'], dtype='<f8')
    wave = np.loadtxt(paths['grid'])
    for values in [original, wave]:
        if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all() or np.any(values <= 0) or np.any(np.diff(values) <= 0):
            raise ValueError('Invalid native frequency grid')
    if len(original) != 39841 or not len(original) <= len(wave) <= 59003 or producer['proposed_points'] != len(wave):
        raise ValueError('Native grid exceeds explicitly scoped refinement')
    index = np.searchsorted(wave, original)
    if np.any(index >= len(wave)) or not np.array_equal(wave[index], original) or not np.array_equal(wave[[0,-1]], original[[0,-1]]):
        raise ValueError('Native grid removed original nodes/support')
    inherited = np.zeros(len(wave), bool); inherited[index] = True
    changed_pair = ~(inherited[:-1] & inherited[1:])
    if np.any(np.diff(wave)[changed_pair] <= float(np.float32(.001))):
        raise ValueError('New native nodes collide under source merge tolerance')
    if control['kind'] == 'baseline' and not np.array_equal(wave, original):
        raise ValueError('Native baseline grid must be exact original')
    if control['kind'] == 'continuum-refinement' and len(wave) <= len(original):
        raise ValueError('Native refinement must add nodes')
    # Recheck all inputs after reading; copy bytes only after every refusal check.
    for key, path in paths.items():
        if fingerprint(path) != control[key]:
            raise ValueError('Native grid changed during validation')
    return paths['grid'], dict(kind=control['kind'], original_points=len(original), points=len(wave),
        identities={str(path):control[key] for key,path in paths.items()},
        all_original_nodes_exact=True, scope='Explicit diagnostic native line_grid_file; physical and resource controls unchanged. No HOP packet or scientific admission.')
