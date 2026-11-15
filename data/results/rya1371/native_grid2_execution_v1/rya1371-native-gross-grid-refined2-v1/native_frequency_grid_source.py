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
    limit = 86048 if control['kind'] == 'continuum-refinement2' else 59003
    if len(original) != 39841 or not len(original) <= len(wave) <= limit or producer['proposed_points'] != len(wave):
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
    predecessor_identities = {}
    if control['kind'] == 'continuum-refinement2':
        parent = producer.get('predecessor', {})
        if set(parent) != {'resource','report','wave'} or parent['resource'] != 'experiments.openmulti3d_oxygen_grid_refined_run':
            raise ValueError('Second refinement requires the first qualified native producer')
        parent_dir = Path(resolver(parent['resource']))
        for name, path in [('report',parent_dir/'report.json'),('wave',parent_dir/'native/save/all_lam.bin')]:
            if fingerprint(path) != parent[name] or producer['inputs'].get(str(path)) != parent[name]:
                raise ValueError('Changed second-refinement predecessor')
            predecessor_identities[str(path)] = parent[name]
        parent_doc = json.loads((parent_dir/'report.json').read_text())
        parent_plan = ROOT/'config/experiments/rya1371_openmulti3d_oxygen_grid_refined.json'
        parent_outputs = [r for r in parent_doc['output_files'] if r['path']=='native/save/all_lam.bin']
        if not parent_doc['computational_completed'] or parent_doc['plan'] != fingerprint(parent_plan) or len(parent_outputs)!=1 or parent_outputs[0]['receipt']!=parent['wave']:
            raise ValueError('Unqualified first-refinement predecessor')
        inherited_parent = np.fromfile(parent_dir/'native/save/all_lam.bin',dtype='<f8')
        j = np.searchsorted(wave,inherited_parent)
        if len(inherited_parent)!=59003 or len(wave)<=59003 or np.any(j>=len(wave)) or not np.array_equal(wave[j],inherited_parent):
            raise ValueError('Second refinement removed predecessor nodes or added none')
    # Recheck all inputs after reading; copy bytes only after every refusal check.
    for key, path in paths.items():
        if fingerprint(path) != control[key]:
            raise ValueError('Native grid changed during validation')
    for path,expected in predecessor_identities.items():
        if fingerprint(path)!=expected:raise ValueError('Predecessor changed during validation')
    return paths['grid'], dict(kind=control['kind'], original_points=len(original), points=len(wave),
        identities={str(path):control[key] for key,path in paths.items()},
        all_original_nodes_exact=True, scope='Explicit diagnostic native line_grid_file; physical and resource controls unchanged. No HOP packet or scientific admission.')
