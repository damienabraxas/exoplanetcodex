"""Compare held native O LTEPOP with cached local LTE at identical T and ne.

No SE, RT, EOS replacement, abundance fit or producer-parity claim.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--comparison', type=Path, required=True)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint, marcs_number_densities
    from pipeline.bride_native_lte import NativeLte
    from pipeline.multi_atom import read_multi_atom
    comparison = json.loads(args.comparison.read_text())
    if fingerprint(args.reference) != comparison['reference']:
        raise ValueError('Changed reference')
    ref = json.loads(args.reference.read_text())
    if ref['atmosphere_file']['sha256'] != 'f4132841a18b0fc62a144d06b1f5ef20e462bf548ca06339a904adae9e84c480':
        raise ValueError('Requires the held native 5750/4.5 atmosphere')
    atom_path = require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(atom_path) != ref['atom']:
        raise ValueError('Changed source atom')
    entry = next(c for c in comparison['cases'] if c['name'] == args.case.name)
    if (fingerprint(args.case / 'populations.npz') != entry['populations']
            or fingerprint(args.case / 'manifest.json') != entry['source_manifest']):
        raise ValueError('Changed cached case')
    manifest = json.loads((args.case / 'manifest.json').read_text())
    if manifest['status'] != 'CONVERGED_EXPERIMENT' or manifest['atmosphere'] != ref['atmosphere_file']:
        raise ValueError('Wrong atmosphere or unconverged case')
    source = read_multi_atom(atom_path)
    structure = {k: np.asarray(v) if isinstance(v, list) else v for k, v in ref['atmosphere'].items()}
    ne, _ = marcs_number_densities(structure)
    data = np.load(args.case / 'populations.npz')
    np.testing.assert_array_equal(data['temperature_K'], structure['temperature_K'])
    native = NativeLte(require_codex_path('experiments.bride_multi_lte_oracle'))
    fractions, constants = native.fractions(structure['temperature_K'], ne,
        [l.E_cm for l in source.levels], [l.g for l in source.levels], [l.stage for l in source.levels])
    local = data['nstar'] / data['nstar'].sum(axis=0)
    ratio = fractions / local
    ion = next(i for i, l in enumerate(source.levels) if l.stage == 2)
    rows = [dict(level_1based=i+1, depth_1based=d+1, stage=l.stage,
        logtau500=float(structure['logtau_500'][d]), temperature_K=float(structure['temperature_K'][d]),
        ne_m3=float(ne[d]), native_LTE_fraction=float(fractions[i,d]),
        local_LTE_fraction=float(local[i,d]), native_over_local=float(ratio[i,d]))
        for i, l in enumerate(source.levels) for d in range(len(ne))]
    result = dict(schema='codex.bride_native_lte_audit/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), wrapper=fingerprint(ROOT / 'pipeline/bride_native_lte.py'),
        oracle=native.receipt, reference=fingerprint(args.reference), comparison=fingerprint(args.comparison),
        atom=ref['atom'], atmosphere=ref['atmosphere_file'], populations=entry['populations'],
        constants=constants, nlevels=len(source.levels), ndepth=len(ne),
        max_abs_native_over_local_minus_one=float(abs(ratio-1).max()),
        ion_ground_1based=ion+1, top_ion_native_over_local=float(ratio[ion,0]),
        ion_native_over_local_range=[float(ratio[ion].min()),float(ratio[ion].max())],
        max_native_conservation_error=float(abs(fractions.sum(axis=0)-1).max()),
        controls='Same full O41 energies/weights/stages, T and ne; native fractions normalized to unit total. No SE/RT updates.',
        limitations='Tests the held LTEPOP formula at matched inputs, not producer ne, total oxygen density, EOS, native opacity or rate parity. Does not identify the departure-grid producer build.')
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / 'fractions.csv').open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    result['fractions'] = fingerprint(args.output / 'fractions.csv')
    (args.output / 'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: result[k] for k in ('max_abs_native_over_local_minus_one', 'top_ion_native_over_local', 'max_native_conservation_error')}))


if __name__ == '__main__':
    main()
