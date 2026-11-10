"""Compact independently checked review of the isolated TS point-opacity runs."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--build-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    run = json.loads(args.run.read_text())
    if run['status'] != 'DIAGNOSTIC_COMPLETE' or fingerprint(args.build_manifest) != run['build']:
        raise ValueError('Incomplete run or changed build manifest')
    original = require_codex_path('experiments.bride_fe_packet')
    original_manifest = json.loads((original/'manifest.json').read_text())
    if fingerprint(original/'packet.npz') != original_manifest['array']:
        raise ValueError('Changed original native Fe packet')
    packet = np.load(original/'packet.npz')
    rows, selected = [], []
    previous = None
    for case in run['cases']:
        path = args.run.parent/case['array']['name']
        if fingerprint(path) != case['array']:
            raise ValueError('Changed point coefficient array')
        data = np.load(path)
        if case['name'] == 'native_identity':
            for new, old in (('FeI_opacity_m2_kg', 'opacity_m2_kg'), ('FeI_opacity_per_neutral_m2', 'opacity_net_per_neutral_Fe_m2'), ('neutral_Fe_per_kg', 'native_neutral_Fe_per_kg')):
                np.testing.assert_array_equal(data[new], packet[old])
        wave = data['wavelength_vac_nm']
        index = int(np.argmin(abs(wave-277.)))
        coefficient = data['FeI_opacity_per_neutral_m2'][index]
        row = {k: v for k, v in case.items() if k not in ('runs', 'array', 'atmosphere_identity')}
        row['selected_wavelength_vac_nm'] = float(wave[index])
        # Native row 29 is a previously declared diagnostic location, not a
        # contribution-function-selected layer. All56 depths are also retained.
        for depth in (0, 28):
            row[f'FeI_per_neutral_m2_depth_{depth+1}'] = float(coefficient[depth])
            row[f'fractional_change_previous_10K_depth_{depth+1}'] = None if previous is None or case['name'] == 'native_identity' else float(coefficient[depth]/previous[depth]-1)
        rows.append(row)
        for depth in range(case['ndep']):
            selected.append(dict(case=case['name'], depth_1based=depth+1, wavelength_vac_nm=float(wave[index]),
                requested_temperature_K=data['requested_temperature_K'][depth], producer_temperature_K=data['producer_temperature_K'][depth],
                requested_electron_pressure_Pa=data['requested_electron_pressure_Pa'][depth], producer_electron_pressure_Pa=data['producer_electron_pressure_Pa'][depth],
                producer_rho_kg_m3=data['producer_rho_kg_m3'][depth], local_rho_kg_m3=data['local_rho_kg_m3'][depth],
                native_neutral_Fe_per_kg=data['neutral_Fe_per_kg'][depth], FeI_opacity_m2_kg=data['FeI_opacity_m2_kg'][index, depth],
                FeI_per_neutral_m2=coefficient[depth]))
        if case['name'] != 'native_identity':
            previous = coefficient.copy()
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_ts_point_review/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), source_run=fingerprint(args.run), source_build=fingerprint(args.build_manifest),
        original_native_packet=fingerprint(original/'manifest.json'), native_packet_array_identity=True,
        method='Exact native array equality; selected diagnostic wavelength closest to277nm, all56 rows; adjacent solar input responses, not cross-section/abundance errors. Producer EOS density remains separate from local input density.', cases=rows)
    for name, data in [('summary.csv', rows), ('selected_points.csv', selected)]:
        with (args.output/name).open('x') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]), lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)
        result[name] = fingerprint(args.output/name)
    for name, source in [('run_report.json', args.run), ('build_manifest.json', args.build_manifest)]:
        shutil.copyfile(source, args.output/name)
        result[name] = fingerprint(args.output/name)
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Verified native identity; retained six cases and336 selected point rows')


if __name__ == '__main__':
    main()
