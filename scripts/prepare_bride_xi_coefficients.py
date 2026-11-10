"""Fresh exact-input producer coefficients for a declared microturbulence slice."""
import argparse
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--reference-run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != 'Darwin':
        raise RuntimeError('Mac-only bounded opacity diagnostic')
    from config.constants import require_codex_path, get_star_params
    from pipeline.bride_reference import fingerprint, read_marcs_structure
    from pipeline.bride_parameter_checks import marcs_filename, interpolate_marcs_structure
    from pipeline.bride_ts_points import point_input
    from pipeline.bride_ts_budget import read_ts_opacity, read_ts_components, native_real32_print
    from pipeline.wavelength_util import air_to_vac

    prior=json.loads((ROOT/'data/results/rya1281/ts_points_v1/run_report.json').read_text())
    if fingerprint(args.build/'manifest.json')!=prior['build'] or fingerprint(args.reference_run)!=prior['reference']:raise ValueError('Changed qualified point producer/reference')
    build = json.loads((args.build/'manifest.json').read_text())
    for entry in build['builds']:
        if fingerprint(args.build/entry['name']/'babsma_lu') != entry['program']:
            raise ValueError('Changed point adapter binary')
    reference = json.loads(args.reference_run.read_text())
    if reference.get('support') != 'nuv200_350' or not reference['reporting_preserves_serialized_coefficients']:
        raise ValueError('Requires verified native broad-NUV baseline/reporting reference')
    sources = []
    native_path = require_codex_path('grids.marcs_standard')/reference['atmosphere']['name']
    if fingerprint(native_path) != reference['atmosphere']:
        raise ValueError('Changed native reference model')
    native = read_marcs_structure(native_path)
    sources.append(fingerprint(native_path))
    from pipeline.bride_xi_sweep import validate_plan,reconstruct,source_input_xi
    plan_path=args.plan;plan=json.loads(plan_path.read_text())
    validate_plan(plan,json.loads((ROOT/'config/experiments/rya1281_mac_solar_bound_control.json').read_text()))
    cases=[('native_identity',native)]
    for case in plan['cases']:
        base,structure,contract,receipts=reconstruct(case)
        sources.extend(receipts);cases.append((case['name'],base))
    original_cases = {}
    for kind in ('baseline', 'reporting'):
        entry = next(x for x in reference['cases'] if x['name'] == 'NUV_source_air_'+kind)
        directory = args.reference_run.parent/entry['name']
        for filename, key in (('continuum.opa', 'output'), ('detail.txt', 'detail'), ('input.txt', 'input')):
            if fingerprint(directory/filename) != entry[key]:
                raise ValueError('Changed native reference output/input')
        original_cases[kind] = (entry, directory)
    args.output.mkdir(parents=True, exist_ok=False)
    data = args.output/'DATA'
    shutil.copytree(args.reference_run.parent/'DATA', data)
    for entry in reference['input_data']:
        if fingerprint(data/entry['relative_path']) != {k: v for k, v in entry.items() if k != 'relative_path'}:
            raise ValueError('Changed native opacity data')
    result = dict(schema='codex.bride_ts_point_run/1', scientific_admission='HOLD', status='RUNNING',
        generator=fingerprint(__file__), adapter=fingerprint(ROOT/'pipeline/bride_ts_points.py'),
        interpolator=fingerprint(ROOT/'pipeline/bride_parameter_checks.py'), plan=fingerprint(plan_path),
        stars_config=fingerprint(ROOT/'config/stars.yaml'), source_models=sources,
        build=fingerprint(args.build/'manifest.json'), reference=fingerprint(args.reference_run),
        interpretation='Pointwise native TS opacity/EOS only. Native geometry/tau are driver placeholders, not changed-input equilibrium fields. No RT/SE, spectrum, O abundance or production packet. Full56 rows, unlike prior55-row solves. Native REAL quantization retained. Changed-input raw per-mass opacity uses TS EOS density; per-neutral coefficient avoids silently applying that density to local atmosphere.', cases=[])

    def save():
        (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')

    save()
    anchor = str(original_cases['baseline'][1].parent/'inputs'/native_path.name)
    for name, structure in cases:
        started = time.perf_counter()
        record = dict(name=name, atmosphere_identity=structure['model_id'],requested_xi_km_s=structure['vturb_m_s']/1000,runs=[])
        outputs = {}
        for kind in ('baseline', 'reporting'):
            directory = args.output/(name+'_'+kind)
            directory.mkdir()
            (directory/'DATA').symlink_to(data, target_is_directory=True)
            points = directory/'point_state.txt'
            points.write_text(point_input(structure['temperature_K'], structure['electron_pressure_Pa']))
            original_text = (original_cases[kind][1]/'input.txt').read_text()
            # Use the fingerprint-verified native file directly, read-only.
            if original_text.count(anchor) != 1:
                raise ValueError('Unknown original model-input layout')
            text = original_text.replace(anchor, str(native_path))
            text=source_input_xi(text,structure['vturb_m_s']/1000)
            (directory/'input.txt').write_text(text)
            run = subprocess.run([str((args.build/kind/'babsma_lu').resolve())], input=text,
                cwd=directory, text=True, capture_output=True, timeout=60)
            (directory/'stdout.log').write_text(run.stdout)
            (directory/'stderr.log').write_text(run.stderr)
            if run.returncode:
                raise RuntimeError('TS point run failed: '+directory.name)
            parsed = read_ts_opacity(directory/'continuum.opa')
            np.testing.assert_array_equal(parsed['temperature_K'], native_real32_print(structure['temperature_K']))
            np.testing.assert_array_equal(parsed['electron_pressure_Pa'], native_real32_print(structure['electron_pressure_Pa']/.1)*.1)
            np.testing.assert_array_equal(parsed['vturb_m_s'],native_real32_print(np.full(structure['ndep'],structure['vturb_m_s']/1000))*1000)
            if parsed['ndep'] != structure['ndep']:
                raise ValueError('Changed point count')
            outputs[kind] = parsed
            record['runs'].append(dict(kind=kind, input=fingerprint(directory/'input.txt'),
                points=fingerprint(points), output=fingerprint(directory/'continuum.opa'), detail=fingerprint(directory/'detail.txt')))
        left, right = [args.output/(name+'_'+kind)/'continuum.opa' for kind in ('baseline', 'reporting')]
        if left.read_bytes() != right.read_bytes():
            raise ValueError('Reporting derivative changes serialized opacity')
        if name == 'native_identity':
            if left.read_bytes() != (original_cases['baseline'][1]/'continuum.opa').read_bytes():
                raise ValueError('Native point override changes reference serialized opacity')
        parsed = outputs['reporting']
        detail = args.output/(name+'_reporting')/'detail.txt'
        comp = read_ts_components(detail, parsed['ndep'], parsed['nwave'])
        fe_index = comp['labels'].index('FeI')
        fe = comp['records'][:, 1:, 2].T*.1*comp['fractions'][:, 1:, fe_index].T
        t, pe, rho, partition, p_fe, p_feii = comp['states'].T
        neutral_per_kg = p_fe/(1.38066e-16*t)*1e6/(rho*1000)
        per_neutral = fe/neutral_per_kg[None, :]
        if not np.isfinite(per_neutral).all() or np.any(per_neutral < 0) or np.any(neutral_per_kg <= 0):
            raise ValueError('Invalid neutral Fe coefficient')
        np.testing.assert_allclose(comp['post'][:, :, 3][:, 1:].T*.1, parsed['absorption_m2_kg'], rtol=6e-5, atol=0)
        array = args.output/(name+'_coefficients.npz')
        np.savez_compressed(array, wavelength_vac_nm=air_to_vac(parsed['wavelength_air_A'])*.1,
            producer_vturb_m_s=parsed['vturb_m_s'],requested_vturb_m_s=np.full(structure['ndep'],structure['vturb_m_s']),
            requested_temperature_K=structure['temperature_K'], requested_electron_pressure_Pa=structure['electron_pressure_Pa'],
            producer_temperature_K=parsed['temperature_K'], producer_electron_pressure_Pa=parsed['electron_pressure_Pa'],
            producer_rho_kg_m3=parsed['rho_kg_m3'], local_rho_kg_m3=structure['rho_kg_m3'],
            neutral_Fe_per_kg=neutral_per_kg, FeI_opacity_m2_kg=fe, FeI_opacity_per_neutral_m2=per_neutral,
            absorption_m2_kg=parsed['absorption_m2_kg'], scattering_m2_kg=parsed['scattering_m2_kg'])
        if name == 'native_identity':
            old_comp = read_ts_components(original_cases['reporting'][1]/'detail.txt', parsed['ndep'], parsed['nwave'])
            for key in ('states', 'records', 'post'):
                np.testing.assert_array_equal(comp[key], old_comp[key])
        record.update(producer_xi_m_s=parsed['vturb_m_s'].tolist(),array=fingerprint(array), ndep=parsed['ndep'], nwave=parsed['nwave'],
            baseline_reporting_byte_identity=True, native_reference_byte_identity=name == 'native_identity',
            max_abs_T_quantization_K=float(abs(parsed['temperature_K']-structure['temperature_K']).max()),
            max_fractional_Pe_quantization=float(abs(parsed['electron_pressure_Pa']/structure['electron_pressure_Pa']-1).max()),
            max_fractional_producer_local_rho_difference=float(abs(parsed['rho_kg_m3']/structure['rho_kg_m3']-1).max()),
            elapsed_seconds=time.perf_counter()-started)
        result['cases'].append(record)
        save()
        print(name+' complete', flush=True)
    result['status'] = 'DIAGNOSTIC_COMPLETE'
    save()


if __name__ == '__main__':
    main()
