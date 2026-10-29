"""Stage one verified MARCS opacity companion from an official download archive.

Retains compressed source bytes and receipts outside the repo. No network, RT,
SE or comparison wavelength conversion is performed by this command.
"""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tarfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
MODEL = 'p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--reader', type=Path, required=True)
    parser.add_argument('--download-response', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint, read_marcs_structure
    from pipeline.bride_marcs_opacity import read_marcs_opacity
    response = json.loads(args.download_response.read_text())
    if response['status'] != 'succes' or response['error']:
        raise ValueError('Official request did not succeed')
    with tarfile.open(args.archive) as archive:
        members = archive.getmembers()
        if (len(members) != 1 or not members[0].isfile()
                or members[0].name != './opa/st/p/'+MODEL+'.opa.gz'
                or members[0].size > 2_000_000):
            raise ValueError('Wrong or oversized archive member')
        payload = archive.extractfile(members[0]).read()
    model_path = require_codex_path('grids.marcs_standard') / (MODEL+'.mod')
    if fingerprint(model_path)['sha256'] != 'f4132841a18b0fc62a144d06b1f5ef20e462bf548ca06339a904adae9e84c480':
        raise ValueError('Changed held atmosphere')
    args.output.mkdir(parents=True, exist_ok=False)
    target = args.output / (MODEL+'.opa.gz')
    target.write_bytes(payload)
    opacity = read_marcs_opacity(target)
    structure = read_marcs_structure(model_path)
    if opacity['ndep'] != structure['ndep'] or opacity['standard_wavelength_source_A'] != 5000.:
        raise ValueError('Incompatible opacity structure')
    checks = {}
    for key in ('temperature_K', 'electron_pressure_Pa', 'rho_kg_m3'):
        checks[key] = dict(max_abs_difference=float(abs(opacity[key]-structure[key]).max()),
            max_fractional_difference=float(abs(opacity[key]/structure[key]-1).max()))
    checks['logtau500_max_abs_difference'] = float(abs(np.log10(opacity['tau_standard'])-structure['logtau_500']).max())
    checks['abundances_max_abs_difference_dex'] = float(abs(opacity['abundances_dex']-structure['abundances_dex']).max())
    checks['vturb_max_abs_difference_m_s'] = float(abs(opacity['vturb_m_s']-structure['vturb_m_s']).max())
    shutil.copy2(args.reader, args.output/'readopa.f')
    result = dict(schema='codex.bride_marcs_opacity_intake/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), parser=fingerprint(ROOT/'pipeline/bride_marcs_opacity.py'),
        source= fingerprint(target), archive=fingerprint(args.archive),
        official_reader=fingerprint(args.reader), atmosphere=fingerprint(model_path),
        official_response=response, official_response_receipt=fingerprint(args.download_response),
        source_urls=['https://marcs.oreme.org/data/',
                     'https://marcs.oreme.org/data/downloadmodels/',
                     'https://marcs.astro.uu.se/documents/auxiliary/readopa.f'],
        request=dict(file_type=['opa'], filenames=[MODEL+'.mod.gz']),
        ndepth=opacity['ndep'], nwave=opacity['nwave'],
        wavelength_range_source_A=[float(opacity['wavelength_source_A'][0]), float(opacity['wavelength_source_A'][-1])],
        structure_checks=checks, wavelength_medium='UNVERIFIED',
        limitations='Independent MARCS continuum absorption/scattering; not MULTI ABSDAT or native producer opacity. Wavelength medium must be established before spectral comparison; source values retained without conversion.')
    (args.output/'manifest.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(checks))


if __name__ == '__main__':
    main()
