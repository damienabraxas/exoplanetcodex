"""Bounded Mac TS baseline/reporting control at the native 277nm comparison node.

Copies source inputs; never calls iSpec's temporary-directory cleanup wrapper.
No spectrum synthesis, NLTE, opacity substitution or abundance fit.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure
    from pipeline.bride_ts_budget import read_ts_opacity,native_real32_print
    from pipeline.wavelength_util import vac_to_air
    build=require_codex_path('experiments.bride_ts_budget_build')
    manifest=json.loads((build/'manifest.json').read_text())
    if manifest['source_commit']!='21269948750e2177fe32c5df38e1d9bb81363374':raise ValueError('Changed build pin')
    for entry in manifest['builds']:
        if fingerprint(build/entry['name']/'babsma_lu')!=entry['program']:raise ValueError('Changed binary')
    intake_path=require_codex_path('experiments.bride_marcs_opacity')/'manifest.json'
    intake=json.loads(intake_path.read_text())
    source_model=require_codex_path('grids.marcs_standard')/intake['atmosphere']['name']
    if fingerprint(source_model)!=intake['atmosphere']:raise ValueError('Changed atmosphere')
    structure=read_marcs_structure(source_model)
    args.output.mkdir(parents=True,exist_ok=False)
    inputs=args.output/'inputs';inputs.mkdir();model=inputs/source_model.name;shutil.copy2(source_model,model)
    original_data=require_codex_path('engines.ispec')/'synthesizer'/'turbospectrum'/'DATA'
    staged=args.output/'DATA';shutil.copytree(original_data,staged)
    data_receipts=[dict(relative_path=str(p.relative_to(original_data)),**fingerprint(p))
                   for p in sorted(original_data.rglob('*')) if p.is_file()]
    if any(fingerprint(original_data/x['relative_path'])!=fingerprint(staged/x['relative_path']) for x in data_receipts):
        raise ValueError('Changed staged TS data')
    report=dict(schema='codex.bride_ts_budget_run/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),parser=fingerprint(ROOT/'pipeline/bride_ts_budget.py'),
        build=fingerprint(build/'manifest.json'),atmosphere=intake['atmosphere'],
        input_data=data_receipts,source_coordinate_A=2769.76,controls=dict(
            seconds_per_case=60,threads=1,NLTE=False,spectrum_synthesis=False,
            global_abundance_scaling_dex=0.,explicit_abundances=structure['abundances_dex']),cases=[])
    for hypothesis,target in [('OPA_is_vacuum',float(vac_to_air(np.array([2769.76]))[0])),
                              ('OPA_is_air',2769.76)]:
        outputs=[]
        for kind in ('baseline','reporting'):
            directory=args.output/(hypothesis+'_'+kind);directory.mkdir();(directory/'DATA').symlink_to(staged,target_is_directory=True)
            text=(f"'LAMBDA_MIN:' '{target-1:.9f}'\n'LAMBDA_MAX:' '{target+1:.9f}'\n"
                  "'LAMBDA_STEP:' '1.0'\n"
                  f"'MODELINPUT:' '{model}'\n'MARCS-FILE:' '.true.'\n"
                  "'MODELOPAC:' 'continuum.opa'\n'LOGFILE:' 'detail.txt'\n"
                  "'METALLICITY:' '0.00'\n'ALPHA/Fe:' '0.00'\n'HELIUM:' '0.00'\n"
                  "'R-PROCESS:' '0.00'\n'S-PROCESS:' '0.00'\n'INDIVIDUAL ABUNDANCES:' '92'\n")
            text+=''.join(f'{i+1} {a:.2f}\n' for i,a in enumerate(structure['abundances_dex']))
            text+=f"'XIFIX:' 'T'\n{structure['vturb_m_s']/1000:.9f}\n"
            (directory/'input.txt').write_text(text)
            start=time.perf_counter();run=subprocess.run([str(build/kind/'babsma_lu')],input=text,
                cwd=directory,text=True,capture_output=True,timeout=60)
            (directory/'stdout.log').write_text(run.stdout);(directory/'stderr.log').write_text(run.stderr)
            if run.returncode:raise RuntimeError('TS failure; outputs retained: '+directory.name)
            output=directory/'continuum.opa';parsed=read_ts_opacity(output)
            if parsed['ndep']!=structure['ndep']:raise ValueError('Changed depth count')
            np.testing.assert_array_equal(parsed['temperature_K'],native_real32_print(structure['temperature_K']))
            np.testing.assert_array_equal(parsed['electron_pressure_Pa'],
                native_real32_print(structure['electron_pressure_Pa']/.1)*.1)
            entry=dict(name=directory.name,elapsed_seconds=time.perf_counter()-start,
                effective_wavelength_air_A_range=[float(parsed['wavelength_air_A'][0]),float(parsed['wavelength_air_A'][-1])],
                effective_nwave=parsed['nwave'],
                max_abs_temperature_quantization_K=float(abs(parsed['temperature_K']-structure['temperature_K']).max()),
                output=fingerprint(output),input=fingerprint(directory/'input.txt'),
                stdout=fingerprint(directory/'stdout.log'),stderr=fingerprint(directory/'stderr.log'),
                detail=fingerprint(directory/'detail.txt'),
                max_fractional_electron_pressure_change=float(abs(parsed['electron_pressure_Pa']/structure['electron_pressure_Pa']-1).max()),
                max_fractional_density_change=float(abs(parsed['rho_kg_m3']/structure['rho_kg_m3']-1).max()),
                recovery_messages=[l for l in run.stdout.splitlines() if any(s in l.lower()for s in ('trying with','giving up','not converged'))])
            report['cases'].append(entry);outputs.append(output)
        if outputs[0].read_bytes()!=outputs[1].read_bytes():raise ValueError('Reporting changed serialized opacity output')
    report['reporting_preserves_serialized_coefficients']=True
    report['limitations']='TS source/EOS control, not native MARCS/MULTI producer. OPA medium remains conditional. Output density changes are reported, not silently overwritten. Original sources/data protected.'
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['cases'],indent=2))


if __name__=='__main__':main()
