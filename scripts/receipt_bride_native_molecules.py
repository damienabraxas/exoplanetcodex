"""Read-only source invariance, native-opacity and preserved-attempt receipts."""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--failed-build',type=Path,required=True);p.add_argument('--failed-run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure
    from pipeline.bride_ts_budget import read_ts_opacity
    run=json.loads((a.source/'report.json').read_text());engine=require_codex_path('engines.ispec');checks=[]
    def check(path,expected):
        if fingerprint(path)!=expected:raise ValueError('Changed protected/declared input: '+str(path))
        checks.append(str(path))
    model=require_codex_path('grids.marcs_standard')/run['model']['name'];check(model,run['model'])
    check(engine/'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv',run['atomic_source'])
    check(engine/'input/isotopes/SPECTRUM.lst',run['isotopes'])
    check(require_codex_path('experiments.bride_current_orch_checkout')/'data/linelists/canonical_gf.csv',run['canonical'])
    check(a.source/'atomic.bsyn',run['atomic_input']);check(a.source/'continuum.opa',run['opacity']);check(a.build/'report.json',run['build'])
    build=json.loads((a.build/'report.json').read_text())
    for item in build['sources']:check(engine/item['relative'],item['receipt'])
    for item in build['builds']:check(a.build/item['target']/item['target'],item['program'])
    for item in run['native_DATA']:check(engine/'synthesizer/turbospectrum/DATA'/item['relative'],item['receipt'])
    inventory=json.loads((ROOT/'data/results/rya1281/missing_opacity_v1/report.json').read_text())
    for item in inventory['files']:check(engine/'input/linelists/turbospectrum/molecules'/item['filename'],item['installed'])
    opacity=read_ts_opacity(a.source/'continuum.opa');structure=read_marcs_structure(model)
    if opacity['ndep']!=structure['ndep']:raise ValueError('Native opacity truncated reference depth count')
    failed_log=a.failed_build/'bsyn_lu/compiler.log'
    if (a.failed_build/'report.json').exists()or 'path=clock'not in failed_log.read_text():raise ValueError('Original build failure not preserved')
    failed_stdout=a.failed_run/'off/stdout.txt'
    if (a.failed_run/'report.json').exists()or(a.failed_run/'off/native_flux.txt').exists()or'input.f; bad keyword: NLTE:'not in failed_stdout.read_text():raise ValueError('Original input failure not preserved')
    if any('NLTE'in line for line in (a.source/'atomic.bsyn').read_text().splitlines()if line.lstrip().startswith("'")):raise ValueError('Accepted atomic export is not LTE')
    result=dict(schema='codex.bride_native_molecule_receipts/1',scientific_admission='HOLD',generator=fingerprint(__file__),source=fingerprint(a.source/'report.json'),
        unchanged_checks=len(checks),protected_inputs_unchanged=True,native_opacity_ndep=opacity['ndep'],native_opacity_nwave=opacity['nwave'],
        native_opacity_air_min=float(opacity['wavelength_air_A'].min()),native_opacity_air_max=float(opacity['wavelength_air_A'].max()),
        failed_build_log=fingerprint(failed_log),failed_run_stdout=fingerprint(failed_stdout),no_accepted_first_attempt=True,
        accepted_atomic_headers_all_LTE=True,scope='Source/build/derived-input invariance and native ABI guards, not physical molecular/EOS validation')
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
