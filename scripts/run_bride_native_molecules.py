"""Matched native LTE molecular-line toggles; retained chemistry and atomic inputs."""
import argparse,csv,importlib.util,json,os,platform,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',type=Path,required=True);p.add_argument('--orch-worktree',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac-only diagnostic')
    plan_path=ROOT/'config/experiments/rya1281_native_molecules.json';plan=json.loads(plan_path.read_text())
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.orch_worktree,text=True).strip()!=plan['orch_commit']:raise ValueError('Wrong current ORCH checkout')
    # Load the unchanged diagnostic MARCS reader; current ORCH owns gf resolution.
    spec=importlib.util.spec_from_file_location('native_diag_reference',ROOT/'pipeline/bride_reference.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    fingerprint=reader.fingerprint
    sys.path.insert(0,str(a.orch_worktree))
    from config.constants import require_codex_path
    from pipeline import gf_resolver
    if Path(gf_resolver.__file__).resolve()!=a.orch_worktree.resolve()/'pipeline/gf_resolver.py':raise ValueError('Wrong gf resolver owner')
    if os.environ.get('CODEX_GF_OVERRIDE'):raise ValueError('No gf override allowed')
    engine=require_codex_path('engines.ispec');sys.path.insert(0,str(engine));import ispec
    build=json.loads((a.build/'report.json').read_text())
    if build['source_commit']!=plan['source_commit']:raise ValueError('Wrong native build')
    for item in build['sources']:
        if fingerprint(engine/item['relative'])!=item['receipt']:raise ValueError('Changed build source')
    programs={}
    for item in build['builds']:
        path=a.build/item['target']/item['target']
        if fingerprint(path)!=item['program']:raise ValueError('Changed native binary')
        programs[item['target']]=path.resolve()
    prior=json.loads((ROOT/'data/results/rya1281/ts_points_v1/run_report.json').read_text())
    receipt=prior['source_models'][0];model=require_codex_path('grids.marcs_standard')/receipt['name']
    if fingerprint(model)!=receipt:raise ValueError('Changed reference atmosphere')
    structure=reader.read_marcs_structure(model)
    if not structure['model_id'].startswith('p5750_g+4.5_') or structure['vturb_m_s']!=1000.:raise ValueError('Not native5750/4.5/xi1 reference')
    abundance=np.array(structure['abundances_dex']);abundance[7]=plan['oxygen_abundance']
    atomic_source=engine/'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv'
    atomic=ispec.read_atomic_linelist(str(atomic_source));gf_report={};gf_resolver.apply_to_synth_array(atomic,report=gf_report)
    lo,hi=plan['range_air_A'];pad=plan['atomic_padding_A']
    selected=atomic[(atomic['wave_A']>=lo-pad)&(atomic['wave_A']<=hi+pad)].copy()
    a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
    (a.output/'DATA').symlink_to(engine/'synthesizer/turbospectrum/DATA',target_is_directory=True)
    atomic_file=a.output/'atomic.bsyn';ispec.write_atomic_linelist(selected,linelist_filename=str(atomic_file),code='turbospectrum')
    isotopes_path=engine/'input/isotopes/SPECTRUM.lst';isotopes=ispec.read_isotope_data(str(isotopes_path))
    inventory=json.loads((ROOT/'data/results/rya1281/missing_opacity_v1/report.json').read_text());molecular={name:[]for name in ('CH','CN','C2','NH','OH')}
    for item in inventory['files']:
        path=engine/'input/linelists/turbospectrum/molecules'/item['filename']
        if fingerprint(path)!=item['installed']:raise ValueError('Changed molecular source')
        molecular[item['molecule']].append(str(path.resolve()))
    result=dict(schema='codex.bride_native_molecules/1',scientific_admission='HOLD',generator=fingerprint(__file__),plan=fingerprint(plan_path),
        build=fingerprint(a.build/'report.json'),model=fingerprint(model),reference_model=structure['model_id'],model_ndep=structure['ndep'],
        composition_dex=abundance.tolist(),header_logg_cgs=structure['logg_cgs'],xi_km_s=plan['xi_km_s'],atomic_source=fingerprint(atomic_source),atomic_input=fingerprint(atomic_file),
        atomic_rows_export_requested=len(selected),atomic_unsupported_rows=int(np.sum(selected['turbospectrum_support']!='T')),gf_resolution=gf_report,current_orch_commit=plan['orch_commit'],gf_resolver=fingerprint(gf_resolver.__file__),
        canonical=fingerprint(a.orch_worktree/'data/linelists/canonical_gf.csv'),isotopes=fingerprint(isotopes_path),inventory=fingerprint(ROOT/'data/results/rya1281/missing_opacity_v1/report.json'),
        scope=plan['scope'],native_EPS=1e-4,normalized_output_decimal_places=5,cases=[])
    result['native_DATA']=[dict(relative=str(path.relative_to(engine/'synthesizer/turbospectrum/DATA')),receipt=fingerprint(path))for path in sorted((engine/'synthesizer/turbospectrum/DATA').rglob('*'))if path.is_file()]
    (a.output/'input.json').write_text(json.dumps(result,indent=2)+'\n')
    common="'METALLICITY:' '0.00'\n'ALPHA/Fe:' '0.00'\n'HELIUM:' '0.00'\n'R-PROCESS:' '0.00'\n'S-PROCESS:' '0.00'\n'INDIVIDUAL ABUNDANCES:' '92'\n"+''.join(f'{i+1} {x:.2f}\n'for i,x in enumerate(abundance))
    def invoke(target,text,directory):
        before=time.perf_counter();(directory/'input.txt').write_text(text)
        (directory/'DATA').symlink_to(engine/'synthesizer/turbospectrum/DATA',target_is_directory=True)
        run=subprocess.run([str(programs[target])],input=text,cwd=directory,text=True,capture_output=True,timeout=plan['seconds_per_program'],env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        (directory/'stdout.txt').write_text(run.stdout);(directory/'stderr.txt').write_text(run.stderr)
        if run.returncode:raise RuntimeError('Native program failed: '+str(directory))
        return time.perf_counter()-before
    opacity=a.output/'continuum.opa';opacity_dir=a.output/'opacity';opacity_dir.mkdir()
    text=f"'LAMBDA_MIN:' '{lo}'\n'LAMBDA_MAX:' '{hi}'\n'LAMBDA_STEP:' '{plan['opacity_step_A']}'\n'MODELINPUT:' '{model}'\n'MARCS-FILE:' '.true.'\n'MODELOPAC:' '{opacity}'\n"+common+f"'XIFIX:' 'T'\n{plan['xi_km_s']}\n"
    result['opacity_seconds']=invoke('babsma_lu',text,opacity_dir)
    # Native formatted opacity dimensions are explicit; no assumed LW depth parity.
    tokens=opacity.read_text().split();result['native_opacity_ndep']=int(tokens[1]);result['native_opacity_nwave']=int(tokens[3]);result['opacity']=fingerprint(opacity)
    isotope_input=f"'ISOTOPES:' '{len(isotopes)}'\n"+''.join(f"{int(r['atomic_code'])}.{int(r['mass_number']):03d} {r['relative_abundance_in_the_solar_system']:.3f}\n"for r in isotopes)
    for case in plan['cases']:
        directory=a.output/case;directory.mkdir();family=case.replace('_refined','')
        files=[]if family=='off'else [p for group in molecular.values()for p in group]if family=='all'else molecular[family]
        step=plan['refined_line_step_A']if case.endswith('_refined')else plan['line_step_A'];raw=directory/'native_flux.txt'
        text=f"'NLTE:' '.false.'\n'LAMBDA_MIN:' '{lo}'\n'LAMBDA_MAX:' '{hi}'\n'LAMBDA_STEP:' '{step}'\n'INTENSITY/FLUX:' 'Flux'\n'COS(THETA):' '1.00'\n'ABFIND:' '.false.'\n'MARCS-FILE:' '.true.'\n'MODELOPAC:' '{opacity}'\n'RESULTFILE:' '{raw}'\n"+common+isotope_input+f"'NFILES:' '{len(files)+2}'\n"+'\n'.join(files+['DATA/Hlinedata',str(atomic_file)])+"\n'SPHERICAL:' 'F'\n30\n300.00\n15\n1.30\n"
        elapsed=invoke('bsyn_lu',text,directory);data=np.loadtxt(raw)
        if data.ndim!=2 or data.shape[1]<3 or not np.isfinite(data).all()or np.any(data[:,1:]<=0)or np.any(np.diff(data[:,0])<=0):raise ValueError('Invalid native spectrum')
        expected=np.round(lo+np.arange(int((hi-lo)/step)+1)*step,3)
        np.testing.assert_array_equal(data[:,0],expected)
        output=directory/'flux.csv';np.savetxt(output,data[:,:3],delimiter=',',header='wavelength_air_A,normalized_flux,native_absolute_flux',comments='')
        record=dict(name=case,step_A=step,molecular_files=files,seconds=elapsed,rows=len(data),flux=fingerprint(output),native_flux=fingerprint(raw),input=fingerprint(directory/'input.txt'),stdout=fingerprint(directory/'stdout.txt'),stderr=fingerprint(directory/'stderr.txt'))
        result['cases'].append(record);(directory/'report.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps(dict(case=case,seconds=elapsed,rows=len(data))),flush=True)
    if fingerprint(opacity)!=result['opacity']:raise ValueError('Native continuum changed between legs')
    for item in result['native_DATA']:
        if fingerprint(engine/'synthesizer/turbospectrum/DATA'/item['relative'])!=item['receipt']:raise ValueError('Protected native DATA changed during run')
    result['seconds']=time.perf_counter()-start;(a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print('Native molecular line controls complete')


if __name__=='__main__':main()
