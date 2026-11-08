"""Reporting-only native CN state and pre-Voigt amplitude with spectrum parity."""
import argparse,csv,json,os,platform,shlex,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
FIELDS=['wavelength_air_A','temperature_K','gas_pressure_dyn_cm2','electron_pressure_dyn_cm2','rho_g_cm3','CN_pressure_dyn_cm2','Q','ntt_per_g','ntot_per_g','lower_population_per_g','stim','doppler_width_Hz','standard_opacity_cm2_g','plez','mass_g','boltz_erg_K','h_erg_s','c_cm_s','e_esu','electron_mass_g','excitation_eV','loggf','gf','calf','constant','C_fraction','N_fraction','g','symmetry_factor','wavelength_cm','voigt_a','pre_voigt_amplitude_over_standard']
ANCHOR='111     continue\n*\n* Start wavelength loop for this line'


def reporting_patch(text,selections):
    if text.count(ANCHOR)!=1:raise ValueError('Ambiguous native reporting insertion')
    clauses=[]
    for row in selections:
        clauses.append(f"     &    (isotope(1).eq.{row['C']}.and.isotope(2).eq.{row['N']}\n     &    .and.abs(xlb-{row['wave']:.3f}).lt.0.0001)")
    condition='\n     &    .or.\n'.join(clauses)
    block="""* BRIDE reporting only: no assignment to native physical variables.
      if (lele.eq.'C N'.and.(
"""+condition+"""
     &    )) then
        do k=1,ntau
          write(*,'(A,4I7,1P,32E25.16)') 'BRIDE_CN_STATE ',
     &      isotope(1),isotope(2),iline,k,
     &      xlb,t(k),pg(k),pe(k),ro(k),partryck(k,molindex),
     &      fpartition(k),ntt(k),ntot(k),n(k),stim(k),dnud(k),
     &      ross(k),plez(k),ma,boltz,h,c,e,m,chie,gfelog,f,
     &      calf,constant,isotopfrac(atom(1),isotope(1)),
     &      isotopfrac(atom(2),isotope(2)),g,symmfactor,xl,a(k),
     &      calf*plez(k)
        enddo
      endif
"""
    # Fixed-form source columns are semantic: refuse an overflowing insertion.
    if any(len(line)>72 for line in block.splitlines()if not line.startswith('*')):raise ValueError('Reporting source exceeds fixed-form width')
    return text.replace(ANCHOR,'111     continue\n'+block+'*\n* Start wavelength loop for this line'),block


def parse_states(stdout):
    rows=[]
    for line in stdout.splitlines():
        if not line.startswith('BRIDE_CN_STATE '):continue
        fields=line.split()[1:]
        if len(fields)!=4+len(FIELDS):raise ValueError('Malformed native CN reporting row')
        ids=list(map(int,fields[:4]));v=np.array([float(x.replace('D','E'))for x in fields[4:]])
        if not np.isfinite(v).all():raise ValueError('Nonfinite native CN state')
        rows.append(dict(zip(['C','N','native_line_index','depth_index']+FIELDS,ids+v.tolist())))
    return rows


def compare_states(rows,partition_record,eps_multiple):
    from scripts.audit_bride_cn_interface import project_partition
    r={k:np.array([row[k]for row in rows])for k in FIELDS};T=r['temperature_K'];k=r['boltz_erg_K']
    if any(np.any(r[x]<=0)for x in FIELDS if x!='loggf'):raise ValueError('Nonpositive CN physical state')
    # Independently evaluate the inspected equations in float64 on native state.
    total=r['CN_pressure_dyn_cm2']/(k*T*r['rho_g_cm3'])*r['C_fraction']*r['N_fraction']*r['symmetry_factor']
    Q=project_partition(T,partition_record['coefficients']).astype(float)
    lower=total*r['g']*np.exp(-r['excitation_eV']*(5040./T)/.43429)/Q
    stim=-np.expm1(-r['h_erg_s']*r['c_cm_s']/(k*r['wavelength_cm']*T))
    gf=10.**r['loggf'];constant=np.sqrt(np.pi)*r['e_esu']**2/(r['electron_mass_g']*r['c_cm_s'])
    plez=lower*stim/r['doppler_width_Hz']/r['standard_opacity_cm2_g']
    derived={'ntt_per_g':total,'ntot_per_g':total,'Q':Q,'lower_population_per_g':lower,'stim':stim,'gf':gf,'constant':constant,'calf':constant*gf,'plez':plez,'pre_voigt_amplitude_over_standard':constant*gf*plez}
    bound=eps_multiple*np.finfo(np.float32).eps;checks={}
    for name,value in derived.items():
        rel=float(np.max(np.abs(value/r[name]-1)));checks[name]=dict(max_relative_difference=rel,prospective_relative_bound=bound,passed=rel<=bound)
    return checks


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac native diagnostic required')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from scripts.build_bride_native_lte import source_list
    from scripts.run_bride_native_molecules import native_status
    plan_path=ROOT/'config/experiments/rya1281_cn_coefficients.json';plan=json.loads(plan_path.read_text())
    source_run=require_codex_path('experiments.bride_native_molecules');source_build=require_codex_path('experiments.bride_native_lte_build');interface=require_codex_path('experiments.bride_cn_interface')
    run=json.loads((source_run/'report.json').read_text());build=json.loads((source_build/'report.json').read_text());contract=json.loads((interface/'report.json').read_text());engine=require_codex_path('engines.ispec');pin=plan['source_commit']
    if build['source_commit']!=pin or fingerprint(source_build/'report.json')!=run['build']:raise ValueError('Wrong native build lineage')
    protected=[]
    def check(path,receipt):
        if fingerprint(path)!=receipt:raise ValueError('Changed protected input '+str(path))
        protected.append((path,receipt))
    for row in build['sources']:check(engine/row['relative'],row['receipt'])
    for row in run['native_DATA']:check(engine/'synthesizer/turbospectrum/DATA'/row['relative'],row['receipt'])
    for row in contract['protected_sources']:check(engine/row['path'],row['receipt'])
    native_program=source_build/'bsyn_lu/bsyn_lu';check(native_program,next(x['program']for x in build['builds']if x['target']=='bsyn_lu'))
    check(source_run/'continuum.opa',run['opacity']);check(source_run/'atomic.bsyn',run['atomic_input'])
    check(require_codex_path('grids.marcs_standard')/run['model']['name'],run['model'])
    selections=[]
    for entry in contract['line_lists']:
        filename=entry['filename'];path=engine/'input/linelists/turbospectrum/molecules'/filename;lines=path.read_text().splitlines();parsed=[]
        for i,line in enumerate(lines[2:],3):
            f=shlex.split(line)
            if f:parsed.append((abs(float(f[0])-plan['target_air_A']),i,float(f[0]),line))
        _,number,wave,line=min(parsed);C,N=map(int,filename.split('C')[0:1]+[filename.split('C')[1].split('N')[0]])
        selections.append(dict(filename=filename,C=C,N=N,source_line=number,wave=wave,raw=line))
    a.output.mkdir(parents=True,exist_ok=False);private=a.output/'source';shutil.copytree(engine/'synthesizer/turbospectrum/source',private)
    original=(private/'bsyn.f').read_text();patched,block=reporting_patch(original,selections);(private/'bsyn.f').write_text(patched);(a.output/'reporting_block.txt').write_text(block)
    recipe=engine/'synthesizer/turbospectrum/exec-gf/Makefile';names=source_list(recipe.read_text(),'bsyn_lu',None)
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Existing compiler required')
    program=private/'bsyn_reporting';command=[compiler,'-fconvert=big-endian','-O3','-fno-automatic','-fbacktrace','-std=legacy','-mcmodel=small','-fallow-argument-mismatch','-o',str(program),*names]
    start=time.perf_counter();r=subprocess.run(command,cwd=private,text=True,capture_output=True,timeout=plan['compile_seconds']);(a.output/'compiler.log').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Reporting build failed; retained compiler log')
    result=dict(schema='codex.bride_cn_coefficients/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),native_run=fingerprint(source_run/'report.json'),native_build=fingerprint(source_build/'report.json'),interface=fingerprint(interface/'report.json'),source_pin=pin,selection=selections,original_source=fingerprint(engine/'synthesizer/turbospectrum/source/bsyn.f'),reporting_source=fingerprint(private/'bsyn.f'),reporting_block=fingerprint(a.output/'reporting_block.txt'),compiler_command=command,compile_seconds=time.perf_counter()-start,program=fingerprint(program),cases=[])
    for name in plan['cases']:
        old=source_run/name;entry=next(c for c in run['cases']if c['name']==name);check(old/'input.txt',entry['input']);check(old/'native_flux.txt',entry['native_flux']);records=[]
        for label,binary in [('original',native_program),('reporting',program)]:
            directory=a.output/name/label;directory.mkdir(parents=True);(directory/'DATA').symlink_to(engine/'synthesizer/turbospectrum/DATA',target_is_directory=True)
            text=(old/'input.txt').read_text();old_name=str(old/'native_flux.txt')
            if text.count(old_name)!=1:raise ValueError('Unknown result path')
            text=text.replace(old_name,str(directory/'native_flux.txt'));(directory/'input.txt').write_text(text);before=time.perf_counter()
            r=subprocess.run([str(binary)],input=text,cwd=directory,text=True,capture_output=True,timeout=plan['seconds_per_program'],env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'));seconds=time.perf_counter()-before
            (directory/'stdout.txt').write_text(r.stdout);(directory/'stderr.txt').write_text(r.stderr);native_status(r.returncode,r.stdout)
            if (directory/'native_flux.txt').read_bytes()!=(old/'native_flux.txt').read_bytes():raise ValueError('Reporting/original spectrum differs from accepted original')
            states=parse_states(r.stdout)
            if label=='original'and states:raise ValueError('Original binary unexpectedly instrumented')
            if label=='reporting':records=states
            result['cases'].append(dict(case=name,program=label,seconds=seconds,spectrum_byte_identical=True,flux=fingerprint(directory/'native_flux.txt'),stdout=fingerprint(directory/'stdout.txt'),input=fingerprint(directory/'input.txt'),state_rows=len(states)))
        expected=0 if name=='off'else 56*len(selections)
        if len(records)!=expected:raise ValueError('Incomplete/duplicate CN reporting rows')
        if records:
            for selection in selections:
                subset=[row for row in records if(row['C'],row['N'])==(selection['C'],selection['N'])]
                if [row['depth_index']for row in subset]!=list(range(1,57)):raise ValueError('Missing or unordered CN depths')
                if any(abs(row['wavelength_air_A']-selection['wave'])>.0005 for row in subset):raise ValueError('Wrong selected native line')
            with (a.output/name/'states.csv').open('w',newline='')as f:
                w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
            result.setdefault('comparisons',{})[name]=compare_states(records,contract['partition_entry'],plan['comparison_float32_eps_multiple'])
    for path,receipt in protected:
        if fingerprint(path)!=receipt:raise ValueError('Protected source changed during reporting run')
    result['protected_input_checks']=len(protected);result['all_equation_checks_pass']=all(x['passed']for case in result['comparisons'].values()for x in case.values())
    result['limitations']=['Same-source equation comparison on exported state, not independent EOS or molecular-primary validation.','The pre-Voigt amplitude is normalized by standard opacity; no profile, spectral integration or abundance claim.','128float32eps is a prospective arithmetic screening bound, not a physical/numerical opacity error budget.','No new per-program peak-memory measurement; previous native repeats were about26MB, not a new measurement.']
    (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result['comparisons']));print(f'Read {len(protected)} protected inputs; wrote private source/build, six native run logs/spectra, two state CSVs and report to {a.output}')
    if not result['all_equation_checks_pass']:raise RuntimeError('CN equation comparison failed; evidence retained without tolerance adjustment')

if __name__=='__main__':main()
