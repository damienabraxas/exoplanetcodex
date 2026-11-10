"""Native full-network hydrogen/density reporting; spectra and inputs unchanged."""
import argparse,csv,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
EQ_ANCHOR='      rho=rho*1.2123e-8/tt\n'
ATM_ANCHOR='        print*,k,pg(k),pgpg,ro(k),rhotsuji\n'
EQ_FIELDS=['temperature_K','electron_pressure_dyn_cm2','input_gas_pressure_dyn_cm2','EOS_gas_pressure_dyn_cm2','EOS_rho_g_cm3','mass_per_H_mass','H_nuclei_pressure_dyn_cm2','H_mean_mass_amu','electron_mass_amu']
ATM_FIELDS=['temperature_K','electron_pressure_dyn_cm2','gas_pressure_dyn_cm2','rho_g_cm3','xi_cm_s','boltz_erg_K','H_mass_g','EOS_rho_g_cm3','mass_per_H_mass']


def reporting_sources(eq,bsyn):
    if eq.count(EQ_ANCHOR)!=1 or eq.count('        implicit none\n')!=1 or bsyn.count(ATM_ANCHOR)!=1:raise ValueError('Ambiguous source anchors')
    block="""* BRIDE read-only EOS inventory. Counters only; no physics assignment.
      write(*,'(A,3I7,1P,9E25.16)') 'BRIDE_BIND_EQ ',
     & kk,nmetal,nmol,tt,pe,pgin,pgas,rho,xmytsuji,fictpres_h,
     & atmass(1),atmass(99)
      do j=1,nmetal
        do m=0,3
          write(*,'(A,3I7,1P,2E25.16)') 'BRIDE_BIND_SPEC ',
     &    kk,j+m*nmetal,merge(1,0,nelemx(j).eq.1),
     &    parptsuji(j+m*nmetal),xmass(j+m*nmetal)
        enddo
        if (kk.eq.1) then
          write(*,'(A,I7,1P,3E25.16)') 'BRIDE_BIND_COMP ',
     &    nelemx(j),ccomp(nelemx(j)),atmass(nelemx(j)),
     &    inputabund(nelemx(j))
        endif
      enddo
      do j=1,nmol
        bride_hcount=0
        do jjj=1,mmax(j)
          if (nelem(jjj,j).eq.1) then
            bride_hcount=bride_hcount+natom(jjj,j)
          endif
        enddo
        write(*,'(A,3I7,1P,2E25.16)') 'BRIDE_BIND_SPEC ',
     &    kk,j+4*nmetal,bride_hcount,parptsuji(j+4*nmetal),
     &    xmass(j+4*nmetal)
      enddo
"""
    atm="""        write(*,'(A,I7,1P,9E25.16)') 'BRIDE_BIND_ATMOS ',
     &  k,t(k),pe(k),pg(k),ro(k),xi(k),boltz,mh,rhotsuji,xmytsuji
"""
    if any(len(line)>72 for line in(block+atm).splitlines()if not line.startswith('*')):raise ValueError('Invalid fixed-form report block')
    return eq.replace('        implicit none\n','        implicit none\n        integer bride_hcount\n').replace(EQ_ANCHOR,EQ_ANCHOR+block),bsyn.replace(ATM_ANCHOR,ATM_ANCHOR+atm),block+atm


def parse_binding(stdout,ndep=56):
    eq={};atmos={};species={};composition=[]
    for line in stdout.splitlines():
        if not line.startswith('BRIDE_BIND_'):continue
        f=line.split();kind=f[0];ids=3 if kind in ('BRIDE_BIND_EQ','BRIDE_BIND_SPEC')else 1
        values=[float(x.replace('D','E'))for x in f[1+ids:]];indices=list(map(int,f[1:1+ids]))
        if not np.isfinite(values).all():raise ValueError('Nonfinite native binding row')
        if kind=='BRIDE_BIND_EQ':
            k,nmetal,nmol=indices
            if k in eq or len(values)!=len(EQ_FIELDS)or min(values)<=0:raise ValueError('Invalid/duplicate EOS header')
            eq[k]=dict(zip(EQ_FIELDS,values),nmetal=nmetal,nmol=nmol)
        elif kind=='BRIDE_BIND_ATMOS':
            k=indices[0]
            if k in atmos or len(values)!=len(ATM_FIELDS)or min(values)<=0:raise ValueError('Invalid/duplicate atmosphere row')
            atmos[k]=dict(zip(ATM_FIELDS,values))
        elif kind=='BRIDE_BIND_SPEC':
            k,i,h=indices
            if len(values)!=2 or values[0]<0 or values[1]<=0 or h<0 or(k,i)in species:raise ValueError('Invalid/duplicate species row')
            species[k,i]=(h,*values)
        elif kind=='BRIDE_BIND_COMP':
            if len(values)!=3 or min(values)<=0:raise ValueError('Invalid native composition')
            composition.append(dict(Z=indices[0],relative_to_H=values[0],mean_mass_amu=values[1],input_relative_abundance=values[2]))
        else:raise ValueError('Unknown native reporting row')
    expected=set(range(1,ndep+1))
    if set(eq)!=expected or set(atmos)!=expected:raise ValueError('Incomplete native binding depths')
    dimensions={(r['nmetal'],r['nmol'])for r in eq.values()}
    if len(dimensions)!=1:raise ValueError('Changed network dimensions')
    nmetal,nmol=dimensions.pop();count=4*nmetal+nmol
    if len(composition)!=nmetal or len({r['Z']for r in composition})!=nmetal:raise ValueError('Incomplete/duplicate elemental inventory')
    if set(species)!={(k,i)for k in expected for i in range(1,count+1)}:raise ValueError('Incomplete species inventory')
    table=np.array([[species[k,i]for i in range(1,count+1)]for k in range(1,ndep+1)])
    return eq,atmos,table,composition


def summarize(eq,atmos,table):
    rows=[]
    for j,k in enumerate(sorted(eq)):
        e,a=eq[k],atmos[k];Hcount,pressure,mass=table[j].T
        Hsum=float(np.sum(Hcount*pressure));masssum=float(np.sum(pressure*mass)+e['electron_pressure_dyn_cm2']*e['electron_mass_amu'])
        rho=masssum*float(np.float32(1.2123e-8))/e['temperature_K'];eps=np.finfo(np.float32).eps
        operations=2*len(pressure)+8;bound=float(operations*eps/(1-operations*eps))
        Hrel=abs(Hsum/e['H_nuclei_pressure_dyn_cm2']-1);rhorel=abs(rho/e['EOS_rho_g_cm3']-1)
        row=dict(depth_index=k,**a,EOS_gas_pressure_dyn_cm2=e['EOS_gas_pressure_dyn_cm2'],H_nuclei_pressure_dyn_cm2=e['H_nuclei_pressure_dyn_cm2'],H_nuclei_sum_dyn_cm2=Hsum,hydrogen_aggregation_relative_difference=Hrel,EOS_density_reconstruction_relative_difference=rhorel,prospective_arithmetic_bound=bound,
            nH_network_m3=e['H_nuclei_pressure_dyn_cm2']/(a['boltz_erg_K']*a['temperature_K'])*1e6,ne_native_m3=a['electron_pressure_dyn_cm2']/(a['boltz_erg_K']*a['temperature_K'])*1e6,
            EOS_over_atmos_rho=e['EOS_rho_g_cm3']/a['rho_g_cm3'],EOS_over_input_gas_pressure=e['EOS_gas_pressure_dyn_cm2']/a['gas_pressure_dyn_cm2'],neutral_H_fraction=float(pressure[0]/Hsum),gas_pressure_sum_over_EOS=float((pressure.sum()+e['electron_pressure_dyn_cm2'])/e['EOS_gas_pressure_dyn_cm2']))
        # The atomic inventory is ordered by native source; locate H via its
        # independently exported stoichiometry instead of assuming first atom.
        atomic_slots=4*e['nmetal'];row['neutral_H_fraction']=float(np.sum(pressure[:e['nmetal']]*Hcount[:e['nmetal']])/Hsum)
        row['H_in_molecules_fraction']=float(np.sum(pressure[atomic_slots:]*Hcount[atomic_slots:])/Hsum)
        if Hrel>bound or rhorel>bound:row['arithmetic_passed']=False
        else:row['arithmetic_passed']=True
        rows.append(row)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac native reporting diagnostic required')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure,marcs_number_densities
    from scripts.build_bride_native_lte import source_list
    from scripts.run_bride_native_molecules import native_status
    plan_path=ROOT/'config/experiments/rya1281_cn_native_binding.json';plan=json.loads(plan_path.read_text());old=require_codex_path('experiments.bride_native_molecules');build_dir=require_codex_path('experiments.bride_native_lte_build');run=json.loads((old/'report.json').read_text());build=json.loads((build_dir/'report.json').read_text());engine=require_codex_path('engines.ispec');checks=[]
    def check(path,receipt):
        if fingerprint(path)!=receipt:raise ValueError('Changed protected source/input '+str(path))
        checks.append((path,receipt))
    if build['source_commit']!=plan['source_commit']or fingerprint(build_dir/'report.json')!=run['build']:raise ValueError('Wrong native lineage')
    for item in build['sources']:check(engine/item['relative'],item['receipt'])
    for item in run['native_DATA']:check(engine/'synthesizer/turbospectrum/DATA'/item['relative'],item['receipt'])
    for item in json.loads((ROOT/'data/results/rya1281/missing_opacity_v1/report.json').read_text())['files']:check(engine/'input/linelists/turbospectrum/molecules'/item['filename'],item['installed'])
    for name,key in [('continuum.opa','opacity'),('atomic.bsyn','atomic_input')]:check(old/name,run[key])
    model=require_codex_path('grids.marcs_standard')/run['model']['name'];check(model,run['model']);structure=read_marcs_structure(model);marcs_ne,marcs_nH=marcs_number_densities(structure)
    original_program=build_dir/'bsyn_lu/bsyn_lu';check(original_program,next(x['program']for x in build['builds']if x['target']=='bsyn_lu'))
    a.output.mkdir(parents=True,exist_ok=False);private=a.output/'source';shutil.copytree(engine/'synthesizer/turbospectrum/source',private)
    eq_path=private/'eqmol_pe_lu.f';bsyn_path=private/'bsyn.f';eq,bsyn,block=reporting_sources(eq_path.read_text(),bsyn_path.read_text());eq_path.write_text(eq);bsyn_path.write_text(bsyn);(a.output/'reporting_block.txt').write_text(block)
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Existing compiler required')
    names=source_list((engine/'synthesizer/turbospectrum/exec-gf/Makefile').read_text(),'bsyn_lu',None);program=private/'bsyn_binding';command=[compiler,'-fconvert=big-endian','-O3','-fno-automatic','-fbacktrace','-std=legacy','-mcmodel=small','-fallow-argument-mismatch','-o',str(program),*names]
    before=time.perf_counter();r=subprocess.run(command,cwd=private,text=True,capture_output=True,timeout=plan['compile_seconds']);(a.output/'compiler.log').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Native binding reporting build failed')
    report=dict(schema='codex.bride_cn_native_binding/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),source_pin=plan['source_commit'],source_run=fingerprint(old/'report.json'),source_build=fingerprint(build_dir/'report.json'),reporting_block=fingerprint(a.output/'reporting_block.txt'),reporting_eq_source=fingerprint(eq_path),reporting_bsyn_source=fingerprint(bsyn_path),program=fingerprint(program),compiler_command=command,compile_seconds=time.perf_counter()-before,cases=[]);reference=None
    for case in plan['cases']:
        entry=next(x for x in run['cases']if x['name']==case);check(old/case/'input.txt',entry['input']);check(old/case/'native_flux.txt',entry['native_flux'])
        for mode,binary in [('original',original_program),('reporting',program)]:
            directory=a.output/case/mode;directory.mkdir(parents=True);(directory/'DATA').symlink_to(engine/'synthesizer/turbospectrum/DATA',target_is_directory=True);text=(old/case/'input.txt').read_text();old_path=str(old/case/'native_flux.txt')
            if text.count(old_path)!=1:raise ValueError('Unknown native result path')
            text=text.replace(old_path,str(directory/'native_flux.txt'));(directory/'input.txt').write_text(text);before=time.perf_counter();r=subprocess.run([str(binary)],input=text,cwd=directory,text=True,capture_output=True,timeout=plan['seconds_per_program'],env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'));seconds=time.perf_counter()-before
            (directory/'stdout.txt').write_text(r.stdout);(directory/'stderr.txt').write_text(r.stderr);native_status(r.returncode,r.stdout)
            if(directory/'native_flux.txt').read_bytes()!=(old/case/'native_flux.txt').read_bytes():raise ValueError('Native binding instrumentation changed spectrum')
            report['cases'].append(dict(case=case,mode=mode,seconds=seconds,spectrum_byte_identical=True,input=fingerprint(directory/'input.txt'),stdout=fingerprint(directory/'stdout.txt'),flux=fingerprint(directory/'native_flux.txt')))
            if mode=='reporting':
                eq,atmos,table,composition=parse_binding(r.stdout);rows=summarize(eq,atmos,table)
                for j,row in enumerate(rows):row.update(MARCS_LW_nH_m3=float(marcs_nH[j]),MARCS_LW_ne_m3=float(marcs_ne[j]),network_nH_over_MARCS_LW=float(row['nH_network_m3']/marcs_nH[j]),native_ne_over_MARCS_LW=float(row['ne_native_m3']/marcs_ne[j]))
                if reference is None:reference=(eq,atmos,table,composition,rows)
                elif eq!=reference[0]or atmos!=reference[1]or not np.array_equal(table,reference[2])or composition!=reference[3]:raise ValueError('Molecular line switch changed native EOS inventory')
    eq,atmos,table,composition,rows=reference
    with(a.output/'state_binding.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    np.savez(a.output/'implemented_species.npz',H_nuclei_per_species=table[:,:,0].astype(np.int32),pressure_dyn_cm2=table[:,:,1],mass_amu=table[:,:,2])
    (a.output/'composition.json').write_text(json.dumps(composition,indent=2,allow_nan=False)+'\n')
    for path,receipt in checks:
        if fingerprint(path)!=receipt:raise ValueError('Protected input changed during probe')
    report.update(protected_input_check_slots=len(checks),molecular_switch_state_identity=True,ndep=56,nmetal=eq[1]['nmetal'],nmol=eq[1]['nmol'],species_per_depth=table.shape[1],arithmetic_passed=all(row['arithmetic_passed']for row in rows),model=run['model'],model_hydrogen_mass_fraction=structure['hydrogen_mass_fraction'],marcs_density_adapter=fingerprint(ROOT/'pipeline/bride_reference.py'),summary={key:[float(min(r[key]for r in rows)),float(max(r[key]for r in rows))]for key in ['EOS_over_atmos_rho','EOS_over_input_gas_pressure','network_nH_over_MARCS_LW','native_ne_over_MARCS_LW','H_in_molecules_fraction','hydrogen_aggregation_relative_difference','EOS_density_reconstruction_relative_difference','gas_pressure_sum_over_EOS']},limitations=['Aggregation on exported native implemented network is not independent chemical/EOS validation or a completeness claim.','Native atmosphere and recomputed EOS densities are kept distinct; no normalization selected or corrected.','MARCS/LW adapter uses its own existing constants/normalization; ratios diagnose source/convention differences, not a new abundance.','No provider admitted or native populations attached to LW; all real spectra remain unchanged.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['nmetal','nmol','species_per_depth','arithmetic_passed','summary']}));print(f'Read {len(checks)} protected input slots; wrote source/build/six replays, complete species array and56depth binding report to {a.output}')
    if not report['arithmetic_passed']:raise RuntimeError('Native aggregation arithmetic check failed; retained without adjustment')

if __name__=='__main__':main()
