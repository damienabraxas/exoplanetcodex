"""Pinned HOP standalone EOS parity via declarations read from its own C header."""
import argparse,ctypes as ct,hashlib,json,os,re,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


from pipeline.hop_abi import bind_header


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import load_engine,fingerprint
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.openmulti3d_output import read_patch_set
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native-run',type=Path,required=True);p.add_argument('--native-analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--trace-selection',action='store_true');p.add_argument('--opacity-plan',type=Path);a=p.parse_args()
    if os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Explicit one-thread environment required')
    if a.opacity_plan:a.opacity_plan=a.opacity_plan.resolve()
    _,app,_,engine=load_engine('openmulti3d');library=app/'HOP/standalone/build/libhop.so';header=app/'HOP/hopacity/libhop.h';files=[app/'input'/name for name in ['absdat/TS_absdat.txt','abund/abund_magg22.txt','isotopes/isotopes_blp.txt','molec/molec_data_BC.txt']];checks={str(x):fingerprint(x)for x in [library,header,*files]}
    native_report=a.native_run/'report.json';analysis_report=a.native_analysis/'report.json';run=json.loads(native_report.read_text());analysis=json.loads(analysis_report.read_text())
    if not run['computational_completed']or run['engine_receipt']!=engine or analysis['native_report']!=fingerprint(native_report):raise ValueError('Native run/analysis/engine identity mismatch')
    checks[str(native_report)]=fingerprint(native_report);checks[str(analysis_report)]=fingerprint(analysis_report)
    for entry in run['output_files']:
        path=a.native_run/entry['path']
        if fingerprint(path)!=entry['receipt']:raise ValueError('Changed native provider reference')
        checks[str(path)]=entry['receipt']
    derived=None
    if a.trace_selection:
        build_file=require_codex_path('engines.hop_trace_provider_build')/'receipt.json';derived=json.loads(build_file.read_text());source=require_codex_path('engines.hop_trace_provider_source');patch=ROOT/'native/hop_trace_provider/trace_continua.patch'
        if not derived['completed']or derived['source']!=str(source.resolve())or derived['upstream_commit']!=engine['source_commit']or fingerprint(patch)!=derived['patch']:raise ValueError('Invalid private provider build identity')
        if subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip()!=derived['upstream_commit']or hashlib.sha256(subprocess.check_output(['git','diff','--binary'],cwd=source)).hexdigest()!=derived['tracked_diff_sha256']:raise ValueError('Changed private provider source')
        base=source/'experiments/OpenMULTI3D/HOP';library=base/'standalone/build/libhop.so';header=base/'hopacity/libhop.h'
        if fingerprint(library)!=derived['library']or fingerprint(header)!=derived['header']:raise ValueError('Changed private provider library/header')
        for path in [build_file,patch,library,header]:checks[str(path)]=fingerprint(path)
    arrays=a.native_analysis/'arrays.npz';checks[str(arrays)]=fingerprint(arrays)
    with np.load(arrays)as d:state=d['consumed_state_cgs'].astype(float);height=d['height_m'].copy()
    T=np.ascontiguousarray(state[:,2]);ne=np.ascontiguousarray(state[:,0]);rho=np.ascontiguousarray(state[:,1]);original_ne=ne.copy()
    ref=json.loads(require_codex_path('experiments.bride_oxygen_reference').read_text());model=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name'];checks[str(model)]=fingerprint(model);s=read_marcs_structure(model);abund=np.ascontiguousarray(np.asarray(s['abundances_dex'],np.float32).astype(float))
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'input').symlink_to(app/'input',target_is_directory=True)
    api=bind_header(ct.CDLL(str(library)),header.read_text())
    def call(name,*args):
        status=api[name](*args)
        if status:
            buf=ct.create_string_buffer(2048);api['hop_last_error'](buf,2048);raise RuntimeError(name+': '+buf.value.decode())
    ptr=lambda x:x.ctypes.data_as(ct.POINTER(ct.c_double))
    selection=None;negative_checks=[];opacity_checks=None
    if a.trace_selection:
        from pipeline.multi_atom import read_multi_atom,CM1_TO_EV
        atom_file=require_codex_path('experiments.bride_oxygen_atom')
        if fingerprint(atom_file)!=run['inputs']['atom']:raise ValueError('Changed explicit trace atom')
        checks[str(atom_file)]=fingerprint(atom_file);atom=read_multi_atom(atom_file)
        stages=np.ascontiguousarray([atom.level(c.i).stage for c in atom.continua],dtype=np.int32)
        exc=np.ascontiguousarray([(atom.level(c.i).E_cm-min(l.E_cm for l in atom.levels if l.stage==atom.level(c.i).stage))*CM1_TO_EV for c in atom.continua])
        def exclude(name=atom.symbol,energy=exc,ions=stages):
            count=ct.c_int(-1);status=api['hop_exclude_model_continua'](name.encode(),atom.weight,len(exc),ions.ctypes.data_as(ct.POINTER(ct.c_int)),ptr(energy),ct.byref(count));return status,count.value
        def refuse(label,**kwargs):
            status,count=exclude(**kwargs)
            if status==0 or count!=0:raise ValueError('Invalid trace operation accepted: '+label)
            negative_checks.append(label)
        refuse('before_init')
    previous=Path.cwd();os.chdir(a.output)
    try:
        call('hop_set_external_abundances',len(abund),ptr(abund));call('hop_set_option',b'abundances',b'O=8.70')
        call('hop_init',*[str(x).encode()for x in files],0.,-99.,0.,0.)
        if a.trace_selection:
            refuse('hydrogen',name='H');refuse('molecular',name='OH');refuse('negative_excitation',energy=np.full_like(exc,-1));refuse('invalid_stage',ions=np.zeros_like(stages))
            status,count=exclude()
            if status or count<=0:raise ValueError('Trace selection failed/selected no records')
            selection=dict(element=atom.symbol,atom=checks[str(atom_file)],stages=stages.tolist(),excitation_eV=exc.tolist(),excluded_records=count)
        species=['O I','O II','CO','OH','H2O','H I','H II']
        for name in species:call('hop_keep_species',name.encode())
        call('hop_load_state',len(T),ptr(T),ptr(ne),ptr(rho),1.,1,1)
        if a.trace_selection:refuse('after_load')
        out={}
        for name in ['temp','ne','rho','chi500']:
            values=np.empty_like(T);call('hop_get_field',name.encode(),len(T),ptr(values));out[name]=values
        for name in species:
            values=np.empty_like(T);call('hop_get_species',name.encode(),len(T),ptr(values));out[name]=values
        if a.opacity_plan:
            from pipeline.wavelength_util import air_to_vac
            plan=json.loads(a.opacity_plan.read_text());checks[str(a.opacity_plan)]=fingerprint(a.opacity_plan)
            if not a.trace_selection or plan['frame']!='air_A'or plan['include_absmet']or plan['include_linelist']:raise ValueError('Require explicit trace-selected continuum-only AIR plan')
            air=np.asarray(plan['wavelengths'],dtype=float)
            if air.ndim!=1 or not 1<=len(air)<=20000 or not np.isfinite(air).all()or np.any(air<2000)or np.any(np.diff(air)<=0):raise ValueError('Invalid bounded continuum grid')
            vac=np.ascontiguousarray(air_to_vac(air));chi=np.zeros((len(T),len(vac)),np.float32);scat=np.zeros_like(chi);separate=np.zeros_like(chi);fptr=lambda v:v.ctypes.data_as(ct.POINTER(ct.c_float))
            call('hop_opacity',len(vac),ptr(vac),fptr(chi),fptr(scat),0,0,100.,0)
            call('hop_scattering',len(vac),ptr(vac),fptr(separate))
            if any(not np.isfinite(v).all()or np.any(v<0)for v in [chi,scat])or np.any(scat>chi):raise ValueError('Invalid total extinction/scattering decomposition')
            np.testing.assert_array_equal(scat,separate)
            for name in species:
                after=np.empty_like(T);call('hop_get_species',name.encode(),len(T),ptr(after));np.testing.assert_array_equal(after,out[name])
            ref_index=np.flatnonzero(air==5000.)
            opacity_checks=dict(wavelengths=len(vac),depths=len(T),scattering_bit_equal=True,species_unchanged_by_opacity=True,canonical5000_max_relative_difference=float(np.max(abs(chi[:,int(ref_index[0])]/state[:,6]-1)))if len(ref_index)else None)
            out.update(wave_air_A=air,wave_vac_A=vac,opacity_chi_cm_inv=chi.T,opacity_scat_cm_inv=scat.T)
        call('hop_free')
    finally:os.chdir(previous)
    errors={}
    for name in species:
        parts=read_patch_set(a.native_run/'native/save',name.replace(' ','_'));native=np.concatenate([parts[k].ravel(order='F')for k in sorted(parts)]).astype(float)
        errors[name]=float(np.max(abs(out[name]/np.maximum(native,np.finfo(float).tiny)-1)))
    np.testing.assert_array_equal(out['temp'],T);np.testing.assert_array_equal(out['ne'],original_ne);np.testing.assert_array_equal(out['rho'],rho)
    if any(not np.isfinite(x).all()or np.any(x<0)for x in out.values()):raise ValueError('Invalid standalone EOS result')
    np.savez_compressed(a.output/'arrays.npz',height_m=height,**out)
    for path,receipt in checks.items():
        if fingerprint(path)!=receipt:raise ValueError('Changed provider input')
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',engine=engine,inputs=checks,derived_build=derived,trace_selection=selection,negative_checks=negative_checks,opacity_checks=opacity_checks,temperature_ne_rho_bit_preserved=True,species_max_relative_differences=errors,chi500_max_relative_difference=float(np.max(abs(out['chi500']/state[:,6]-1))),limitations=['Standalone EOS preflight, not a Bride background/SE provider or model admission.','C signatures parsed from pinned committed header, no CFFI/environment install.','Native atom-overwrite atomic weight versus standalone table weights require explicit audit; no exact species-parity claim.','Original native species are saved float32; standalone returns float64.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['temperature_ne_rho_bit_preserved','species_max_relative_differences','chi500_max_relative_difference']}))


if __name__=='__main__':main()
