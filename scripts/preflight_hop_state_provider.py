"""Pinned HOP standalone EOS parity via declarations read from its own C header."""
import argparse,ctypes as ct,json,os,re,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def bind_header(library,header):
    types={'int':ct.c_int,'float':ct.c_float,'double':ct.c_double,'char *':ct.c_char_p,'const char *':ct.c_char_p}
    for name,base in [('int',ct.c_int),('float',ct.c_float),('double',ct.c_double)]:
        types[name+' *']=ct.POINTER(base);types['const '+name+' *']=ct.POINTER(base)
    bound={}
    for result,name,args in re.findall(r'(int|double)\s+(hop_\w+)\s*\(([^)]*)\);',header):
        params=[]
        if args.strip():
            for arg in args.split(','):
                kind=re.sub(r'\b[A-Za-z_]\w*\s*$','',arg).strip();kind=re.sub(r'\s+',' ',kind)
                if kind not in types:raise ValueError('Unsupported native ABI type '+kind)
                params.append(types[kind])
        fun=getattr(library,name);fun.argtypes=params;fun.restype=types[result];bound[name]=fun
    return bound


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import load_engine,fingerprint
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.openmulti3d_output import read_patch_set
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native-run',type=Path,required=True);p.add_argument('--native-analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Explicit one-thread environment required')
    _,app,_,engine=load_engine('openmulti3d');library=app/'HOP/standalone/build/libhop.so';header=app/'HOP/hopacity/libhop.h';files=[app/'input'/name for name in ['absdat/TS_absdat.txt','abund/abund_magg22.txt','isotopes/isotopes_blp.txt','molec/molec_data_BC.txt']];checks={str(x):fingerprint(x)for x in [library,header,*files]}
    native_report=a.native_run/'report.json';analysis_report=a.native_analysis/'report.json';run=json.loads(native_report.read_text());analysis=json.loads(analysis_report.read_text())
    if not run['computational_completed']or run['engine_receipt']!=engine or analysis['native_report']!=fingerprint(native_report):raise ValueError('Native run/analysis/engine identity mismatch')
    checks[str(native_report)]=fingerprint(native_report);checks[str(analysis_report)]=fingerprint(analysis_report)
    for entry in run['output_files']:
        path=a.native_run/entry['path']
        if fingerprint(path)!=entry['receipt']:raise ValueError('Changed native provider reference')
        checks[str(path)]=entry['receipt']
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
    previous=Path.cwd();os.chdir(a.output)
    try:
        call('hop_set_external_abundances',len(abund),ptr(abund));call('hop_set_option',b'abundances',b'O=8.70')
        call('hop_init',*[str(x).encode()for x in files],0.,-99.,0.,0.)
        species=['O I','O II','CO','OH','H2O','H I','H II']
        for name in species:call('hop_keep_species',name.encode())
        call('hop_load_state',len(T),ptr(T),ptr(ne),ptr(rho),1.,1,1)
        out={}
        for name in ['temp','ne','rho','chi500']:
            values=np.empty_like(T);call('hop_get_field',name.encode(),len(T),ptr(values));out[name]=values
        for name in species:
            values=np.empty_like(T);call('hop_get_species',name.encode(),len(T),ptr(values));out[name]=values
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
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',engine=engine,inputs=checks,temperature_ne_rho_bit_preserved=True,species_max_relative_differences=errors,chi500_max_relative_difference=float(np.max(abs(out['chi500']/state[:,6]-1))),limitations=['Standalone EOS preflight, not a Bride background/SE provider or model admission.','C signatures parsed from pinned committed header, no CFFI/environment install.','Native atom-overwrite atomic weight versus standalone table weights require explicit audit; no exact species-parity claim.','Original native species are saved float32; standalone returns float64.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['temperature_ne_rho_bit_preserved','species_max_relative_differences','chi500_max_relative_difference']}))


if __name__=='__main__':main()
