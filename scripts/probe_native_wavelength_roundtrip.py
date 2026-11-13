"""Measure native internal table roundtrip through actual typed library functions."""
import argparse,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint,load_engine
    from pipeline.multi_atom import read_multi_atom
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin'or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Require one-thread Mac probe')
    _,app,_,engine=load_engine('openmulti3d');build=app/'HOP/standalone/build';lib=build/'libhop.so';module=build/'handyfuncs_mod.mod';driver=ROOT/'native/airvac_profile_probe/probe.f90';compiler=Path(shutil.which('gfortran'));atom_file=require_codex_path('experiments.bride_oxygen_atom');table=require_codex_path('experiments.openmulti3d_oxygen_linear_run')/'native/atoms/O.txt';atom=read_multi_atom(atom_file);text=table.read_text().splitlines();start=next(i for i,x in enumerate(text)if x.split()[:2]==['i','j']and'lam_air'in x)+1;grid=[]
    for index,line in enumerate(atom.lines):
        f=text[start+index].split()
        if tuple(map(int,f[:3]))!=(index,line.i-1,line.j-1):raise ValueError('Native table line identity mismatch')
        grid.append([1e8/(atom.levels[line.j-1].E_cm-atom.levels[line.i-1].E_cm),float(f[3])])
    inputs={str(x):fingerprint(x)for x in [lib,module,driver,compiler,atom_file,table,app/'HOP/handyfuncs_mod.f90']};a.output.mkdir(parents=True,exist_ok=False);exe=a.output.resolve()/'airvac_probe';start_time=time.perf_counter();manifest=dict(generator=fingerprint(__file__),engine=engine,inputs=inputs,threads=1,seconds_limit=600,scientific_admission='HOLD');(a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    command=[str(compiler),'-O2','-I',str(build),str(driver),str(lib),'-Wl,-rpath,'+str(build),'-o',str(exe)]
    with(a.output/'build.log').open('w')as f:subprocess.run(command,env=dict(os.environ,OMP_NUM_THREADS='1'),stdout=f,stderr=subprocess.STDOUT,timeout=600,check=True)
    np.savetxt(a.output/'input.txt',grid,fmt='%.17e')
    with(a.output/'input.txt').open()as fi,(a.output/'output.txt').open('w')as fo,(a.output/'run.log').open('w')as log:subprocess.run([str(exe)],cwd=build.parent,stdin=fi,stdout=fo,stderr=log,env=dict(os.environ,OMP_NUM_THREADS='1'),timeout=600-(time.perf_counter()-start_time),check=True)
    d=np.loadtxt(a.output/'output.txt');np.testing.assert_allclose(d[:,:2],grid,rtol=1e-15,atol=1e-15)
    if d.shape!=(256,5)or not np.isfinite(d).all():raise ValueError('Invalid native wavelength probe')
    rows={}
    for name,mask in [('FUV_vacuum_canonical',d[:,0]<2000),('AIR_canonical',d[:,0]>=2000)]:
        rows[name]=dict(n=int(mask.sum()),max_native_internal_unrounded_roundtrip_abs_A=float(np.max(abs(d[mask,4]-d[mask,0]))),max_consumed_rounded_table_to_source_vac_abs_A=float(np.max(abs(d[mask,2]-d[mask,0]))),max_consumed_relative_vac_difference=float(np.max(abs(d[mask,2]/d[mask,0]-1))))
    for path,expected in inputs.items():
        if fingerprint(path)!=expected:raise ValueError('Changed roundtrip source')
    _,_,_,after=load_engine('openmulti3d')
    if after!=engine:raise ValueError('Original engine changed')
    report=dict(**manifest,seconds=time.perf_counter()-start_time,executable=fingerprint(exe),output=fingerprint(a.output/'output.txt'),results=rows,limitations=['Native internal extrapolated air-table encoding is not a valid atmospheric-air coordinate below2000A; project canonical observation coordinates there remain vacuum.','Roundtrip/line-centre test only; no FUV opacity, transition/rate, continuum or all-band science admission.','Actual source photon frequency agreement must be distinguished from large internal table-coordinate differences. No source, atom or original binary change.']);(a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(rows))
if __name__=='__main__':main()
