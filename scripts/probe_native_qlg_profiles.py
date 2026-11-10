"""Bounded Mac-only independent comparison of actual pinned QLG-bot profile kernel."""
import argparse,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint,load_engine
    from scipy.special import wofz
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--damping-samples',type=Path);a=p.parse_args()
    if platform.system()!='Darwin'or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Require one-thread Mac probe')
    _,app,_,reference=load_engine('openmulti3d');source=app/'HOP/qlg_mod.f90';driver=ROOT/'native/qlg_profile_probe/probe.f90';compiler=Path(shutil.which('gfortran'));inputs={str(x):fingerprint(x)for x in [source,driver,compiler]};a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter();exe=a.output.resolve()/'profile_probe'
    damp=np.logspace(-5,1,61);v=np.linspace(-20,20,801)
    if a.damping_samples is not None:
        inputs[str(a.damping_samples.resolve())]=fingerprint(a.damping_samples)
        with np.load(a.damping_samples)as d:damp=np.unique(d['a_damp'])
        if len(damp)>560 or not np.isfinite(damp).all()or np.any(damp<=0):raise ValueError('Invalid/budget-exceeding damping support')
        v=np.linspace(-20,20,401)
    manifest=dict(generator=fingerprint(__file__),reference=reference,inputs=inputs,compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],seconds_limit=600,threads=1,scientific_admission='HOLD')
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');cmd=[str(compiler),'-O2','-fno-fast-math','-J',str(a.output.resolve()),'-o',str(exe),str(source),str(driver)]
    with(a.output/'build.log').open('w')as f:subprocess.run(cmd,cwd=a.output,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=600,check=True)
    aa,vv=np.meshgrid(damp,v,indexing='ij');grid=np.column_stack([aa.ravel(),vv.ravel()]);np.savetxt(a.output/'input.txt',grid,fmt='%.17e')
    with(a.output/'input.txt').open()as fi,(a.output/'profiles.txt').open('w')as fo:subprocess.run([str(exe)],cwd=a.output,env=env,stdin=fi,stdout=fo,timeout=600-(time.perf_counter()-start),check=True)
    data=np.loadtxt(a.output/'profiles.txt');np.testing.assert_allclose(data[:,:2],grid,rtol=1e-15,atol=1e-15);truth=wofz(data[:,1]+1j*data[:,0]).real;bot=data[:,2];native_ref=data[:,3];mask=truth>1e-8
    if not np.isfinite(data).all()or np.any(bot<0):raise ValueError('Invalid native kernel output')
    bot_area=np.trapezoid(bot.reshape(aa.shape),v,axis=1);ref_area=np.trapezoid(truth.reshape(aa.shape),v,axis=1)
    np.savez_compressed(a.output/'profiles.npz',grid_a_v=grid,native_bot=bot,native_fp64_reference=native_ref,scipy_reference=truth)
    report=dict(**manifest,finite_v_integral_max_relative_difference=float(np.max(abs(bot_area/ref_area-1))),relative_error_floor_H=1e-8,seconds=time.perf_counter()-start,executable=fingerprint(exe),profiles=fingerprint(a.output/'profiles.txt'),shape=list(aa.shape),a_range=[float(damp.min()),float(damp.max())],damping_samples=None if a.damping_samples is None else str(a.damping_samples.resolve()),v_range=[-20,20],bot_max_abs_error=float(np.max(abs(bot-truth))),bot_max_relative_above_floor=float(np.max(abs(bot[mask]/truth[mask]-1))),native_fp64_ref_max_abs_error=float(np.max(abs(native_ref-truth))),native_fp64_ref_max_relative_above_floor=float(np.max(abs(native_ref[mask]/truth[mask]-1))),limitations=[('Kernel-only generic diagnostic grid, not an emergent-intensity bound.'if a.damping_samples is None else'Kernel sampled on actual matched LW damping support; native stored damping parity and emergent-intensity error remain unqualified.'),'Actual upstream QLG-bot fp32 and upstream fp64 reference compared with independent SciPy wofz. RT frequency Jacobians, profile normalization and broadening inputs require separate accounting.','No source, native binary, profile default, atom, population or production change.'])
    for path,f in inputs.items():
        if fingerprint(path)!=f:raise ValueError('Changed probe input')
    _,_,_,after=load_engine('openmulti3d')
    if after!=reference:raise ValueError('Reference source/binary changed')
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['seconds','bot_max_abs_error','bot_max_relative_above_floor','native_fp64_ref_max_abs_error']}))
if __name__=='__main__':main()
