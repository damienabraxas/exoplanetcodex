"""Link a read-only continuum-interpolation caller against pinned native HOP objects."""
import argparse,json,os,platform,subprocess,sys,time
import numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from config.constants import require_codex_path
    from pipeline.hop_packet import checked_runtime
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin'or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Require one-thread Mac')
    _,_,checks,_,_,_=checked_runtime('collider')
    base=require_codex_path('engines.hop_collider_provider_source')/'experiments/OpenMULTI3D/HOP';build=base/'standalone/build';driver=ROOT/'native/hop_continuum_oracle/driver.f90';objects=sorted(build.glob('*.o'))
    if not objects or any(p.name in ['hop_prune.o','atom_convert.o']for p in objects):raise ValueError('Unexpected object set')
    for path in [driver,base/'standalone/Makefile',*sorted(base.glob('*.f90')),*objects,*sorted(build.glob('*.mod'))]:checks[str(path)]=fingerprint(path)
    a.output.mkdir(parents=True,exist_ok=False);exe=a.output/'continuum_probe'
    command=['gfortran','-O3','-fopenmp','-I'+str(build),str(driver),*[str(p)for p in objects],'-o',str(exe.resolve())]
    report=dict(generator=fingerprint(__file__),inputs=checks,command=command,threads=1,seconds_limit=600,completed=False);start=time.perf_counter()
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    with(a.output/'compiler.log').open('w')as log:result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=600,env=dict(os.environ,OMP_NUM_THREADS='1'))
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed native oracle build input')
    report.update(completed=result.returncode==0,seconds=time.perf_counter()-start)
    if report['completed']:
        report['executable']=fingerprint(exe)
        wave=np.linspace(1000.,2000.,5);query=np.unique(np.r_[999.,np.linspace(1000.,2000.,17),2001.])
        controls=[]
        for label,alpha in [('constant',np.full(5,2e-18)),('linear_frequency',(2+2.99792458e18/wave/1e15)*1e-18)]:
            request=f'{len(wave)} {len(query)}\n'+'\n'.join(f'{w:.17e} {v:.17e}'for w,v in zip(wave,alpha))+'\n'+'\n'.join(f'{w:.17e}'for w in query)+'\n'
            (a.output/(label+'.in')).write_text(request)
            probe=subprocess.run([str(exe.resolve())],input=request,text=True,capture_output=True,timeout=30,check=True,env=dict(os.environ,OMP_NUM_THREADS='1'))
            (a.output/(label+'.out')).write_text(probe.stdout)
            actual=np.fromstring(probe.stdout,sep=' ')
            expected=np.full(len(query),2e-18)if label=='constant'else(2+2.99792458e18/query/1e15)*1e-18
            expected[(query<1000)|(query>2000)]=0
            np.testing.assert_allclose(actual,expected,rtol=1e-7,atol=0)
            controls.append(dict(case=label,queries=len(query),max_relative_nonzero_error=float(np.max(abs(actual[expected>0]/expected[expected>0]-1))),outside_support_exact_zero=bool(np.all(actual[expected==0]==0))))
        report['analytic_controls']=controls
    (a.output/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['completed']:raise RuntimeError('Native driver build failed; evidence retained')
    print(json.dumps(dict(seconds=report['seconds'],executable=report['executable'])))

if __name__=='__main__':main()
