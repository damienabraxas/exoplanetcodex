"""Link a read-only line-frequency-grid caller against pinned native HOP objects."""
import argparse,json,os,platform,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from config.constants import require_codex_path
    from pipeline.hop_packet import checked_runtime
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin'or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Require one-thread Mac')
    _,_,checks,_,_,_=checked_runtime('collider')
    base=require_codex_path('engines.hop_collider_provider_source')/'experiments/OpenMULTI3D/HOP';build=base/'standalone/build';driver=ROOT/'native/hop_line_grid_oracle/driver.f90';objects=sorted(build.glob('*.o'))
    if not objects or any(p.name in ['hop_prune.o','atom_convert.o']for p in objects):raise ValueError('Unexpected object set')
    for path in [driver,base/'standalone/Makefile',*sorted(base.glob('*.f90')),*objects,*sorted(build.glob('*.mod'))]:checks[str(path)]=fingerprint(path)
    a.output.mkdir(parents=True,exist_ok=False);exe=a.output/'line_grid_probe'
    command=['gfortran','-O3','-fopenmp','-I'+str(build),str(driver),*[str(p)for p in objects],'-o',str(exe.resolve())]
    report=dict(generator=fingerprint(__file__),inputs=checks,command=command,threads=1,seconds_limit=600,completed=False);start=time.perf_counter()
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    with(a.output/'compiler.log').open('w')as log:result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=600,env=dict(os.environ,OMP_NUM_THREADS='1'))
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed native oracle build input')
    report.update(completed=result.returncode==0,seconds=time.perf_counter()-start)
    if report['completed']:report['executable']=fingerprint(exe)
    (a.output/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['completed']:raise RuntimeError('Native driver build failed; evidence retained')
    print(json.dumps(dict(seconds=report['seconds'],executable=report['executable'])))

if __name__=='__main__':main()
