"""Fresh bounded private HOP derivative build; original reference remains immutable."""
import argparse,hashlib,json,os,platform,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import load_engine,fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--patch',type=Path,default=ROOT/'native/hop_trace_provider/trace_continua.patch');a=p.parse_args()
    if platform.system()!='Darwin':raise ValueError('Authorized Mac-only build')
    spec,app,_,reference=load_engine('openmulti3d');checkout=app.parents[1];patch=a.patch.resolve();patch_receipt=fingerprint(patch)
    if a.source.exists()or a.source.resolve().is_relative_to(checkout.resolve())or a.output.resolve().is_relative_to(checkout.resolve()):raise ValueError('Require fresh source and output outside reference checkout')
    a.output.mkdir(parents=True,exist_ok=False);started=time.perf_counter();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');report=dict(generator=fingerprint(__file__),upstream_commit=spec['source_commit'],patch=patch_receipt,reference=reference,source=str(a.source.resolve()),completed=False,seconds_limit=600,threads=1)
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    commands=[(['git','worktree','add','--detach',str(a.source.resolve()),spec['source_commit']],checkout),(['git','apply','--check',str(patch)],a.source),(['git','apply',str(patch)],a.source),(['make','-j1','bindings','build/libhop.so'],a.source/spec['application_rel']/'HOP/standalone')]
    try:
        with(a.output/'build.log').open('w')as log:
            for command,cwd in commands:
                elapsed=time.perf_counter()-started
                if elapsed>=600:raise subprocess.TimeoutExpired(command,600)
                result=subprocess.run(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600-elapsed)
                if result.returncode:raise subprocess.CalledProcessError(result.returncode,command)
        diff=subprocess.check_output(['git','diff','--binary'],cwd=a.source)
        if diff!=patch.read_bytes():raise ValueError('Unexpected derived tracked source changes')
        _,_,_,after=load_engine('openmulti3d')
        if after!=reference or fingerprint(patch)!=patch_receipt:raise ValueError('Original reference/patch changed during build')
        base=a.source/spec['application_rel']/'HOP'
        report.update(completed=True,tracked_diff_sha256=hashlib.sha256(diff).hexdigest(),library=fingerprint(base/'standalone/build/libhop.so'),header=fingerprint(base/'hopacity/libhop.h'),compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0])
    except (subprocess.CalledProcessError,subprocess.TimeoutExpired,ValueError)as error:
        report['error']=str(error)
    report['seconds']=time.perf_counter()-started;(a.output/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    if not report['completed']:raise RuntimeError('Private build failed; output preserved')


if __name__=='__main__':main()
