"""Build fresh isolated active-profile derivative; preserve original reference."""
import argparse,hashlib,json,os,platform,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import load_engine,fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin':raise ValueError('Authorized Mac-only derivative')
    spec,app,_,original=load_engine('openmulti3d');checkout=app.parents[1];patch=ROOT/'native/qlg_profile_reference/active_profile.patch';patch_receipt=fingerprint(patch)
    if a.source.exists()or a.source.resolve().is_relative_to(checkout.resolve())or a.output.resolve().is_relative_to(checkout.resolve()):raise ValueError('Require fresh private source/output outside reference')
    if not patch.read_bytes().startswith(b'diff --git '):raise ValueError('Require canonical Git patch before creating any worktree')
    a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');report=dict(generator=fingerprint(__file__),source_commit=spec['source_commit'],original_reference=original,patch=patch_receipt,source=str(a.source.resolve()),completed=False,threads=1,seconds_limit=600,scientific_admission='HOLD')
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    commands=[(['git','worktree','add','--detach',str(a.source.resolve()),spec['source_commit']],checkout),(['git','apply','--check',str(patch)],a.source),(['git','apply',str(patch)],a.source),(['make','-j1'],a.source/spec['application_rel'])]
    try:
        with(a.output/'build.log').open('w')as log:
            for command,cwd in commands:
                remaining=600-(time.perf_counter()-start)
                if remaining<=0:raise subprocess.TimeoutExpired(command,600)
                subprocess.run(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=remaining,check=True)
        diff=subprocess.check_output(['git','diff','--binary'],cwd=a.source)
        if diff!=patch.read_bytes():raise ValueError('Unexpected derived tracked diff')
        _,_,_,after=load_engine('openmulti3d')
        if after!=original or fingerprint(patch)!=patch_receipt:raise ValueError('Original reference/patch changed')
        derived=a.source/spec['application_rel'];report.update(completed=True,tracked_diff_sha256=hashlib.sha256(diff).hexdigest(),executable=fingerprint(derived/'dispatch.x'),libraries=dict(hop=fingerprint(derived/'HOP/standalone/build/libhop.so')),compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0])
    except (subprocess.CalledProcessError,subprocess.TimeoutExpired,ValueError)as error:report['error']=str(error)
    report['seconds']=time.perf_counter()-start;(a.output/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report.get(k)for k in ['completed','seconds','error']}))
    if not report['completed']:raise RuntimeError('Isolated profile build failed; output preserved')
if __name__=='__main__':main()
