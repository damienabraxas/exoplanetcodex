"""Run the isolated, serial native physical-rate observer; no science admission."""
import argparse,json,os,platform,subprocess,sys,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def verify(source,build):
    from pipeline.reference_engines import fingerprint,load_plan_engine
    spec,_,_,original=load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'))
    receipt=json.loads((build/'receipt.json').read_text())
    if not receipt['completed'] or receipt['source_commit']!=spec['source_commit']:raise ValueError('Incomplete/unpinned build')
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip()!=spec['source_commit']:raise ValueError('Wrong source pin')
    diff=subprocess.check_output(['git','diff','--binary'],cwd=source)
    if hashlib.sha256(diff).hexdigest()!=receipt['tracked_diff_sha256']:raise ValueError('Changed derivative source')
    for key,relative in [('patch','native/native_gross_rates/gross_rates.patch'),('profile_patch','native/qlg_profile_reference/active_profile.patch')]:
        if fingerprint(ROOT/relative)!=receipt[key]:raise ValueError('Changed patch')
    for relative,expected in receipt['changed_sources'].items():
        if fingerprint(source/relative)!=expected:raise ValueError('Changed diagnostic implementation')
    app=source/spec['application_rel'];exe=app/'dispatch.x'
    if fingerprint(exe)!=receipt['executable']or fingerprint(app/'HOP/standalone/build/libhop.so')!=receipt['libraries']['hop']:raise ValueError('Changed executable/HOP')
    return app,exe,dict(build=fingerprint(build/'receipt.json'),original_accurate_engine=original,executable=fingerprint(exe))

def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import prepare_run,fingerprint
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','build','plan','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--execute',action='store_true');a=p.parse_args()
    if platform.system()!='Darwin'or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Mac-only one-thread observer')
    plan=json.loads(a.plan.read_text());app,exe,engine=verify(a.source,a.build)
    from pipeline.native_reference_controls import validate_observer_plan
    if plan['controls'].get('decouple_continuum')is False:
        from pipeline.native_reference_controls import validate_shared_radiation_plan
        validate_shared_radiation_plan(plan,json.loads((ROOT/'config/experiments/rya1371_openmulti3d_oxygen_frequency4_accurate.json').read_text()))
    else:validate_observer_plan(plan,json.loads((ROOT/'config/experiments/rya1371_openmulti3d_oxygen_accurate_profile.json').read_text()))
    inputs=prepare_run(plan,a.output,app,resolver=require_codex_path)
    report=dict(generator=fingerprint(__file__),plan=fingerprint(a.plan),engine=engine,inputs=inputs,purpose='stage-matched-physical-rate-observer',scientific_admission='HOLD',computational_completed=False,dry_run=not a.execute)
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    start=time.perf_counter()
    try:
        if a.execute:
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_STACKSIZE='256M')
            with (a.output/'native.log').open('w')as log:r=subprocess.run([str(exe),'input.nml'],cwd=a.output,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600,check=False)
            report.update(returncode=r.returncode,computational_completed=r.returncode==0)
        if verify(a.source,a.build)[2]!=engine:raise ValueError('Observer identities changed')
    except (ValueError,subprocess.TimeoutExpired)as error:report.update(computational_completed=False,error=str(error))
    finally:
        report['seconds']=time.perf_counter()-start
        report['output_files']=[dict(path=str(q.relative_to(a.output)),receipt=fingerprint(q))for q in sorted(a.output.rglob('*'))if q.is_file()and not q.is_symlink()and'input'not in q.relative_to(a.output).parts and q.name not in ['manifest.json','report.json']]
        (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['computational_completed','seconds','error']if k in report}))
    if a.execute and not report['computational_completed']:raise RuntimeError('Observer failed; preserved evidence')
if __name__=='__main__':main()
