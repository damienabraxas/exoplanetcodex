"""Prepare or execute a pinned external reference engine in a fresh owned directory."""
import argparse,json,os,platform,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import load_plan_engine,prepare_run,fingerprint,validate_plan
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('This initial reference integration is authorized on the Mac')
    plan=json.loads(a.plan.read_text());validate_plan(plan);spec,app,exe,engine=load_plan_engine(plan);inputs=prepare_run(plan,a.output,app,resolver=require_codex_path)
    report=dict(schema='codex.reference_engine_result/1',scientific_admission='HOLD',engine=plan['engine'],mode=plan['mode'],profile_variant=plan.get('profile_variant','original'),scope=plan['scope'],generator=fingerprint(__file__),adapter=fingerprint(ROOT/'pipeline/reference_engines.py'),plan=fingerprint(a.plan),engine_receipt=engine,inputs=inputs,dry_run=not a.execute,computational_completed=False)
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    if not a.execute:
        (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('Dry preparation complete; no native solve');return
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_STACKSIZE='256M');start=time.perf_counter();error=None
    try:
        with(a.output/'native.log').open('w')as log:
            run=subprocess.run([str(exe),'input.nml'],cwd=a.output,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=plan['controls']['seconds'],check=False)
        report.update(returncode=run.returncode,computational_completed=run.returncode==0)
    except subprocess.TimeoutExpired as e:error='TimeoutExpired';report.update(error=error,computational_completed=False)
    report['seconds']=time.perf_counter()-start
    try:
        _,_,_,after=load_plan_engine(plan)
        if after!=engine or fingerprint(a.plan)!=report['plan']:raise ValueError('Engine/plan changed during run')
        _,app2,_,_=load_plan_engine(plan)
        for key,expected in [('atom',inputs['atom']),('atmosphere',inputs['atmosphere'])]:
            item=plan[key];path=require_codex_path(item['resource'])if key=='atom'and'resource'in item else app2/'input/atoms'/item['supplied_file']if key=='atom'else require_codex_path(item['resource'])/item['filename']
            if fingerprint(path)!=expected:raise ValueError('Physical input changed during reference solve')
    except Exception as e:report.update(computational_completed=False,protection_error=str(e));error=str(e)
    report['output_files']=[dict(path=str(f.relative_to(a.output)),receipt=fingerprint(f))for f in sorted(a.output.rglob('*'))if f.is_file()and not f.is_symlink()and 'input'not in f.relative_to(a.output).parts and f.name not in ['manifest.json','report.json']]
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(computational_completed=report['computational_completed'],seconds=report['seconds'],scientific_admission='HOLD')))
    if not report['computational_completed']:raise RuntimeError(error or 'Reference native execution failed; inspect preserved output')

if __name__=='__main__':main()
