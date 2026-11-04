"""Run both predeclared Mac Ng controls serially, preserving individual failures."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from importlib.metadata import version
    if platform.system()!='Darwin' or version('lightweaver')!='0.16.1':raise ValueError('Mac/LW0.16.1 only')
    if args.output.resolve()!=require_codex_path('experiments.bride_ng_controls').resolve():raise ValueError('Use registered fresh experiment output (override resource for replay)')
    plan_path=ROOT/'config/experiments/rya1281_ng_control.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_ng_experiment/1' or plan['modes']!=['record','off']:raise ValueError('Changed declared two-control experiment')
    store=require_codex_path('experiments.bride_mac_grid');coeff=require_codex_path('experiments.bride_thermal_coefficients')
    args.output.mkdir(parents=True,exist_ok=False);codes=[]
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    with (args.output/'attempts.jsonl').open('x') as receipt:
        for mode in plan['modes']:
            command=[sys.executable,str(ROOT/'scripts/run_bride_thermal_sweep.py'),
                '--plan',str(ROOT/plan['physical_plan']),'--case',plan['case'],
                '--refine-numerics','--warm-start-run',str(store/'rya1281-thermal-coupled-v1'/plan['case']),
                '--ng-control',mode,'--coefficients',str(coeff),
                '--case-directory',str(store/'rya1281-solar-bound-v1/solar5772_bound_Fe_sampled'),
                '--lines',str(store/'rya1281-trace-fe-lines-v2/lines.json'),
                '--pressure-source',str(store/'rya1281-fe-pressure-source-v2/report.json'),
                '--stark-source',str(store/'rya1281-fe-stark-source-v2/report.json'),
                '--output',str(args.output/mode)]
            start=time.perf_counter()
            with (args.output/(mode+'.log')).open('x') as log:
                completed=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            row=dict(mode=mode,returncode=completed.returncode,wall_seconds=time.perf_counter()-start,
                command=command,code_commit=commit,generator=fingerprint(__file__),plan=fingerprint(plan_path),
                log=fingerprint(args.output/(mode+'.log')),scientific_admission='HOLD')
            receipt.write(json.dumps(row,allow_nan=False)+'\n');receipt.flush();codes.append(completed.returncode)
            print(json.dumps(dict(mode=mode,returncode=completed.returncode)),flush=True)
    raise SystemExit(1 if any(codes) else 0)


if __name__=='__main__':main()
