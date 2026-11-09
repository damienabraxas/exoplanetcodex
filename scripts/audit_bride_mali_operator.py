"""Bounded fixed-population operator audit; emits no admitted spectrum."""
import argparse,json, runpy, sys, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from pipeline import bride_coupled_mali_control as control
PLUGIN=None

def audit(ctx,atom,original,output,start):
 from config.constants import require_codex_path
 from pipeline.bride_native_spline import NativeTautsp,install_native_collision_control
 from pipeline.bride_charge_contract import install_author_mn_pair
 from pipeline.bride_top_cap import converge_cap_radiation
 from pipeline.bride_numerical_validation import physical_rate_matrix
 from pipeline.bride_reference import fingerprint
 plugin=None
 if PLUGIN:
  from pipeline.bride_overlap_plugin import load_overlap_plugin
  scheme,plugin=load_overlap_plugin(PLUGIN);ctx.set_fs_iter_scheme(scheme)
 oracle=NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
 install_native_collision_control(ctx,atom,oracle,source_directions=True,reversible_ch0=False)
 install_author_mn_pair(ctx,atom,require_codex_path('experiments.bride_barklem_oxygen'),oracle)
 active=ctx.activeAtoms[0]
 # Use the already completed native-baseline populations, not a failed seed.
 with np.load(ROOT/'data/results/rya1281/solar_matched_v1/run/baseline/populations.npz') as d:active.n[:]=d['populations']
 n=active.n.copy()
 converge_cap_radiation(ctx,rtol=1e-8,maxiter=100,deadline=start+120)
 arrays={'populations':n};rows=[]
 for mode in [True,False,True]:
  ctx.formal_sol_gamma_matrices(lambdaIterate=mode)
  gamma=active.Gamma.copy();physical=physical_rate_matrix(active)
  flow=lambda m:np.einsum('ijd,jd->id',m,n)
  gross=np.einsum('ijd,jd->id',abs(physical),n)
  key=str(len(rows));arrays['gamma_'+key]=gamma;arrays['physical_'+key]=physical
  rows.append(dict(lambda_iterate=mode,physical_residual=float(np.max(abs(flow(physical))/gross)),operator_action_discrepancy=float(np.max(abs(flow(gamma)-flow(physical))/gross)),operator_residual=float(np.max(abs(flow(gamma))/gross))))
 np.savez_compressed(output/'operator_audit.npz',**arrays)
 (output/'operator_audit.json').write_text(json.dumps(dict(scientific_admission='HOLD',generator=fingerprint(__file__),plugin=plugin,rows=rows),indent=2)+'\n')
 raise RuntimeError('Intentional stop after bounded fixed-population operator audit')

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('--plugin',type=Path);args=p.parse_args();PLUGIN=args.plugin
 control.solve_coupled_mali=audit
 sys.argv=[str(ROOT/'scripts/run_bride_solar_matched_mali.py'),'--case','baseline','--output',str(args.output)]
 runpy.run_path(sys.argv[0],run_name='__main__')
