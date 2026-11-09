"""Frozen-matrix compiled versus relative-population linear-solver comparison."""
import argparse,json,runpy,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from pipeline import bride_coupled_mali_control as control


class AuditComplete(Exception):
    pass


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 from pipeline.bride_reference import fingerprint
 source=fingerprint(a.source)
 def audit(ctx,atom,original,output,start):
  with np.load(a.source)as d:n=d['populations'].copy();matrix=d['gamma_1'].copy()
  active=ctx.activeAtoms[0]
  np.testing.assert_allclose(n.sum(axis=0),active.nTotal,rtol=64*np.finfo(float).eps,atol=0)
  active.n[:]=n;scale=np.max(abs(matrix),axis=1,keepdims=True);active.Gamma[:]=matrix/scale
  update=ctx.stat_equil()
  if any(update.ngAccelerated):raise ValueError('Linear comparison accidentally accelerated')
  compiled=active.n.copy();row_results=[];relative_results=[];rows=[]
  for k in range(n.shape[1]):
   row=int(n[:,k].argmax());total=float(active.nTotal[k]);m=matrix[:,:,k].copy();m/=np.max(abs(m),axis=1)[:,None];m[row]=1.;rhs=np.zeros(len(n));rhs[row]=total
   row_solution=np.linalg.solve(m,rhs);row_results.append(row_solution)
   q=matrix[:,:,k]*n[:,k][None,:];q/=np.max(abs(q),axis=1)[:,None];q[row]=n[:,k]/total;rhs[row]=1.
   relative_solution=np.linalg.solve(q,rhs)*n[:,k];relative_results.append(relative_solution)
   residual=lambda x:float(np.max(abs(matrix[:,:,k]@x)/(abs(matrix[:,:,k])@abs(x))))
   rows.append(dict(depth=k,compiled_vs_relative_fraction=float(np.max(abs(compiled[:,k]/relative_solution-1))),numpy_row_vs_relative_fraction=float(np.max(abs(row_solution/relative_solution-1))),compiled_backward_residual=residual(compiled[:,k]),relative_backward_residual=residual(relative_solution)))
  np.savez_compressed(output/'linear_solutions.npz',initial=n,compiled=compiled,numpy_row=np.array(row_results).T,numpy_relative=np.array(relative_results).T)
  if fingerprint(a.source)!=source:raise ValueError('Input matrix changed during comparison')
  result=dict(schema='codex.bride_population_linear_audit/1',scientific_admission='HOLD',generator=fingerprint(__file__),source=source,rows=rows,max_compiled_vs_relative_fraction=max(r['compiled_vs_relative_fraction']for r in rows),max_numpy_row_vs_relative_fraction=max(r['numpy_row_vs_relative_fraction']for r in rows),scope='Same saved approximate matrix and O total; one unaccelerated compiled SE solve versus two NumPy coordinate systems. No RT, physical rate or atmosphere/model admission.')
  (output/'linear_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k.startswith('max_')}),flush=True)
  raise AuditComplete('Intentional end of frozen linear-algebra control')
 control.solve_coupled_mali=audit
 sys.argv=[str(ROOT/'scripts/run_bride_solar_matched_mali.py'),'--case','baseline','--output',str(a.output)]
 try:runpy.run_path(sys.argv[0],run_name='__main__')
 except AuditComplete:pass

if __name__=='__main__':main()
