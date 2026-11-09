"""Offline conditioning/stationary iteration diagnostics, never nonlinear proof."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from pipeline.bride_reference import fingerprint


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 source=fingerprint(a.source)
 with np.load(a.source)as d:n=d['populations'].copy();physical=d['physical_1'].copy();mali=d['gamma_1'].copy()
 a.output.mkdir(parents=True,exist_ok=False);rows=[]
 for blend in [0.,.01,.05,.1,.25,.5,.75,1.]:
  for depth in range(n.shape[1]):
   nn=n[:,depth];r=physical[:,:,depth];op=(1-blend)*r+blend*mali[:,:,depth];row=int(nn.argmax())
   original=op/np.max(abs(op),axis=1)[:,None];original[row]=1.
   # Relative-population coordinates and row equilibration, same equations.
   matrix=op*nn[None,:];rhs=r*nn[None,:];scale=np.max(abs(matrix),axis=1)
   matrix/=scale[:,None];rhs/=scale[:,None]
   matrix[row]=nn/nn.sum();rhs[row]=matrix[row]
   error_map=np.eye(len(nn))-np.linalg.solve(matrix,rhs)
   radius=float(np.max(abs(np.linalg.eigvals(error_map))))
   rows.append(dict(mali_fraction=blend,depth=depth,conservation_row=row,row_only_condition=float(np.linalg.cond(original)),relative_population_condition=float(np.linalg.cond(matrix)),frozen_rate_error_spectral_radius=radius))
 with (a.output/'conditioning.csv').open('x')as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 summary=[]
 for blend in [0.,.01,.05,.1,.25,.5,.75,1.]:
  group=[r for r in rows if r['mali_fraction']==blend]
  summary.append(dict(mali_fraction=blend,unstable_layers=sum(r['frozen_rate_error_spectral_radius']>=1 for r in group),maximum_spectral_radius=max(r['frozen_rate_error_spectral_radius']for r in group),maximum_row_only_condition=max(r['row_only_condition']for r in group),maximum_relative_population_condition=max(r['relative_population_condition']for r in group)))
 if fingerprint(a.source)!=source:raise ValueError('Source changed during audit')
 report=dict(schema='codex.bride_solver_conditioning/1',scientific_admission='HOLD',generator=fingerprint(__file__),source=source,layers=n.shape[1],rows=len(rows),summary=summary,limitations=['Frozen radiation/rates and saved baseline only; not the derivative or spectral radius of the full coupled RT/SE map.','MALI is intended to precondition radiation feedback; frozen-rate instability alone does not prove full nonlinear instability.','Relative-population scaling changes numerical coordinates, not physical rate equations or source inputs.','No blend is selected for use and no new radiative-transfer solve or spectrum is produced.'],table=fingerprint(a.output/'conditioning.csv'))
 (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 print(json.dumps(report['summary']))

if __name__=='__main__':main()
