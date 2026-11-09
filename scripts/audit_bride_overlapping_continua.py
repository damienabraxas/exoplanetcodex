"""Minimal compiled-backend thermal control for overlapping ion continua.

Synthetic three-level oxygen-labelled atom, no scientific oxygen interpretation,
custom collisions, finite cap, NLTE seed, source atom, interpolation or Ng.
"""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 import lightweaver as lw
 import lightweaver.LwCompiled as compiled
 from lightweaver.atomic_model import AtomicModel,AtomicLevel,HydrogenicContinuum
 from lightweaver.atomic_table import Element
 from lightweaver.collisional_rates import CI,CE
 from lightweaver.rh_atoms import H_6_atom
 from pipeline.bride_numerical_validation import physical_rate_matrix
 from pipeline.bride_reference import fingerprint
 a.output.mkdir(parents=True,exist_ok=False);rows=[];arrays={};started=time.perf_counter()
 for variant,upper in [('single_ground',[1]),('single_excited',[2]),('overlapping',[1,2])]:
  atmosphere=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=np.linspace(1e5,0,9),temperature=np.full(9,6000.),ne=np.full(9,1e16),nHTot=np.full(9,1e20),vlos=np.zeros(9),vturb=np.zeros(9),logG=2.,lowerBc=lw.ThermalisedRadiation(),upperBc=lw.ThermalisedRadiation())
  atmosphere.quadrature(2)
  levels=[AtomicLevel(E=0.,g=1.,label='synthetic neutral',stage=0),AtomicLevel(E=9e4,g=2.,label='synthetic ion ground',stage=1),AtomicLevel(E=1.1e5,g=4.,label='synthetic ion excited',stage=1)]
  model=AtomicModel(element=Element(8),levels=levels,lines=[],continua=[HydrogenicContinuum(j=j,i=0,NlambdaGen=41,alpha0=1e-20,minWavelength=50.)for j in upper],collisions=[CI(j=1,i=0,temperature=[3000.,10000.],rates=[1e-15]*2),CE(j=2,i=1,temperature=[3000.,10000.],rates=[1e-15]*2)])
  rad=lw.RadiativeSet([H_6_atom(),model]);rad.set_active('O');eq=rad.compute_eq_pops(atmosphere)
  spect=rad.compute_wavelength_grid(extraWavelengths=np.array([500.,600.]));ctx=lw.Context(atmosphere,spect,eq,Nthreads=1,formalSolver='piecewise_besser_1d',conserveCharge=False,ngOptions=None)
  active=ctx.activeAtoms[0];n=active.nStar.copy();active.n[:]=n
  for iteration in range(100):
   if time.perf_counter()-started>60:raise TimeoutError('Synthetic overlap audit budget exceeded')
   update=ctx.formal_sol_gamma_matrices(lambdaIterate=True)
   if iteration>=2 and update.dJMax<1e-12:break
  else:raise RuntimeError('Fixed LTE radiation failed convergence')
  for mode in [True,False,True]:
   ctx.formal_sol_gamma_matrices(lambdaIterate=mode)
   physical=physical_rate_matrix(active);operator=active.Gamma.copy();flow=lambda m:np.einsum('ijd,jd->id',m,n);gross=np.einsum('ijd,jd->id',abs(physical),n)
   key=str(len(rows));arrays['populations_'+key]=n;arrays['physical_'+key]=physical;arrays['operator_'+key]=operator
   rows.append(dict(variant=variant,lambda_iterate=mode,radiation_iterations=iteration+1,physical_residual=float(np.max(abs(flow(physical))/gross)),operator_action_discrepancy=float(np.max(abs(flow(operator)-flow(physical))/gross)),population_change=float(np.max(abs(active.n/n-1)))))
   np.testing.assert_array_equal(active.n,n)
 np.savez_compressed(a.output/'operators.npz',**arrays)
 report=dict(schema='codex.bride_overlap_thermal_control/1',scientific_admission='HOLD',generator=fingerprint(__file__),backend=fingerprint(compiled.__file__),rows=rows,arrays=fingerprint(a.output/'operators.npz'),seconds=time.perf_counter()-started,scope=__doc__)
 (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(rows))

if __name__=='__main__':main()
