"""Actual compiled cap-response diagnostic versus independent derivative checks."""
import argparse,json,sys
from pathlib import Path
from types import SimpleNamespace as S
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plugin',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 import lightweaver as lw
 from lightweaver.LwCompiled import BackgroundProvider
 from scipy.special import exprel
 from pipeline.bride_top_cap import FiniteColumnCap,cap_intensity
 from pipeline.bride_overlap_plugin import load_overlap_plugin
 from pipeline.bride_reference import fingerprint
 scheme,plugin=load_overlap_plugin(a.plugin)
 if scheme!='bride_full_coupling_cap_scalar':raise ValueError('Require explicit cap-response plugin')
 a.output.mkdir(parents=True,exist_ok=False);chi_values=np.array([1e-8,.1,2.,100.]);eta_values=3*chi_values;rows=[];arrays={}
 class Bottom(lw.BoundaryCondition):
  def compute_bc(self,atmos,spect):return np.full((len(spect.wavelength),len(self.muz),1),9.)
 class Coefficients(BackgroundProvider):
  def compute_background(self,atmos,chi,eta,sca):chi[:]=chi_values[:,None];eta[:]=eta_values[:,None];sca[:]=0.
 for length in [0.,.1,3.,30.]:
  cap=FiniteColumnCap(length);atmos=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=np.linspace(1.,0.,9),temperature=np.full(9,6000.),ne=np.full(9,1e16),nHTot=np.full(9,1e20),vturb=np.zeros(9),vlos=np.zeros(9),logG=2.,lowerBc=Bottom(),upperBc=cap);atmos.quadrature(2)
  spect=S(wavelength=np.array([400.,500.,600.,700.]),radSet=S(activeAtoms=[],detailedAtoms=[]));ctx=lw.Context(atmos,spect,{},backgroundProvider=Coefficients,Nthreads=1,formalSolver='piecewise_besser_1d',fsIterScheme=scheme);ctx.depthData.fill=True;cap=ctx.atmos.pyAtmos.zUpperBc
  for _ in range(3):
   top=np.array(ctx.spect.J)[:,0].copy();ctx.formal_sol_gamma_matrices(lambdaIterate=True);cap.update_from_context(ctx,top)
  old=np.array(ctx.depthData.I).copy();response=np.zeros((4,2));ctx.formal_sol_gamma_matrices(lambdaIterate=False,extraParams={'BrideCapLengthM':length,'BrideCapResponse':response})
  np.testing.assert_array_equal(np.array(ctx.depthData.I),old)
  mu=np.asarray(atmos.muz);distance=length/mu;analytic=exprel(-chi_values[:,None]*distance[None,:])*distance[None,:]
  step=2.**-20;finite=(cap_intensity(chi_values,eta_values*(1+step),length,mu)-cap_intensity(chi_values,eta_values*(1-step),length,mu))/(2*step*eta_values[:,None])
  np.testing.assert_allclose(response,analytic,rtol=64*np.finfo(float).eps,atol=0.)
  np.testing.assert_allclose(response,finite,rtol=128*np.finfo(float).eps/step,atol=0.)
  np.testing.assert_allclose(np.array(ctx.depthData.I)[:,:,0,0],cap_intensity(chi_values,eta_values,length,mu),rtol=64*np.finfo(float).eps,atol=0.)
  key=str(len(rows));arrays['compiled_'+key]=response;arrays['analytic_'+key]=analytic;arrays['finite_difference_'+key]=finite
  rows.append(dict(length_m=length,max_absolute_analytic_error=float(np.max(abs(response-analytic))),max_absolute_finite_difference_error=float(np.max(abs(response-finite))),physical_intensity_exact=True))
 np.savez(a.output/'responses.npz',**arrays)
 result=dict(schema='codex.bride_cap_response/1',scientific_admission='HOLD',generator=fingerprint(__file__),plugin=plugin,chi_m_inv=chi_values.tolist(),rows=rows,scope='Fixed-opacity derivative of actual incoming slab intensity with respect to cap emissivity; zero/thin/thick and angular controls. No abundance, atmosphere or physical-model admission.')
 (a.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(rows))

if __name__=='__main__':main()
