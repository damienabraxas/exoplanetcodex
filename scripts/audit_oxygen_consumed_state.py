"""Recompute Bride LTE/background at native consumed T/ne/rho/height, no RT/SE fit."""
import argparse,csv,json,sys
from dataclasses import asdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.bride_reference import read_marcs_structure,marcs_number_densities
    from pipeline.bride_binding_controls import BindingControl
    from pipeline.bride_H_nlte_control import make_context
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac
    import lightweaver.constants as lc
    from lightweaver.atomic_table import PeriodicTable
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    report=json.loads((a.native/'report.json').read_text());file=a.native/'arrays.npz';before=fingerprint(file)
    if report['stored_rate_max_residual']>1e-6:raise ValueError('Native stored rate failure')
    with np.load(file)as z:fields=z['consumed_state_cgs'][::-1].astype(float);height=z['height_m'][::-1];foreign=z['native_nstar_cm3'][:,::-1]*1e6
    ref=json.loads(require_codex_path('experiments.bride_oxygen_reference').read_text());structure=read_marcs_structure(require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name'])
    structure['temperature_K']=fields[:,2];structure['rho_kg_m3']=fields[:,1]*1e3;structure['electron_pressure_Pa']=fields[:,0]*1e6*lc.KBoltzmann*fields[:,2]
    ne,nH=marcs_number_densities(structure)
    state=BindingControl('model_mass_H',tuple(height),tuple(fields[:,2]),tuple(structure['electron_pressure_Pa']),tuple(ne),tuple(nH),tuple(fields[:,7]*1e3),'Native consumed ne (cm-3 to m-3), retained through Pe/T identity','Native consumed rho with original MARCS hydrogen mass fraction and pinned LW H mass')
    sample=float(air_to_vac(np.array([5000.]))[0]/10)
    ctx,_=make_context(state,structure,read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')),rays=5,qscale=1,extra_vac_nm=np.array([sample,845.]),cap_length_m=0.,ng=False)
    local=ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].nStar.copy();index=int(np.argmin(abs(np.asarray(ctx.spect.wavelength)-sample)))
    if ctx.spect.wavelength[index]!=sample:raise ValueError('Requested continuum sample absent')
    # FastBackground.compute_background already adds scattering to chi
    # (pinned Source/LwMiddleLayer.pyx); do not count it a second time.
    chi=np.asarray(ctx.background.chi)[index].copy();native_chi=fields[:,6]*100
    ratio=(local/local.sum(axis=0))/(foreign/foreign.sum(axis=0))
    rows=[dict(depth_top_to_bottom=k,height_above_source_bottom_m=height[k],temperature_K=fields[k,2],ne_m3=ne[k],rho_kg_m3=structure['rho_kg_m3'][k],model_mass_H_m3=nH[k],native_chi500_m_inv=native_chi[k],Bride_chi500_m_inv=chi[k],Bride_over_native_chi500=chi[k]/native_chi[k],Bride_over_native_total_O=local[:,k].sum()/foreign[:,k].sum())for k in range(len(height))]
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'state.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    np.savez_compressed(a.output/'arrays.npz',local_nstar_m3=local,native_nstar_m3=foreign,normalized_fraction_ratio=ratio,local_chi500_m_inv=chi,native_chi500_m_inv=native_chi)
    if fingerprint(file)!=before:raise ValueError('Changed native state')
    result=dict(generator=fingerprint(__file__),scientific_admission='HOLD',native_arrays=before,state=asdict(state),triplet_fraction_ratio_range=[float(np.min(ratio[[5,7,8,9]])),float(np.max(ratio[[5,7,8,9]]))],total_O_ratio_range=[float(np.min(local.sum(axis=0)/foreign.sum(axis=0))),float(np.max(local.sum(axis=0)/foreign.sum(axis=0)))],chi500_ratio_range=[float(np.min(chi/native_chi)),float(np.max(chi/native_chi))],limitations=['Diagnostic at native consumed thermodynamics; not a new SE or RT solution.','H normalization remains source MARCS mass convention, not a native EOS packet transplant.','Native reference opacity uses its internal AIR converter; Bride uses canonical converter. No wavelength-shift fit or exact wavelength-parity claim.','Background source tables and partition functions differ; fraction normalization is shown separately from raw total differences.'])
    (a.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k.endswith('range')}))


if __name__=='__main__':main()
