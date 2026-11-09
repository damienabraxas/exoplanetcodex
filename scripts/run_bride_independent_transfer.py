"""Fixed-population vertical spectra for the exact held O benchmark, no new SE fit."""
import argparse,json,platform,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--plan',type=Path,default=ROOT/'config/experiments/rya1371_bride_independent_transfer.json');a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure
    from pipeline.bride_binding_controls import BindingControl
    from pipeline.bride_H_nlte_control import make_context,scale_initial_fractions
    from pipeline.bride_top_cap import converge_cap_radiation,emergent_rays
    from pipeline.bride_observables import flux_quadrature,integrate_flux
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac
    import lightweaver.LwCompiled as compiled
    from lightweaver.atomic_table import PeriodicTable
    if platform.system()!='Darwin':raise RuntimeError('Mac-only diagnostic')
    plan_path=a.plan;plan=json.loads(plan_path.read_text());source=ROOT/plan['source'];r=json.loads((source/'report.json').read_text());ref_path=require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(ref_path.read_text());native=Path(require_codex_path(plan.get('native_analysis_resource','experiments.openmulti3d_oxygen_analysis')))/'spectra.npz';checks=[]
    def check(path,receipt):
        if fingerprint(path)!=receipt:raise ValueError('Changed declared input '+str(path))
        checks.append((path,receipt))
    check(plan_path,fingerprint(plan_path));check(source/'report.json',plan['source_report']);check(source/'populations.npz',r['coupled']['populations']);check(source/'flux.csv',r['flux']);check(ref_path,plan['reference']);check(native,plan['native_spectra']);check(Path(compiled.__file__),r['backend']);check(ROOT/'pipeline/bride_H_nlte_control.py',r['context_helper'])
    atom_path=require_codex_path('experiments.bride_oxygen_atom');check(atom_path,r['atom']);model_path=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name'];check(model_path,ref['atmosphere_file']);structure=read_marcs_structure(model_path);state=BindingControl(**r['state']);controls=r['controls']
    with np.load(source/'populations.npz')as d:n=d['populations'].copy();s=d['nstar'].copy()
    fields=dict(local=n,held_multi_b=scale_initial_fractions(np.array(ref['reference']['departures'])*s,n.sum(axis=0)),LTE=s.copy())
    with np.load(native)as d:
        vacA=d['wave_vac_A'].copy();air=d['wave_air_A'].copy()
    mask=np.zeros(len(air),bool)
    for lo,hi in plan['windows_air_A']:mask|=(air>=lo)&(air<=hi)
    vac=vacA[mask]/10;wave=air[mask]
    if np.any(np.diff(vac)<=0):raise ValueError('Independent output wavelengths must increase')
    old=np.genfromtxt(source/'flux.csv',delimiter=',',names=True);oldvac=air_to_vac(old['wavelength_air_A'])/10;extra=np.unique([vac[0],vac[-1],oldvac[0],oldvac[-1]])
    a.output.mkdir(parents=True,exist_ok=False);started=time.perf_counter();records=[]
    for name,pop in fields.items():
        ctx,audits=make_context(state,structure,read_multi_atom(atom_path),rays=controls['rays'],qscale=controls['qscale'],extra_vac_nm=extra,cap_length_m=r['cap_length_m'],ng=False,fixed=pop)
        continuum,_=make_context(state,structure,read_multi_atom(atom_path),rays=controls['rays'],qscale=controls['qscale'],extra_vac_nm=extra,cap_length_m=r['cap_length_m'],ng=False,fixed=pop,continuum_only=True)
        rad=converge_cap_radiation(ctx,rtol=1e-8,maxiter=1000,deadline=started+plan['seconds']);crad=converge_cap_radiation(continuum,rtol=1e-8,maxiter=1000,deadline=started+plan['seconds'])
        _,I=emergent_rays(ctx,vac,np.array([1.]),np.array([1.]));_,Ic=emergent_rays(continuum,vac,np.array([1.]),np.array([1.]));I=np.asarray(I).squeeze();Ic=np.asarray(Ic).squeeze()
        if I.shape!=vac.shape or Ic.shape!=vac.shape:raise ValueError('Vertical intensity must have exactly one value per wavelength')
        if any(not np.isfinite(x).all()or np.any(x<=0)for x in [I,Ic]):raise ValueError('Invalid vertical intensity')
        np.savez_compressed(a.output/(name+'.npz'),wave_vac_A=vacA[mask],wave_air_A=wave,Inu_SI=I,continuum_Inu_SI=Ic,normalized=I/Ic,populations=pop)
        endpoint=None
        if name=='local':
            mu,weights=flux_quadrature(controls['flux_rays']);_,If=emergent_rays(ctx,oldvac,mu,weights);_,Icf=emergent_rays(continuum,oldvac,mu,weights);F=integrate_flux(If,mu,weights);Fc=integrate_flux(Icf,mu,weights);endpoint=dict(max_normalized_flux_difference=float(np.max(abs(F/Fc-old['normalized_flux']))),max_relative_continuum_difference=float(np.max(abs(Fc/old['continuum_Fnu_W_m2_Hz']-1))))
            if max(endpoint.values())>plan['endpoint_replay_screen']:raise ValueError('Original baseline endpoint replay failed')
        for context in [ctx,continuum]:np.testing.assert_array_equal(context.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,pop)
        records.append(dict(case=name,radiation=rad,continuum_radiation=crad,endpoint_replay=endpoint,SE_updated=False,populations=fingerprint(a.output/(name+'.npz'))));print(name,flush=True)
    for path,receipt in checks:
        if fingerprint(path)!=receipt:raise ValueError('Input changed during fixed transfer')
    doc=dict(schema='codex.bride_independent_transfer/1',scientific_admission='HOLD',generator=fingerprint(__file__),plan=fingerprint(plan_path),cases=records,seconds=time.perf_counter()-started,scope=plan['scope'],source_report=fingerprint(source/'report.json'),reference=fingerprint(ref_path),native_spectra=fingerprint(native))
    (a.output/'report.json').write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
