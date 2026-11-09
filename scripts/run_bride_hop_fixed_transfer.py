"""Complete native-state fixed transfer with trace-selected HOP background; no SE."""
import argparse,json,os,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_engine
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_hop_transfer import make_context
    from pipeline.bride_top_cap import converge_cap_radiation,emergent_rays
    from lightweaver.atomic_table import PeriodicTable
    import lightweaver.LwCompiled as compiled
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    plan=json.loads(a.plan.read_text());checks={str(a.plan.resolve()):fingerprint(a.plan)}
    if fingerprint(Path(compiled.__file__))!=plan['backend']:raise ValueError('Changed qualified LW backend')
    options_file=ROOT/'config/experiments/rya1371_hop_native_options.json'
    if fingerprint(options_file)!=plan['native_options']:raise ValueError('Changed predeclared native module options')
    checks[str(options_file)]=plan['native_options']
    for relative in ['pipeline/bride_hop_transfer.py','pipeline/hop_packet.py','pipeline/hop_abi.py','pipeline/bride_top_cap.py','pipeline/multi_atom.py','pipeline/bride_reference.py']:
        path=ROOT/relative;checks[str(path)]=fingerprint(path)
    if os.environ.get('OMP_NUM_THREADS')!='1'or not 0<plan['seconds']<=600 or plan['purpose']!='fixed-native-oxygen-transfer':raise ValueError('Require explicit bounded one-thread fixed-transfer plan')
    native=Path(require_codex_path(plan['native_analysis_resource']));arrays=native/'arrays.npz';spectra=native/'spectra.npz';report_file=native/'report.json';r=json.loads(report_file.read_text());_,_,_,reference=load_engine('openmulti3d')
    if r['engine_receipt']!=reference or r['stored_rate_max_residual']>1e-6:raise ValueError('Invalid native reference state')
    for path,expected in [(arrays,plan['native_arrays']),(spectra,plan['native_spectra']),(report_file,plan['native_report'])]:
        if fingerprint(path)!=expected:raise ValueError('Changed native state identity')
        checks[str(path)]=expected
    with np.load(arrays)as d:f=d['consumed_state_cgs'][::-1].astype(float);height=d['height_m'][::-1].copy();n=d['native_populations_cm3'][:,::-1].copy()*1e6;s=d['native_nstar_cm3'][:,::-1].copy()*1e6
    with np.load(spectra)as d:vac=d['wave_vac_A'].copy();air=d['wave_air_A'].copy()
    if len(vac)>10000 or np.any(np.diff(vac)<=0):raise ValueError('Unbounded/invalid reference sampling')
    state=dict(height_m=height,temperature_K=f[:,2],ne_m3=f[:,0]*1e6,rho_kg_m3=f[:,1]*1e3,vturb_m_s=f[:,7]*1e3)
    atom_file=require_codex_path('experiments.bride_oxygen_atom');ref_file=require_codex_path('experiments.bride_oxygen_reference');held=json.loads(ref_file.read_text());model_file=require_codex_path('grids.marcs_standard')/held['atmosphere_file']['name']
    if r['reference']!=fingerprint(ref_file):raise ValueError('Changed declared native/held reference relation')
    if fingerprint(atom_file)!=held['atom']or fingerprint(model_file)!=held['atmosphere_file']:raise ValueError('Changed physical input')
    for path in [atom_file,ref_file,model_file,Path(compiled.__file__)]:checks[str(path)]=fingerprint(path)
    structure=read_marcs_structure(model_file);extra=np.unique([vac[0]/10,vac[-1]/10]);a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter();records=[]
    manifest=dict(generator=fingerprint(__file__),plan=fingerprint(a.plan),scientific_admission='HOLD',inputs=checks,scope=plan['scope'],dry_run=a.dry_run)
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for name,pop in [('native_NLTE',n),('native_LTE',s)]:
        ctx,audit=make_context(state,structure,read_multi_atom(atom_file),pop,s,extra_vac_nm=extra,formal_solver=plan.get('formal_solver','piecewise_besser_1d'))
        continuum,caudit=make_context(state,structure,read_multi_atom(atom_file),pop,s,extra_vac_nm=extra,continuum_only=True,formal_solver=plan.get('formal_solver','piecewise_besser_1d'))
        if a.dry_run:
            records.append(dict(case=name,audit=audit,continuum_audit=caudit,constructed_only=True));continue
        rad=converge_cap_radiation(ctx,rtol=1e-8,maxiter=1000,deadline=start+plan['seconds']);crad=converge_cap_radiation(continuum,rtol=1e-8,maxiter=1000,deadline=start+plan['seconds'])
        _,I=emergent_rays(ctx,vac/10,np.array([1.]),np.array([1.]));_,Ic=emergent_rays(continuum,vac/10,np.array([1.]),np.array([1.]));I=I[:,0];Ic=Ic[:,0]
        if I.shape!=vac.shape or Ic.shape!=vac.shape or any(not np.isfinite(x).all()or np.any(x<=0)for x in [I,Ic]):raise ValueError('Invalid fixed-transfer output')
        for context in [ctx,continuum]:
            atom=context.kwargs['eqPops'].atomicPops[PeriodicTable['O']];np.testing.assert_array_equal(atom.n,pop);np.testing.assert_array_equal(atom.nStar,s);np.testing.assert_array_equal(atom.nTotal,s.sum(axis=0))
        np.savez_compressed(a.output/(name+'.npz'),wave_vac_A=vac,wave_air_A=air,Inu_SI=I,continuum_Inu_SI=Ic,normalized=I/Ic,populations=pop,nstar=s)
        records.append(dict(case=name,radiation=rad,continuum_radiation=crad,audit=audit,continuum_audit=caudit,output=fingerprint(a.output/(name+'.npz')),SE_updated=False));print(name,flush=True)
    for path,expected in checks.items():
        if fingerprint(path)!=expected:raise ValueError('Changed fixed-transfer input')
    manifest.update(records=records,seconds=time.perf_counter()-start,limitations=['No population or collision/SE solve. Complete native O state, source geometry/thermodynamics, EOS and selected background are explicitly bound.','Native perturber densities with original LW broadening formulas; native versus LW atomic/profile/constants and lower-ghost/formal-solver conventions remain distinct.','Native historical-constant LTE source function retained; no empirical adjustment or physical-model/observation admission.'])
    (a.output/'report.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
