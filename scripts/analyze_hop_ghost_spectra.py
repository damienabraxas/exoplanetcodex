"""Forward coupled O populations on the exact native optical grid, without overrides."""
import argparse,csv,json,os,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.multi_width_law import load_width_law
    from pipeline.bride_hop_transfer import make_context
    from pipeline.native_ghost_boundary import NativeLinearGhost,converge_ghost_radiation as converge_cap_radiation
    from pipeline.bride_observables import window_integral
    from lightweaver.atomic_table import PeriodicTable
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['plan','coupled-run','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());rp=a.coupled_run/'report.json';run=json.loads(rp.read_text())
    if os.environ.get('OMP_NUM_THREADS')!='1'or not run['computational_completed']or run['dry_run']or run['fresh_physical_rate_residual']>=run['controls']['rtol']or 'table_control'in plan:raise ValueError('Require completed full-source coupled diagnostic')
    checks={str(rp):fingerprint(rp),str(a.plan.resolve()):fingerprint(a.plan)}
    for path,fp in run['inputs'].items():
        if fingerprint(path)!=fp:raise ValueError('Changed coupled source '+path)
        checks[path]=fp
    native=require_codex_path(plan['native_analysis_resource']);arrays=native/'arrays.npz';spectra=native/'spectra.npz';pp=a.coupled_run/'populations.npz'
    for path,expected in [(arrays,plan['native_arrays']),(spectra,plan['native_spectra'])]:
        if fingerprint(path)!=expected:raise ValueError('Changed native comparison grid/state')
        checks[str(path)]=expected
    checks[str(pp)]=fingerprint(pp);checks[str(ROOT/'pipeline/native_ghost_boundary.py')]=fingerprint(ROOT/'pipeline/native_ghost_boundary.py')
    with np.load(arrays)as d:f=d['consumed_state_cgs'][::-1].astype(float);state=dict(height_m=d['height_m'][::-1].copy(),temperature_K=f[:,2],ne_m3=f[:,0]*1e6,rho_kg_m3=f[:,1]*1e3,vturb_m_s=f[:,7]*1e3);native_star=d['native_nstar_cm3'][:,::-1]*1e6
    with np.load(pp)as d:pop=d['populations'].copy();star=d['nstar'].copy();np.testing.assert_array_equal(star,native_star)
    with np.load(spectra)as d:vac=d['wave_vac_A'].copy();air=d['wave_air_A'].copy();reference=d['Inu_cgs']/d['continuum_Inu_cgs']
    atom_file=require_codex_path('experiments.bride_oxygen_atom');ref_file=require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(ref_file.read_text());model_file=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name']
    if fingerprint(atom_file)!=ref['atom']or fingerprint(model_file)!=ref['atmosphere_file']:raise ValueError('Changed source atom/atmosphere')
    atom=read_multi_atom(atom_file);law,more=load_width_law(plan['original_width_law'],atom_file,atom);checks.update(more);structure=read_marcs_structure(model_file)
    a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter();deadline=start+600
    contexts=[];audits=[]
    for continuum in [False,True]:
        ctx,audit=make_context(state,structure,atom,pop,star,extra_vac_nm=np.array([vac[0]/10,vac[-1]/10]),formal_solver=plan['formal_solver'],helium_broadening=False,vdw_broadening=True,original_width_law=law,continuum_only=continuum,top_boundary=NativeLinearGhost(run['native_ghost_step_m']))
        rad=converge_cap_radiation(ctx,rtol=1e-8,maxiter=1000,deadline=deadline);contexts.append(ctx);audits.append(dict(context=audit,radiation=rad))
    values=[]
    for ctx in contexts:
        # Diagnostic rays keep the parent full-angle J for scattering; no one-ray
        # angular-equilibrium replacement. Recompute ghost coefficients on the
        # actual dense grid instead of interpolating the sharp line opacities.
        ctx.atmos.pyAtmos.zUpperBc.wave=None
        _,rayctx=ctx.compute_rays(wavelengths=vac/10,mus=np.array([1.]),returnCtx=True,upOnly=False)
        rayctx.atmos.pyAtmos.wmu[:]=1.;rayctx.depthData.fill=True
        fixedJ=np.array(rayctx.spect.J).copy();rayctx.formal_sol_gamma_matrices(lambdaIterate=True)
        rayctx.spect.J[:]=fixedJ;rayctx.atmos.pyAtmos.zUpperBc.update_from_context(rayctx,fixedJ[:,0])
        rayctx.formal_sol_gamma_matrices(lambdaIterate=True)
        out=np.asarray(rayctx.spect.I).reshape(len(vac),1);values.append(out[:,0].copy())
        np.testing.assert_array_equal(rayctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,pop)
        np.testing.assert_array_equal(ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,pop)
    I,Ic=values
    if any(not np.isfinite(v).all()or np.any(v<=0)for v in values):raise ValueError('Invalid coupled forward intensities')
    normalized=I/Ic;rows=[]
    for centre in [6158.17,6300.304,6363.776,7771.944,7774.166,7775.388,8446.4]:
        lo,hi=centre-.5,centre+.5;mask=abs(air-centre)<=.5
        n=1000*window_integral(air,reference,lo,hi,max_gap=.0051);own=1000*window_integral(air,normalized,lo,hi,max_gap=.0051)
        rows.append(dict(centre_air_A=centre,scope='absent-transition control'if centre==6363.776 else'finite O-only vertical intensity',native_area_mA=n,coupled_area_mA=own,difference_mA=own-n,max_normalized_difference=float(np.max(abs(normalized[mask]-reference[mask])))))
    np.savez_compressed(a.output/'spectra.npz',wave_vac_A=vac,wave_air_A=air,Inu_SI=I,continuum_Inu_SI=Ic,normalized=normalized,native_normalized=reference)
    with(a.output/'windows.csv').open('w',newline='')as stream:w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for residual in [False,True]:
        fig,axes=plt.subplots(3,2,figsize=(12,10),constrained_layout=True)
        for ax,centre in zip(axes.flat,[6158.17,6300.304,7771.944,7774.166,7775.388,8446.4]):
            mask=abs(air-centre)<=.5
            if residual:ax.plot(air[mask]-centre,(normalized-reference)[mask]);ax.axhline(0,color='0.6',lw=.7);ax.set_ylabel('Coupled minus native I/Ic')
            else:ax.plot(air[mask]-centre,reference[mask],label='Native reference');ax.plot(air[mask]-centre,normalized[mask],'--',label='Coupled source O41');ax.set_ylabel('I/Ic');ax.legend(fontsize=8)
            ax.set_title(f'{centre:.3f} AIR A');ax.set_xlabel('AIR offset / A')
        fig.suptitle('Full source atom, coupled populations; no fitted shifts/scales or optical overrides\nNative zero-ghost boundary; intensity at the sampled top node\nTheoretical disk-centre only; physical/observational admission HOLD')
        fig.savefig(a.output/('residuals.png'if residual else'comparison.png'),dpi=160);plt.close(fig)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed spectrum input')
    report=dict(generator=fingerprint(__file__),inputs=checks,scientific_admission='HOLD',seconds=time.perf_counter()-start,audits=audits,windows=rows,scope='Forward fixed completed coupled populations. No SE update during this postprocess; coupled solve provenance is separate. Finite±0.5A O-only theoretical intensity, no full EW/abundance/observational fit.6363 absent; source/native table rounding, constants and outer layers remain distinct.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
