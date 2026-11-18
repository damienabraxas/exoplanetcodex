"""Bounded full-source O41/HOP coupled diagnostic; foreign LTE conventions explicit."""
import argparse,json,os,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def source_atom_audit(atom,star,T):
    from pipeline.multi_atom import to_lightweaver
    from lightweaver.atomic_model import VoigtLine
    import lightweaver.constants as c
    m=to_lightweaver(atom);ratios=[];centres=[];continua=[]
    for raw,level in zip(atom.levels,m.levels):
        if (raw.E_cm,raw.g,raw.stage-1)!=(level.E,level.g,level.stage):raise ValueError('Source level identity changed')
    for raw,line in zip(atom.lines,m.lines):
        if type(line)is not VoigtLine or (line.i,line.j,line.f)!=(raw.i-1,raw.j-1,raw.f):raise ValueError('Transfer override in SE atom')
        gap=m.levels[line.j].E_SI-m.levels[line.i].E_SI
        centres.append(abs(line.lambda0_m/(c.HPlanck*c.CLight/gap)-1))
        ratio=m.levels[line.j].g/m.levels[line.i].g*np.exp(-gap/(c.KBoltzmann*T))
        ratios.append(float(np.max(abs((star[line.j]/star[line.i])/ratio-1))))
    for raw,bf in zip(atom.continua,m.continua):
        order=np.argsort(raw.wavelength);np.testing.assert_array_equal(bf.wavelengthGrid,raw.wavelength[order]*.1);np.testing.assert_array_equal(bf.alphaGrid,raw.alpha[order]*1e-4)
        kept=np.asarray(bf.wavelengthGrid)<=bf.lambdaEdge
        if np.count_nonzero(kept)<2:raise ValueError('Empty source BF support at physical edge')
        continua.append(dict(i=bf.i,j=bf.j,edge_vac_nm=bf.lambdaEdge,source_nodes=len(kept),nodes_not_redward=int(np.count_nonzero(kept))))
    if max(centres)>64*np.finfo(float).eps:raise ValueError('Source energy/photon inconsistency')
    return dict(levels=len(m.levels),lines=len(m.lines),continua=continua,source_energy_photon_max_relative=max(centres),foreign_nstar_vs_source_Boltzmann_max_relative=max(ratios),scope='Source level/f/wavelength/BF units and support checked; imported native LTE and historical-constant HOP background remain explicitly different thermodynamic conventions. Not a self-consistent LTE or observed-model admission.')

def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_collision_binding import load_collision_binding
    from pipeline.multi_width_law import load_width_law
    from pipeline.bride_hop_transfer import make_context
    from pipeline.native_ghost_boundary import converge_ghost_radiation as converge_cap_radiation
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.bride_parameter_checks import positive_stat_equil
    import lightweaver.LwCompiled as compiled
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--dry-run',action='store_true');a=p.parse_args();plan=json.loads(a.plan.read_text())
    if plan['purpose']!='full-source-O41-HOP-native-ghost-coupled-diagnostic'or 'table_control'in plan or os.environ.get('OMP_NUM_THREADS')!='1'or plan['seconds']!=600 or fingerprint(Path(compiled.__file__))!=plan['backend']:raise ValueError('Require explicit full-source one-thread600s diagnostic')
    control_file=ROOT/plan['controls_source'];control_report=json.loads(control_file.read_text())
    if fingerprint(control_file)!=plan['controls_receipt']:raise ValueError('Changed predeclared solver controls')
    from pipeline.bride_coupled_contract import validate_operator_contract
    validate_operator_contract(plan,control_report)
    controls=control_report['controls']
    if controls['seconds_per_case']!=plan['seconds']or controls['rtol']!=1e-6 or controls['rays']!=5 or controls['qscale']!=2 or not controls['ng']:raise ValueError('Unqualified original solver controls')
    native=require_codex_path(plan['native_analysis_resource']);checks={str(a.plan.resolve()):fingerprint(a.plan),str(control_file):fingerprint(control_file)}
    for filename,key in [('arrays.npz','native_arrays'),('report.json','native_report')]:
        q=native/filename
        if fingerprint(q)!=plan[key]:raise ValueError('Changed native state')
        checks[str(q)]=plan[key]
    nr=json.loads((native/'report.json').read_text());_,engine_app,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant=nr['profile_variant']))
    if nr['engine_receipt']!=engine:raise ValueError('Changed native engine')
    with np.load(native/'arrays.npz')as d:f=d['consumed_state_cgs'][::-1].astype(float);state=dict(height_m=d['height_m'][::-1].copy(),temperature_K=f[:,2],ne_m3=f[:,0]*1e6,rho_kg_m3=f[:,1]*1e3,vturb_m_s=f[:,7]*1e3);pop=d['native_populations_cm3'][:,::-1].copy()*1e6;star=d['native_nstar_cm3'][:,::-1].copy()*1e6
    atom_file=require_codex_path('experiments.bride_oxygen_atom');ref_file=require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(ref_file.read_text());model_file=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name']
    if fingerprint(atom_file)!=ref['atom']or fingerprint(model_file)!=ref['atmosphere_file']:raise ValueError('Changed source atom/atmosphere')
    native_seed=pop.copy()
    if 'restart'in plan:
        item=plan['restart'];parent=require_codex_path(item['resource']);rp=parent/'report.json';pp=parent/'populations.npz'
        if fingerprint(rp)!=item['report_receipt']or fingerprint(pp)!=item['populations_receipt']:raise ValueError('Changed restart parent')
        prior=json.loads(rp.read_text())
        if not prior['computational_completed']or prior['dry_run']or prior['fresh_physical_rate_residual']>=controls['rtol']or prior['iteration_plugin']!=control_report['iteration_plugin']:raise ValueError('Unqualified restart parent')
        with np.load(pp)as d:
            np.testing.assert_array_equal(d['nstar'],star);pop=d['populations'].copy()
        if pop.shape!=star.shape or not np.isfinite(pop).all()or np.any(pop<=0):raise ValueError('Invalid restart populations')
        np.testing.assert_allclose(pop.sum(axis=0),star.sum(axis=0),rtol=1e-10,atol=0)
        checks[str(rp)]=item['report_receipt'];checks[str(pp)]=item['populations_receipt']
    atom=read_multi_atom(atom_file);law,more=load_width_law(plan['original_width_law'],atom_file,atom);checks.update(more);binding,more=load_collision_binding(plan['collision_binding'],atom_file,atom,state,star);checks.update(more)
    for q in [atom_file,ref_file,model_file,Path(compiled.__file__),ROOT/'pipeline/bride_hop_transfer.py',ROOT/'pipeline/bride_collision_binding.py',ROOT/'pipeline/hop_packet.py',ROOT/'pipeline/bride_numerical_validation.py',ROOT/'pipeline/bride_top_cap.py',ROOT/'pipeline/bride_parameter_checks.py',ROOT/'pipeline/bride_coupled_contract.py',ROOT/'pipeline/multi_atom.py',ROOT/'pipeline/multi_width_law.py',ROOT/'pipeline/native_oxygen_collisions.py',ROOT/'pipeline/bride_reference.py',ROOT/'pipeline/hop_abi.py']:checks[str(q)]=fingerprint(q)
    import lightweaver.constants as lw_constants
    import lightweaver.atomic_model as lw_atomic_model
    import lightweaver.broadening as lw_broadening
    for module in [lw_constants,lw_atomic_model,lw_broadening]:checks[str(Path(module.__file__))]=fingerprint(module.__file__)
    if fingerprint(lw_constants.__file__)!=control_report['LW_constants']:raise ValueError('Inherited LW constants changed')
    preflight=source_atom_audit(atom,star,state['temperature_K']);a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter();deadline=start+plan['seconds']
    checks[str(Path(__file__).resolve())]=fingerprint(__file__)
    scheme=None;plugin=None
    if 'iteration_plugin'in plan:
        from pipeline.bride_overlap_plugin import load_overlap_plugin
        directory=require_codex_path(plan['iteration_plugin']['resource'])
        if fingerprint(directory/'receipt.json')!=plan['iteration_plugin']['build_receipt']:raise ValueError('Changed predeclared iteration plugin')
        scheme,plugin=load_overlap_plugin(directory)
        if scheme!='bride_full_coupling_cap_scalar':raise ValueError('Unqualified numerical operator')
        for q in [directory/'receipt.json',directory/'bride_overlap_guard.dylib',ROOT/'pipeline/bride_overlap_plugin.py']:checks[str(q)]=fingerprint(q)
    report=dict(generator=fingerprint(__file__),inputs=checks,controls=controls,scientific_admission='HOLD',computational_completed=False,dry_run=a.dry_run,iteration_plugin=plugin,source_atom=preflight,scope=plan['scope'])
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    native_run=require_codex_path(plan['native_run_resource']);mesh=native_run/'native/atmos_mesh.txt';native_run_report=native_run/'report.json'
    if fingerprint(native_run_report)!=nr['native_report']:raise ValueError('Wrong native boundary producer')
    rows=mesh.read_text().splitlines();steps=np.concatenate([np.full(int(rows[i].split()[2]),float(rows[i+1].split()[2]))for i in range(0,len(rows),2)])
    np.testing.assert_array_equal(np.cumsum(steps)[::-1]*1e6,state['height_m'])
    from pipeline.native_ghost_boundary import NativeLinearGhost
    boundary=NativeLinearGhost(float(steps[-1]*1e6))
    for path in [native_run_report,mesh,ROOT/'pipeline/native_ghost_boundary.py']:checks[str(path)]=fingerprint(path)
    native_doc=json.loads(native_run_report.read_text())
    for entry in native_doc['output_files']:
        path=native_run/entry['path']
        if fingerprint(path)!=entry['receipt']:raise ValueError('Changed native boundary input')
        checks[str(path)]=entry['receipt']
    import re
    params=(native_run/'native/params.nml').read_text();group=params[params.index('&RT_ITERATION_PARAMS'):].split('/')[0]
    if re.search(r"SOLVER='linear +?'",group)is None:raise ValueError('Ghost interval requires actual native linear solver')
    for path in [engine_app/'SRC/rt_solvers_mod.f90',engine_app/'SRC/rt_patch_mod.f90',engine_app/'SRC/atom_mod.f90']:checks[str(path)]=fingerprint(path)
    report['native_ghost_step_m']=boundary.step_m
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    ctx=None
    try:
        ctx,audit=make_context(state,read_marcs_structure(model_file),atom,pop,star,extra_vac_nm=np.array([500.,845.]),formal_solver=plan['formal_solver'],helium_broadening=False,vdw_broadening=True,original_width_law=law,collision_binding=binding,coupled_diagnostic=True,top_boundary=boundary)
        if len(ctx.activeAtoms)!=1 or ctx.detailedAtoms:raise ValueError('Unexpected active model set')
        if scheme is not None:ctx.set_fs_iter_scheme(scheme);ctx.kwargs['fsIterScheme']=scheme
        audit['scope']='Source-defined native linear zero-ghost boundary coupled diagnostic; full source atom. Incoming boundary emissivity feedback is iterated externally; MALI boundary derivative left zero, physical action/fresh rate identity required. No physical admission.'
        active=ctx.activeAtoms[0];active.compute_collisions(fillDiagonal=False);np.testing.assert_array_equal(active.C,binding.matrix);report['context_audit']=audit;report['wavelengths']=len(ctx.spect.wavelength)
        if not a.dry_run:
            with(a.output/'history.jsonl').open('x')as stream:
                for iteration in range(controls['maxiter']):
                    if time.perf_counter()>deadline:raise TimeoutError('Original600s coupled budget exceeded')
                    old=active.n.copy();oldJ=np.array(ctx.spect.J).copy();rad=converge_cap_radiation(ctx,rtol=controls['cap_radiation_rtol'],maxiter=controls['cap_radiation_maxiter'],deadline=deadline)
                    refresh=ctx.formal_sol_gamma_matrices(lambdaIterate=False,extraParams={'BrideCapLengthM':float(ctx.atmos.pyAtmos.zUpperBc.length_m)})
                    physical=physical_rate_matrix(active);gross_old=np.einsum('ijd,jd->id',abs(physical),old);difference=np.einsum('ijd,jd->id',active.Gamma-physical,old);action_error=float(np.max(abs(difference)/np.maximum(gross_old,np.finfo(float).tiny)))
                    if scheme is not None and action_error>256*np.finfo(float).eps:raise ValueError('Declared rate-preserving operator failed physical action identity')
                    scale=np.max(abs(active.Gamma),axis=1,keepdims=True)
                    if np.any(scale==0):raise ValueError('Disconnected SE equation')
                    active.Gamma[:]/=scale;rejected=positive_stat_equil(ctx,active,old,reject_invalid_ng=controls['reject_invalid_ng'])
                    rate=physical_rate_matrix(active);net=np.einsum('ijd,jd->id',rate,active.n);gross=np.einsum('ijd,jd->id',abs(rate),active.n)
                    row=dict(iteration=iteration+1,operator_action_discrepancy=action_error,population_change=float(np.max(abs(active.n-old)/active.n)),rate_residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny))),J_change=float(np.max(abs(np.array(ctx.spect.J)-oldJ)/np.maximum(abs(np.array(ctx.spect.J)),np.finfo(float).tiny))),conservation=float(np.max(abs(active.n.sum(axis=0)/active.nTotal-1))),rejected_invalid_ng=bool(rejected),fixed_radiation=rad,refresh_dJ=float(refresh.dJMax),seconds=time.perf_counter()-start)
                    stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
                    if row['conservation']>1e-10:raise ValueError('Conservation failed')
                    if iteration>=2 and max(row['population_change'],row['rate_residual'],row['J_change'])<controls['rtol']:break
                else:raise RuntimeError('Original maxiter exhausted')
            final=active.n.copy();fresh=converge_cap_radiation(ctx,rtol=controls['cap_radiation_rtol'],maxiter=controls['cap_radiation_maxiter'],deadline=deadline);active.compute_collisions(fillDiagonal=False);rate=physical_rate_matrix(active);net=np.einsum('ijd,jd->id',rate,final);gross=np.einsum('ijd,jd->id',abs(rate),final);res=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
            np.testing.assert_array_equal(active.n,final);np.testing.assert_array_equal(active.nStar,star)
            report.update(final_iteration=row,fresh_radiation=fresh,fresh_physical_rate_residual=res,max_relative_population_response=float(np.max(abs(final/pop-1))),max_relative_native_population_response=float(np.max(abs(final/native_seed-1))))
            if res>=controls['rtol']:raise ValueError('Fresh physical rate residual failed')
            np.savez_compressed(a.output/'rates.npz',physical_R_plus_C=rate,collision_offdiag=active.C.copy(),radiative_up_s_inv=np.array([t.Rij for t in active.trans]),radiative_down_s_inv=np.array([t.Rji for t in active.trans]),transition_i=np.array([t.i for t in active.trans]),transition_j=np.array([t.j for t in active.trans]),net=net,gross=gross,wave_vac_nm=np.array(ctx.spect.wavelength),Jnu=np.array(ctx.spect.J))
        report['computational_completed']=True
    except (ValueError,RuntimeError,TimeoutError)as error:report['error']=str(error)
    finally:
        report['seconds']=time.perf_counter()-start;report['collision_calls']=binding.calls
        if ctx is not None and ctx.activeAtoms:np.savez_compressed(a.output/'populations.npz',populations=ctx.activeAtoms[0].n.copy(),nstar=ctx.activeAtoms[0].nStar.copy(),initial_populations=pop)
        for path,fp in checks.items():
            if fingerprint(path)!=fp:report['computational_completed']=False;report['input_changed']=path
        (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['computational_completed','dry_run','seconds','collision_calls','error','fresh_physical_rate_residual']if k in report}))
    if not report['computational_completed']:raise RuntimeError('Coupled diagnostic refused/failed; evidence retained')

if __name__=='__main__':main()
