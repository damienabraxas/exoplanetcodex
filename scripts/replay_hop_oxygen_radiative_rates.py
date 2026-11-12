"""Matched-population full-source O radiative-rate replay; no SE update."""
import argparse,json,os,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

from scripts.run_hop_oxygen_coupled_diagnostic import source_atom_audit

def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_collision_binding import load_collision_binding
    from pipeline.multi_width_law import load_width_law
    from pipeline.bride_hop_transfer import make_context
    from pipeline.bride_top_cap import converge_cap_radiation
    import lightweaver.LwCompiled as compiled
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text())
    if plan['purpose']!='full-source-O41-HOP-fixed-population-rates'or 'table_control'in plan or os.environ.get('OMP_NUM_THREADS')!='1'or plan['seconds']!=600 or fingerprint(Path(compiled.__file__))!=plan['backend']:raise ValueError('Require explicit full-source one-thread600s diagnostic')
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
    if 'restart'in plan:raise ValueError('Matched-state replay forbids restart seed')
    atom=read_multi_atom(atom_file);law,more=load_width_law(plan['original_width_law'],atom_file,atom);checks.update(more);binding,more=load_collision_binding(plan['collision_binding'],atom_file,atom,state,star);checks.update(more)
    for q in [atom_file,ref_file,model_file,Path(compiled.__file__),ROOT/'pipeline/bride_hop_transfer.py',ROOT/'pipeline/bride_collision_binding.py',ROOT/'pipeline/hop_packet.py',ROOT/'pipeline/bride_numerical_validation.py',ROOT/'pipeline/bride_top_cap.py',ROOT/'pipeline/bride_parameter_checks.py',ROOT/'pipeline/bride_coupled_contract.py',ROOT/'pipeline/multi_atom.py',ROOT/'pipeline/multi_width_law.py',ROOT/'pipeline/native_oxygen_collisions.py',ROOT/'pipeline/bride_reference.py',ROOT/'pipeline/hop_abi.py',ROOT/'scripts/run_hop_oxygen_coupled_diagnostic.py',ROOT/'pipeline/openmulti3d_output.py']:checks[str(q)]=fingerprint(q)
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
    from pipeline.openmulti3d_output import read_rate_pair
    native_run=require_codex_path(plan['native_run_resource']);native_report=native_run/'report.json'
    if fingerprint(native_report)!=nr['native_report']:raise ValueError('Wrong native producer for matched-state replay')
    native_doc=json.loads(native_report.read_text())
    for entry in native_doc['output_files']:
        path=native_run/entry['path']
        if fingerprint(path)!=entry['receipt']:raise ValueError('Changed native output')
        checks[str(path)]=entry['receipt']
    checks[str(native_report)]=fingerprint(native_report)
    boundary_mode=plan.get('boundary_control','top-node-zero')
    if boundary_mode not in ['top-node-zero','native-zero-ghost-linear']:raise ValueError('Unknown explicit top-boundary control')
    top_boundary=None
    if boundary_mode=='native-zero-ghost-linear':
        from pipeline.native_ghost_boundary import NativeLinearGhost
        mesh=native_run/'native/atmos_mesh.txt';lines=mesh.read_text().splitlines();steps=np.concatenate([np.full(int(lines[i].split()[2]),float(lines[i+1].split()[2]))for i in range(0,len(lines),2)])
        np.testing.assert_array_equal(np.cumsum(steps)[::-1]*1e6,state['height_m'])
        params=(native_run/'native/params.nml').read_text();group=params[params.index('&RT_ITERATION_PARAMS'):].split('/')[0]
        if __import__('re').search(r"SOLVER='linear +?'",group)is None:raise ValueError('Native ghost control requires the linear solver')
        top_boundary=NativeLinearGhost(float(steps[-1]*1e6))
        for path in [mesh,native_run/'native/params.nml',engine_app/'SRC/rt_solvers_mod.f90',engine_app/'SRC/rt_patch_mod.f90',engine_app/'SRC/atom_mod.f90',ROOT/'pipeline/native_ghost_boundary.py']:checks[str(path)]=fingerprint(path)
    parts=[]
    for path in sorted((native_run/'native/save').glob('rates_O_*.bin')):parts.append(read_rate_pair(path,path.with_name(path.name.replace('rates_','colrates_'))))
    if len(parts)!=8:raise ValueError('Missing native rate slabs')
    native_n=np.concatenate([x['populations_cm3']for x in parts],axis=1)[:,::-1]*1e6;native_s=np.concatenate([x['nstar_cm3']for x in parts],axis=1)[:,::-1]*1e6
    np.testing.assert_array_equal(native_n,pop);np.testing.assert_array_equal(native_s,star)
    foreign_R=np.concatenate([x['radiative_rates_s_inv']for x in parts],axis=2)[:,:,::-1];foreign_C=np.concatenate([x['collision_rates_s_inv']for x in parts],axis=2)[:,:,::-1]
    physical_control=plan.get('native_physical_control')
    if physical_control is not None:
        directory=require_codex_path(physical_control['resource'])
        for name in ['report.json','rates.npz']:
            path=directory/name
            if fingerprint(path)!=physical_control[name]:raise ValueError('Changed native physical-rate control')
            checks[str(path)]=physical_control[name]
        physical_report=json.loads((directory/'report.json').read_text())
        if not physical_report['computational_completed']or physical_report['original_binary_files_bit_identical']!=145:raise ValueError('Unqualified native observer')
        if physical_report['inputs'].get(str(native_report))!=fingerprint(native_report):raise ValueError('Observer has wrong native baseline')
        with np.load(directory/'rates.npz')as d:
            np.testing.assert_array_equal(d['post'][:,::-1]*1e6,pop)
            np.testing.assert_array_equal(d['star'][:,::-1]*1e6,star)
            pop=d['n'][:,::-1].copy()*1e6
            foreign_R=d['R'][:,:,::-1].copy();foreign_C=d['C'][:,:,::-1].copy()
        native_seed=pop.copy()
    report=dict(generator=fingerprint(__file__),inputs=checks,scientific_admission='HOLD',computational_completed=False,scope=plan['scope'],source_atom=preflight,boundary_control=boundary_mode,native_rate_semantics='Saved native coefficients contain SE preconditioning; individual up/down values are not gross Einstein rates. Net-flow interpretation additionally requires the native linearization population. Fresh unpreconditioned native-field replay remains owed.')
    if physical_control is not None:report['native_rate_semantics']='Unpreconditioned physical rate integrals at exact PRE-update radiation populations; source FP32 contributions summed in FP64, normalized with source profile/LTE factors. Not a post-update fresh native formal solve.'
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');ctx=None
    try:
        ctx,audit=make_context(state,read_marcs_structure(model_file),atom,pop,star,extra_vac_nm=np.array([500.,845.]),formal_solver=plan['formal_solver'],helium_broadening=False,vdw_broadening=True,original_width_law=law,collision_binding=binding,coupled_diagnostic=True,top_boundary=top_boundary)
        ctx.set_fs_iter_scheme(scheme);ctx.kwargs['fsIterScheme']=scheme;active=ctx.activeAtoms[0];active.compute_collisions(fillDiagonal=False);np.testing.assert_array_equal(active.C,binding.matrix)
        if top_boundary is None:radiation=converge_cap_radiation(ctx,rtol=controls['cap_radiation_rtol'],maxiter=controls['cap_radiation_maxiter'],deadline=deadline)
        else:
            from pipeline.native_ghost_boundary import converge_ghost_radiation
            radiation=converge_ghost_radiation(ctx,rtol=controls['cap_radiation_rtol'],maxiter=controls['cap_radiation_maxiter'],deadline=deadline)
        np.testing.assert_array_equal(active.n,pop);np.testing.assert_array_equal(active.nStar,star)
        own_R=np.zeros_like(foreign_R)
        for transition in active.trans:own_R[transition.j,transition.i]+=transition.Rij;own_R[transition.i,transition.j]+=transition.Rji
        native_pairs=np.any(foreign_R>0,axis=2);own_pairs=np.any(own_R>0,axis=2)
        if not np.array_equal(native_pairs,own_pairs):raise ValueError('Source/native radiative graph differs')
        def residual(R,C):
            off=R+C;incoming=np.sum(off*pop[None,:,:],axis=1);outgoing=np.sum(off,axis=0)*pop;gross=incoming+outgoing;net=incoming-outgoing
            return float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny))),net,gross
        native_res,native_net,native_gross=residual(foreign_R,foreign_C);own_res,own_net,own_gross=residual(own_R,active.C)
        focus=plan['focus_level_1based']-1;depth=plan['focus_depth_index'];changed=own_R-foreign_R
        contributions=changed[focus,:,depth]*pop[:,depth]-changed[:,focus,depth]*pop[focus,depth]
        ranking=[]
        for partner in np.argsort(abs(contributions))[::-1]:
            if partner==focus or contributions[partner]==0:continue
            lo=min(focus,int(partner));hi=max(focus,int(partner));gap=atom.level(hi+1).E_cm-atom.level(lo+1).E_cm
            ranking.append(dict(partner_1based=int(partner+1),kind='BB'if atom.level(lo+1).stage==atom.level(hi+1).stage else'BF',energy_gap_vac_A=1e8/gap,net_flow_difference_m3_s=float(contributions[partner]),native_effective_up_s_inv=float(foreign_R[hi,lo,depth]),replay_physical_up_s_inv=float(own_R[hi,lo,depth]),native_effective_down_s_inv=float(foreign_R[lo,hi,depth]),replay_physical_down_s_inv=float(own_R[lo,hi,depth])))
        from pipeline.rate_semantics import gross_rate_relative_difference
        checks[str(ROOT/'pipeline/rate_semantics.py')]=fingerprint(ROOT/'pipeline/rate_semantics.py')
        if physical_control is not None:
            comparison=gross_rate_relative_difference(foreign_R,own_R,reference_semantics='physical-gross',replay_semantics='physical-gross')
            comparison_refusal=None
            for row in ranking:
                for direction in ['up','down']:row['native_physical_'+direction+'_s_inv']=row.pop('native_effective_'+direction+'_s_inv')
        else:
            comparison=None
            try:gross_rate_relative_difference(foreign_R,own_R,reference_semantics='native-preconditioned-effective',replay_semantics='physical-gross')
            except ValueError as error:comparison_refusal=str(error)
            else:raise ValueError('Mixed rate semantics were incorrectly accepted')
        report.update(computational_completed=True,radiation=radiation,context_audit=audit,native_saved_effective_rate_residual=native_res,matched_population_replay_physical_residual=own_res,positive_directional_pairs=int(native_pairs.sum()),gross_coefficient_comparison_refused=comparison_refusal,focus=dict(level_1based=focus+1,depth=depth,temperature_K=float(state['temperature_K'][depth]),native_net_m3_s=float(native_net[focus,depth]),replay_net_m3_s=float(own_net[focus,depth]),ranked_partner_flow_changes=ranking))
        if physical_control is not None:report['native_physical_rate_residual']=report.pop('native_saved_effective_rate_residual');report['gross_rate_relative_difference']=comparison
        np.savez_compressed(a.output/'rates.npz',native_R_s_inv=foreign_R,replay_R_s_inv=own_R,native_C_s_inv=foreign_C,replay_C_s_inv=active.C.copy(),populations=pop,nstar=star,wave_vac_nm=np.array(ctx.spect.wavelength),Jnu=np.array(ctx.spect.J))
    except (ValueError,RuntimeError,TimeoutError)as error:report['error']=str(error)
    finally:
        report['seconds']=time.perf_counter()-start
        for path,fp in checks.items():
            if fingerprint(path)!=fp:report['computational_completed']=False;report['input_changed']=path
        (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['computational_completed','seconds','native_saved_effective_rate_residual','matched_population_replay_physical_residual','error']if k in report}))
    if not report['computational_completed']:raise RuntimeError('Matched-state rate replay failed; evidence retained')

if __name__=='__main__':main()
