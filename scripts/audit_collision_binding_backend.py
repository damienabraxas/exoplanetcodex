"""Verify actual compiled width and whole-O collision routing; fixed populations only."""
import argparse,json,os,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_hop_transfer import make_context
    from pipeline.multi_width_law import load_width_law
    from pipeline.native_vacuum_centres import load_centres
    from lightweaver.atomic_table import PeriodicTable
    import lightweaver.constants as const
    import lightweaver.LwCompiled as compiled
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--oracle-comparison',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());checks={str(a.plan.resolve()):fingerprint(a.plan)}
    if plan['purpose']!='fixed-native-oxygen-collision-backend-audit':raise ValueError('Require explicit collision backend-audit plan')
    if os.environ.get('OMP_NUM_THREADS')!='1'or fingerprint(Path(compiled.__file__))!=plan['backend']:raise ValueError('Require qualified one-thread backend')
    checks[str(Path(compiled.__file__))]=plan['backend'];native=require_codex_path(plan['native_analysis_resource'])
    for name,key in [('arrays.npz','native_arrays'),('spectra.npz','native_spectra'),('report.json','native_report')]:
        path=native/name
        if fingerprint(path)!=plan[key]:raise ValueError('Changed native source state')
        checks[str(path)]=plan[key]
    nr=json.loads((native/'report.json').read_text());_,_,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant=nr.get('profile_variant','original')))
    if nr['engine_receipt']!=engine:raise ValueError('Changed native reference')
    with np.load(native/'arrays.npz')as d:f=d['consumed_state_cgs'][::-1].astype(float);height=d['height_m'][::-1].copy();pop=d['native_populations_cm3'][:,::-1].copy()*1e6;star=d['native_nstar_cm3'][:,::-1].copy()*1e6
    with np.load(native/'spectra.npz')as d:vac=d['wave_vac_A'].copy()
    state=dict(height_m=height,temperature_K=f[:,2],ne_m3=f[:,0]*1e6,rho_kg_m3=f[:,1]*1e3,vturb_m_s=f[:,7]*1e3)
    atom_file=require_codex_path('experiments.bride_oxygen_atom');ref_file=require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(ref_file.read_text());model_file=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name']
    if fingerprint(atom_file)!=ref['atom']or fingerprint(model_file)!=ref['atmosphere_file']:raise ValueError('Changed original atom/atmosphere')
    atom=read_multi_atom(atom_file);law,more=load_width_law(plan['original_width_law'],atom_file,atom);checks.update(more);item=plan['table_control'];table=require_codex_path(item['resource'])/'native/atoms/O.txt'
    if fingerprint(table)!=item['receipt']:raise ValueError('Changed native table')
    checks[str(table)]=item['receipt'];centres,more=load_centres(item['native_vacuum_packet'],table,atom);checks.update(more)
    law_file=require_codex_path(plan['original_width_law']['resource'])/'law.json';law_data=json.loads(law_file.read_text());op=a.oracle_comparison/'widths.npz'
    if law_data['source_inputs'].get(str(op))!=fingerprint(op):raise ValueError('Wrong native oracle comparison for width law')
    checks[str(op)]=fingerprint(op)
    with np.load(op)as d:original=d['native_widths'].copy();T=d['temperature_K'].copy();HI=d['native_HI_m3'].copy()
    np.testing.assert_array_equal(T,state['temperature_K'])
    for path in [atom_file,ref_file,model_file,Path(const.__file__),ROOT/'pipeline/multi_width_law.py',ROOT/'pipeline/bride_hop_transfer.py']:checks[str(path)]=fingerprint(path)
    from pipeline.bride_collision_binding import load_collision_binding
    binding,more=load_collision_binding(plan['collision_binding'],atom_file,atom,state,star);checks.update(more);checks[str(ROOT/'pipeline/bride_collision_binding.py')]=fingerprint(ROOT/'pipeline/bride_collision_binding.py')
    ctx,audit=make_context(state,read_marcs_structure(model_file),atom,pop,star,extra_vac_nm=np.unique([vac[0]/10,vac[-1]/10]),formal_solver=plan['formal_solver'],helium_broadening=plan['helium_broadening'],vdw_broadening=plan['vdw_broadening'],table_control=(table,item['mode'],centres),original_width_law=law,collision_binding=binding)
    eq=ctx.kwargs['eqPops'];oxygen=eq.atomicPops[PeriodicTable['O']];model=oxygen.model;comp=ctx.detailedAtoms[0]
    if ctx.activeAtoms or len(model.lines)!=256:raise ValueError('Unexpected coupled/source model')
    actual_gamma=[];actual_damp=[];lam=[]
    for line,transition in zip(model.lines,comp.trans[:256]):
        if transition.type!='Line':raise ValueError('Wrong computational transition type')
        np.testing.assert_allclose(transition.lambda0,line.lambda0,rtol=1e-15);actual_gamma.append(transition.Qelast.copy());actual_damp.append(transition.aDamp.copy());lam.append(line.lambda0_m)
    g=np.asarray(actual_gamma);ad=np.asarray(actual_damp);speed=np.sqrt(2*const.KBoltzmann/(const.Amu*float(model.element.mass))*T+state['vturb_m_s']**2);expected=original[:,:,0]*np.array(lam)[:,None]/(4*np.pi*speed[None,:]);bound=128*np.finfo(float).eps
    np.testing.assert_allclose(g,original[:,:,2],rtol=bound,atol=0);np.testing.assert_allclose(ad,expected,rtol=bound,atol=0);np.testing.assert_array_equal(oxygen.n,pop);np.testing.assert_array_equal(oxygen.nStar,star)
    # Detailed-static atoms deliberately allocate no C/Gamma. Construct a detached
    # computational atom with real active-mode storage, without registering it in
    # the Context or permitting any HOP coupled solve. Constructor runs no SE.
    collision_atom=compiled.LwAtom(model,ctx.atmos,eq,ctx.kwargs['spect'],ctx.background,detailed=False)
    comp=collision_atom
    comp.compute_collisions(fillDiagonal=False);actual_C=comp.C.copy();np.testing.assert_array_equal(actual_C,binding.matrix)
    refusals=[];atmos=ctx.atmos.pyAtmos
    for name in ['temperature','ne','nHTot','vturb','height']:
        value=getattr(atmos,name);old=value.copy();value[0]*=1.01
        try:
            try:comp.compute_collisions(fillDiagonal=False)
            except ValueError:refusals.append(name)
            else:raise AssertionError('Stale state reached compiled collision matrix '+name)
        finally:value[:]=old
    old=oxygen.nStar.copy();oxygen.nStar[0,0]*=1.01
    try:
        try:comp.compute_collisions(fillDiagonal=False)
        except ValueError:refusals.append('O nStar')
        else:raise AssertionError('Changed O nStar accepted')
    finally:oxygen.nStar[:]=old
    comp.compute_collisions(fillDiagonal=False);np.testing.assert_array_equal(comp.C,actual_C);np.testing.assert_array_equal(oxygen.n,pop);np.testing.assert_array_equal(oxygen.nStar,star)
    if len(refusals)!=6 or binding.calls<2:raise ValueError('Missing actual compiled refusal/replay')
    comp.compute_collisions(fillDiagonal=True);diagonalized=comp.C.copy()
    scale=np.maximum(np.sum(abs(diagonalized),axis=0),np.finfo(float).tiny)
    column_residual=float(np.max(abs(diagonalized.sum(axis=0))/scale))
    if column_residual>256*np.finfo(float).eps:raise ValueError('Compiled collision conservation diagonal failed')
    from pipeline.bride_rate_audit import process_flow
    net,gross=process_flow(actual_C,star)
    lte_residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
    if lte_residual>256*np.finfo(float).eps:raise ValueError('Source contract LTE closure failed')
    a.output.mkdir(parents=True,exist_ok=False);np.savez_compressed(a.output/'backend.npz',compiled_Qelast_s_inv=g,compiled_adamp=ad,expected_adamp=expected,native_gamma_vdw_s_inv=original[:,:,2],compiled_collisions=actual_C,source_collisions=binding.matrix);report=dict(generator=fingerprint(__file__),collision_calls=binding.calls,compiled_conservation_column_residual=column_residual,source_contract_LTE_residual=lte_residual,collision_matrix_bit_identical=True,collision_refusals=refusals,collision_inputs=binding.receipt,inputs=checks,audit=audit,compiled_gamma_max_relative_error=float(np.max(abs(g/original[:,:,2]-1))),compiled_adamp_max_relative_error=float(np.max(abs(ad/expected-1))),all_256x56_VDW_widths_positive=bool(np.all(g>0)),source_mass_amu=atom.weight,LW_mass_amu=float(model.element.mass),scientific_admission='HOLD',scope='Actual compiled fixed-population width routing and detached active-storage collision atom; HOP Context stays detailed-static with no active atoms. LTE inverse closure is by contract, not independent reverse-rate evidence. No coupled SE or observed solar admission.')
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed actual-backend audit source')
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['compiled_gamma_max_relative_error','compiled_adamp_max_relative_error','all_256x56_VDW_widths_positive','collision_matrix_bit_identical','collision_calls','collision_refusals','compiled_conservation_column_residual','source_contract_LTE_residual']}))
if __name__=='__main__':main()
