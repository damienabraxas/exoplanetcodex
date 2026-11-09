"""Labeled MALI population-preconditioning control; physical admission unchanged."""
import json,time
from pathlib import Path
import numpy as np


def solve_coupled_cap_mali(ctx,atom,original,output,start):
    from config.constants import require_codex_path
    from pipeline.bride_native_spline import NativeTautsp,install_native_collision_control
    from pipeline.bride_charge_contract import install_author_mn_pair
    from pipeline.bride_top_cap import converge_cap_radiation,cap_intensity
    from pipeline.bride_parameter_checks import positive_stat_equil
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.bride_reference import fingerprint
    controls=original['controls']
    if original['experimental_collision_convention']!='barklem2018_MN_LTE_pair':raise ValueError('Unsupported source collision contract')
    oracle=NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
    spline=install_native_collision_control(ctx,atom,oracle,source_directions=True,reversible_ch0=False)
    mn=install_author_mn_pair(ctx,atom,require_codex_path('experiments.bride_barklem_oxygen'),oracle)
    if spline!=original['native_spline_control']or mn!=original['author_charge_control']:raise ValueError('Collision source receipt changed')
    active=ctx.activeAtoms[0];initial=active.n.copy();deadline=start+controls['seconds_per_case']
    with (output/'history.jsonl').open('x') as stream:
        for iteration in range(controls['maxiter']):
            if time.perf_counter()>deadline:raise TimeoutError('Original coupled budget exceeded')
            old=active.n.copy()
            previous_outer_J=np.array(ctx.spect.J).copy();cap=ctx.atmos.pyAtmos.zUpperBc;mu=abs(np.asarray(ctx.atmos.pyAtmos.muz));previous_outer_boundary=np.zeros((len(ctx.spect.wavelength),len(mu)))if cap.wave is None else cap_intensity(cap.chi,cap.eta,cap.length_m,mu)
            radiation=converge_cap_radiation(ctx,rtol=controls['cap_radiation_rtol'],maxiter=controls['cap_radiation_maxiter'],deadline=deadline)
            previous_top=np.array(ctx.spect.J)[:,0].copy();refresh=ctx.formal_sol_gamma_matrices(lambdaIterate=False,extraParams={"BrideCapLengthM":float(cap.length_m)});boundary_refresh=cap.update_from_context(ctx,previous_top)
            if not np.isfinite(refresh.dJMax)or max(refresh.dJMax,boundary_refresh)>controls['cap_radiation_rtol']:raise ValueError('MALI refresh disturbed fixed-population radiation')
            relative=lambda new,old:float(np.max(abs(new-old)/np.maximum(abs(new),np.finfo(float).tiny)))
            radiation['outer_J_change']=relative(np.array(ctx.spect.J),previous_outer_J);radiation['outer_boundary_change']=relative(cap_intensity(cap.chi,cap.eta,cap.length_m,mu),previous_outer_boundary)
            radiation['mali_refresh_J_change']=float(refresh.dJMax);radiation['mali_refresh_boundary_change']=boundary_refresh
            matrix=active.Gamma;scale=np.max(abs(matrix),axis=1,keepdims=True)
            if np.any(scale==0):raise ValueError('Disconnected SE equation')
            matrix[:]/=scale
            rejected=positive_stat_equil(ctx,active,old,reject_invalid_ng=controls['reject_invalid_ng'])
            n=active.n.copy()
            if not np.isfinite(n).all()or np.any(n<=0):raise ValueError('Invalid coupled populations')
            rates=physical_rate_matrix(active);net=np.einsum('ijd,jd->id',rates,n);gross=np.einsum('ijd,jd->id',abs(rates),n)
            rate=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
            change=float(np.max(abs(n-old)/n));conservation=float(np.max(abs(n.sum(axis=0)/active.nTotal-1)))
            row=dict(iteration=iteration+1,population_change=change,physical_rate_residual=rate,
                J_change=radiation['outer_J_change'],cap_boundary_change=radiation['outer_boundary_change'],
                particle_conservation_error=conservation,rejected_invalid_ng=bool(rejected),fixed_population_radiation=radiation)
            stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
            if conservation>1e-10:raise ValueError('Particle conservation failed')
            if iteration>=2 and max(change,rate,row['J_change'],row['cap_boundary_change'])<controls['rtol']:break
        else:raise RuntimeError('Joint coupled convergence not reached')
    final=active.n.copy()
    fresh=converge_cap_radiation(ctx,rtol=1e-8,maxiter=1000,deadline=deadline)
    rates=physical_rate_matrix(active);net=np.einsum('ijd,jd->id',rates,final);gross=np.einsum('ijd,jd->id',abs(rates),final)
    residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
    np.testing.assert_array_equal(final,active.n)
    if residual>=controls['rtol']:raise ValueError('Fresh coupled rate residual exceeds original outer tolerance')
    collisions=active.C.copy()
    for i in range(len(collisions)):collisions[i,i]=0.
    radiative=np.zeros_like(collisions)
    for transition in active.trans:
        radiative[transition.j,transition.i]+=transition.Rij
        radiative[transition.i,transition.j]+=transition.Rji
    np.savez_compressed(output/'rates.npz',total_rates_s_inv=rates,collision_rates_offdiag_s_inv=collisions,
        radiative_rates_offdiag_s_inv=radiative,net_flow_m3_s=net,gross_flow_m3_s=gross,
        wavelength_vac_nm=np.array(ctx.spect.wavelength),Jnu_W_m2_sr_Hz=np.array(ctx.spect.J))
    np.savez_compressed(output/'populations.npz',populations=final,nstar=active.nStar,initial_populations=initial)
    return final,dict(controls=controls,iterations=iteration+1,final_diagnostics=row,fresh_physical_rate_residual=residual,
        fresh_radiation=fresh,max_relative_population_change_from_initial=float(np.max(abs(final/initial-1))),
        source_spline=spline,source_mn=mn,history=fingerprint(output/'history.jsonl'),populations=fingerprint(output/'populations.npz'),
        rates=fingerprint(output/'rates.npz'),helper=fingerprint(__file__),operator='Complete weighted MALI with incoming finite-cap emissivity derivative; unchanged physical rates',legacy_helper=fingerprint(Path(__file__).with_name('bride_coupled_control.py')))
