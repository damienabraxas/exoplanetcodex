"""Bounded inverse-secant solve of the corrected physical MALI fixed-point map."""
import json,time
from pathlib import Path
import numpy as np


class InvalidMapPopulation(ValueError):
    pass


def solve_coupled_secant(ctx,atom,original,output,start):
    from config.constants import require_codex_path
    from pipeline.bride_native_spline import NativeTautsp,install_native_collision_control
    from pipeline.bride_charge_contract import install_author_mn_pair
    from pipeline.bride_parameter_checks import positive_stat_equil
    from pipeline.bride_top_cap import converge_cap_radiation,cap_intensity
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_secant import InverseSecant,log_coordinates,normalized_log_populations,populations_from_coordinates
    controls=original['controls'];parameters=original['numerical_algorithm']['parameters']
    if controls['ng']or not 0<parameters['radiation_rtol']<=controls['cap_radiation_rtol']:raise ValueError('Require explicit Ng off and no weaker radiation tolerance')
    if original['experimental_collision_convention']!='barklem2018_MN_LTE_pair':raise ValueError('Unsupported collision contract')
    oracle=NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
    spline=install_native_collision_control(ctx,atom,oracle,source_directions=True,reversible_ch0=False)
    mn=install_author_mn_pair(ctx,atom,require_codex_path('experiments.bride_barklem_oxygen'),oracle)
    if spline!=original['native_spline_control']or mn!=original['author_charge_control']:raise ValueError('Collision source receipt changed')
    active=ctx.activeAtoms[0];initial=active.n.copy();total=active.nTotal.copy();deadline=start+controls['seconds_per_case'];cap=ctx.atmos.pyAtmos.zUpperBc;mu=abs(np.asarray(ctx.atmos.pyAtmos.muz));evaluations=0
    inverse=InverseSecant(memory=parameters['memory'],denominator_rtol=parameters['denominator_rtol'])

    def evaluate(coordinates):
        nonlocal evaluations
        if time.perf_counter()>deadline:raise TimeoutError('Original coupled budget exceeded')
        n=populations_from_coordinates(coordinates,total);active.n[:]=n
        radiation=converge_cap_radiation(ctx,rtol=parameters['radiation_rtol'],maxiter=controls['cap_radiation_maxiter'],deadline=deadline)
        top=np.array(ctx.spect.J)[:,0].copy();refresh=ctx.formal_sol_gamma_matrices(lambdaIterate=False);boundary=cap.update_from_context(ctx,top)
        if not np.isfinite(refresh.dJMax)or max(refresh.dJMax,boundary)>parameters['radiation_rtol']:raise ValueError('Secant MALI refresh disturbed radiation')
        rates=physical_rate_matrix(active);net=np.einsum('ijd,jd->id',rates,n);gross=np.einsum('ijd,jd->id',abs(rates),n);physical=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
        matrix=active.Gamma;scale=np.max(abs(matrix),axis=1,keepdims=True)
        if np.any(scale==0):raise ValueError('Disconnected SE equation')
        matrix[:]/=scale
        try:positive_stat_equil(ctx,active,n,reject_invalid_ng=controls['reject_invalid_ng'])
        except ValueError as error:
            active.n[:]=n
            if str(error)!='Invalid populations':raise
            raise InvalidMapPopulation('Trial lies outside the positive MALI-map domain')from error
        mapped=active.n.copy();active.n[:]=n
        residual=coordinates-log_coordinates(mapped)
        if not np.isfinite(residual).all()or not np.isfinite(physical):raise ValueError('Nonfinite coupled residual')
        evaluations+=1
        return dict(populations=n,mapped=mapped,residual=residual,merit=float(np.linalg.norm(residual)),map_change=float(np.max(abs(mapped/n-1))),physical=physical,radiation=radiation)

    coordinates=log_coordinates(initial);current=evaluate(coordinates)
    with(output/'history.jsonl').open('x')as stream:
        for iteration in range(controls['maxiter']):
            old=active.n.copy();old_J=np.array(ctx.spect.J).copy();old_boundary=cap_intensity(cap.chi,cap.eta,cap.length_m,mu);old_coordinates=coordinates.copy();old_residual=current['residual'].copy();rejections=0
            if current['map_change']<controls['rtol']:
                coordinates=log_coordinates(current['mapped']);chosen=evaluate(coordinates);direction_kind='direct_verification';alpha=1.;update=dict(updated=False,restarted=False,rank=len(inverse.u))
            else:
                directions=[('secant',-inverse.apply(old_residual.ravel()).reshape(coordinates.shape)),('direct_fallback',-old_residual)]
                chosen=None;logs=normalized_log_populations(old_coordinates)
                for direction_kind,direction in directions:
                    if direction_kind=='direct_fallback':inverse.u=[];inverse.v=[]
                    alpha=1.
                    for _ in range(parameters['line_search_halvings']):
                        trial=old_coordinates+alpha*direction
                        if not np.isfinite(trial).all()or np.max(abs(normalized_log_populations(trial)-logs))>parameters['max_log_step']:
                            alpha*=.5;continue
                        try:candidate=evaluate(trial)
                        except InvalidMapPopulation:
                            rejections+=1;alpha*=.5;continue
                        if candidate['merit']<(1-parameters['armijo']*alpha)*current['merit']:
                            coordinates=trial;chosen=candidate;break
                        rejections+=1;alpha*=.5
                    if chosen is not None:break
                if chosen is None:raise RuntimeError('Secant and direct-map line searches failed; no tolerance relaxation')
                update=inverse.update((coordinates-old_coordinates).ravel(),(chosen['residual']-old_residual).ravel())
            current=chosen;n=active.n.copy();np.testing.assert_array_equal(n,current['populations'])
            relative=lambda new,old:float(np.max(abs(new-old)/np.maximum(abs(new),np.finfo(float).tiny)))
            conservation=float(np.max(abs(n.sum(axis=0)/total-1)))
            row=dict(iteration=iteration+1,population_change=relative(n,old),physical_rate_residual=current['physical'],J_change=relative(np.array(ctx.spect.J),old_J),cap_boundary_change=relative(cap_intensity(cap.chi,cap.eta,cap.length_m,mu),old_boundary),particle_conservation_error=conservation,map_change=current['map_change'],map_merit=current['merit'],direction=direction_kind,alpha=alpha,secant=update,line_search_rejections=rejections,map_evaluations=evaluations,fixed_population_radiation=current['radiation'])
            stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
            if conservation>1e-10:raise ValueError('Particle conservation failed')
            if iteration>=2 and max(row[k]for k in ['population_change','physical_rate_residual','J_change','cap_boundary_change','map_change'])<controls['rtol']:break
        else:raise RuntimeError('Original coupled iteration limit exceeded')
    final=active.n.copy();fresh=converge_cap_radiation(ctx,rtol=controls['cap_radiation_rtol'],maxiter=1000,deadline=deadline)
    rates=physical_rate_matrix(active);net=np.einsum('ijd,jd->id',rates,final);gross=np.einsum('ijd,jd->id',abs(rates),final);residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
    np.testing.assert_array_equal(final,active.n)
    if residual>=controls['rtol']:raise ValueError('Fresh physical admission failed')
    collision=active.C.copy();i=np.arange(len(collision));collision[i,i]=0.;radiative=np.zeros_like(collision)
    for t in active.trans:radiative[t.j,t.i]+=t.Rij;radiative[t.i,t.j]+=t.Rji
    np.savez_compressed(output/'rates.npz',total_rates_s_inv=rates,collision_rates_offdiag_s_inv=collision,radiative_rates_offdiag_s_inv=radiative,net_flow_m3_s=net,gross_flow_m3_s=gross,wavelength_vac_nm=np.array(ctx.spect.wavelength),Jnu_W_m2_sr_Hz=np.array(ctx.spect.J))
    np.savez_compressed(output/'populations.npz',populations=final,nstar=active.nStar,initial_populations=initial)
    return final,dict(controls=controls,iterations=iteration+1,final_diagnostics=row,fresh_physical_rate_residual=residual,fresh_radiation=fresh,map_evaluations=evaluations,source_spline=spline,source_mn=mn,history=fingerprint(output/'history.jsonl'),populations=fingerprint(output/'populations.npz'),rates=fingerprint(output/'rates.npz'),helper=fingerprint(__file__),secant_helper=fingerprint(__import__('pipeline.bride_secant',fromlist=['']).__file__),numerical_algorithm=original['numerical_algorithm'])
