"""Mac-only serial oxygen parameter experiment authorized by Ryan; no fitting."""
from __future__ import annotations
import argparse
from dataclasses import replace
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cases',nargs='+',help='Explicit named subset; no silent case substitution')
    args=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('This experimental Mac authorization does not enable other hosts')
    if version('lightweaver')!='0.16.1':raise ValueError('Pinned LW0.16.1 required')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.bride_parameter_checks import controlled_structure,assert_collision_coverage,marcs_filename,interpolate_marcs_structure,diagnostic_collision_scaling,positive_stat_equil
    from pipeline.bride_observables import flux_quadrature,integrate_flux
    from pipeline.wavelength_util import air_to_vac
    plan=json.loads(args.plan.read_text());controls=plan['controls']
    cases=plan['cases']
    if args.cases:
        if set(args.cases)-{c['name']for c in cases}:raise ValueError('Unknown requested case')
        cases=[c for c in cases if c['name']in args.cases]
    args.output.mkdir(parents=True,exist_ok=False)
    source_path=require_codex_path('experiments.bride_oxygen_atom')
    source=read_multi_atom(source_path)
    if fingerprint(source_path)['sha256']!='d058b75bdec80c5317523a6d7f185a64ca59144d3d494a17eeabee5f10d04b57':
        raise ValueError('Changed held oxygen atom')
    lo,hi=plan['window_air_A'];wave=np.arange(lo,hi+plan['step_A']/2,plan['step_A']);vac=air_to_vac(wave)*.1
    summary=dict(schema='codex.bride_mac_parameter_results/1',status='RUNNING',
        scientific_admission='HOLD: diagnostic input responses, not a solar fit or uncertainty budget',
        host=dict(system=platform.system(),machine=platform.machine()),
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),
        plan=fingerprint(args.plan),atom=fingerprint(source_path),cases=[])
    def save():(args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    save()
    for case in cases:
        name=case['name'];out=args.output/name;out.mkdir();start=time.perf_counter()
        m=dict(case=case,status='PREPARING',controls=controls,threads=1,
            interpretation='xi changes formal-transfer/SE microturbulence at fixed MARCS t01 structure; Teff/logg use actual independent atmosphere files')
        def save_case():(out/'manifest.json').write_text(json.dumps(m,indent=2,allow_nan=False)+'\n')
        save_case()
        try:
            z=case.get('feh',0.);alpha=case.get('alpha',0.);store=require_codex_path('grids.marcs_standard')
            if case.get('interpolate',False):
                nodes=[];sources=[]
                for teff in case['bracket_teff']:
                    for logg in case['bracket_logg']:
                        path=store/marcs_filename(teff,logg,z,alpha)
                        nodes.append(((teff,logg),read_marcs_structure(path)));sources.append(fingerprint(path))
                original,interpolation=interpolate_marcs_structure(nodes,(case['teff'],case['logg']),
                    thermodynamic_scale=case.get('thermodynamic_scale','linear'))
                provenance=dict(kind='EXPERIMENTAL_INTERPOLATED_MARCS',sources=sources,interpolation=interpolation)
                m['interpretation']='Full thermodynamic-field interpolation on common native tauRoss layers; not a fresh atmosphere-equilibrium solution'
            else:
                path=store/marcs_filename(case['teff'],case['logg'],z,alpha)
                original=read_marcs_structure(path);provenance=fingerprint(path)
            structure=controlled_structure(original,xi=case['xi'],drop_bottom=case['drop_bottom'])
            scales=case.get('diagnostic_collision_scales',{})
            atom=replace(diagnostic_collision_scaling(source,scales),abundance=case['abundance'])
            m['diagnostic_collision_scales']=scales
            if scales:m['interpretation']+='; explicit collision-class scaling is a causal sensitivity control, not adopted atomic data'
            assert_collision_coverage(atom,structure['temperature_K'])
            m.update(atmosphere=provenance,original_depths=original['ndep'],retained_depths=structure['ndep'],
                temperature_range_K=[float(structure['temperature_K'].min()),float(structure['temperature_K'].max())],
                bottom_logtau500=float(structure['logtau_500'][-1]),status='SOLVING')
            hydrogen_lines=case.get('lte_h_background_lines',False)
            uv_lines=case.get('lte_uv_background_lines',False)
            top_cap=case.get('finite_column_top_cap',False)
            if top_cap and case.get('interpolate',False):
                raise ValueError('Finite-column cap requires a native atmosphere column, not interpolated metadata')
            save_case();ctx=make_reference_context(structure,atom,rays=controls['rays'],qscale=controls['qscale'],ng=controls.get('ng',True),
                lte_h_background_lines=hydrogen_lines,lte_uv_background_lines=uv_lines,
                finite_column_top_cap=top_cap)
            if top_cap:
                from pipeline import bride_top_cap
                m['finite_column_top_cap']=dict(adapter=fingerprint(bride_top_cap.__file__),
                    column_mass_top_kg_m2=float(structure['column_mass_kg_m2'][0]),
                    density_top_kg_m3=float(structure['rho_kg_m3'][0]),
                    effective_length_m=ctx.atmos.pyAtmos.zUpperBc.length_m,
                    interpretation='Native top mass column with uniform top opacity/source; no fitted thickness or external illumination; not exact producer proof')
            if uv_lines:
                from lightweaver import rh_atoms
                m['lte_uv_background_lines']=dict(source=fingerprint(rh_atoms.__file__),
                    species=[dict(element=str(a.atomicModel.element),lines=len(a.atomicModel.lines),
                                  continua=len(a.atomicModel.continua))for a in ctx.detailedAtoms],
                    interpretation='Native passive atoms with 100-300 nm line centers moved to fixed LTE/CRD; native continua and spectral sampling also change; incomplete line blanketing')
            if hydrogen_lines:
                from lightweaver import rh_atoms
                m['lte_h_background_lines']=dict(source=fingerprint(rh_atoms.__file__),
                    interpretation='Native H_6_atom fixed LTE populations, CRD lines; not an H NLTE/PRD model')
            convention=case.get('experimental_collision_convention')
            if convention:
                if convention in ('native_tautsp_rh_directions','native_tautsp_source_directions','native_tautsp_lte_balanced_CH0','barklem2018_MN_LTE_pair'):
                    from pipeline.bride_native_spline import NativeTautsp,install_native_collision_control
                    oracle=NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
                    m['native_spline_control']=install_native_collision_control(ctx,atom,oracle,
                        source_directions=convention!='native_tautsp_rh_directions',
                        reversible_ch0=convention=='native_tautsp_lte_balanced_CH0')
                    if convention=='barklem2018_MN_LTE_pair':
                        from pipeline.bride_charge_contract import install_author_mn_pair
                        m['author_charge_control']=install_author_mn_pair(ctx,atom,
                            require_codex_path('experiments.bride_barklem_oxygen'),oracle)
                elif convention in ('directed_CH0_lograte','gencol_CH_down_directed_CH0_lograte'):
                    from pipeline.bride_gencol_controls import install_gencol_control
                    install_gencol_control(ctx,atom,change_ch=convention.startswith('gencol_CH_down'))
                else:
                    raise ValueError('Unknown experimental collision convention')
                m['experimental_collision_convention']=convention
                m['interpretation']+='; non-adopted GENCOL convention experiment, exact reference-build parity unproven'
            active=ctx.activeAtoms[0]
            with (out/'history.jsonl').open('x')as history:
                for iteration in range(controls['maxiter']):
                    if time.perf_counter()-start>controls['seconds_per_case']:raise TimeoutError('Case budget exceeded')
                    old=active.n.copy()
                    old_J_top=np.array(ctx.spect.J)[:,0].copy() if top_cap else None
                    update=ctx.formal_sol_gamma_matrices(lambdaIterate=True)
                    cap_change=ctx.atmos.pyAtmos.zUpperBc.update_from_context(ctx,old_J_top) if top_cap else 0.
                    matrix=active.Gamma;scale=np.max(abs(matrix),axis=1,keepdims=True)
                    if np.any(scale==0):raise ValueError('Disconnected SE equation')
                    matrix[:]/=scale
                    rejected=positive_stat_equil(ctx,active,old,reject_invalid_ng=controls.get('reject_invalid_ng',False))
                    n=active.n.copy()
                    if not np.isfinite(n).all()or np.any(n<=0):raise ValueError('Invalid populations')
                    rates=physical_rate_matrix(active)
                    balance=np.einsum('ijd,jd->id',rates,n);flow=np.einsum('ijd,jd->id',abs(rates),n)
                    rate=float(np.max(abs(balance)/np.maximum(flow,np.finfo(float).tiny)))
                    delta=float(np.max(abs(n-old)/n));conservation=float(np.max(abs(n.sum(axis=0)/active.nTotal-1)))
                    entry=dict(iteration=iteration+1,population_change=delta,J_change=float(update.dJMax),
                        physical_rate_residual=rate,particle_conservation_error=conservation,rejected_invalid_ng=bool(rejected))
                    if top_cap:entry['cap_boundary_change']=cap_change
                    history.write(json.dumps(entry,allow_nan=False)+'\n');history.flush()
                    if conservation>1e-10:raise ValueError('Particle conservation failed')
                    if iteration>=2 and max(delta,rate,update.dJMax,cap_change)<controls['rtol']:break
                else:raise RuntimeError('Joint convergence not reached')
            m.update(status='FORMING_FLUX',iterations=iteration+1,final_diagnostics=entry);save_case()
            continuum=make_reference_context(structure,atom,rays=controls['rays'],qscale=controls['qscale'],
                fixed_populations=n,continuum_only=True,lte_h_background_lines=hydrogen_lines,
                lte_uv_background_lines=uv_lines,finite_column_top_cap=top_cap)
            for j in range(1000):
                if time.perf_counter()-start>controls['seconds_per_case']:raise TimeoutError('Continuum budget exceeded')
                old_J_top=np.array(continuum.spect.J)[:,0].copy() if top_cap else None
                u=continuum.formal_sol_gamma_matrices()
                cap_change=continuum.atmos.pyAtmos.zUpperBc.update_from_context(continuum,old_J_top) if top_cap else 0.
                if max(u.dJMax,cap_change)<1e-8:break
            else:raise RuntimeError('Continuum J convergence failed')
            mu,w=flux_quadrature(controls['flux_rays'])
            if top_cap:
                from pipeline.bride_top_cap import emergent_rays
                line_grid,line_out=emergent_rays(ctx,vac,mu)
                cont_grid,cont_out=emergent_rays(continuum,vac,mu)
                line=integrate_flux(line_out,mu,w);cont=integrate_flux(cont_out,mu,w)
                np.savetxt(out/'gridpoint_flux.csv',np.column_stack((wave,integrate_flux(line_grid,mu,w),integrate_flux(cont_grid,mu,w))),
                    delimiter=',',header='wavelength_air_A,gridpoint_Fnu,gridpoint_continuum_Fnu',comments='')
                cap=ctx.atmos.pyAtmos.zUpperBc
                np.savez_compressed(out/'top_cap.npz',wavelength_vac_nm=cap.wave,chi_m_inv=cap.chi,
                    eta_including_scattering=cap.eta,effective_length_m=cap.length_m)
                m['flux_boundary']='Outer cap face; separate gridpoint flux retained'
            else:
                line=integrate_flux(ctx.compute_rays(wavelengths=vac,mus=mu),mu,w)
                cont=integrate_flux(continuum.compute_rays(wavelengths=vac,mus=mu),mu,w)
            np.savetxt(out/'flux.csv',np.column_stack((wave,line,cont,line/cont)),delimiter=',',
                header='wavelength_air_A,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,normalized_flux',comments='')
            np.savez_compressed(out/'populations.npz',populations=n,nstar=active.nStar,temperature_K=structure['temperature_K'])
            m.update(status='CONVERGED_EXPERIMENT',elapsed_seconds=time.perf_counter()-start,
                flux=fingerprint(out/'flux.csv'),continuum_iterations=j+1);save_case()
        except Exception as error:
            m.update(status='FAILED',error=f'{type(error).__name__}: {error}',elapsed_seconds=time.perf_counter()-start);save_case()
        summary['cases'].append(dict(name=name,status=m['status'],manifest=fingerprint(out/'manifest.json')));save()
        print(name,m['status'],m.get('elapsed_seconds'),m.get('error',''),flush=True)
    summary['status']='EXPERIMENTS_FINISHED';save()
    if any(c['status']!='CONVERGED_EXPERIMENT' for c in summary['cases']):raise SystemExit(1)


if __name__=='__main__':main()
