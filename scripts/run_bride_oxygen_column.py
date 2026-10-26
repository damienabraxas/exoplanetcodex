"""Bounded Sirius-only O41f/MARCS validation experiment; no abundance fitting.

Default prepares an immutable manifest only. --run executes the model-under-test
against an independently produced exact-node MULTI population reference. It does
not admit a solar abundance, assert paper replication, publish, or run a cube.
"""
from __future__ import annotations
import argparse
from importlib.metadata import version
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--run',action='store_true')
    p.add_argument('--rays',type=int,default=3)
    p.add_argument('--qscale',type=int,default=1)
    p.add_argument('--rtol',type=float,default=1e-6)
    p.add_argument('--maxiter',type=int,default=300)
    p.add_argument('--seconds',type=float,default=900.)
    p.add_argument('--initial',choices=['lte','tilted'],default='lte')
    p.add_argument('--lambda-iterate',action='store_true',help='Diagnostic unaccelerated transfer/SE iteration')
    p.add_argument('--ng',action='store_true',help='Native third-order Ng extrapolation every five iterations after ten')
    args=p.parse_args()
    from config.constants import assert_on_sirius,require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.wavelength_util import air_to_vac,vac_to_air
    assert_on_sirius('Bride O41f validation experiment')
    if version('lightweaver')!='0.16.1':raise ValueError('Only LW0.16.1 is validated')
    if not(args.rays>=2 and args.qscale>=1 and 0<args.rtol<1 and args.maxiter>0 and args.seconds>0):
        raise ValueError('Invalid numerical controls')
    ref=json.loads(args.reference.read_text())
    if ref['schema']!='codex.bride_oxygen_reference/1':raise ValueError('Wrong reference schema')
    held=require_codex_path('grids.gerber_ts')
    source=held/ref['atom']['name']
    atmosphere_path=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name']
    for path,record in [(source,ref['atom']),(atmosphere_path,ref['atmosphere_file'])]:
        if fingerprint(path)!=record:raise ValueError(f'Changed input: {path}')
    structure=read_marcs_structure(atmosphere_path)
    source_atom=read_multi_atom(source)
    if source_atom.abundance!=ref['reference']['abundance']:raise ValueError('Reference abundance mismatch')
    if structure['model_id']!=ref['reference']['model_id']:raise ValueError('Reference model identity mismatch')
    manifest=dict(schema='codex.bride_oxygen_experiment/1',kind='exact-node model validation',
        scientific_admission='HOLD: no solar abundance or full-3D claim',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),converter=fingerprint(ROOT/'pipeline/multi_atom.py'),
        continuum_adapter=fingerprint(ROOT/'pipeline/multi_continuum.py'),
        reference_file=fingerprint(args.reference),reference_record_sha256=ref['reference']['record_sha256'],
        atom=fingerprint(source),atmosphere=fingerprint(atmosphere_path),
        atmosphere_key='grids.marcs_standard',atom_key='grids.gerber_ts',
        fixed_abundance=source_atom.abundance,abundance_origin='unchanged atom header and exact reference grid axis',
        background_abundances_origin='all 92 MARCS header values; O overridden by reference trace abundance',
        electron_density='fixed source Pe/(pinned-engine k*T)',hydrogen_total='source rho*X/(pinned-engine amu*H mass)',
        gravity='source MARCS benchmark cgs converted to SI; this node is not canonical solar STAR_PARAMS',
        boundaries=dict(lower='ThermalisedRadiation',upper='ZeroRadiation'),
        velocities='zero static reference atmosphere; microturbulence from MARCS header',
        geometry='1D plane parallel; emergent vertical intensity, not flux',
        output_spectrum=dict(wavelength_air_A=[7770.5,7776.8],points=1261,intensity_unit='W m^-2 sr^-1 Hz^-1'),
        background='LW default continuum with explicitly inventoried passive native atoms; molecular opacity not included',
        known_reference_differences=['LW vs MULTI background/EOS/numerical implementations not proven identical',
            'Held 41-level reduction not proven identical to cited 122-level published atom',
            'No Amarsi2018 atom or 3D solar abundance replication claimed'],
        controls=dict(rays=args.rays,qscale=args.qscale,rtol=args.rtol,maxiter=args.maxiter,
            seconds=args.seconds,initial=args.initial,lambda_iterate=args.lambda_iterate,
            ng=dict(order=3,period=5,delay=10,threshold=.05,lower_threshold=0.) if args.ng else None,
            threads=1,formal_solver='piecewise_besser_1d',
            row_equilibration='linear equations scaled by max absolute row coefficient; physical rates unchanged'),
        environment={name:version(name)for name in ('lightweaver','numpy','scipy','numba')},
        host_load=list(__import__('os').getloadavg()),status='PREPARING',science_executed=False,solver_started=False)
    args.output.mkdir(parents=True,exist_ok=False)
    def save():
        (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    save()
    import lightweaver as lw
    try:
        preview=to_lightweaver(source_atom)
        preview_set=lw.RadiativeSet([preview]);preview_set.set_active('O')
        preview_spect=preview_set.compute_wavelength_grid()
        for continuum in preview.continua:
            alpha=continuum.alpha(continuum.wavelength())
            if not np.isfinite(alpha).all()or np.any(alpha<0):
                raise ValueError('Invalid tabulated continuum cross-sections')
        manifest.update(status='PREPARED',base_frequency_points=len(preview_spect.wavelength),
            continuum_semantics='MULTI explicit absolute wavelength tables retained; energy normalization unchanged',
            continuum_edges=[dict(lower=c.i+1,upper=c.j+1,energy_edge_nm=c.energy_edge_nm,
                table_edge_nm=c.lambdaEdge,points=len(c.wavelength()))for c in preview.continua])
        save()
    except Exception as error:
        manifest.update(status='INPUT_REJECTED',error=f'{type(error).__name__}: {error}')
        save();raise
    if not args.run:
        print('Prepared only:',args.output/'manifest.json');return
    import lightweaver.constants as constants
    from lightweaver.atomic_table import PeriodicTable
    from lightweaver import rh_atoms
    start_wall=time.perf_counter();start_cpu=time.process_time()
    manifest.update(status='SETTING_UP');save()
    with (args.output/'history.jsonl').open('x')as log:
        try:
            ne=structure['electron_pressure_Pa']/(constants.KBoltzmann*structure['temperature_K'])
            nH=structure['rho_kg_m3']*structure['hydrogen_mass_fraction']/(constants.Amu*PeriodicTable['H'].mass)
            abundance_data={PeriodicTable[i+1]:a for i,a in enumerate(structure['abundances_dex'])}
            abundance_data[PeriodicTable['O']]=source_atom.abundance
            abund=lw.AtomicAbundance(abundance_data)
            passive_names=['H_6_atom','He_9_atom','C_atom','N_atom','MgI_66_atom',
                'Si_atom','Al_atom','Fe_atom','CaII_atom','Na_atom','S_atom']
            def setup(lte):
                n=structure['ndep']
                atmosphere=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,
                    depthScale=structure['height_m'].copy(),temperature=structure['temperature_K'].copy(),
                    ne=ne.copy(),nHTot=nH.copy(),vturb=np.full(n,structure['vturb_m_s']),
                    vlos=np.zeros(n),logG=structure['logg_cgs']-np.log10(100.),
                    lowerBc=lw.ThermalisedRadiation(),upperBc=lw.ZeroRadiation())
                atmosphere.quadrature(args.rays)
                model=to_lightweaver(source_atom)
                if args.qscale!=1:
                    for line in model.lines:line.quadrature.Nlambda=args.qscale*(line.quadrature.Nlambda-1)+1
                radset=lw.RadiativeSet([getattr(rh_atoms,name)()for name in passive_names]+[model],abundance=abund)
                if lte:radset.set_detailed_static('O')
                else:radset.set_active('O')
                eq=radset.compute_eq_pops(atmosphere)
                spect=radset.compute_wavelength_grid()
                ng=lw.NgOptions(Norder=3,Nperiod=5,Ndelay=10,threshold=.05,lowerThreshold=0.) if args.ng and not lte else None
                ctx=lw.Context(atmosphere,spect,eq,Nthreads=1,formalSolver='piecewise_besser_1d',conserveCharge=False,ngOptions=ng)
                return ctx
            ctx=setup(False);atom=ctx.activeAtoms[0]
            manifest.update(passive_atoms=passive_names,wavelength_points=len(ctx.spect.wavelength),
                setup_seconds=time.perf_counter()-start_wall)
            if args.initial=='tilted':
                atom.n[:]*=np.geomspace(.5,2.,atom.n.shape[0])[:,None]
                atom.n[:]*=atom.nTotal[None,:]/atom.n.sum(axis=0)[None,:]
            save()
            manifest.update(status='RUNNING',science_executed=True,solver_started=True);save()
            converged=False
            for iteration in range(args.maxiter):
                if time.perf_counter()-start_wall>args.seconds:raise TimeoutError('Bounded oxygen solver time budget exceeded')
                old=atom.n.copy()
                update=ctx.formal_sol_gamma_matrices(lambdaIterate=args.lambda_iterate)
                matrix=atom.Gamma
                row_scale=np.max(np.abs(matrix),axis=1,keepdims=True)
                if np.any(row_scale==0):raise ValueError('Disconnected oxygen SE equation')
                matrix[:]/=row_scale
                ctx.stat_equil()
                actual=atom.n.copy()
                if not np.isfinite(actual).all()or np.any(actual<=0):raise ValueError('Invalid oxygen populations')
                rates=physical_rate_matrix(atom)
                balance=np.einsum('ijd,jd->id',rates,actual)
                flow=np.einsum('ijd,jd->id',np.abs(rates),actual)
                residual=float(np.max(np.abs(balance)/np.maximum(flow,np.finfo(float).tiny)))
                delta=float(np.max(np.abs(actual-old)/actual))
                conservation=float(np.max(np.abs(actual.sum(axis=0)/atom.nTotal-1)))
                entry=dict(iteration=iteration+1,population_change=delta,J_change=float(update.dJMax),
                    physical_rate_residual=residual,particle_conservation_error=conservation,
                    elapsed_seconds=time.perf_counter()-start_wall)
                log.write(json.dumps(entry,allow_nan=False)+'\n');log.flush()
                if conservation>1e-10:raise ValueError('Particle conservation failed')
                if iteration>=2 and max(delta,residual,update.dJMax)<args.rtol:
                    converged=True;break
            if not converged:
                np.savez_compressed(args.output/'unconverged_diagnostics.npz',
                    populations=actual,nstar=atom.nStar,physical_rate_matrix=rates,
                    scaled_mali_matrix=matrix,collisions=atom.C,
                    radiative_up=np.asarray([t.Rij for t in atom.trans]),
                    radiative_down=np.asarray([t.Rji for t in atom.trans]),
                    transition_indices=np.asarray([(t.i,t.j) for t in atom.trans]))
                raise RuntimeError('Oxygen did not meet joint population/J/physical-rate convergence')
            b=atom.n/atom.nStar
            reference_b=np.asarray(ref['reference']['departures'])
            if b.shape!=reference_b.shape:raise ValueError('Reference population layout mismatch')
            np.savez_compressed(args.output/'populations.npz',populations=atom.n,nstar=atom.nStar,
                departures=b,reference_departures=reference_b,physical_rate_matrix=rates,
                temperature_K=structure['temperature_K'],height_m=structure['height_m'],ne_m3=ne,
                nH_total_m3=nH,reference_tau500=np.asarray(ref['reference']['tau']))
            wave_air=np.linspace(7770.5,7776.8,1261)
            wave_nm=air_to_vac(wave_air)*.1
            intensity=ctx.compute_rays(wavelengths=wave_nm,mus=[1.])
            lte_ctx=setup(True)
            # Fixed LTE populations still require continuum-scattering J convergence.
            lte_converged=False
            for _ in range(args.maxiter):
                if time.perf_counter()-start_wall>args.seconds:raise TimeoutError('LTE comparator time budget exceeded')
                u=lte_ctx.formal_sol_gamma_matrices()
                if u.dJMax<args.rtol:lte_converged=True;break
            if not lte_converged:raise RuntimeError('LTE comparator radiation did not converge')
            intensity_lte=lte_ctx.compute_rays(wavelengths=wave_nm,mus=[1.])
            if not(np.isfinite(intensity).all()and np.isfinite(intensity_lte).all()and np.all(intensity_lte>0)and np.all(intensity>0)):
                raise ValueError('Invalid emergent intensity')
            np.savetxt(args.output/'triplet_profiles.csv',np.column_stack((wave_air,wave_nm,intensity,intensity_lte,intensity/intensity_lte)),
                delimiter=',',header='wavelength_air_A,wavelength_vacuum_nm,I_nlte_W_m2_sr_Hz,I_lte_W_m2_sr_Hz,I_nlte_over_I_lte',comments='')
            triplet_levels=sorted({index-1 for line in source_atom.lines
                if 776 < 1e7/(source_atom.level(line.j).E_cm-source_atom.level(line.i).E_cm) < 779
                for index in (line.i,line.j)})
            comparison={str(level+1):dict(min_b=float(b[level].min()),max_b=float(b[level].max()),
                max_absolute_difference_to_reference=float(np.max(np.abs(b[level]-reference_b[level]))))for level in triplet_levels}
            manifest.update(status='CONVERGED_EXPERIMENT',iterations=iteration+1,
                final_diagnostics=entry,reference_comparison=comparison,
                comparison_verdict='DIAGNOSTIC_ONLY: model/background differences not yet adjudicated',
                total_wall_seconds=time.perf_counter()-start_wall,total_cpu_seconds=time.process_time()-start_cpu,
                peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            save();print(json.dumps({k:manifest[k]for k in ('status','iterations','final_diagnostics','reference_comparison','total_wall_seconds')},indent=2))
        except Exception as error:
            manifest.update(status='FAILED',error=f'{type(error).__name__}: {error}',
                            elapsed_seconds=time.perf_counter()-start_wall)
            save();raise

if __name__=='__main__':main()
