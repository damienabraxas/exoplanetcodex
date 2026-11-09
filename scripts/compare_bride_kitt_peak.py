"""Fixed-abundance Kitt Peak diagnostic: replay O populations as disk flux.

Not a solar abundance fit. No observed continuum, RV, width or physics fitting.
The benchmark atmosphere, missing blends and unbroadened line-shape limitations
remain explicit. All inputs are read-only; output directory must be new.
"""
from __future__ import annotations
import argparse
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--experiment',type=Path,required=True)
    p.add_argument('--comparison',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seconds',type=float,default=900.)
    args=p.parse_args()
    from config.constants import assert_on_sirius,require_codex_path,get_star_params
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.multi_atom import read_multi_atom
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.wavelength_util import air_to_vac
    from pipeline.bride_observables import flux_quadrature,integrate_flux,window_integral
    assert_on_sirius('Bride Kitt Peak fixed-input diagnostic')
    if version('lightweaver')!='0.16.1':raise ValueError('Pinned LW0.16.1 required')
    original=json.loads((args.experiment/'manifest.json').read_text())
    if original['status']!='CONVERGED_EXPERIMENT':raise ValueError('Unconverged source')
    receipt=json.loads(args.comparison.read_text())
    hits=[r for r in receipt['experiments'] if r['name']==args.experiment.name]
    if len(hits)!=1:raise ValueError('Missing/ambiguous experiment receipt')
    for name in ('manifest.json','populations.npz','triplet_profiles.csv','history.jsonl'):
        if fingerprint(args.experiment/name)!=hits[0]['inputs'][name]:raise ValueError(f'Changed experiment {name}')
    ref=json.loads(args.reference.read_text())
    if fingerprint(args.reference)!=original['reference_file']:raise ValueError('Changed reference')
    atom_path=require_codex_path('grids.gerber_ts')/ref['atom']['name']
    atmosphere_path=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name']
    for path,expected in [(atom_path,original['atom']),(atmosphere_path,original['atmosphere'])]:
        if fingerprint(path)!=expected:raise ValueError('Changed physical input')
    observed_path=require_codex_path('data.kurucz2005_residual')
    with observed_path.open() as f:header=''.join(next(f)for _ in range(5))
    if 'KURUCZ RESIDUAL IRRADIANCE ATLAS 2005' not in header or 'vacuum wavelength including gravitational red shift' not in header:
        raise ValueError('Unexpected Kitt Peak product/frame')
    args.output.mkdir(parents=True,exist_ok=False)
    m=dict(schema='codex.bride_kitt_peak_diagnostic/1',status='PREPARING',
        scientific_admission='HOLD: fixed benchmark-node comparison, no solar abundance',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),
        source_manifest=fingerprint(args.experiment/'manifest.json'),comparison=fingerprint(args.comparison),
        populations=fingerprint(args.experiment/'populations.npz'),
        observed_holding='solar_kpno_kurucz2005_corrected',observed_key='data.kurucz2005_residual',
        observed=fingerprint(observed_path),observed_header=header,
        observed_normalization='Published residual irradiance; no continuum refit',
        observed_frame='vacuum nm including gravitational redshift -> canonical air Angstrom; no velocity shift fitted',
        model=original['atmosphere'],fixed_abundance=original['fixed_abundance'],canonical_solar=get_star_params('solar'),
        observable='disk-integrated plane-parallel F_nu=2pi integral I_nu(mu)*mu dmu',
        continuum='Separate bound-bound-free context retaining fixed NLTE bound-free populations',
        broadening='Intrinsic model only: source microturbulence; no macro/rotation/FTS convolution',
        caveats=['5750/4.5 benchmark node differs from canonical Sun',
                 'Missing observed blends/molecular opacity; original MULTI mismatch remains',
                 'No observational uncertainty vector or local telluric residual budget; no chi-square/admission',
                 'Delivered gravitational shift and convective shifts not fitted; shape residuals are diagnostic'],
        controls=dict(fixed_population_replay=True,no_SE_update=True,J_threshold=1e-8,
                      maxiter=1000,flux_rays=[5,10],seconds=args.seconds,threads=1))
    def save():(args.output/'manifest.json').write_text(json.dumps(m,indent=2,allow_nan=False)+'\n')
    save();start=time.perf_counter()
    try:
        structure=read_marcs_structure(atmosphere_path);source=read_multi_atom(atom_path)
        a=np.load(args.experiment/'populations.npz');pop=a['populations']
        c=original['controls'];contexts=[];iterations=[]
        for continuum in (False,True):
            ctx=make_reference_context(structure,source,rays=c['rays'],qscale=c['qscale'],
                fixed_populations=pop,continuum_only=continuum)
            from lightweaver.atomic_table import PeriodicTable
            np.testing.assert_allclose(ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].nStar,
                                       a['nstar'],rtol=64*np.finfo(float).eps,atol=0)
            for i in range(1000):
                if time.perf_counter()-start>args.seconds:raise TimeoutError('Replay budget exceeded')
                update=ctx.formal_sol_gamma_matrices()
                if update.dJMax<1e-8:break
            else:raise RuntimeError('Fixed-population radiation failed convergence')
            contexts.append(ctx);iterations.append(i+1)
        old=np.genfromtxt(args.experiment/'triplet_profiles.csv',delimiter=',',names=True)
        wave=old['wavelength_air_A'];vac=air_to_vac(wave)*.1
        diskcentre=contexts[0].compute_rays(wavelengths=vac,mus=[1.])
        m['diskcentre_replay_max_relative_difference']=float(np.max(abs(diskcentre/old['I_nlte_W_m2_sr_Hz']-1)))
        fluxes=[]
        for rays in (5,10):
            mu,weight=flux_quadrature(rays);pair=[]
            for ctx in contexts:
                intensity=ctx.compute_rays(wavelengths=vac,mus=mu)
                if intensity.shape!=(len(wave),rays):raise ValueError('Unexpected emergent array geometry')
                pair.append(integrate_flux(intensity,mu,weight))
            fluxes.append(pair)
        norm=fluxes[1][0]/fluxes[1][1]
        wobs,fobs=read_kurucz2005_residual(observed_path,wave[0],wave[-1])
        if len(wobs)<10 or wobs[0]>wave[0] or wobs[-1]<wave[-1]:raise ValueError('Observed coverage missing')
        if not np.isfinite(fobs).all() or not np.all(np.diff(wobs)>0):raise ValueError('Invalid observation')
        ongrid=np.interp(wave,wobs,fobs)
        # Immutable line identifiers, already audited in the parent O777 decision.
        centres=[7771.944,7774.166,7775.388]
        rows=[]
        for centre in centres:
            # Explicit narrow fixed-aperture diagnostic, not a fitted total line EW.
            lo,hi=centre-.5,centre+.5;mask=(wave>=lo)&(wave<=hi)
            rows.append(dict(centre_air_A=centre,window_air_A=[lo,hi],
                model_window_area_mA=1000*window_integral(wave,norm,lo,hi),
                observed_window_area_mA=1000*window_integral(wobs,fobs,lo,hi),
                model_min_normalized_flux=float(norm[mask].min()),
                observed_min_normalized_flux=float(fobs[(wobs>=lo)&(wobs<=hi)].min()),
                unweighted_profile_RMS=float(np.sqrt(np.mean((norm[mask]-ongrid[mask])**2)))))
        m.update(status='DIAGNOSTIC_COMPLETE',radiation_iterations=iterations,
            flux_quadrature_max_relative_difference=float(np.max(abs((fluxes[0][0]/fluxes[0][1])/norm-1))),
            observed_samples_in_window=int(((wobs>=wave[0])&(wobs<=wave[-1])).sum()),
            lines=rows,elapsed_seconds=time.perf_counter()-start)
        np.savetxt(args.output/'profiles.csv',np.column_stack((wave,vac,fluxes[1][0],fluxes[1][1],norm,ongrid,norm-ongrid)),
            delimiter=',',header='wavelength_air_A,wavelength_vacuum_nm,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,model_normalized_flux,kitt_peak_residual_flux,model_minus_observed',comments='')
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,1,sharex=True,figsize=(10,6),height_ratios=[3,1])
        axes[0].plot(wave,ongrid,label='Kitt Peak 2005 corrected residual flux',color='black')
        axes[0].plot(wave,norm,label='Bride fixed A(O)=8.70, intrinsic flux',color='tab:orange')
        axes[0].set_ylabel('Continuum-normalized flux');axes[0].legend(fontsize=9)
        axes[0].set_title('Diagnostic only: MARCS5750/4.5; no abundance/RV/width fit')
        axes[1].plot(wave,norm-ongrid,color='tab:blue');axes[1].axhline(0,color='grey',lw=.5)
        axes[1].set_ylabel('Model − observed');axes[1].set_xlabel('Air wavelength (Å)')
        fig.tight_layout();fig.savefig(args.output/'comparison.png',dpi=160);plt.close(fig)
        save();print(json.dumps({k:m[k]for k in ['status','diskcentre_replay_max_relative_difference',
             'flux_quadrature_max_relative_difference','lines','elapsed_seconds']},indent=2))
    except Exception as error:
        m.update(status='FAILED',error=f'{type(error).__name__}: {error}');save();raise


if __name__=='__main__':main()
