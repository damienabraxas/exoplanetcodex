"""Fixed-population Mac rate/radiation audit at the exact native reference node."""
from pathlib import Path
import argparse
import csv
import json
import platform
import subprocess
import sys
import time
from importlib.metadata import version
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if platform.system()!='Darwin' or version('lightweaver')!='0.16.1':
        raise RuntimeError('Requires authorized Mac and LW0.16.1')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.bride_rate_audit import process_flow,collision_classes
    from pipeline.multi_atom import read_multi_atom
    from scripts.audit_bride_oxygen_populations import conserve_departure_replay
    import lightweaver.constants as const
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if fingerprint(args.reference)['sha256']!='c3113efdb21e3f40d79c1836b0ddb0ebb41dbceecaf04d622af6b7a8d8c080ec':
        raise ValueError('Changed reference')
    ref=json.loads(args.reference.read_text())
    atom_path=require_codex_path('experiments.bride_oxygen_atom')
    path=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name']
    for file,receipt in ((atom_path,ref['atom']),(path,ref['atmosphere_file'])):
        if fingerprint(file)!=receipt:raise ValueError('Changed native physical input')
    origin=require_codex_path('experiments.bride_mac_grid')/'rya1281-v1/baseline'
    archive=ROOT/'data/results/rya1281/mac_parameter_checks_v1/analysis.json'
    receipt=next(c for c in json.loads(archive.read_text())['cases'] if c['name']=='baseline')
    if fingerprint(origin/'populations.npz')!=receipt['original_population_file']:
        raise ValueError('Changed cached populations')
    data=np.load(origin/'populations.npz');nstar=data['nstar'];native=data['populations']
    external,factors=conserve_departure_replay(np.array(ref['reference']['departures']),nstar)
    source=read_multi_atom(atom_path);structure=read_marcs_structure(path)
    np.testing.assert_array_equal(structure['temperature_K'],data['temperature_K'])
    args.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema='codex.bride_rate_audit/1',status='RUNNING',scientific_admission='HOLD',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),
        accounting=fingerprint(ROOT/'pipeline/bride_rate_audit.py'),reference=fingerprint(args.reference),
        atom=ref['atom'],atmosphere=ref['atmosphere_file'],population_receipt=fingerprint(archive),
        populations=fingerprint(origin/'populations.npz'),threads=1,
        controls=dict(rays=3,qscale=1,J_rtol=1e-8,maxiter=1000,seconds_per_mode=180,SE_updates=False),
        caveat='External b on local LTE is a shape intervention, not native MULTI radiation/rates. No native MULTI rate dump is held here.',modes=[])
    def save():(args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save();rows=[]
    for mode,pop in [('local_nlte',native),('local_lte',nstar),('external_b_local_lte',external)]:
        start=time.perf_counter();ctx=make_reference_context(structure,source,rays=3,qscale=1)
        atom=ctx.activeAtoms[0];np.testing.assert_allclose(atom.nStar,nstar,rtol=64*np.finfo(float).eps,atol=0)
        atom.n[:]=pop
        for iteration in range(1000):
            if time.perf_counter()-start>180:raise TimeoutError('Rate replay budget exceeded')
            update=ctx.formal_sol_gamma_matrices(lambdaIterate=True)
            if update.dJMax<1e-8:break
        else:raise RuntimeError('Radiation did not converge')
        np.testing.assert_array_equal(atom.n,pop)
        eq=ctx.kwargs['eqPops'];pyatmos=ctx.atmos.pyAtmos
        groups=collision_classes(atom.atomicModel,pyatmos,eq)
        atom.compute_collisions(fillDiagonal=False)
        summed=sum(groups.values())
        np.testing.assert_allclose(summed,atom.C,rtol=256*np.finfo(float).eps,atol=0)
        radiative=np.zeros_like(summed);processes=dict(groups)
        identities={}
        for index,transition in enumerate(atom.trans):
            key=f'R{index}'
            r=np.zeros_like(summed)
            r[transition.j,transition.i]=transition.Rij
            r[transition.i,transition.j]=transition.Rji
            radiative+=r;processes[key]=r
            identities[key]=dict(lower_0based=transition.i,upper_0based=transition.j,
                kind=type(transition.transModel).__name__,wavelength_vac_nm=float(transition.lambda0))
        physical=physical_rate_matrix(atom)
        total=summed+radiative
        # Independently reconstruct the same destination/source operator.
        for i in range(len(pop)):total[i,i]=-total[:,i].sum(axis=0)
        np.testing.assert_allclose(total,physical,rtol=256*np.finfo(float).eps,atol=0)
        net,gross=process_flow(total,pop)
        diagnostic=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
        net_groups={key:process_flow(matrix,pop)[0] for key,matrix in processes.items()}
        np.testing.assert_allclose(sum(net_groups.values()),net,
            atol=512*np.finfo(float).eps*float(gross.max()),rtol=0)
        selected=[]
        for tau in (-3.,-2.,-1.,0.):
            d=int(np.argmin(abs(structure['logtau_500']-tau)))
            for level in (5,7,8,9,35):
                rank=sorted(net_groups,key=lambda k:abs(net_groups[k][level,d]),reverse=True)
                selected.append(dict(logtau500=float(structure['logtau_500'][d]),level_0based=level,
                    level_label=source.levels[level].label,net_over_gross=float(net[level,d]/gross[level,d]),
                    strongest=[dict(process=k,net_over_total_gross=float(net_groups[k][level,d]/gross[level,d]),
                                    **identities.get(k,{})) for k in rank[:10]]))
                for key in processes:
                    rows.append(dict(mode=mode,logtau500=structure['logtau_500'][d],level_0based=level,
                        process=key,net_m3_s=net_groups[key][level,d],total_gross_m3_s=gross[level,d]))
        wave=np.array(ctx.spect.wavelength);J=np.array(ctx.spect.J)
        nu=const.CLight/(wave[:,None]*1e-9)
        B=2*const.HPlanck*nu**3/const.CLight**2/np.expm1(const.HPlanck*nu/(const.KBoltzmann*structure['temperature_K']))
        np.savez_compressed(args.output/f'{mode}_rates.npz',populations=pop,nstar=nstar,
            logtau500=structure['logtau_500'],temperature_K=structure['temperature_K'],
            physical_rates_s_inv=total,radiative_rates_s_inv=radiative,
            **{f'C_{key}_s_inv':value for key,value in groups.items()},
            wavelength_vac_nm=wave,Jnu_W_m2_sr_Hz=J,Bnu_W_m2_sr_Hz=B)
        report['modes'].append(dict(name=mode,iterations=iteration+1,elapsed_seconds=time.perf_counter()-start,
            physical_imbalance_max=diagnostic,selected_depth_budgets=selected,
            process_identities=identities,array=fingerprint(args.output/f'{mode}_rates.npz')));save()
        fig,axes=plt.subplots(1,2,figsize=(12,4))
        for level in (5,9,35):
            for key,matrix in groups.items():
                ng=process_flow(matrix,pop)[0]
                axes[0].plot(structure['logtau_500'],ng[level]/gross[level],label=f'{key} level{level}')
        axes[0].set(xlabel='log tau500',ylabel='Signed collision flow / total gross',xlim=(-4,1))
        for w in (130.,277.,777.4):
            i=int(np.argmin(abs(wave-w)));axes[1].plot(structure['logtau_500'],J[i]/B[i],label=f'{wave[i]:.2f} nm vacuum')
        axes[1].set(xlabel='log tau500',ylabel='Jnu / Bnu',xlim=(-4,1))
        for ax in axes:ax.legend(fontsize=6)
        fig.suptitle(f'{mode}: local LW rates and radiation, fixed populations');fig.tight_layout()
        fig.savefig(args.output/f'{mode}_budgets.png',dpi=150);plt.close(fig)
        print(mode,diagnostic,flush=True)
    with (args.output/'selected_flows.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    report['status']='DIAGNOSTIC_COMPLETE';save()


if __name__=='__main__':main()
