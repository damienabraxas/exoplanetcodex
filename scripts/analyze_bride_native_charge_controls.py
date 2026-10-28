"""Native interpolation, thermal closure and reference-population diagnostics."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',nargs='+',required=True,type=Path)
    p.add_argument('--reference',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_native_spline import NativeTautsp,install_native_collision_control
    from pipeline.bride_rate_audit import process_flow
    from pipeline.bride_observables import window_integral
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if fingerprint(args.reference)['sha256']!='c3113efdb21e3f40d79c1836b0ddb0ebb41dbceecaf04d622af6b7a8d8c080ec':raise ValueError('Changed reference')
    reference=json.loads(args.reference.read_text())
    path=require_codex_path('grids.marcs_standard')/reference['atmosphere_file']['name']
    atom_path=require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(path)!=reference['atmosphere_file'] or fingerprint(atom_path)!=reference['atom']:raise ValueError('Changed physical source')
    structure=read_marcs_structure(path);source=read_multi_atom(atom_path)
    tau=structure['logtau_500'];ref=np.array(reference['reference']['departures'])
    photosphere=(tau>=-3)&(tau<=0);sample=int(np.argmin(abs(tau+1.2621)))
    oracle=NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
    args.output.mkdir(parents=True,exist_ok=False)
    result=dict(schema='codex.bride_native_charge_diagnostics/1',scientific_admission='HOLD',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),native_adapter=fingerprint(ROOT/'pipeline/bride_native_spline.py'),
        context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),reference=fingerprint(args.reference),
        oracle=oracle.receipt,source_atom=reference['atom'],source_atmosphere=reference['atmosphere_file'],
        interpretation='Native coefficients do not establish producer identity. Reverse CH0 is an explicit LTE-reservoir hypothesis, not an adopted rate change.',
        population_comparison_scope='All levels/depths reported; logtau[-3,0] is a declared diagnostic subset, not a contribution function',
        thermal_controls=[],cases=[])
    fig,axes=plt.subplots(1,3,figsize=(15,4))
    for ax,level,label in zip(axes,(5,9,39),('Triplet lower','Triplet upper','O II ground')):
        ax.plot(tau,ref[level],color='black',label='Independent reference');ax.set(xlabel='log tau500',ylabel='Departure b',title=label)
        ax.axvspan(-3,0,alpha=.06,color='grey')
        if level==39:ax.set_yscale('log')
    baseline=require_codex_path('experiments.bride_mac_grid')/'rya1281-v1/baseline/populations.npz'
    receipt_path=ROOT/'data/results/rya1281/mac_parameter_checks_v1/analysis.json'
    receipt=next(c for c in json.loads(receipt_path.read_text())['cases'] if c['name']=='baseline')
    if fingerprint(baseline)!=receipt['original_population_file']:raise ValueError('Changed baseline population')
    base=np.load(baseline);base_b=base['populations']/base['nstar']
    result['baseline_populations']=fingerprint(baseline)
    for ax,level in zip(axes,(5,9,39)):ax.plot(tau,base_b[level],label='Original RH convention',ls=':')
    arrays={}
    for name,directions,reverse in [('original_rh',None,False),('native_rh',False,False),('native_source',True,False),('native_source_reversible_CH0',True,True)]:
        ctx=make_reference_context(structure,source,rays=3,qscale=1)
        if directions is not None:install_native_collision_control(ctx,source,oracle,source_directions=directions,reversible_ch0=reverse)
        atom=ctx.activeAtoms[0];atom.compute_collisions(fillDiagonal=False)
        C=np.array(atom.C);lte=np.array(atom.nStar)
        net,gross=process_flow(C,lte)
        residual=np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny))
        if reverse and residual>1024*np.finfo(float).eps:raise ValueError('Paired full-atom collision LTE closure failed')
        result['thermal_controls'].append(dict(name=name,max_relative_LTE_imbalance=float(residual),
            collider='Native passive LTE H at fixed ne/nH; no radiation included'))
        arrays[name+'_LTE_net_over_gross']=net/np.maximum(gross,np.finfo(float).tiny)
    for directory in args.runs:
        summary=json.loads((directory/'summary.json').read_text())
        if summary['status']!='EXPERIMENTS_FINISHED' or summary['atom']!=reference['atom']:raise ValueError('Unfinished/changed run')
        for entry in summary['cases']:
            root=directory/entry['name'];m=json.loads((root/'manifest.json').read_text())
            if fingerprint(root/'manifest.json')!=entry['manifest']:raise ValueError('Changed run manifest')
            if m['status']!='CONVERGED_EXPERIMENT':raise ValueError('Unconverged source control')
            if m['atmosphere']!=reference['atmosphere_file']:raise ValueError('Wrong reference node')
            a=np.load(root/'populations.npz');b=a['populations']/a['nstar']
            np.testing.assert_allclose(a['nstar'],base['nstar'],rtol=64*np.finfo(float).eps,atol=0)
            if fingerprint(root/'flux.csv')!=m['flux']:raise ValueError('Changed flux')
            curve=np.genfromtxt(root/'flux.csv',delimiter=',',names=True)
            history=[json.loads(line) for line in (root/'history.jsonl').read_text().splitlines()]
            where=np.unravel_index(np.argmax(abs(b-ref)),b.shape)
            result['cases'].append(dict(name=entry['name'],source_run=directory.name,
                source_inputs=m['case'],lte_h_background_lines=m.get('lte_h_background_lines'),
                source_manifest=entry['manifest'],source_summary=fingerprint(directory/'summary.json'),
                populations=fingerprint(root/'populations.npz'),flux=m['flux'],history=fingerprint(root/'history.jsonl'),
                controls=m['controls'],native_spline_control=m['native_spline_control'],
                iterations=m['iterations'],elapsed_seconds=m['elapsed_seconds'],final_diagnostics=m['final_diagnostics'],
                rejected_invalid_ng_count=sum(x.get('rejected_invalid_ng',False) for x in history),
                window_area_mA=[1000*window_integral(curve['wavelength_air_A'],curve['normalized_flux'],c-.5,c+.5) for c in (7771.944,7774.166,7775.388)],
                selected_logtau500=float(tau[sample]),b_lower=float(b[5,sample]),b_upper=float(b[9,sample]),b_OII_ground=float(b[39,sample]),
                ref_b_lower=float(ref[5,sample]),ref_b_upper=float(ref[9,sample]),ref_b_OII_ground=float(ref[39,sample]),
                upper_max_absolute_difference_all_depths=float(np.max(abs(b[9]-ref[9]))),
                upper_max_absolute_difference_logtau_minus3_to0=float(np.max(abs(b[9,photosphere]-ref[9,photosphere]))),
                all_level_depth_max_absolute_difference=float(np.max(abs(b-ref))),
                largest_difference_location=dict(level_0based=int(where[0]),depth_index=int(where[1]),logtau500=float(tau[where[1]]))))
            np.savetxt(args.output/f'{entry["name"]}_departures.csv',np.column_stack((tau,b[5],b[9],b[39],ref[5],ref[9],ref[39])),
                delimiter=',',header='logtau500,b_lower,b_upper,b_OII_ground,ref_b_lower,ref_b_upper,ref_b_OII_ground',comments='')
            for ax,level in zip(axes,(5,9,39)):ax.plot(tau,b[level],label=entry['name'])
    arrays['logtau500']=tau
    np.savez_compressed(args.output/'thermal_closure.npz',**arrays)
    for ax in axes:ax.legend(fontsize=5)
    fig.suptitle('Native coefficient / direction / reverse-path controls; no abundance tuning')
    fig.tight_layout();fig.savefig(args.output/'departures_all_depths.png',dpi=150);plt.close(fig)
    result['status']='DIAGNOSTIC_COMPLETE'
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(thermal_controls=result['thermal_controls'],cases=result['cases']),indent=2))


if __name__=='__main__':main()
