"""Compare every declared xi case with fresh O SE and independently saved rates."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_observables import window_integral
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.telluric_policy import gate_holding
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    import inspect,matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    from pipeline.bride_xi_sweep import validate_plan
    plan=json.loads(args.plan.read_text());validate_plan(plan,json.loads((ROOT/'config/experiments/rya1281_mac_solar_bound_control.json').read_text()))
    receipts=[json.loads((d/'report.json').read_text())for d in args.runs]
    if [r['case'] for r in receipts]!=plan['cases']:raise ValueError('Require every declared case in plan order; no best-case filter')
    for d,r in zip(args.runs,receipts):
        if r['schema']!='codex.bride_xi_flux/1' or r['input_manifest']!=fingerprint(d/'input_manifest.json') or r['plan']!=fingerprint(args.plan):raise ValueError('Wrong/changed xi control input')
    names=[r['case']['name']for r in receipts]
    for key in ('generator','stark_source','stark_kernel','pressure_source','pressure_provider','provider','line_coefficients','lines','partition','source_manifest','source_populations','source_binding','wide_baseline','sweep_adapter'):
        if any(r[key]!=receipts[0][key]for r in receipts):raise ValueError('Mixed sweep inputs: '+key)
    fields=[]
    for r in receipts:
        import copy
        b=copy.deepcopy(r['local_binding']['local_binding']);b['scalars'].pop('vturb_m_s');fields.append(b)
    if any(f!=fields[0]for f in fields):raise ValueError('Sweep changes fields beyond xi')
    baseline=next(i for i,r in enumerate(receipts)if r['case']['xi']==1.)
    rate_checks=[]
    for directory,r in zip(args.runs,receipts):
        for name,key in [('populations.npz','populations'),('rates.npz','rates'),('history.jsonl','history')]:
            if fingerprint(directory/name)!=r['coupled'][key]:raise ValueError('Changed coupled evidence: '+name)
        with np.load(directory/'populations.npz')as p:pop=p['populations']
        with np.load(directory/'rates.npz')as p:
            off=p['radiative_rates_offdiag_s_inv']+p['collision_rates_offdiag_s_inv']
            np.testing.assert_array_equal(off[np.arange(len(off)),np.arange(len(off))],0)
            total=off.copy();outgoing=off.sum(axis=0)
            for i in range(len(total)):total[i,i]=-outgoing[i]
            np.testing.assert_allclose(total,p['total_rates_s_inv'],rtol=128*np.finfo(float).eps,atol=0)
            net=np.einsum('ijd,jd->id',total,pop);gross=np.einsum('ijd,jd->id',abs(total),pop)
            residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
            if residual>=plan['controls']['rtol']:raise ValueError('Fresh independent rate residual exceeds original tolerance')
            rate_checks.append(dict(case=r['case']['name'],residual=residual,stored_residual=r['coupled']['fresh_physical_rate_residual']))
    arrays=[]
    for d,r in zip(args.runs,receipts):
        if fingerprint(d/'flux.csv')!=r['flux']:raise ValueError('Changed control flux')
        arrays.append(np.genfromtxt(d/'flux.csv',delimiter=',',names=True))
    wave=arrays[0]['wavelength_air_A']
    for array in arrays[1:]:np.testing.assert_array_equal(wave,array['wavelength_air_A'])
    previous=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json';prior=json.loads(previous.read_text())
    allowed,reason=gate_holding('solar_kpno_kurucz2005_corrected','kpno_solar_atlas')
    if not allowed:raise ValueError(reason)
    path=require_codex_path('data.kurucz2005_residual')
    if fingerprint(path)!=prior['observation']:raise ValueError('Changed Kitt observation')
    ow,obs=read_kurucz2005_residual(path,7763,7783)
    b=prior['broadening'];segments=np.array([(wave[0]/10,wave[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
    result=dict(schema='codex.bride_xi_review/1',plan=fingerprint(args.plan),independent_rate_checks=rate_checks,scientific_admission='HOLD',generator=fingerprint(__file__),
        controls=[fingerprint(d/'report.json')for d in args.runs],observation=fingerprint(path),policy=reason,
        broadening=b,ispec_effects=fingerprint(inspect.getsourcefile(ispec.apply_post_fundamental_effects)),cases=[])
    profiles=[wave,np.interp(wave,ow,obs)];labels=['wavelength_air_A','shipped_observation']
    args.output.mkdir(parents=True,exist_ok=False);fig,axes=plt.subplots(3,1,figsize=(12,9))
    for ax in axes:ax.plot(ow,obs,color='black',lw=.9,label='Kitt shipped continuum')
    for name,data in zip(names,arrays):
        case=dict(mode=name,intrinsic=[],variants=[])
        for c in (7771.944,7774.166,7775.388):case['intrinsic'].append(dict(centre_air_A=c,area_mA=1000*window_integral(wave,data['normalized_flux'],c-.5,c+.5)))
        for label,rv in [('canonical_broadened',0.),('plus_gravitational_shift',SOLAR_GRAVITATIONAL_REDSHIFT)]:
            f=ispec.apply_post_fundamental_effects(wave/10,data['normalized_flux'].copy(),segments,
                macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(rv,))
            profiles.append(f);labels.append(name+'_'+label)
            metrics=[]
            for c in (7771.944,7774.166,7775.388):metrics.append(dict(centre_air_A=c,model_area_mA=1000*window_integral(wave,f,c-.5,c+.5),observed_area_mA=1000*window_integral(ow,obs,c-.5,c+.5)))
            case['variants'].append(dict(name=label,rv_km_s=rv,metrics=metrics))
            if rv!=0:
                for ax in axes:ax.plot(wave,f,lw=.9,label=name.replace('_',' '))
        result['cases'].append(case)
    result['area_changes_vs_nominal_mA']={c['mode']:[v['area_mA']-base['area_mA']for v,base in zip(c['intrinsic'],result['cases'][baseline]['intrinsic'])]for c in result['cases']}
    result['max_profile_change_vs_nominal']={name:float(np.max(abs(array['normalized_flux']-arrays[baseline]['normalized_flux'])))for name,array in zip(names,arrays)}
    ranges=[(7770.5,7776.5), (7777.8,7781.2), (7766.9,7770.7)]
    for ax,limits in zip(axes,ranges):
        ax.set_xlim(*limits)
        inside=(wave>=limits[0])&(wave<=limits[1])
        floor=min(float(np.min(a[inside]))for a in profiles[1:])
        ax.set_ylim(max(0,floor-.02),1.005)
        ax.legend(fontsize=8);ax.grid(alpha=.2);ax.set_ylabel('Normalized flux')
    axes[-1].set_xlabel('AIR wavelength (Angstrom)')
    fig.suptitle('Solar microturbulence slice: fresh coupled O at fixed thermal structure\n41 source Fe I pressure lines, Cowley Stark; H2 omitted; no fit')
    fig.tight_layout();figure=args.output/'xi_control.png';fig.savefig(figure,dpi=150);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(8,8),sharex=True)
    xis=np.array([r['case']['xi']for r in receipts]);order=np.argsort(xis)
    for j,(ax,centre)in enumerate(zip(axes,(7771.944,7774.166,7775.388))):
        values=np.array([result['area_changes_vs_nominal_mA'][name][j]for name in names])
        ax.plot(xis[order],values[order],'o-');ax.axhline(0,color='black',lw=.5)
        ax.set_ylabel('Area change (mA)');ax.set_title(f'{centre:.3f} AIR A');ax.grid(alpha=.2)
    axes[-1].set_xlabel('Microturbulence (km/s)')
    fig.suptitle('First local xi slice: change from nominal1.0km/s; fresh coupled O\nFixed thermal structure and O8.70; no abundance fit or uncertainty bound')
    fig.tight_layout();response=args.output/'xi_response.png';fig.savefig(response,dpi=150);plt.close(fig)
    result['response_figure']=fingerprint(response)
    path=args.output/'profiles.csv';np.savetxt(path,np.column_stack(profiles),delimiter=',',header=','.join(labels),comments='')
    result.update(profiles=fingerprint(path),figure=fingerprint(figure),caveats='First local xi slice of broader program, not full sweep. Fresh coupled O, fixed solar thermal structure; no self-consistent atmosphere response, H2 EOS, abundance or full3D admission. Both RV variants retained, no optimization.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Areas vs nominal',result['area_changes_vs_nominal_mA']);print('Independent rates',rate_checks)


if __name__=='__main__':main()
