"""Compare natural versus source-width FeII, both references and fixed RV controls."""
import argparse,csv,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stark-mode',required=True,choices=('cowley','provided'));p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from scripts.audit_bride_kitt_residuals import integral
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    plan=json.loads((ROOT/'config/experiments/rya1281_feii_width_transfer.json').read_text())
    frozen=json.loads((ROOT/'config/experiments/rya1281_continuum_anchor.json').read_text())
    screen=json.loads((ROOT/'data/results/rya1281/lte_opacity_screen_v1/report.json').read_text())
    comparison_path=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json';comparison=json.loads(comparison_path.read_text());b=comparison['broadening']
    args.output.mkdir(parents=True,exist_ok=False)
    result=dict(schema='codex.bride_feii_width_analysis/1',scientific_admission='HOLD',stark_mode=args.stark_mode,leg_labels={'zero':'natural','feii':'W+'+args.stark_mode},generator=fingerprint(__file__),
        source_plan=fingerprint(ROOT/'config/experiments/rya1281_feii_width_transfer.json'),comparison=fingerprint(comparison_path),
        caveats='Fixed O populations; natural FeII versus W pressure plus explicit ionic Stark. Direct model response and observed continuum-reference response are separate. Historical3percent estimator retained; main15percent rule unexercised. No abundance or preferred reference/RV adopted.',cases=[])
    rows=[]
    for name in plan['cases']:
        origin=args.source/name;run=json.loads((origin/'report.json').read_text())
        if run['atom']!=next(c for c in screen['cases'] if c['name']==name)['atom']:raise ValueError('Changed frozen oxygen atom')
        case=next(c for c in frozen['cases'] if c['name']==name)
        if run['stark_mode']!=args.stark_mode:raise ValueError('Changed Stark alternative')
        raw={};audits={};broadened={}
        for mode in ('zero','feii'):
            leg=origin/mode;receipt=json.loads((leg/'report.json').read_text())
            if fingerprint(leg/'flux.csv')!=receipt['flux'] or receipt!=next(r for r in run['legs'] if r['mode']==mode):raise ValueError('Changed transfer output')
            raw[mode]=np.genfromtxt(leg/'flux.csv',delimiter=',',names=True);w=raw[mode]['wavelength_air_A']
            segment=np.array([(w[0]/10,w[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
            broadened[mode]={}
            for rvlabel,rv in [('canonical_broadened',0.),('plus_gravitational_shift',SOLAR_GRAVITATIONAL_REDSHIFT)]:
                broadened[mode][rvlabel]=ispec.apply_post_fundamental_effects(w/10,raw[mode]['normalized_flux'].copy(),segment,
                    macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(rv,))
            audits[mode]={}
            for reference in ('synthesis','iag'):
                destination=args.output/name/mode/reference
                command=[sys.executable,str(ROOT/'scripts/audit_bride_kitt_residuals.py'),'--comparison',str(comparison_path),
                    '--wide-spectrum',str(leg),'--continuum-reference',reference,'--output',str(destination)]
                completed=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,timeout=60)
                log=args.output/(name+'_'+mode+'_'+reference+'.log');log.write_text(completed.stdout+completed.stderr)
                if completed.returncode:raise RuntimeError('Continuum audit failed: '+str(log))
                audit=json.loads((destination/'report.json').read_text());audits[mode][reference]=audit
                if mode=='zero':
                    prior=json.loads((ROOT/'data/results/rya1281/feii_transfer_v1/analysis'/name/'feii'/reference/'report.json').read_text())
                    if [r['model_guided'] for r in audit['lines']]!=[r['model_guided'] for r in prior['lines']]:raise ValueError('Baseline continuum records differ')
        np.testing.assert_array_equal(raw['zero']['wavelength_air_A'],raw['feii']['wavelength_air_A'])
        np.testing.assert_array_equal(raw['zero']['continuum_Fnu_W_m2_Hz'],raw['feii']['continuum_Fnu_W_m2_Hz'])
        item=dict(name=name,source=fingerprint(origin/'report.json'),direct=[],audit_receipts={mode:{ref:fingerprint(args.output/name/mode/ref/'report.json') for ref in audits[mode]} for mode in audits})
        fig,axes=plt.subplots(3,1,figsize=(11,10))
        view=(w>7776.7)&(w<7777.6)
        for mode in ('zero','feii'):axes[0].plot(w[view],broadened[mode]['plus_gravitational_shift'][view],label=('natural' if mode=='zero' else 'W+'+args.stark_mode))
        axes[0].set(title='Natural versus W+source Stark FeII neighborhood; fixed O, canonical broadening + gravitational shift',xlabel='AIR A',ylabel='Normalized model flux');axes[0].legend()
        for li,centre in enumerate((7771.944,7774.166,7775.388)):
            intrinsic={mode:1000*integral(w,1-raw[mode]['normalized_flux'],centre-.5,centre+.5) for mode in raw}
            direct=dict(centre_air_A=centre,intrinsic_area_delta_mA=intrinsic['feii']-intrinsic['zero'],
                max_intrinsic_O_window_flux_change=float(np.max(abs(raw['feii']['normalized_flux'][abs(w-centre)<=.5]-raw['zero']['normalized_flux'][abs(w-centre)<=.5]))),rv_variants={})
            for rv in broadened['zero']:
                pair={mode:1000*integral(w,1-broadened[mode][rv],centre-.5,centre+.5) for mode in raw}
                direct['rv_variants'][rv]=dict(zero_area_mA=pair['zero'],feii_area_mA=pair['feii'],delta_mA=pair['feii']-pair['zero'])
            item['direct'].append(direct)
            for reference in ('synthesis','iag'):
                pair={}
                for mode in raw:
                    record=audits[mode][reference]['lines'][li]
                    guided=[v for v in record['variants'] if v['name']=='model_guided_'+reference+'_diagnostic']
                    pair[mode]=dict(applied=record['model_guided']['applied'],level=record['model_guided']['level_at_centre'],
                        area_mA=guided[0]['observed_area_mA'] if guided else None)
                delta=pair['feii']['area_mA']-pair['zero']['area_mA'] if all(p['area_mA']is not None for p in pair.values()) else None
                rows.append(dict(case=name,centre_air_A=centre,reference=reference,zero_applied=pair['zero']['applied'],feii_applied=pair['feii']['applied'],
                    zero_level=pair['zero']['level'],feii_level=pair['feii']['level'],zero_observed_area_mA=pair['zero']['area_mA'],
                    feii_observed_area_mA=pair['feii']['area_mA'],observed_area_delta_mA=delta,
                    direct_model_area_delta_mA=direct['rv_variants']['plus_gravitational_shift']['delta_mA']))
        x=np.arange(3)
        for rv in broadened['zero']:axes[1].plot(x,[d['rv_variants'][rv]['delta_mA'] for d in item['direct']],'o-',label=rv)
        axes[1].set(title='Direct model O-window response',ylabel='Area change (mA)');axes[1].legend(fontsize=8)
        for reference in ('synthesis','iag'):
            rr=[r for r in rows if r['case']==name and r['reference']==reference]
            axes[2].plot(x,[r['observed_area_delta_mA'] for r in rr],'o-',label=reference)
        axes[2].set(title='Observed-area response from each continuum-reference leg',ylabel='Area change (mA)');axes[2].legend()
        for ax in axes[1:]:ax.set_xticks(x,['7771.944','7774.166','7775.388']);ax.set_xlabel('O centre, AIR A');ax.axhline(0,color='grey',lw=.6)
        for ax in axes:ax.grid(alpha=.2)
        fig.suptitle(name+': fixed-O natural versus source-width FeII differential, HOLD');fig.tight_layout();fig.savefig(args.output/(name+'_response.png'),dpi=140);plt.close(fig)
        result['cases'].append(item)
    with (args.output/'summary.csv').open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows))


if __name__=='__main__':main()
