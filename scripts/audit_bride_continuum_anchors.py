"""Decompose fixed-spectrum continuum selection/reference effects; no adoption."""
import argparse,csv,json,os,sys
from dataclasses import asdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_continuum_audit import selection,trace_fit,two_factor_changes
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.telluric_policy import gate_holding
    from pipeline.prenormalised_guard import assert_not_renormalising
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    # Enables only the existing atlas reader; direct estimator reference labels
    # below remain explicit and both denominators are always retained.
    os.environ['CODEX_CONT_REF']='iag'
    from pipeline import local_continuum as lc
    if (lc.MODEL_Q,lc.DEGREE,lc.CLIP_MAD_MODEL,lc.MIN_MODEL_PIX,lc.MAX_SHIFT,lc.ENV_HALF_WIDTH_A)!=(90.,1,3.,20,.03,5.):raise ValueError('Require unchanged nominal continuum settings')
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    from scripts.audit_bride_kitt_residuals import integral
    plan_path=ROOT/'config/experiments/rya1281_continuum_anchor.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_continuum_anchor_plan/1' or plan['reference_modes']!=['synthesis','iag']:raise ValueError('Changed declared reference accounting')
    comparison_path=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json';comparison=json.loads(comparison_path.read_text())
    holding='solar_kpno_kurucz2005_corrected';allowed,policy=gate_holding(holding,'kpno_solar_atlas')
    if not allowed:raise ValueError(policy)
    assert_not_renormalising(holding,pre_normalised=True,local_window_envelope=True)
    allowed,reference_policy=gate_holding('solar_iag','iag_fts_solar_atlas')
    if not allowed:raise ValueError(reference_policy)
    observation=require_codex_path('data.kurucz2005_residual')
    if fingerprint(observation)!=comparison['observation']:raise ValueError('Changed Kitt observation')
    w,obs=read_kurucz2005_residual(observation,7765,7782)
    if not np.isfinite(obs).all() or np.any(obs<=0):raise ValueError('Invalid observed flux in fixed audit segment')
    source=next(c for c in comparison['cases'] if c['name']=='solar5772_bound_Fe_sampled')
    b=comparison['broadening']
    result=dict(schema='codex.bride_continuum_anchor_review/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        adapter=fingerprint(ROOT/'pipeline/bride_continuum_audit.py'),native=fingerprint(lc.__file__),plan=fingerprint(plan_path),
        comparison=fingerprint(comparison_path),observation=fingerprint(observation),holding=holding,policy=policy,
        reference_policy=reference_policy,reference_atlas=fingerprint(require_codex_path('data.solar_iag_baker2020')),
        broadening=b,reference_rv_km_s=SOLAR_GRAVITATIONAL_REDSHIFT,cases=[],
        caveats='Frozen spectra; native nominal rule unchanged. Crossed/frozen/side/postclip fits are diagnostics, not adopted estimators or uncertainty terms. No raw external flux arrays archived; no abundance/RV/LSF/telluric fit or preferred reference.')
    payloads={};pixel_rows=[]
    for case in plan['cases']:
        d=require_codex_path(case['resource'])/case['relative'];run=json.loads((d/'report.json').read_text())
        if (run['source_manifest']!=source['source_manifest'] or run['source_flux']!=source['flux']
                or run['flux']!=fingerprint(d/'flux.csv')):raise ValueError('Changed fixed model identity')
        raw=np.genfromtxt(d/'flux.csv',delimiter=',',names=True)
        segments=np.array([(raw['wavelength_air_A'][0]/10,raw['wavelength_air_A'][-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        flux=ispec.apply_post_fundamental_effects(raw['wavelength_air_A']/10,raw['normalized_flux'].copy(),segments,
            macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(SOLAR_GRAVITATIONAL_REDSHIFT,))
        prior={m:json.loads((ROOT/case['prior_'+m]).read_text()) for m in plan['reference_modes']}
        for m,r in prior.items():
            if r['wide_spectrum']!=run or r['continuum_code']!=fingerprint(lc.__file__) or r['observation']!=fingerprint(observation):raise ValueError('Prior diagnostic/source identity differs')
        record=dict(name=case['name'],physical_case=run['case'],source_report=fingerprint(d/'report.json'),prior_reports={m:fingerprint(ROOT/case['prior_'+m]) for m in prior},lines=[])
        fig,axes=plt.subplots(3,1,figsize=(12,11))
        for i,centre in enumerate((7771.944,7774.166,7775.388)):
            references={};masks={};native={};corners={};contours={}
            for m in plan['reference_modes']:
                if m=='synthesis':mw,mf=raw['wavelength_air_A'],flux
                else:
                    ref=lc.solar_atlas_reference(centre-5,centre+5,b['R'])
                    if ref is None:raise ValueError('Missing IAG reference, no synthesis fallback')
                    mw,mf=ref
                keep=(mw>=centre-5.01)&(mw<=centre+5.01)
                if mw[keep].min()>centre-5 or mw[keep].max()<centre+5:raise ValueError('Truncated reference neighborhood')
                references[m]=np.interp(w,mw[keep],mf[keep]);masks[m],_=selection(w,obs,centre,references[m])
                rec,c=lc.fit_model_guided(w,obs,centre,mw[keep],mf[keep],exclude_half_width_A=.5,reference=m)
                if asdict(rec)!=prior[m]['lines'][i]['model_guided']:raise ValueError('Native replay differs from stored continuum record')
                native[m]=(rec,c)
            for selector in plan['reference_modes']:
                for denominator in plan['reference_modes']:
                    key=selector+'_'+denominator
                    trace,cont=trace_fit(w,obs,centre,references[denominator],masks[selector]);corners[key]=trace;contours[key]=cont
                    if cont is not None:trace['observed_area_mA']=1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                    if selector==denominator:
                        rec,c=native[selector]
                        if (trace['level']!=rec.level_at_centre or trace['slope_per_A']!=rec.slope_per_A or trace['reported_n_bins']!=rec.n_bins or trace['applied']!=rec.applied):raise ValueError('Exact trace differs from native coefficients')
                        if c is not None:np.testing.assert_array_equal(c,cont)
            line=dict(centre_air_A=centre,corners=corners,level_accounting=two_factor_changes(corners,'level'),
                slope_accounting=two_factor_changes(corners,'slope_per_A'),frozen_support={},side_controls={},postclip_refit={})
            if all('observed_area_mA' in t for t in corners.values()):line['area_accounting_mA']=two_factor_changes(corners,'observed_area_mA')
            else:line['area_accounting_withheld']='At least one crossed correction refused; no replacement area'
            eligible=abs(w-centre)<=5
            for m in plan['reference_modes']:
                nominal=corners[m+'_'+m];used=np.array(nominal['fit_indices']);fixed=np.zeros(len(w),bool);fixed[used]=True
                line['frozen_support'][m]={}
                for denominator in plan['reference_modes']:
                    trace,cont=trace_fit(w,obs,centre,references[denominator],fixed,clip=False)
                    if cont is not None:trace['observed_area_mA']=1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                    line['frozen_support'][m][denominator]=trace
                for side,side_mask in (('blue',w<centre),('red',w>centre)):
                    trace,cont=trace_fit(w,obs,centre,references[m],masks[m]&side_mask)
                    if cont is not None:trace['observed_area_mA']=1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                    line['side_controls'][m+'_'+side]=trace
                if nominal['postfit_survivors_differ']:
                    mask=np.zeros(len(w),bool);mask[nominal['surviving_indices']]=True
                    trace,cont=trace_fit(w,obs,centre,references[m],mask,clip=False)
                    if cont is not None:trace['observed_area_mA']=1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                    line['postclip_refit'][m]=trace
                ratio=obs/references[m];ax=axes[i]
                ax.scatter(w[used]-centre,ratio[used],s=10,label=m+' actual coefficient support')
                coefficients=nominal['passes'][-1]['coefficients']
                ax.plot(w[eligible]-centre,np.polyval(coefficients,w[eligible]-centre),label=m+' returned continuum')
                weights=dict(zip(nominal['fit_indices'],nominal['geometry']['centre_weights']))
                for index in np.flatnonzero(eligible):
                    pixel_rows.append(dict(case=case['name'],centre_air_A=centre,reference=m,offset_air_A=float(w[index]-centre),
                        initially_selected=bool(masks[m][index]),coefficient_support=bool(fixed[index]),centre_weight=weights.get(int(index),0.),
                        relative_ratio=float(ratio[index]),fit_residual=float(ratio[index]-np.polyval(coefficients,w[index]-centre))))
            ax=axes[i];ax.axvspan(-.5,.5,color='grey',alpha=.12,label='excluded O aperture');ax.axvline(0,color='grey',ls=':')
            ax.set(xlim=(-5,5),ylabel='Observed/reference ratio',title=f'{centre:.3f} AIR A: fixed spectrum, native support and extrapolation')
            ax.grid(alpha=.2);ax.legend(fontsize=8)
            record['lines'].append(line)
        axes[-1].set_xlabel('Offset from O centre, AIR Angstrom');fig.suptitle(case['name']+': nominal native rule, diagnostics only; scientific HOLD');fig.tight_layout()
        payloads[case['name']]=fig;result['cases'].append(record)
    args.output.mkdir(parents=True,exist_ok=False)
    for name,fig in payloads.items():
        path=args.output/(name+'_anchors.png');fig.savefig(path,dpi=140);plt.close(fig)
        next(c for c in result['cases'] if c['name']==name)['figure']=fingerprint(path)
    path=args.output/'pixel_diagnostics.csv'
    with path.open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(pixel_rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(pixel_rows)
    result['derived_pixels']=fingerprint(path)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps([dict(case=c['name'],lines=[dict(centre=l['centre_air_A'],level=l['level_accounting'],area=l.get('area_accounting_mA'),postclip=list(l['postclip_refit'])) for l in c['lines']]) for c in result['cases']]))


if __name__=='__main__':main()
