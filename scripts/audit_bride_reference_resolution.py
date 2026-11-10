"""Frozen O777 reference-resolution sensitivity; no FTS calibration or adoption."""
import argparse,csv,json,os,sys
from dataclasses import asdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_resolution_audit import broaden_reference
    from pipeline.bride_continuum_audit import selection,trace_fit,two_factor_changes
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.telluric_policy import gate_holding
    from pipeline.prenormalised_guard import assert_not_renormalising
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    from scripts.audit_bride_kitt_residuals import integral
    os.environ['CODEX_CONT_REF']='iag'
    from pipeline import local_continuum as lc
    from scripts.measure_band_ew import iag_atlas
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    plan_path=ROOT/'config/experiments/rya1281_resolution_audit.json';plan=json.loads(plan_path.read_text())
    previous_path=ROOT/'data/results/rya1281/continuum_support_v1/report.json';previous=json.loads(previous_path.read_text())
    frozen=json.loads((ROOT/'config/experiments/rya1281_continuum_anchor.json').read_text())
    if plan['schema']!='codex.bride_resolution_plan/1' or previous['native']!=fingerprint(lc.__file__):raise ValueError('Changed rule or plan')
    if (lc.MODEL_Q,lc.DEGREE,lc.CLIP_MAD_MODEL,lc.MIN_MODEL_PIX,lc.MAX_SHIFT,lc.ENV_HALF_WIDTH_A)!=(90.,1,3.,20,.03,5.):raise ValueError('Changed nominal continuum settings')
    b=previous['broadening']
    if b['R']!=plan['target_R']:raise ValueError('Changed target R')
    policies={}
    for holding,instrument in [('solar_kpno_kurucz2005_corrected','kpno_solar_atlas'),('solar_iag','iag_fts_solar_atlas')]:
        allowed,policy=gate_holding(holding,instrument)
        if not allowed:raise ValueError(policy)
        policies[holding]=policy
    assert_not_renormalising('solar_kpno_kurucz2005_corrected',pre_normalised=True,local_window_envelope=True)
    observation=require_codex_path('data.kurucz2005_residual');atlas=require_codex_path('data.solar_iag_baker2020')
    if fingerprint(observation)!=previous['observation'] or fingerprint(atlas)!=previous['reference_atlas']:raise ValueError('Changed observation/reference')
    w,obs=read_kurucz2005_residual(observation,7765,7782);iw,ifl=iag_atlas()
    result=dict(schema='codex.bride_reference_resolution/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),adapter=fingerprint(ROOT/'pipeline/bride_resolution_audit.py'),
        continuum_adapter=fingerprint(ROOT/'pipeline/bride_continuum_audit.py'),native=fingerprint(lc.__file__),
        plan=fingerprint(plan_path),previous=fingerprint(previous_path),observation=fingerprint(observation),
        reference_atlas=fingerprint(atlas),policies=policies,broadening=b,
        caveats='Frozen spectra, original AIR holding/internal RT vacuum. Gaussian approximations only; no measured FTS kernel, fitted R, preferred reference, abundance or physical-error inference. No second telluric correction. Both earlier RV residual legs remain unchanged.',cases=[])
    rows=[];figures={}
    for case,oldcase in zip(frozen['cases'],previous['cases']):
        if case['name']!=oldcase['name']:raise ValueError('Changed case order')
        d=require_codex_path(case['resource'])/case['relative'];run=json.loads((d/'report.json').read_text())
        if fingerprint(d/'report.json')!=oldcase['source_report'] or fingerprint(d/'flux.csv')!=run['flux']:raise ValueError('Changed frozen spectrum')
        raw=np.genfromtxt(d/'flux.csv',delimiter=',',names=True)
        segments=np.array([(raw['wavelength_air_A'][0]/10,raw['wavelength_air_A'][-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        synth=ispec.apply_post_fundamental_effects(raw['wavelength_air_A']/10,raw['normalized_flux'].copy(),segments,
            macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(SOLAR_GRAVITATIONAL_REDSHIFT,))
        record=dict(name=case['name'],physical_case=run['case'],source_report=fingerprint(d/'report.json'),lines=[])
        fig,axes=plt.subplots(3,1,figsize=(11,10))
        for li,oldline in enumerate(oldcase['lines']):
            centre=oldline['centre_air_A'];k=(iw>=centre-6)&(iw<=centre+6);nw,nf=iw[k],ifl[k]
            native=lc.solar_atlas_reference(centre-5,centre+5,b['R'])
            if native is None:raise ValueError('No IAG reference; no fallback')
            np.testing.assert_array_equal(nw,native[0])
            keep=(nw>=centre-5.01)&(nw<=centre+5.01)
            references={};masks={};kernels={}
            for variant in plan['variants']:
                f,kernel=broaden_reference(nw,nf,centre,b['R'],native_R=variant['native_R'],no_extra=variant['no_extra'])
                if variant['name']=='shipped_gaussian':np.testing.assert_array_equal(f,native[1])
                name=variant['name'];references[name]=np.interp(w,nw[keep],f[keep]);masks[name],_=selection(w,obs,centre,references[name]);kernels[name]=kernel
            base='shipped_gaussian';ref0=references[base];mask0=masks[base]
            synref=np.interp(w,raw['wavelength_air_A'],synth);synmask,_=selection(w,obs,centre,synref)
            syntrace,syncont=trace_fit(w,obs,centre,synref,synmask)
            syntrace['observed_area_mA']=1000*integral(w,1-obs/syncont,centre-.5,centre+.5)
            if syntrace!=oldline['corners']['synthesis_synthesis']:raise ValueError('Changed synthesis baseline')
            baseline,_=trace_fit(w,obs,centre,ref0,mask0)
            baseline['observed_area_mA']=1000*integral(w,1-obs/np.polyval(baseline['passes'][-1]['coefficients'],w-centre),centre-.5,centre+.5)
            if baseline!=oldline['corners']['iag_iag']:raise ValueError('Changed IAG baseline')
            line=dict(centre_air_A=centre,synthesis_baseline=syntrace,iag_baseline=baseline,variants=[])
            for variant in plan['variants']:
                name=variant['name'];ref=references[name];mask=masks[name];corners={}
                for s,sm in [('baseline',mask0),('variant',mask)]:
                    for t,denom in [('baseline',ref0),('variant',ref)]:
                        trace,cont=trace_fit(w,obs,centre,denom,sm)
                        if cont is not None:trace['observed_area_mA']=1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                        corners[s+'_'+t]=trace
                guided=corners['variant_variant']
                # Direct native fit uses original native samples, avoiding a second interpolation.
                vf,_=broaden_reference(nw,nf,centre,b['R'],native_R=variant['native_R'],no_extra=variant['no_extra'])
                rec,native_cont=lc.fit_model_guided(w,obs,centre,nw[keep],vf[keep],exclude_half_width_A=.5,reference='iag')
                if (guided['level'],guided['slope_per_A'],guided['reported_n_bins'],guided['applied'])!=(rec.level_at_centre,rec.slope_per_A,rec.n_bins,rec.applied):raise ValueError('Native variant trace mismatch')
                if native_cont is not None:np.testing.assert_array_equal(native_cont,np.polyval(guided['passes'][-1]['coefficients'],w-centre))
                fixed={}
                for label,indices in [('iag',baseline['fit_indices']),('synthesis',syntrace['fit_indices'])]:
                    fm=np.zeros(len(w),bool);fm[indices]=True;fixed[label]={}
                    for t,denom in [('baseline',ref0),('variant',ref)]:
                        tr,cont=trace_fit(w,obs,centre,denom,fm,clip=False)
                        if cont is not None:tr['observed_area_mA']=1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                        fixed[label][t]=tr
                mapped={a:corners[z] for a,z in zip(('synthesis_synthesis','synthesis_iag','iag_synthesis','iag_iag'),('baseline_baseline','baseline_variant','variant_baseline','variant_variant'))}
                accounting=two_factor_changes(mapped,'observed_area_mA') if all('observed_area_mA' in t for t in corners.values()) else None
                item=dict(name=name,kernel=kernels[name],native_record=asdict(rec),corners=corners,frozen_support=fixed,
                    area_accounting_mA=accounting,factor_labels='selector and denominator: baseline shipped IAG versus variant IAG; synthesis labels in utility are only positional',
                    selection_changed_pixels=int(np.count_nonzero(mask!=mask0)),
                    coefficient_support_changed_pixels=len(set(guided['fit_indices'])^set(baseline['fit_indices'])),
                    max_reference_flux_change=float(np.max(abs(ref[abs(w-centre)<=5]-ref0[abs(w-centre)<=5]))))
                line['variants'].append(item)
                row=dict(case=case['name'],centre_air_A=centre,variant=name,applied=guided['applied'],
                    selected_pixels_changed=item['selection_changed_pixels'],fit_pixels_changed=item['coefficient_support_changed_pixels'],
                    guided_area_delta_mA=None if accounting is None else accounting['total'])
                for label in fixed:
                    pair=fixed[label];row[label+'_fixed_delta_mA']=(pair['variant']['observed_area_mA']-pair['baseline']['observed_area_mA']) if all('observed_area_mA' in t for t in pair.values()) else None
                rows.append(row)
            rr=rows[-len(plan['variants']):];x=np.arange(len(rr));ax=axes[li]
            for j,key in enumerate(('guided_area_delta_mA','iag_fixed_delta_mA','synthesis_fixed_delta_mA')):
                values=[r[key] for r in rr];ax.plot(x,values,'o-',label=key)
            ax.axhline(0,color='grey',lw=.8);ax.set_xticks(x,[r['variant'] for r in rr],fontsize=8)
            ax.set(title=f'{centre:.3f} AIR A; fixed spectra, change from shipped IAG',ylabel='Aperture area change (mA)');ax.grid(alpha=.2);ax.legend(fontsize=8)
            record['lines'].append(line)
        fig.suptitle(case['name']+': Gaussian sensitivity only, scientific HOLD');fig.tight_layout();figures[case['name']]=fig;result['cases'].append(record)
    args.output.mkdir(parents=True,exist_ok=False)
    for name,fig in figures.items():
        path=args.output/(name+'_resolution.png');fig.savefig(path,dpi=140);plt.close(fig)
    with (args.output/'summary.csv').open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows,allow_nan=False))


if __name__=='__main__':main()
