"""Compare immutable ORCH continuum vintages on qualified O columns, without new RT."""
import argparse,csv,importlib.util,json,os,subprocess,sys,types
from dataclasses import asdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def load_git_module(commit,name):
    """Execute the exact committed standalone estimator; retain its blob as evidence."""
    source=subprocess.check_output(['git','show',commit+':pipeline/local_continuum.py'],cwd=ROOT)
    module=types.ModuleType(name);module.__file__=str(ROOT/'pipeline/local_continuum.py')
    sys.modules[name]=module;exec(compile(source,module.__file__,'exec'),module.__dict__)
    return module,source


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--plan',type=Path,default=ROOT/'config/experiments/rya1281_orch_continuum_parity.json');a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.prenormalised_guard import assert_not_renormalising
    from pipeline.telluric_policy import gate_holding
    from pipeline.bride_observables import window_integral
    from scripts.analyze_bride_H_kitt_orch import validate_processed_flux
    from pipeline.bride_continuum_audit import selection,trace_fit
    plan=json.loads(a.plan.read_text());checks=[]
    def check(path,receipt):
        if fingerprint(path)!=receipt:raise ValueError('Changed source '+str(path))
        checks.append((path,receipt))
    check(a.plan,fingerprint(a.plan));check(Path(__file__),fingerprint(__file__))
    os.environ['CODEX_CONT_REF']='iag'
    if any(os.environ.get(k) for k in ['CODEX_CONT_Q','CODEX_CONT_PCT']):raise ValueError('Unexpected continuum overrides')
    modules={};sources={}
    for label,commit in plan['orch_commits'].items():
        modules[label],sources[label]=load_git_module(commit,'bride_continuum_'+label)
        if modules[label].MAX_SHIFT!=.03 or modules[label].MAX_SHIFT_MODEL!=.15 or modules[label].MODEL_Q!=90:raise ValueError('Changed estimator policy')
    current=modules['current']
    if(current.NEIGH_A,current.NEIGH_STEP_A,current.MIN_NEIGH,current.SMOOTH_TOL)!=(50.,10.,3,.02):raise ValueError('Changed neighbourhood policy')
    parent_path=ROOT/plan['parent_report'];parent=json.loads(parent_path.read_text());check(parent_path,plan['parent_receipt'])
    observation=require_codex_path('data.kurucz2005_residual');check(observation,parent['observation'])
    iag=require_codex_path('data.solar_iag_baker2020');check(iag,parent['iag'])
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec,inspect
    check(Path(inspect.getfile(ispec.apply_post_fundamental_effects)),parent['ispec_effects'])
    w,obs=read_kurucz2005_residual(observation,min(plan['centres_air_A'])-56,max(plan['centres_air_A'])+56)
    sys.path.insert(0,str(ROOT/'scripts'));import measure_band_ew as mbe
    iw,iff=mbe.iag_atlas()
    if mbe.IAG_FITS.resolve()!=iag.resolve():raise ValueError('Wrong IAG holding')
    observations={'kitt_shipped':(w,obs),'iag_native':(iw,iff)}
    holdings={'kitt_shipped':('solar_kpno_kurucz2005_corrected','kpno_solar_atlas'),'iag_native':('solar_iag','iag_fts_solar_atlas')}
    # Both source holdings remain pre-normalised; only the declared local operation is permitted.
    for label,(holding,instrument) in holdings.items():
        allowed,reason=gate_holding(holding,instrument)
        if not allowed:raise ValueError(reason)
        assert_not_renormalising(holding,pre_normalised=True,local_window_envelope=True)
    models={};rows=[];arrays={};source_receipts=[];b=parent['broadening']
    for case in plan['cases']:
        folder=ROOT/plan['source']/case;r=json.loads((folder/'report.json').read_text());check(folder/'report.json',plan['source_reports'][case]);check(folder/'flux.csv',r['flux'])
        if r['case']!=case or r.get('dry_run') or r['independent']['physical_rate_residual']>=r['controls']['rtol']:raise ValueError('Unqualified source solve')
        raw=np.genfromtxt(folder/'flux.csv',delimiter=',',names=True);mw=raw['wavelength_air_A'];seg=np.array([(mw[0]/10,mw[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        source_receipts.append(dict(case=case,report=fingerprint(folder/'report.json'),flux=r['flux']))
        for rv in plan['RVs']:
            model=ispec.apply_post_fundamental_effects(mw/10,raw['normalized_flux'].copy(),seg,macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(parent['RV_km_s'][rv],))
            validate_processed_flux(mw,model,plan['centres_air_A']);models[case,rv]=model
            beta=parent['RV_km_s'][rv]*1000/299792458.
            for ci,rest in enumerate(plan['centres_air_A']):
                centre=float(rest)if plan['continuum_centre']=='rest_air'else float(rest*np.sqrt((1+beta)/(1-beta)))
                for label,(ow,of) in observations.items():
                    cut=abs(ow-centre)<=current.ENV_HALF_WIDTH_A;cw,cf=ow[cut],of[cut]
                    if cw[0]>centre-4.99 or cw[-1]<centre+4.99 or np.any(cf<=0)or not np.isfinite(cf).all():raise ValueError('Missing/nonpositive observation support')
                    for reference in plan['references']:
                        if reference=='synthesis':rw,rf=mw,model
                        else:
                            atlas=current.solar_atlas_reference(centre-56,centre+56,b['R'])
                            if atlas is None:raise ValueError('Missing admitted IAG reference')
                            rw,rf=atlas
                        if rw[0]>centre-5 or rw[-1]<centre+5:raise ValueError('Incomplete local reference')
                        callback=current.neighbourhood_check(ow,of,centre,rw,rf,exclude_half_width_A=plan['half_width_A'])
                        old,oldcont=modules['previous'].fit_model_guided(cw,cf,centre,rw,rf,exclude_half_width_A=plan['half_width_A'],reference=reference)
                        new,newcont=current.fit_model_guided(cw,cf,centre,rw,rf,exclude_half_width_A=plan['half_width_A'],reference=reference,corroborate=callback)
                        if not np.isfinite(new.level_at_centre):raise ValueError('Unconstrained actual O window')
                        if abs(new.level_at_centre-1)<=current.MAX_SHIFT:
                            if asdict(old)!=asdict(new):raise ValueError('Unchanged-policy record differs')
                            np.testing.assert_array_equal(oldcont,newcont)
                        ref_on_pixels=np.interp(cw,rw,rf);selected,_=selection(cw,cf,centre,ref_on_pixels);trace,_=trace_fit(cw,cf,centre,ref_on_pixels,selected)
                        if abs(trace['level']-new.level_at_centre)>512*np.finfo(float).eps*trace['geometry']['design_condition_2']:raise ValueError('Independent fit trace differs')
                        geometry=trace['geometry']
                        prefix=f'{case}_{rv}_{ci}_{label}_{reference}';arrays[prefix+'_wave_air_A']=cw;arrays[prefix+'_shipped_flux']=cf;arrays[prefix+'_model']=np.interp(cw,mw,model)
                        corrected={}
                        for version,rec,cont in [('previous',old,oldcont),('current',new,newcont)]:
                            corrected[version]=cf.copy()if cont is None else cf/cont
                            if not np.isfinite(corrected[version]).all()or np.any(corrected[version]<=0):raise ValueError('Invalid corrected spectrum')
                            arrays[prefix+'_'+version]=corrected[version]
                        support=abs(cw-centre)<=plan['half_width_A'];mf=np.interp(cw,mw,model)
                        area=lambda f:1000*window_integral(cw,f,centre-plan['half_width_A'],centre+plan['half_width_A'])
                        ok,note=callback(new.level_at_centre)
                        rows.append(dict(case=case,rv=rv,rest_air_A=rest,continuum_centre_air_A=centre,observation=label,reference=reference,centre_outside_anchor_hull_A=geometry['centre_outside_hull_A'],anchor_blue_count=geometry['blue_count'],anchor_red_count=geometry['red_count'],design_condition_2=geometry['design_condition_2'],centre_weight_l1=geometry['centre_weight_l1'],coefficient_support_count=trace['coefficient_support_count'],clipping_limit_reached=trace['clipping_limit_reached'],level=new.level_at_centre,slope_per_A=new.slope_per_A,reported_n_bins=new.n_bins,previous_applied=old.applied,current_applied=new.applied,current_reason=new.reason,neighbourhood_required=abs(new.level_at_centre-1)>.03,neighbourhood_ok=ok,neighbourhood_note=note,reference_full_neighbourhood=bool(rw[0]<=centre-55 and rw[-1]>=centre+55),max_corrected_flux_change=float(np.max(abs(corrected['current']-corrected['previous']))),shipped_area_residual_mA=area(mf)-area(cf),previous_area_residual_mA=area(mf)-area(corrected['previous']),current_area_residual_mA=area(mf)-area(corrected['current']),current_core_rms=float(np.sqrt(np.mean((mf[support]-corrected['current'][support])**2)))))
    if len(rows)!=len(plan['cases'])*len(plan['RVs'])*len(plan['centres_air_A'])*len(observations)*len(plan['references']):raise ValueError('Incomplete comparison')
    for path,receipt in checks:
        if fingerprint(path)!=receipt:raise ValueError('Source changed during replay')
    a.output.mkdir(parents=True,exist_ok=False)
    for label,source in sources.items():(a.output/('local_continuum_'+label+'.py')).write_bytes(source)
    with(a.output/'comparison.csv').open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    np.savez_compressed(a.output/'profiles.npz',**arrays)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for oi,label in enumerate(observations):
        for ci,rest in enumerate(plan['centres_air_A']):
            ax=axes[oi,ci];prefix=f'solar_plus_gravity_{ci}_{label}_synthesis';x=arrays[prefix+'_wave_air_A']-rest
            ax.plot(x,arrays[prefix+'_model'],label='solar model');ax.plot(x,arrays[prefix+'_shipped_flux'],label='shipped observation')
            for reference in plan['references']:
                pre=f'solar_plus_gravity_{ci}_{label}_{reference}'
                for version in modules:ax.plot(x,arrays[pre+'_'+version],label=version+' / '+reference,alpha=.7)
            ax.set_xlim(-.9,.9);ax.set_title(f'{label} {rest} AIR A');ax.set_xlabel('AIR offset (A)');ax.legend(fontsize=7)
    fig.suptitle('Pinned ORCH continuum policies; fixed source abundance, no inversion');fig.tight_layout();fig.savefig(a.output/'continuum_parity.png',dpi=140);plt.close(fig)
    report=dict(schema='codex.bride_orch_continuum_parity/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(a.plan),orch_commits=plan['orch_commits'],continuum_centre=plan['continuum_centre'],source_receipts=source_receipts,protected_inputs=[dict(path=str(p),receipt=r)for p,r in checks],rows=len(rows),required_neighbourhood=sum(r['neighbourhood_required']for r in rows),changed_rows=sum(r['max_corrected_flux_change']>0 for r in rows),level_range=[min(r['level']for r in rows),max(r['level']for r in rows)],max_corrected_flux_change=max(r['max_corrected_flux_change']for r in rows),limitations=['Actual estimator replay, not complete ORCH synthesis or abundance inversion.','Synthesis covers only7763-7783AIR A; no invented neighbourhood opacity. Missing wide synthesis would refuse a beyond3percent correction.','No extra telluric operation, observation shift or source renormalisation.','Native IAG and Kitt holding identities retained; IAG reference leg uses established resolving-power Gaussian.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['rows','required_neighbourhood','changed_rows','level_range','max_corrected_flux_change']}))

if __name__=='__main__':main()
