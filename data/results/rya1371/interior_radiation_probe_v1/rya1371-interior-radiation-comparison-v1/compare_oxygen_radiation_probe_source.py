"""Compare actual interior J to parent interpolation; sampled-panel shadow only."""
import argparse,csv,json,os,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.saved_field_quadrature import interpolate_saved_field
    from pipeline.radiation_probe_grid import gaussian_probe_grid,METHOD
    from pipeline.native_continuum_control import load_native_continuum_control
    from pipeline.multi_atom import read_multi_atom
    from config.constants import require_codex_path
    import lightweaver.constants as c
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['plan','grid','baseline','probe','audit','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();start=time.perf_counter();checks={};docs={};arrays={}
    if os.environ.get('OMP_NUM_THREADS')!='1'or os.environ.get('OPENBLAS_NUM_THREADS')!='1':raise ValueError('Require single-thread bounded probe comparison')
    def pin(path,expected=None):
        path=Path(path).resolve();fp=fingerprint(path)
        if expected is not None and fp!=expected:raise ValueError('Changed probe input '+str(path))
        checks[str(path)]=fp;return path
    plan=json.loads(pin(a.plan).read_text());spec=plan['wavelength_grid_control']
    if spec['method']!=METHOD:raise ValueError('Wrong probe plan')
    grid_doc=json.loads(pin(a.grid/'report.json',spec['producer']).read_text());packet=pin(a.grid/'probe_grid.npz',spec['packet'])
    with np.load(packet)as d:grid={k:d[k].copy()for k in d.files}
    expanded,query,weights,panels=gaussian_probe_grid(grid['inherited_wave_vac_nm'],grid['panel_scores'],32)
    for actual,expected in [(grid['wave_vac_nm'],expanded),(grid['probe_wave_vac_nm'],query),(grid['probe_weights_nm'],weights),(grid['selected_panels'],panels)]:np.testing.assert_array_equal(actual,expected)
    for key,directory in [('baseline',a.baseline),('probe',a.probe),('parent',Path(grid_doc['fine_rates']).parent)]:
        docs[key]=json.loads(pin(directory/'report.json').read_text());pin(directory/'rates.npz')
        if not docs[key]['computational_completed']:raise ValueError('Incomplete radiation source')
        with np.load(directory/'rates.npz')as d:arrays[key]={k:d[k].copy()for k in d.files}
    if set(arrays['parent'])!=set(arrays['baseline']):raise ValueError('Baseline inventory changed')
    for key in arrays['parent']:np.testing.assert_array_equal(arrays['parent'][key],arrays['baseline'][key])
    base=arrays['baseline'];probe=arrays['probe'];np.testing.assert_array_equal(base['wave_vac_nm'],grid['inherited_wave_vac_nm']);np.testing.assert_array_equal(probe['wave_vac_nm'],expanded)
    for key in ['populations','nstar','replay_C_s_inv','native_R_s_inv','native_C_s_inv']:np.testing.assert_array_equal(base[key],probe[key])
    for path in [a.grid/'report.json',packet]:
        if docs['probe']['inputs'].get(str(path.resolve()))!=fingerprint(path):raise ValueError('Probe did not consume declared producer')
    inherited=np.searchsorted(expanded,base['wave_vac_nm']);indices=np.searchsorted(expanded,query);np.testing.assert_array_equal(expanded[indices],query)
    predicted=interpolate_saved_field(base['wave_vac_nm'],base['Jnu'],query);actual=probe['Jnu'][indices]
    if np.any(predicted<=0)or not np.isfinite(actual).all()or np.any(actual<=0):raise ValueError('Invalid positive probe field')
    relative=actual/predicted-1;worst=np.unravel_index(np.argmax(abs(relative)),relative.shape)
    audit=json.loads(pin(a.audit/'report.json').read_text())
    if not audit['computational_completed']or not audit['absolute_replay']['within_FP64_bound']or audit['coefficient_replay']['counts']!={'BB':256,'BF':22}or audit['inputs'].get(str((a.probe/'rates.npz').resolve()))!=fingerprint(a.probe/'rates.npz'):raise ValueError('Wrong/failed independent probe audit')
    packets={};coeffs={}
    for key,directory in [('baseline',a.baseline),('probe',a.probe),('parent',Path(grid_doc['fine_rates']).parent)]:
        for filename in ['emission_profile_arrays.npz','stationary_coefficients.npz']:pin(directory/filename)
        with np.load(directory/'emission_profile_arrays.npz')as d:packets[key]={k:d[k].copy()for k in d.files}
        with np.load(directory/'stationary_coefficients.npz')as d:coeffs[key]={k:d[k].copy()for k in d.files}
    for key in coeffs['parent']:np.testing.assert_array_equal(coeffs['parent'][key],coeffs['baseline'][key])
    for key in packets['parent']:np.testing.assert_array_equal(packets['parent'][key],packets['baseline'][key])
    def maximum(new,old):
        mask=old>0
        if not np.isfinite(new).all()or not np.isfinite(old).all()or np.any(new[~mask]!=old[~mask]):raise ValueError('Invalid changed positive support')
        return float(np.max(abs(new[mask]/old[mask]-1)))
    inherited_metrics=dict(J_exact=bool(np.array_equal(probe['Jnu'][inherited],base['Jnu'])),max_relative_J=maximum(probe['Jnu'][inherited],base['Jnu']),max_relative_total_chi=maximum(coeffs['probe']['total_chi_m_inv'][inherited],coeffs['baseline']['total_chi_m_inv']),max_relative_total_thermal_eta=maximum(coeffs['probe']['total_thermal_eta'][inherited],coeffs['baseline']['total_thermal_eta']),background_arrays_exact={key:bool(np.array_equal(packets['probe'][key][inherited],packets['baseline'][key]))for key in ['background_chi','background_eta','background_sca']})
    native=require_codex_path(plan['native_analysis_resource']);pin(native/'arrays.npz',plan['native_arrays'])
    with np.load(native/'arrays.npz')as d:T=d['consumed_state_cgs'][::-1,2].astype(float)
    atom_path=pin(require_codex_path('experiments.bride_oxygen_atom'),plan['continuum_interpolation_control']['atom']);atom=read_multi_atom(atom_path)
    pin(c.__file__,docs['probe']['inputs'][str(Path(c.__file__))]);pin(__file__)
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'compare_oxygen_radiation_probe_source.py').write_bytes(Path(__file__).read_bytes())
    control,more=load_native_continuum_control(plan['continuum_interpolation_control'],atom_path,atom,resolver=require_codex_path,output=a.output/'native_alpha_queries',deadline=start+600);checks.update(more)
    delta=np.zeros_like(base['replay_R_s_inv'])
    for t in atom.continua:
        pair=(t.i-1,t.j-1);alpha=control.evaluate(pair,query);frequency=c.CLight/(query[:,None]*1e-9);beta=np.exp(-c.HPlanck*frequency/(c.KBoltzmann*T[None,:]));pref=(4*np.pi/c.HPlanck*weights/query*alpha)[:,None];difference=actual-predicted
        delta[pair[1],pair[0]]+=np.sum(pref*difference,axis=0);delta[pair[0],pair[1]]+=base['nstar'][pair[0]]/base['nstar'][pair[1]]*np.sum(pref*beta*difference,axis=0)
    shadow=base['replay_R_s_inv']+delta
    if np.any(shadow<0)or not np.array_equal(shadow>0,base['replay_R_s_inv']>0):raise ValueError('Invalid sampled-panel shadow rate graph')
    def residual(R):
        rates=R+base['replay_C_s_inv'];n=base['populations'];inc=np.sum(rates*n[None,:,:],axis=1);out=np.sum(rates,axis=0)*n;local=abs(inc-out)/(inc+out);worst=np.unravel_index(np.argmax(local),local.shape)
        return dict(max_relative_local=float(local[worst]),level_1based=int(worst[0]+1),depth_top_to_bottom=int(worst[1]))
    mask=base['replay_R_s_inv']>0
    rows=[]
    for point in range(len(query)):
        for depth in range(len(T)):rows.append(dict(probe_index=point,selected_panel_rank=point//4,parent_interval_index=int(panels[point//4]),wave_vac_nm=float(query[point]),depth_top_to_bottom=depth,temperature_K=T[depth],linear_Jnu=predicted[point,depth],actual_Jnu=actual[point,depth],actual_over_linear_minus_one=relative[point,depth]))
    report=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=True,scientific_admission='HOLD',baseline_exact_arrays=list(base),nodes=dict(parent=len(base['wave_vac_nm']),probe=len(expanded),new=len(query),remaining_original_cap=20000-len(expanded)),probe_field=dict(max_relative_interpolation_error=float(np.max(abs(relative))),rms_relative_interpolation_error=float(np.sqrt(np.mean(relative**2))),worst_probe_index=int(worst[0]),worst_wave_vac_nm=float(query[worst[0]]),worst_depth_top_to_bottom=int(worst[1]),signed_relative_at_worst=float(relative[worst]),panel_max_relative=[float(np.max(abs(relative[i:i+4])))for i in range(0,len(query),4)]),inherited=inherited_metrics,compiled_rate_refinement_max_relative=maximum(probe['replay_R_s_inv'],base['replay_R_s_inv']),cases=dict(baseline=dict(seconds=docs['baseline']['seconds'],radiation=docs['baseline']['radiation'],residual=residual(base['replay_R_s_inv'])),probe=dict(seconds=docs['probe']['seconds'],radiation=docs['probe']['radiation'],residual=residual(probe['replay_R_s_inv']))),independent_audit=dict(absolute_rates=audit['absolute_replay'],coefficients=audit['coefficient_replay']),sampled_panel_field_response=dict(max_abs_delta_over_original_positive_rate=float(np.max(abs(delta[mask]/base['replay_R_s_inv'][mask]))),whole_graph_shadow_residual=residual(shadow),scope='Gauss4weighted(actualJ-linearJ) on ONLY8selected parent panels,22BF/56depths; spontaneous term cancels, BB/C untouched. Not a global error bound, converged integral, replacement solver operator or unique attribution of compiled-grid changes.'),scope='Actual32interior samples versus parent interpolation under identical physical state; all19963ancestors retained within20000cap. Selected-source panels are not globally representative or independent observed UV validation. No physical/numerical model admission.',seconds=time.perf_counter()-start)
    control.verify_inputs()
    for path,expected in checks.items():
        if fingerprint(path)!=expected:raise ValueError('Input changed during probe comparison')
    # Verify/snapshot actual producing components after-run; no private tables.
    for doc in [docs['probe'],audit,grid_doc]:
        for name,expected in doc['inputs'].items():
            source=Path(name)
            if source.is_relative_to(ROOT)and source.suffix=='.py':
                if fingerprint(source)!=expected:raise ValueError('Producing source changed before snapshot')
                target=a.output/(source.stem+'_source.py')
                if target.exists()and target.read_bytes()!=source.read_bytes():raise ValueError('Component snapshot collision')
                target.write_bytes(source.read_bytes())
    with(a.output/'new_radiation_samples.csv').open('w',newline='')as out:writer=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    np.savez_compressed(a.output/'sampled_field_response.npz',probe_wave_vac_nm=query,actual_Jnu=actual,linear_Jnu=predicted,delta_R_s_inv=delta,shadow_R_s_inv=shadow)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='inputs'},indent=2))


if __name__=='__main__':main()
