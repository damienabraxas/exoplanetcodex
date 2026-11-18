"""Read-only joint LTE/inverse diagnostic checks; not reference-fit selection."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.thermal_profile_contract import collision_summary
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['prior','baseline','joint','coarse','baseline-audit','joint-audit','coarse-audit','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();checks={};reports={};arrays={};audit_reports={}
    def pin(path):checks[str(path.resolve())]=fingerprint(path)
    for key in ['prior','baseline','joint','coarse']:
        directory=getattr(a,key)
        for name in ['report.json','rates.npz']:pin(directory/name)
        reports[key]=json.loads((directory/'report.json').read_text())
        if not reports[key]['computational_completed']:raise ValueError('Incomplete replay')
        with np.load(directory/'rates.npz')as d:arrays[key]={k:d[k].copy()for k in d.files}
    if set(arrays['prior'])!=set(arrays['baseline']):raise ValueError('Baseline inventory differs')
    for key in arrays['prior']:np.testing.assert_array_equal(arrays['prior'][key],arrays['baseline'][key])
    for key in ['populations','wave_vac_nm','native_R_s_inv','native_C_s_inv']:np.testing.assert_array_equal(arrays['baseline'][key],arrays['joint'][key])
    np.testing.assert_allclose(arrays['joint']['nstar'].sum(axis=0),arrays['baseline']['nstar'].sum(axis=0),rtol=1e-12,atol=0)
    np.testing.assert_array_equal(arrays['joint']['replay_C_s_inv']>0,arrays['baseline']['replay_C_s_inv']>0)
    channel=reports['joint']['joint_lte_control']['channels']
    if {k:v['records']for k,v in channel.items()}!={'CE':741,'CI':39,'CH':302,'CH0':24}or not all(v['source_forward_exact']for v in channel.values()):raise ValueError('Incomplete source-direction preservation')
    packets={};coefficients={}
    for key in ['baseline','joint','coarse']:
        directory=getattr(a,key)
        for name in ['emission_profile_arrays.npz','stationary_coefficients.npz']:pin(directory/name)
        with np.load(directory/'emission_profile_arrays.npz')as d:packets[key]={k:d[k].copy()for k in d.files}
        with np.load(directory/'stationary_coefficients.npz')as d:coefficients[key]={k:d[k].copy()for k in d.files}
        audit_path=getattr(a,key+'_audit')/'report.json';pin(audit_path);doc=json.loads(audit_path.read_text())
        if not doc['computational_completed']or not doc['absolute_replay']['within_FP64_bound']or doc['coefficient_replay']['counts']!={'BB':256,'BF':22}:raise ValueError('Failed/incomplete independent audit')
        if doc['inputs'].get(str((directory/'rates.npz').resolve()))!=fingerprint(directory/'rates.npz')or doc['inputs'].get(str((directory/'stationary_coefficients.npz').resolve()))!=fingerprint(directory/'stationary_coefficients.npz'):raise ValueError('Audit has wrong saved state/coefficient packet')
        audit_reports[key]=doc
    for key in ['background_chi','background_eta','background_sca']:np.testing.assert_array_equal(packets['baseline'][key],packets['joint'][key])
    def relative(new,old):
        if not np.isfinite(new).all()or not np.isfinite(old).all():raise ValueError('Nonfinite comparison')
        mask=old>0
        if np.any(new[~mask]!=old[~mask]):raise ValueError('Changed zero support')
        return float(np.max(abs(new[mask]/old[mask]-1)))
    metrics=dict(max_relative_nstar_change=relative(arrays['joint']['nstar'],arrays['baseline']['nstar']),max_relative_J_change=relative(arrays['joint']['Jnu'],arrays['baseline']['Jnu']),max_relative_R_change=relative(arrays['joint']['replay_R_s_inv'],arrays['baseline']['replay_R_s_inv']),max_relative_C_change=relative(arrays['joint']['replay_C_s_inv'],arrays['baseline']['replay_C_s_inv']),max_relative_total_chi_change=relative(coefficients['joint']['total_chi_m_inv'],coefficients['baseline']['total_chi_m_inv']),max_relative_total_thermal_eta_change=relative(coefficients['joint']['total_thermal_eta'],coefficients['baseline']['total_thermal_eta']),collision_equilibrium_at_declared_star=collision_summary(arrays['joint']['replay_C_s_inv'],arrays['joint']['nstar']))
    fine=arrays['joint'];coarse=arrays['coarse']
    for key in ['populations','nstar','native_R_s_inv','native_C_s_inv','replay_C_s_inv']:np.testing.assert_array_equal(fine[key],coarse[key])
    idx=np.searchsorted(fine['wave_vac_nm'],coarse['wave_vac_nm'])
    np.testing.assert_array_equal(fine['wave_vac_nm'][idx],coarse['wave_vac_nm'])
    numerical=dict(coarse_nodes=len(coarse['wave_vac_nm']),fine_nodes=len(fine['wave_vac_nm']),inherited_J_exact=bool(np.array_equal(fine['Jnu'][idx],coarse['Jnu'])),max_relative_inherited_J_change=relative(fine['Jnu'][idx],coarse['Jnu']),max_relative_rate_change=relative(fine['replay_R_s_inv'],coarse['replay_R_s_inv']),scope='Matched-state inherited17263->19963nodes; both retained, numerical response only, not a convergence certificate or mesh selected by native closeness.')
    conservation={}
    for key in ['baseline','joint','coarse']:
        values=arrays[key];rates=values['replay_R_s_inv']+values['replay_C_s_inv'];n=values['populations']
        incoming=np.sum(rates*n[None,:,:],axis=1);outgoing=np.sum(rates,axis=0)*n
        if np.any(incoming+outgoing<=0):raise ValueError('Missing whole-graph rate coverage')
        local=abs(incoming-outgoing)/(incoming+outgoing);worst=np.unravel_index(np.argmax(local),local.shape)
        error=float(np.max(abs(np.sum(incoming-outgoing,axis=0))/np.sum(incoming+outgoing,axis=0)))
        bound=256*np.finfo(float).eps*len(n)
        if error>bound:raise ValueError('Whole-graph particle balance arithmetic failure')
        conservation[key]=dict(max_local_relative_balance=float(local[worst]),worst_level_1based=int(worst[0]+1),worst_depth_top_to_bottom=int(worst[1]),max_relative_global_net_sum=error,arithmetic_screen=bound,levels=len(n),depths=n.shape[1])
    report=dict(generator=fingerprint(__file__),inputs=checks,scientific_admission='HOLD',computational_completed=True,baseline_exact_arrays=list(arrays['baseline']),fixed_actual_population_and_reference_arrays=['populations','wave_vac_nm','native_R_s_inv','native_C_s_inv'],background_arrays_exact=['chi','eta','sca'],source_channels=channel,metrics=metrics,numerical=numerical,global_particle_balance=conservation,source_snapshot_timing='After-run component snapshots checked against actual producing input receipts; primary report and script snapshots retained.',cases={key:dict(seconds=reports[key]['seconds'],radiation=reports[key]['radiation'],physical_residual=reports[key]['matched_population_replay_physical_residual'],independent_rates=audit_reports[key]['absolute_replay'],independent_coefficients=audit_reports[key]['coefficient_replay'],integrated_thermal_canaries=audit_reports[key]['max_absolute_integrated_defects'])for key in ['baseline','joint','coarse']},scope='Coordinated finite-level LTE/source-channel inverse fixed-population diagnostic; compiled opacity/emissivity/rates jointly consume LTE. Both outcomes retained without residual selection. Absolute EOS/partition, coupled/preconditioner/restart/numerical/observed-solar admission open.')
    a.output.mkdir(parents=True,exist_ok=False)
    # After-run source snapshots are verified against each actual producing report.
    for key in ['baseline','joint','coarse']:
        for path,expected in reports[key]['inputs'].items():
            source=Path(path)
            if source.is_relative_to(ROOT)and source.suffix=='.py':
                if fingerprint(source)!=expected:raise ValueError('Producing source changed before snapshot')
                target=a.output/(source.stem+'_source.py')
                if target.exists()and target.read_bytes()!=source.read_bytes():raise ValueError('Snapshot collision')
                target.write_bytes(source.read_bytes())
    for doc in audit_reports.values():
        for name,expected in doc['inputs'].items():
            source=Path(name)
            if source.is_relative_to(ROOT)and source.suffix=='.py':
                if fingerprint(source)!=expected:raise ValueError('Audit source changed before snapshot')
                target=a.output/(source.stem+'_source.py')
                if target.exists()and target.read_bytes()!=source.read_bytes():raise ValueError('Audit snapshot collision')
                target.write_bytes(source.read_bytes())
    for path in [Path(__file__),ROOT/'scripts/audit_oxygen_integrated_thermal_rates.py',ROOT/'pipeline/stationary_oxygen_coefficients.py',ROOT/'scripts/archive_oxygen_reference_checks.py']:(a.output/(path.stem+'_source.py')).write_bytes(path.read_bytes())
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(metrics=metrics,numerical=numerical,global_particle_balance=conservation,source_snapshot_timing='After-run component snapshots checked against actual producing input receipts; primary report and script snapshots retained.',cases=report['cases']),indent=2))


if __name__=='__main__':main()
