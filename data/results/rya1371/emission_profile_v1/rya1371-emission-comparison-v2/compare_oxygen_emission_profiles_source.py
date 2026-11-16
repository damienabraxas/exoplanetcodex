"""Read-only matched-state emission/boundary comparison, not model admission."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['prior','default','lw','native','lw-audit','native-audit','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();checks={};reports={};arrays={}
    for name in ['prior','default','lw','native']:
        directory=getattr(a,name)
        for filename in ['report.json','rates.npz']:checks[str(directory/filename)]=fingerprint(directory/filename)
        reports[name]=json.loads((directory/'report.json').read_text())
        if not reports[name]['computational_completed']:raise ValueError('Incomplete source run')
        with np.load(directory/'rates.npz')as d:arrays[name]={k:d[k].copy()for k in d.files}
    if set(arrays['prior'])!=set(arrays['default']):raise ValueError('Default inventory differs')
    for k,v in arrays['prior'].items():np.testing.assert_array_equal(v,arrays['default'][k])
    for name in ['lw','native']:
        for k in ['native_R_s_inv','native_C_s_inv','replay_C_s_inv','populations','nstar','wave_vac_nm']:np.testing.assert_array_equal(arrays['default'][k],arrays[name][k])
        expected='held-LW'if name=='lw'else'source-native'
        if reports[name]['context_audit']['emission_profile_control']['profile']!=expected:raise ValueError('Wrong declared profile')
    packets={}
    for name in ['lw','native']:
        path=getattr(a,name)/'emission_profile_arrays.npz';checks[str(path)]=fingerprint(path)
        with np.load(path)as d:packets[name]={k:d[k].copy()for k in d.files}
    for k in ['background_chi','background_sca','bottom_total_chi']:np.testing.assert_array_equal(packets['lw'][k],packets['native'][k])
    def relative(x,y):
        if not np.isfinite(x).all()or not np.isfinite(y).all():raise ValueError('Nonfinite comparison')
        positive=y>0
        if np.any((~positive)&(x!=y)):raise ValueError('New positive support in comparison')
        return float(np.max(abs(x[positive]/y[positive]-1)))
    result={}
    for name in ['default','lw','native']:
        report=reports[name];result[name]=dict(seconds=report['seconds'],radiation=report['radiation'],physical_residual=report['matched_population_replay_physical_residual'])
        if name!='default':
            result[name].update(max_relative_J_change=relative(arrays[name]['Jnu'],arrays['default']['Jnu']),max_relative_positive_rate_change=relative(arrays[name]['replay_R_s_inv'],arrays['default']['replay_R_s_inv']))
            audit_dir=getattr(a,name+'_audit');path=audit_dir/'report.json';checks[str(path)]=fingerprint(path)
            doc=json.loads(path.read_text())
            # Keep actual independent report metrics; no residual-based choice.
            result[name]['independent_audit']=doc.get('absolute_replay')
            if not doc['computational_completed']or not doc['absolute_replay']['within_FP64_bound']:raise ValueError('Independent audit incomplete or failed')
            if doc['inputs'].get(str((getattr(a,name)/'rates.npz').resolve()))!=fingerprint(getattr(a,name)/'rates.npz'):raise ValueError('Independent audit has wrong replay')
    report=dict(generator=fingerprint(__file__),inputs=checks,scientific_admission='HOLD',computational_completed=True,default_exact_arrays=list(arrays['default']),profile_identity_arrays=['native_R_s_inv','native_C_s_inv','replay_C_s_inv','populations','nstar','wave_vac_nm','background_chi','background_sca','bottom_total_chi'],cases=result,between_profiles=dict(max_relative_background_eta=relative(packets['lw']['background_eta'],packets['native']['background_eta']),max_relative_bottom_intensity=relative(packets['lw']['bottom_intensity'],packets['native']['bottom_intensity'])),scope='Background thermal emission and lower diffusion profile axis only. Active BB/BF/LTE/collision conventions unchanged; no preferred profile, coupled closure, observational or scientific admission.')
    a.output.mkdir(parents=True,exist_ok=False)
    snapshots={}
    for path,expected in reports['lw']['inputs'].items():
        source=Path(path)
        if source.is_relative_to(ROOT)and source.suffix=='.py':
            if fingerprint(source)!=expected:raise ValueError('Producing source changed before after-run snapshot')
            target=a.output/(source.stem+'_source.py');target.write_bytes(source.read_bytes())
            snapshots[path]=dict(producing_receipt=expected,snapshot=fingerprint(target),timing='after-run, producing input receipt verified')
    for source in [Path(__file__),ROOT/'scripts/archive_oxygen_reference_checks.py',ROOT/'scripts/audit_oxygen_integrated_thermal_rates.py']:
        target=a.output/(source.stem+'_source.py')
        if target.exists()and target.read_bytes()!=source.read_bytes():raise ValueError('Source snapshot name collision')
        target.write_bytes(source.read_bytes())
    report['source_snapshots']=snapshots
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(report['cases'],indent=2))


if __name__=='__main__':main()
