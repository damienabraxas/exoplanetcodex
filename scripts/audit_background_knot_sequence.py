"""Audit the independent source-background-knot numerical branch."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--refined-native-state',action='store_true');a=p.parse_args()
    groups=[('hop-physical-shared16-v1',None),('background-knots-v1',None),('background-midpoint-v1','background-adaptive-plan-v1'),('background-adaptive-v1','background-final-audit-v1')]
    if a.refined_native_state:groups=[('hop-physical-grid-refined-v1',None),('hop-grid-refined-background-knots-v1',None),('hop-grid-refined-midpoint-v1','grid-refined-midpoint-audit-v1'),('hop-grid-refined-adaptive-v1','grid-refined-final-audit-v1')]
    checks={str(Path(__file__)):fingerprint(__file__)};rows=[];old=None;old_parab=None
    atom=read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'));mask=np.zeros((41,41),bool)
    for t in atom.continua:mask[t.i-1,t.j-1]=mask[t.j-1,t.i-1]=True
    for name,audit in groups:
        directory=a.raw_root/('rya1371-'+name)
        for file in ['report.json','rates.npz']:checks[str(directory/file)]=fingerprint(directory/file)
        doc=json.loads((directory/'report.json').read_text())
        if not doc['computational_completed']:raise ValueError('Incomplete background-grid case')
        with np.load(directory/'rates.npz')as x:d={k:x[k].copy()for k in x.files}
        row=dict(case=name,points=len(d['wave_vac_nm']),max_physical_residual=doc['matched_population_replay_physical_residual'],seconds=doc['seconds'],radiation_iterations=doc['radiation']['iterations'],max_inherited_J_relative=None,max_BF_coefficient_relative_step=None,parabolic_BF_max_residual=None,max_parabolic_BF_coefficient_relative_step=None)
        q=d['replay_R_s_inv']+d['replay_C_s_inv'];n=d['populations'];inc=np.sum(q*n[None,:,:],axis=1);out=np.sum(q,axis=0)*n
        residual=abs(inc-out)/(inc+out)
        if float(np.max(residual))!=row['max_physical_residual']:raise ValueError('Physical residual mismatch')
        loc=np.unravel_index(np.argmax(residual),residual.shape);row['worst_level_1based']=int(loc[0]+1);row['worst_depth_top_to_bottom']=int(loc[1])
        if old is not None:
            for key in ['populations','nstar','replay_C_s_inv']:np.testing.assert_array_equal(d[key],old[key])
            np.testing.assert_array_equal(d['replay_R_s_inv']>0,old['replay_R_s_inv']>0)
            idx=np.searchsorted(d['wave_vac_nm'],old['wave_vac_nm']);np.testing.assert_array_equal(d['wave_vac_nm'][idx],old['wave_vac_nm'])
            row['max_inherited_J_relative']=float(np.max(abs(d['Jnu'][idx]/old['Jnu']-1)))
            good=mask[:,:,None]&(old['replay_R_s_inv']>0)
            row['max_BF_coefficient_relative_step']=float(np.max(abs(d['replay_R_s_inv'][good]/old['replay_R_s_inv'][good]-1)))
        if audit:
            path=a.raw_root/('rya1371-'+audit)
            for file in ['report.json','parabolic_rates.npz']:checks[str(path/file)]=fingerprint(path/file)
            report=json.loads((path/'report.json').read_text())
            if not report['computational_completed']or report['inputs'].get(str(directory/'rates.npz'))!=fingerprint(directory/'rates.npz'):raise ValueError('Unmatched parabolic audit')
            row['parabolic_BF_max_residual']=report['adaptive_continuum']['independent_parabolic_BF_residual']
            with np.load(path/'parabolic_rates.npz')as x:parab=x['R_s_inv'].copy()
            if old_parab is not None:
                good=mask[:,:,None]&(old_parab>0);row['max_parabolic_BF_coefficient_relative_step']=float(np.max(abs(parab[good]/old_parab[good]-1)))
            old_parab=parab
        rows.append(row);old=d
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed background sequence input')
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'sequence.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    scope='Independent source-feature branch retains original O nodes and all native opacity knots. Initial inherited-J difference is1.46e-12 with radiation iterations8->9; subsequent grids retain inherited J exactly. No physical data change or certified numerical convergence. The residual at native reference populations is not an abundance error; full coupled/solar/flux/Fe qualification remains owed.'
    if a.refined_native_state:scope='Same source-feature quadrature branch at exact PRE-update59003-node native-reference populations. Original own nodes/physical graph and all source opacity knots retained, n/nStar/C identical across own cases. Inherited J differences measured in each row, not assumed zero. No SE update, error bound, numerical/physical/model/solar/Fe admission. Native refinement and own-grid refinement are distinct numerical axes.'
    report=dict(generator=fingerprint(__file__),inputs=checks,cases=rows,computational_completed=True,scientific_admission='HOLD',scope=scope)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(rows))


if __name__=='__main__':main()
