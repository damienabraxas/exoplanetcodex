"""Compare fresh fixed-population alpha controls; never select adopted physics."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline','default','joint','audit','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();checks={str(Path(__file__)):fingerprint(__file__)}
    for directory in [a.baseline,a.default,a.joint,a.audit]:
        for name in ['report.json','rates.npz'] if directory!=a.audit else ['report.json']:
            path=directory/name;checks[str(path)]=fingerprint(path)
    reports=[json.loads((d/'report.json').read_text())for d in [a.default,a.joint]]
    if not all(r['computational_completed']for r in reports):raise ValueError('Incomplete replay')
    audit=json.loads((a.audit/'report.json').read_text())
    if not audit['absolute_replay']['within_FP64_bound']:raise ValueError('Unqualified absolute replay')
    control=reports[1]['context_audit']['continuum_interpolation_control']
    if control['source_tables_reproduced_at_knots']!=22 or len(control['compiled_alpha_checks'])!=44:raise ValueError('Incomplete native alpha qualification')
    with np.load(a.baseline/'rates.npz')as old,np.load(a.default/'rates.npz')as x,np.load(a.joint/'rates.npz')as y:
        if set(old.files)!=set(x.files):raise ValueError('Default packet inventory changed')
        for key in old.files:np.testing.assert_array_equal(old[key],x[key])
        held=['populations','nstar','wave_vac_nm','native_R_s_inv','native_C_s_inv','replay_C_s_inv']
        for key in held:np.testing.assert_array_equal(x[key],y[key])
        n=x['populations'];C=x['replay_C_s_inv'];R=[x['replay_R_s_inv'],y['replay_R_s_inv']]
        def balance(matrix):
            off=matrix+C;incoming=np.sum(off*n[None,:,:],axis=1);outgoing=np.sum(off,axis=0)*n
            return (incoming-outgoing)/np.maximum(incoming+outgoing,np.finfo(float).tiny)
        b=[balance(r)for r in R];idx=np.unravel_index(np.argmax(abs(b[1])),b[1].shape)
        positive=R[0]>0
        delta_J=y['Jnu']-x['Jnu'];nonzero=x['Jnu']>0
        result=dict(generator=fingerprint(__file__),inputs=checks,scientific_admission='HOLD',computational_completed=True,
            default_reproduces_prior_arrays=sorted(old.files),held_arrays_exact=held,wavelength_points=len(x['wave_vac_nm']),source_continua=22,compiled_alpha_checks=44,
            default_residual=float(abs(b[0]).max()),joint_residual=float(abs(b[1]).max()),
            worst_joint_level_1based=int(idx[0]+1),worst_joint_depth_top_to_bottom=int(idx[1]),
            level41_top_signed=[float(v[40,0])for v in b],
            max_rate_relative_change=float(np.max(abs(R[1][positive]/R[0][positive]-1))),
            max_J_absolute_change=float(abs(delta_J).max()),max_positive_J_relative_change=float(np.max(abs(delta_J[nonzero]/x['Jnu'][nonzero]))),
            J_changed=not np.array_equal(x['Jnu'],y['Jnu']),absolute_replay=audit['absolute_replay'],
            scope='Raw-source native frequency-Hermite alpha shared by active opacity/emissivity and rates; fresh radiation at identical fixed populations. Native FP32 alpha output and exact canonical source-knot coordinates are part of this control. Held LW thermal profile, BB, collisions, source grid and background provider unchanged. No SE, full thermal consistency, adopted model, solar abundance or observed UV/VIS/IR qualification.')
    a.output.mkdir(parents=True,exist_ok=False)
    (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
