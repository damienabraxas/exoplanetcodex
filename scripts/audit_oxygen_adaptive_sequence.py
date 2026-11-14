"""Independent state/field/rate accounting for bounded continuum adaptation."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    atom_path=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_path);model=to_lightweaver(atom)
    checks={str(Path(__file__)):fingerprint(__file__),str(atom_path):fingerprint(atom_path)};bf=np.zeros((41,41),bool)
    for t in atom.continua:bf[t.j-1,t.i-1]=bf[t.i-1,t.j-1]=True
    groups=[('shared16-midpoint-v1',4),('adaptive-continuum1-v1',5),('adaptive-continuum2-v1',6),('adaptive-continuum3-v1',7)]
    rows=[];previous=None;old_alpha=None
    for group,audit_index in groups:
        path=a.raw_root/('rya1371-'+group);audit=a.raw_root/f'rya1371-adaptive-continuum-audit-v{audit_index}'
        for directory,files in [(path,['report.json','rates.npz']),(audit,['report.json','parabolic_rates.npz'])]:
            for name in files:checks[str(directory/name)]=fingerprint(directory/name)
        doc=json.loads((path/'report.json').read_text());ad=json.loads((audit/'report.json').read_text())
        if not doc['computational_completed']or not ad['computational_completed']:raise ValueError('Incomplete adaptive case')
        with np.load(path/'rates.npz')as x:d={k:x[k].copy()for k in x.files}
        with np.load(audit/'parabolic_rates.npz')as x:parabolic=x['R_s_inv'].copy()
        wave=d['wave_vac_nm'];R=d['replay_R_s_inv'];n=d['populations'];C=d['replay_C_s_inv']
        if len(wave)>MAX_PACKET_WAVELENGTHS:raise ValueError('Adaptive packet exceeds unchanged size cap')
        def balance(matrix):
            q=matrix+C;inc=np.sum(q*n[None,:,:],axis=1);out=q.sum(axis=0)*n;r=abs(inc-out)/(inc+out);loc=np.unravel_index(np.argmax(r),r.shape)
            return float(np.max(r)),int(loc[0]+1),int(loc[1])
        maximum,level,depth=balance(R);higher,_,_=balance(parabolic)
        if maximum!=doc['matched_population_replay_physical_residual']:raise ValueError('Adaptive producer residual mismatch')
        row=dict(case=group,points=len(wave),trapezoid_max_residual=maximum,parabolic_BF_max_residual=higher,worst_level_1based=level,worst_depth_top_to_bottom=depth,seconds=doc['seconds'],previous_J_bit_identical=None,max_BF_coefficient_relative_step=None,max_parabolic_BF_coefficient_relative_step=None)
        if previous is not None:
            old,old_parabolic=previous
            for key in ['populations','nstar','replay_C_s_inv']:np.testing.assert_array_equal(d[key],old[key])
            idx=np.searchsorted(wave,old['wave_vac_nm']);np.testing.assert_array_equal(wave[idx],old['wave_vac_nm']);np.testing.assert_array_equal(d['Jnu'][idx],old['Jnu'])
            row['previous_J_bit_identical']=True
            mask=bf[:,:,None]&(old['replay_R_s_inv']>0)
            row['max_BF_coefficient_relative_step']=float(np.max(abs(R[mask]/old['replay_R_s_inv'][mask]-1)))
            row['max_parabolic_BF_coefficient_relative_step']=float(np.max(abs(parabolic[mask]/old_parabolic[mask]-1)))
        alpha=[]
        for t in model.continua:
            local=wave[(wave>=t.minLambda)&(wave<=t.lambdaEdge)];alpha.append((local,t.alpha(local)))
        if old_alpha is not None:
            for (ow,oa),(w,values)in zip(old_alpha,alpha):
                idx=np.searchsorted(w,ow);np.testing.assert_array_equal(w[idx],ow);np.testing.assert_array_equal(values[idx],oa)
        old_alpha=alpha;previous=(d,parabolic);rows.append(row)
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'sequence.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    x=[r['points']for r in rows]
    for key,label in [('trapezoid_max_residual','Compiled trapezoid'),('parabolic_BF_max_residual','Independent parabolic BF')]:axes[0].plot(x,[r[key]for r in rows],'o-',label=label)
    axes[0].set_ylabel('Maximum physical rate-balance residual');axes[0].set_yscale('log');axes[0].legend(fontsize=8)
    for key,label in [('max_BF_coefficient_relative_step','Trapezoid BF'),('max_parabolic_BF_coefficient_relative_step','Parabolic BF')]:axes[1].plot(x[1:],[r[key]for r in rows[1:]],'o-',label=label)
    axes[1].set_ylabel('Maximum BF rate change from preceding grid');axes[1].set_yscale('log');axes[1].legend(fontsize=8)
    for ax in axes:ax.set_xlabel('All integration wavelengths retained');ax.axvline(MAX_PACKET_WAVELENGTHS,color='gray',ls=':',label='Existing cap')
    fig.suptitle('Bounded adaptive continuum checks — numerical convergence remains OPEN');fig.savefig(a.output/'adaptive_sequence.png',dpi=160);plt.close(fig)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed adaptive audit input')
    report=dict(generator=fingerprint(__file__),inputs=checks,cases=rows,computational_completed=True,scientific_admission='HOLD',original_node_cap=MAX_PACKET_WAVELENGTHS,inherited_alpha_bit_identical=True,scope='Every prior wavelength/population/LTE/collision state retained. J and continuum cross sections at inherited nodes are bit-identical. All22continua/all56depths contribute error indicators; source transitions retained. Absolute integral audits validate each sampled-field compiled rate packet before estimating errors. Neither nonmonotonic trapezoid nor same-node positive parabolic BF sequence establishes convergence. No model/flux/observed/all-band/3D admission. Final19989-node RT case is distinct from the unexecuted8-point proposal in audit-v7.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(rows))


if __name__=='__main__':main()
