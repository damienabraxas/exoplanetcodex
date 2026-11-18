"""Compare completed fixed-population field/strength/quadrature controls."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    groups=['hop-physical-shared16-v1','frozen-feedback-baseline-v1','bb-thermal-feedback-v1','bf-thermal-feedback-v1','thermal-feedback-v2','native-strength-feedback-v1','thermal-strength-feedback-v1','shared16-midpoint-v1','shared16-local1-v1','shared16-local2-v1']
    checks={str(Path(__file__)):fingerprint(__file__)};rows=[];base=None
    for group in groups:
        directory=a.raw_root/('rya1371-'+group);doc=json.loads((directory/'report.json').read_text())
        if not doc['computational_completed']:raise ValueError('Incomplete feedback case')
        for name in ['report.json','rates.npz']:checks[str(directory/name)]=fingerprint(directory/name)
        with np.load(directory/'rates.npz')as data:d={k:data[k].copy()for k in data.files}
        if base is None:base=d
        for key in ['populations','nstar','replay_C_s_inv']:np.testing.assert_array_equal(d[key],base[key])
        idx=np.searchsorted(d['wave_vac_nm'],base['wave_vac_nm']);np.testing.assert_array_equal(d['wave_vac_nm'][idx],base['wave_vac_nm'])
        R=d['replay_R_s_inv'];n=d['populations'];off=R+d['replay_C_s_inv'];inc=np.sum(off*n[None,:,:],axis=1);out=off.sum(axis=0)*n;res=abs(inc-out)/(inc+out)
        if not np.isfinite(res).all()or not np.array_equal(R>0,base['replay_R_s_inv']>0):raise ValueError('Invalid residual/source graph')
        loc=np.unravel_index(np.argmax(res),res.shape)
        maximum=float(np.max(res))
        if maximum!=doc['matched_population_replay_physical_residual']:raise ValueError('Producer residual does not reconstruct')
        rows.append(dict(case=group,points=len(d['wave_vac_nm']),seconds=doc['seconds'],max_physical_residual=maximum,worst_level_1based=int(loc[0]+1),worst_depth_top_to_bottom=int(loc[1]),level13_top_residual=float(res[12,0]),level40_top_residual=float(res[39,0]),max_common_node_J_relative=float(np.max(abs(d['Jnu'][idx]/base['Jnu']-1))),common_node_J_bit_identical=bool(np.array_equal(d['Jnu'][idx],base['Jnu']))))
    baseline=rows[0];frozen=rows[1]
    if frozen['max_common_node_J_relative']>1e-12:raise ValueError('Frozen representation baseline failed')
    refinement=rows[-3:]
    if not all(x['common_node_J_bit_identical']for x in refinement):raise ValueError('Quadrature controls changed inherited-node field')
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'controls.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    labels=['Source','Frozen','BB thermal','BF thermal','BB+BF','loggf','Both']
    axes[0].bar(np.arange(7),[x['max_physical_residual']for x in rows[:7]])
    axes[0].set_xticks(np.arange(7),labels,rotation=35,ha='right');axes[0].set_ylabel('Maximum physical balance residual');axes[0].set_title('Fixed-population field controls')
    selected=[rows[0]]+refinement
    for key,label in [('max_physical_residual','All levels: maximum'),('level13_top_residual','O I level13: top'),('level40_top_residual','O II ground: top')]:
        axes[1].plot([x['points']for x in selected],[x[key]for x in selected],'o-',label=label)
    axes[1].set_yscale('log');axes[1].set_xlabel('Integration wavelengths (all source transitions retained)');axes[1].set_ylabel('Physical balance residual');axes[1].legend(fontsize=8);axes[1].set_title('Grid sensitivity; scientific HOLD')
    fig.savefig(a.output/'feedback_controls.png',dpi=160);plt.close(fig)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed feedback-audit input')
    report=dict(generator=fingerprint(__file__),inputs=checks,cases=rows,computational_completed=True,scientific_admission='HOLD',findings=['Frozen complete coefficients reproduce original J within9e-16; populations/LTE/C exact.','Combined thermal field feedback improves coarse-grid maximum only about2percent; consumed-loggf control does not improve it.','Global node doubling leaves inherited-node J bit-identical, improves the former OI level13 discrepancy, but exposes a larger OII ground/BF balance error.','Local8446 refinements leave the OII maximum essentially unchanged. Continuum quadrature convergence remains unresolved; a closer coarse result is not validated.'],scope='Same matched PRE-update populations and full source graph. Frozen coefficient representation, thermal full-frequency sensitivity, explicit all256 consumed-loggf precision and numerical integration controls remain distinct. Historical strength-case scope prose saying strengths unchanged is superseded by its explicit strength-factor packet and this interpretation; source atom stays immutable, effective absorption strengths do change in those cases. No full native thermal-split reproduction, coupled SE, observational solar, flux, all-band or3D admission. Next: bounded adaptive continuum quadrature under original20000-point cap, then coherent coupled/observed checks.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(rows))


if __name__=='__main__':main()
