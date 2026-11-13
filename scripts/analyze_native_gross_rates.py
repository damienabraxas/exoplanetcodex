"""Qualify non-invasive stage-matched native gross-rate exports, independently."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.native_gross_rates import read_state,read_integrals,physical_matrix
    from pipeline.openmulti3d_output import read_rate_pair
    from pipeline.reference_engines import fingerprint
    from pipeline.bride_rate_audit import process_flow
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','baseline','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();checks={}
    for directory in [a.source,a.baseline]:
        doc=json.loads((directory/'report.json').read_text())
        if not doc['computational_completed']:raise ValueError('Incomplete native producer')
        for row in doc['output_files']:
            q=directory/row['path']
            if fingerprint(q)!=row['receipt']:raise ValueError('Changed producer output')
            checks[str(q)]=row['receipt']
        checks[str(directory/'report.json')]=fingerprint(directory/'report.json')
    if json.loads((a.source/'report.json').read_text())['engine']['original_accurate_engine']!=json.loads((a.baseline/'report.json').read_text())['engine_receipt']:raise ValueError('Wrong baseline engine')
    old=a.baseline/'native/save';new=a.source/'native/save';original_files=sorted(old.glob('*.bin'))
    for q in original_files:
        if q.read_bytes()!=(new/q.name).read_bytes():raise ValueError('Observer changed native output '+q.name)
    table=(a.source/'native/atoms/O.txt').read_text().splitlines();i=next(i for i,x in enumerate(table)if x.split()==['ev','g','ion','label']);weights=np.array([float(x.split()[2])for x in table[i+1:i+42]])
    parts=[];iterations=[];counts=0
    for q in sorted(new.glob('rates_O_*.bin')):
        native=read_rate_pair(q,q.with_name(q.name.replace('rates_','colrates_')));it=native['iteration'];zs=int(q.stem.split('_')[-1]);state=read_state(new/f'audit_state_{it:05d}_{zs:05d}.bin')
        sums=read_integrals(new/f'audit_gross_{it-1:05d}_{zs:05d}.bin',state,radiation_iteration=it-1);R=physical_matrix(state,sums,weights);C=native['collision_rates_s_inv'];n=state['populations_cm3'];star=state['nstar_cm3']
        np.testing.assert_array_equal(star,native['nstar_cm3']);net,gross=process_flow(R+C,n);res=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)));change=float(np.max(abs(native['populations_cm3']/n-1)))
        parts.append(dict(state=state,R=R,C=C,n=n,star=star,post=native['populations_cm3'],net=net,gross=gross));iterations.append(dict(zs=zs,iteration=it,physical_residual_at_actual_radiation_population=res,post_update_max_population_change=change,integrals=sums[2]));counts+=sums[2]['records']
    if len(parts)!=8 or len({x['iteration']for x in iterations})!=1:raise ValueError('Not one complete final global stage')
    arrays={key:np.concatenate([x[key]for x in parts],axis=2 if key in ['R','C']else 1)for key in ['R','C','n','star','post','net','gross']}
    if arrays['n'].shape!=(41,56):raise ValueError('Wrong observer mesh')
    report=dict(generator=fingerprint(__file__),reader=fingerprint(ROOT/'pipeline/native_gross_rates.py'),inputs=checks,scientific_admission='HOLD',computational_completed=True,original_binary_files_bit_identical=len(original_files),final_stage=iterations,records=counts,max_physical_residual=float(np.max(abs(arrays['net'])/np.maximum(arrays['gross'],np.finfo(float).tiny))),max_update_relative=float(np.max(abs(arrays['post']/arrays['n']-1))),semantics='Unpreconditioned rate integrations at exact PRE-update radiation populations; not a new formal solve at post-update populations. FP32 native contributions accumulated in FP64; physical conventions and source data remain independently unvalidated.')
    a.output.mkdir(parents=True,exist_ok=False);np.savez_compressed(a.output/'rates.npz',**arrays);(a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['original_binary_files_bit_identical','records','max_physical_residual','max_update_relative']}))
if __name__=='__main__':main()
