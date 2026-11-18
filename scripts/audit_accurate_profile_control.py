"""Separate isolated native kernel SE response from matched-state transfer residuals."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint,load_plan_engine,load_engine
    from pipeline.openmulti3d_output import read_rate_pair
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={};arrays=[]
    for name in ['rya1371-oxygen-linear-no-stark-analysis-v1','rya1371-oxygen-accurate-profile-analysis-v1']:
        path=a.raw_root/name/'arrays.npz';inputs[str(path)]=fingerprint(path)
        with np.load(path)as d:arrays.append({k:d[k].copy()for k in d.files})
    old,new=arrays
    for key in ['height_m','native_nstar_cm3','consumed_state_cgs','source_logtau500_on_native_height']:np.testing.assert_array_equal(old[key],new[key])
    count=0
    for name in sorted((a.raw_root/'rya1371-oxygen-linear-no-stark-v1/native/save').glob('rates_O_*.bin')):
        pairs=[]
        for folder in ['rya1371-oxygen-linear-no-stark-v1','rya1371-oxygen-accurate-profile-v1']:
            rate=a.raw_root/folder/'native/save'/name.name;col=rate.parent/name.name.replace('rates_','colrates_')
            for path in [rate,col]:inputs[str(path)]=fingerprint(path)
            pairs.append(read_rate_pair(rate,col))
        np.testing.assert_array_equal(pairs[0]['collision_rates_s_inv'],pairs[1]['collision_rates_s_inv']);count+=1
    if count!=8:raise ValueError('Incomplete native rate slabs')
    tables=[a.raw_root/n/'native/atoms/O.txt'for n in ['rya1371-oxygen-linear-no-stark-v1','rya1371-oxygen-accurate-profile-v1']]
    for path in tables:inputs[str(path)]=fingerprint(path)
    if tables[0].read_bytes()!=tables[1].read_bytes():raise ValueError('Native consumed atom changed beyond kernel')
    for name in ['rya1371-oxygen-accurate-profile-analysis-v1','rya1371-hop-accurate-profile-transfer-v1']:
        path=a.raw_root/name/'report.json';inputs[str(path)]=fingerprint(path);r=json.loads(path.read_text())
        if name.startswith('rya1371-oxygen'):
            if r['profile_variant']!='accurate-kernel-control'or r['stored_rate_max_residual']>1e-6 or r['max_conservation']>1e-12 or len(r['distinct_final_iterations'])!=1:raise ValueError('Unqualified native solve/rates')
        else:
            if r['dry_run']or any(x['SE_updated']for x in r['records']):raise ValueError('Unqualified fixed pair')
    _,_,_,original=load_engine('openmulti3d');_,_,_,variant=load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'))
    a.output.mkdir(parents=True,exist_ok=False);report=dict(generator=fingerprint(__file__),inputs=inputs,original=original,variant=variant,scientific_admission='HOLD',mesh_thermodynamics_nstar_atom_table_collision_rates_bit_identical=True,native_max_relative_n_response=float(np.max(abs(new['native_populations_cm3']/old['native_populations_cm3']-1))),limitations=['Native stored physical-rate balance independently reconstructed, not a fresh fixed-field replay.','Accurate native kernel is the existing fp64 reference retaining fp32 input/output and frequency Jacobian; no production/default adoption.','Own-state fixed pair consumes complete new n/nStar; small triplet transfer residuals do not qualify collider/EOS/SE or observed solar physics.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(n_response=report['native_max_relative_n_response'],same_collisions_all_depths=True)))
if __name__=='__main__':main()
