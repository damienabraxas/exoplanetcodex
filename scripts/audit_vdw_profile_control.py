"""Verify native VDW switch, unchanged physical state/rates and actual fixed-pair routing."""
import argparse,json,re,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.openmulti3d_output import read_rate_pair
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={};arrays=[]
    for name in ['rya1371-oxygen-accurate-profile-analysis-v1','rya1371-oxygen-no-vdw-analysis-v1']:
        path=a.raw_root/name/'arrays.npz';inputs[str(path)]=fingerprint(path)
        with np.load(path)as d:arrays.append({k:d[k].copy()for k in d.files})
    for key in ['height_m','native_nstar_cm3','consumed_state_cgs']:np.testing.assert_array_equal(arrays[0][key],arrays[1][key])
    params=a.raw_root/'rya1371-oxygen-no-vdw-v1/native/params.nml';inputs[str(params)]=fingerprint(params);text=params.read_bytes().decode('latin1')
    if not re.search(r'\bDISABLE_VDW\s*=\s*T\b',text):raise ValueError('Native saved VDW control not true')
    for path in sorted((a.raw_root/'rya1371-oxygen-accurate-profile-v1/native/save').glob('rates_O_*.bin')):
        pairs=[]
        for name in ['rya1371-oxygen-accurate-profile-v1','rya1371-oxygen-no-vdw-v1']:
            rate=a.raw_root/name/'native/save'/path.name;col=rate.parent/path.name.replace('rates_','colrates_')
            for f in [rate,col]:inputs[str(f)]=fingerprint(f)
            pairs.append(read_rate_pair(rate,col))
        np.testing.assert_array_equal(pairs[0]['collision_rates_s_inv'],pairs[1]['collision_rates_s_inv'])
    path=a.raw_root/'rya1371-hop-no-vdw-transfer-v1/report.json';inputs[str(path)]=fingerprint(path);r=json.loads(path.read_text())
    if r['dry_run']or len(r['records'])!=2:raise ValueError('Incomplete fixed pair')
    for record in r['records']:
        if record['SE_updated']:raise ValueError('Unexpected SE update')
        for key in ['audit','continuum_audit']:
            audit=record[key]
            if audit['vdw_broadening']or audit['elastic_terms_remaining']or audit['perturbers']:raise ValueError('VDW-off switch did not reach actual line model')
    path=a.raw_root/'rya1371-oxygen-no-vdw-analysis-v1/report.json';inputs[str(path)]=fingerprint(path);n=json.loads(path.read_text())
    if n['stored_rate_max_residual']>1e-6 or n['max_conservation']>1e-12 or len(n['distinct_final_iterations'])!=1:raise ValueError('Unqualified native rate balance')
    a.output.mkdir(parents=True,exist_ok=False);report=dict(generator=fingerprint(__file__),inputs=inputs,scientific_admission='HOLD',mesh_thermodynamics_nstar_collisions_bit_identical=True,native_max_relative_n_response=float(np.max(abs(arrays[1]['native_populations_cm3']/arrays[0]['native_populations_cm3']-1))),native_stored_rate_residual=n['stored_rate_max_residual'],limitations=['Explicit VDW-off diagnostic only, not physical omission or width/default adoption. Natural widths retained.','Own-state transfer and stored-rate balance do not qualify fresh-field/collider/SE/observed solar physics.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['native_max_relative_n_response','native_stored_rate_residual']}))
if __name__=='__main__':main()
