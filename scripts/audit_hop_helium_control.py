"""Check He broadening-only control changed no populations, continuum or wavelength."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={};cases=[]
    for case in ['native_NLTE','native_LTE']:
        arrays=[]
        for run in ['rya1371-hop-no-stark-transfer-v1','rya1371-hop-no-he-transfer-v2']:
            path=a.raw_root/run/(case+'.npz');inputs[str(path)]=fingerprint(path)
            with np.load(path)as d:arrays.append({k:d[k].copy()for k in d.files})
        for key in ['populations','nstar','continuum_Inu_SI','wave_vac_A','wave_air_A']:np.testing.assert_array_equal(arrays[0][key],arrays[1][key])
        cases.append(dict(case=case,unchanged_arrays=['populations','nstar','continuum_Inu_SI','wave_vac_A','wave_air_A']))
    path=a.raw_root/'rya1371-hop-no-he-transfer-v2/report.json';inputs[str(path)]=fingerprint(path);r=json.loads(path.read_text())
    for item in r['records']:
        audit=item['audit']
        if audit['helium_broadening']or not any('explicit He broadening-off'in s for s in audit['perturbers'])or any('effective HeI/HI'in s for s in audit['perturbers']):raise ValueError('He-off diagnostic did not reach actual broadening interface')
    a.output.mkdir(parents=True,exist_ok=False);report=dict(generator=fingerprint(__file__),inputs=inputs,cases=cases,actual_He_off_reads_verified=True,scientific_admission='HOLD',invalid_prior='rya1371-hop-no-he-transfer-v1: plan flag not wired to Perturbers(seed); actual reads retained He. No scientific use.',limitations=['Diagnostic only: native HI-only damping contrasted with LW H+He. No physical preference or abundance/EOS/default change.','Small integrated differences do not certify pointwise profiles, coupled SE, observed solar agreement or all-band applicability.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('He interface, unchanged populations/continuum/wavelength verified')
if __name__=='__main__':main()
