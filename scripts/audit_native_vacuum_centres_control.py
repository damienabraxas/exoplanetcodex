"""Verify photon-centre packet control preserves native state and continuum exactly."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={}
    packet=a.raw_root/'rya1371-native-wavelength-roundtrip-v2/output.txt';inputs[str(packet)]=fingerprint(packet);centres=np.loadtxt(packet)[:,2]
    for case in ['native_NLTE','native_LTE']:
        ds=[]
        for name in ['rya1371-hop-no-vdw-transfer-v1','rya1371-hop-native-vacuum-centres-v1']:
            path=a.raw_root/name/(case+'.npz');inputs[str(path)]=fingerprint(path)
            with np.load(path)as d:ds.append({k:d[k].copy()for k in d.files})
        for key in ['populations','nstar','continuum_Inu_SI','wave_air_A','wave_vac_A']:np.testing.assert_array_equal(ds[0][key],ds[1][key])
    path=a.raw_root/'rya1371-hop-native-vacuum-centres-v1/report.json';inputs[str(path)]=fingerprint(path);report=json.loads(path.read_text())
    if report['dry_run']or any(r['SE_updated']for r in report['records']):raise ValueError('Invalid fixed-state transfer')
    for record in report['records']:
        audit=record['audit']['table_control']
        if 'no conversion'in audit['vacuum_centre_source']and len(audit['changes'])==10:
            for line in audit['changes']:np.testing.assert_allclose(line['transfer_vac_nm'],centres[line['line_index']]/10,rtol=1e-15)
        else:raise ValueError('Wrong wavelength source or transition support')
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'report.json').write_text(json.dumps(dict(generator=fingerprint(__file__),inputs=inputs,native_state_continuum_axes_bit_identical=True,actual_transfer_centres_match_native_packet=True,scientific_admission='HOLD',limitations=['Fixed-population diagnostic only. No global converter, source atom, observation axis, physical widths or production default changed.']),indent=2)+'\n');print('Exact state/continuum/axes and actual native photon-centre use verified')
if __name__=='__main__':main()
