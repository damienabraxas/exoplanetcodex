"""Verify runtime He guard preserves all numerical arrays of qualified transfer."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={}
    for case in ['native_NLTE','native_LTE']:
        ds=[]
        for v in [2,3]:
            path=a.raw_root/f'rya1371-hop-no-he-transfer-v{v}'/(case+'.npz');inputs[str(path)]=dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            with np.load(path)as d:ds.append({k:d[k].copy()for k in d.files})
        for key in ds[0]:np.testing.assert_array_equal(ds[0][key],ds[1][key])
    a.output.mkdir(exist_ok=False);(a.output/'report.json').write_text(json.dumps(dict(inputs=inputs,all_arrays_bit_identical=True,scientific_admission='HOLD',scope='Runtime He routing guard preserves every completed v2 output array; v3 refuses miswired controls at context construction.'),indent=2)+'\n')
if __name__=='__main__':main()
