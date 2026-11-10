"""Check opt-in variant plumbing preserves every original native output analysis array."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['first','second','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();inputs={};checked={}
    for name in ['arrays.npz','spectra.npz']:
        ds=[]
        for folder in [a.first,a.second]:
            path=folder/name;inputs[str(path)]=fingerprint(path)
            with np.load(path)as d:ds.append({k:d[k].copy()for k in d.files})
        if set(ds[0])!=set(ds[1]):raise ValueError('Changed native export schema')
        for key in ds[0]:np.testing.assert_array_equal(ds[0][key],ds[1][key])
        checked[name]=list(ds[0])
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'report.json').write_text(json.dumps(dict(generator=fingerprint(__file__),inputs=inputs,checked=checked,all_arrays_bit_identical=True,scientific_admission='HOLD'),indent=2)+'\n');print('Original native arrays and spectra bit-identical')
if __name__=='__main__':main()
