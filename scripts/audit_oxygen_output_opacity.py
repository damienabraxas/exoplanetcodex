"""Paired native output-opacity control, explicitly checking unchanged SE state."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.bride_observables import window_integral
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline','bound-free','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();inputs={};data=[]
    for folder in [a.baseline,a.bound_free]:
        d={}
        for name in ['arrays','spectra']:
            file=folder/(name+'.npz');inputs[str(file)]=fingerprint(file)
            with np.load(file)as z:d.update({k:z[k].copy()for k in z.files})
        data.append(d)
    b,d=data
    for key in ['native_populations_cm3','native_nstar_cm3','height_m','wave_vac_A']:np.testing.assert_array_equal(b[key],d[key])
    wave=b['wave_air_A'];rows=[]
    for c in [6158.17,6300.304,7771.944,7774.166,7775.388,8446.4]:
        row=dict(centre_air_A=c)
        for label,key in [('NLTE','Inu_cgs'),('LTE','LTE_Inu_cgs')]:
            yb=b[key]/b['continuum_Inu_cgs'];yd=d[key]/d['continuum_Inu_cgs']
            row[label+'_area_change_mA']=1000*(window_integral(wave,yd,c-.5,c+.5,max_gap=.0051)-window_integral(wave,yb,c-.5,c+.5,max_gap=.0051))
        rows.append(row)
    for path,receipt in inputs.items():
        if fingerprint(path)!=receipt:raise ValueError('Changed output-control input')
    a.output.mkdir(parents=True,exist_ok=False)
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',inputs=inputs,SE_populations_bit_identical=True,windows=rows,limitations=['Native output bound-free switch only; not complete opacity parity or observational validation.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(rows))


if __name__=='__main__':main()
