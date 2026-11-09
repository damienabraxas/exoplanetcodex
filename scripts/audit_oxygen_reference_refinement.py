"""Compare declared native one-axis numerical controls without selecting physics."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.bride_observables import window_integral
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--frequency',type=Path,required=True);p.add_argument('--angle',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    inputs={};loaded={}
    for name in ['baseline','frequency','angle']:
        folder=getattr(a,name);r=json.loads((folder/'report.json').read_text())
        if r['stored_rate_max_residual']>1e-6 or len(r['distinct_final_iterations'])!=1:raise ValueError('Reference numerical audit failed or depths stopped at different iterations')
        for file in [folder/'report.json',folder/'spectra.npz',folder/'arrays.npz']:inputs[str(file)]=fingerprint(file)
        with np.load(folder/'spectra.npz')as d:loaded[name]={k:d[k].copy()for k in d.files}
    base=loaded['baseline'];wave=base['wave_air_A'];records=[]
    for name,d in loaded.items():
        np.testing.assert_array_equal(d['wave_vac_A'],base['wave_vac_A'])
        y=d['Inu_cgs']/d['continuum_Inu_cgs'];b=base['Inu_cgs']/base['continuum_Inu_cgs'];changes=[]
        for c in [6158.17,6300.304,7771.944,7774.166,7775.388,8446.4]:
            diff=1000*(window_integral(wave,y,c-.5,c+.5,max_gap=.0051)-window_integral(wave,b,c-.5,c+.5,max_gap=.0051));changes.append(dict(centre_air_A=c,area_change_mA=diff))
        maximum=max(abs(x['area_change_mA'])for x in changes)
        records.append(dict(control=name,windows=changes,max_area_change_mA=maximum,diagnostic_0p1mA_screen=maximum<=.1,max_normalized_profile_change=float(np.max(abs(y-b)))))
    for path,receipt in inputs.items():
        if fingerprint(path)!=receipt:raise ValueError('Changed refinement input')
    a.output.mkdir(parents=True,exist_ok=False)
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',inputs=inputs,records=records,limitations=['One-axis native frequency/angle tests only, not a bound on all numerical or physical errors.','No selected default, abundance fit, observational validation or production admission.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(records,indent=2))


if __name__=='__main__':main()
