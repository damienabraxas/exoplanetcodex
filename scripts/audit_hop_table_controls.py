"""Account separately for source-determined optical centre and gf rounding controls."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={};results={}
    runs=dict(baseline='rya1371-hop-no-he-transfer-v3',centres='rya1371-hop-table-centres-v2',strengths='rya1371-hop-table-strengths-v1',combined='rya1371-hop-table-combined-v1')
    for case in ['native_NLTE','native_LTE']:
        baseline=None
        for mode,run in runs.items():
            path=a.raw_root/run/(case+'.npz');inputs[str(path)]=fingerprint(path)
            with np.load(path)as d:arrays={k:d[k].copy()for k in d.files}
            if baseline is None:baseline=arrays
            else:
                for key in ['populations','nstar','wave_vac_A','wave_air_A','continuum_Inu_SI']:np.testing.assert_array_equal(baseline[key],arrays[key])
    names=dict(baseline='rya1371-hop-no-he-comparison-v2',centres='rya1371-hop-table-centres-comparison-v2',strengths='rya1371-hop-table-strengths-comparison-v1',combined='rya1371-hop-table-combined-comparison-v1')
    for mode,run in names.items():
        path=a.raw_root/run/'windows.csv';inputs[str(path)]=fingerprint(path)
        with path.open()as f:results[mode]=list(csv.DictReader(f))
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'report.json').write_text(json.dumps(dict(generator=fingerprint(__file__),inputs=inputs,scientific_admission='HOLD',results=results,populations_continuum_wavelengths_bit_identical=True,limitations=['Centre and gf overrides follow actual consumed source table, not a fitted adjustment. Only ten optical transitions in predeclared windows.','All original levels, continua, collision data and other transitions retained; transfer wavelength overrides are not a self-consistent SE atom.','Strength rounding does not uniformly improve area agreement; no recipe/default selected by agreement. Pointwise shape residuals remain.']),indent=2)+'\n');print('All four controls preserve native n/nStar, continuum and wavelengths exactly')
if __name__=='__main__':main()
