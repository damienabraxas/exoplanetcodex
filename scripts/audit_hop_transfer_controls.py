"""Audit matched linear and Stark semantic controls; no physical default adoption."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint,load_engine
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={}
    def arrays(name):
        path=a.raw_root/name/'arrays.npz';inputs[str(path)]=fingerprint(path)
        with np.load(path)as d:return {k:d[k].copy()for k in d.files}
    baseline=arrays('rya1371-oxygen-dense-analysis-v3');linear=arrays('rya1371-oxygen-linear-analysis-v1');off=arrays('rya1371-oxygen-linear-no-stark-analysis-v1')
    if set(baseline)!=set(linear)or set(baseline)!=set(off):raise ValueError('Changed native export schema')
    for key in baseline:np.testing.assert_array_equal(baseline[key],linear[key])
    for key in ['height_m','consumed_state_cgs','native_nstar_cm3','source_logtau500_on_native_height']:np.testing.assert_array_equal(baseline[key],off[key])
    reports={}
    for tag,name in [('baseline','rya1371-hop-fixed-transfer-analysis-v1'),('linear','rya1371-hop-linear-comparison-v1'),('stark_off','rya1371-hop-no-stark-comparison-v1')]:
        path=a.raw_root/name/'report.json';inputs[str(path)]=fingerprint(path);reports[tag]=json.loads(path.read_text())
    _,app,_,_=load_engine('openmulti3d')
    sources={str(path):fingerprint(path)for path in [app/'HOP/bbtrans_mod.f90',app/'SRC/rt_params_mod.f90',ROOT/'pipeline/multi_atom.py']}
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',inputs=inputs,sources=sources,linear_native_all_arrays_bit_identical=True,stark_off_thermodynamics_mesh_nstar_bit_identical=True,stark_off_max_relative_n_change=float(np.max(abs(off['native_populations_cm3']/baseline['native_populations_cm3']-1))),comparisons=reports,interpretation=['Matched linear formal transfer reduces continuum differences from0.109-0.142percent to below0.002percent, with identical native populations.','Triplet finite-window intensity areas still differ0.56-0.82mA under matched linear transfer.','Pinned HOP gstk0 requests automatic Cowley estimate; Bride gStark0 disables Stark. Native Stark-off control changes populations, so paired transfer consumes its own complete native population state.','No Stark/He/model/default/source-atom adoption. Remaining line profile and damping conventions require explicit qualification; absent6363 is still excluded.'])
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(linear_arrays_identical=True,stark_off_n_response=report['stark_off_max_relative_n_change'])))
if __name__=='__main__':main()
