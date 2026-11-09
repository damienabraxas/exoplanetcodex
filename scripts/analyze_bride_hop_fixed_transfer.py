"""Compare fixed native-state HOP/LW transfer against independent native intensities."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.bride_observables import window_integral
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['native','transfer','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();inputs={}
    def read(path):
        inputs[str(path)]=fingerprint(path)
        with np.load(path)as f:return {k:f[k].copy()for k in f.files}
    native=read(a.native/'spectra.npz');report_path=a.transfer/'report.json';r=json.loads(report_path.read_text());inputs[str(report_path)]=fingerprint(report_path)
    if r['dry_run']or len(r['records'])!=2 or any(x['SE_updated']for x in r['records']):raise ValueError('Not the declared completed fixed-transfer pair')
    for path,receipt in r['inputs'].items():
        if fingerprint(path)!=receipt:raise ValueError('Changed transfer input')
    cases={}
    for record in r['records']:
        path=a.transfer/(record['case']+'.npz')
        if fingerprint(path)!=record['output']:raise ValueError('Changed transfer output')
        cases[record['case']]=read(path)
    wave=native['wave_air_A'];rows=[]
    for name,d in cases.items():
        np.testing.assert_array_equal(native['wave_vac_A'],d['wave_vac_A']);np.testing.assert_array_equal(wave,d['wave_air_A'])
    for centre in [6158.17,6300.304,6363.776,7771.944,7774.166,7775.388,8446.4]:
        lo,hi=centre-.5,centre+.5;mask=(wave>=lo)&(wave<=hi)
        row=dict(centre_air_A=centre,scope='missing-transition control'if centre==6363.776 else 'finite O-only vertical-intensity window')
        for name,key in [('native_NLTE','Inu_cgs'),('native_LTE','LTE_Inu_cgs')]:
            d=cases[name];foreign=native[key]/native['continuum_Inu_cgs'];local=d['normalized'];ref=1000*window_integral(wave,foreign,lo,hi,max_gap=.0051);current=1000*window_integral(wave,local,lo,hi,max_gap=.0051)
            row[name+'_reference_area_mA']=ref;row[name+'_transfer_area_mA']=current;row[name+'_difference_mA']=current-ref;row[name+'_max_normalized_difference']=float(np.max(abs(local[mask]-foreign[mask])));row[name+'_max_relative_continuum_difference']=float(np.max(abs(d['continuum_Inu_SI'][mask]/(native['continuum_Inu_cgs'][mask]*1e-3)-1)))
        rows.append(row)
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'windows.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,figsize=(12,10))
    for ax,centre in zip(axes.flat,[6158.17,6300.304,7771.944,7774.166,7775.388,8446.4]):
        mask=abs(wave-centre)<=.5
        for name,key in [('native_NLTE','Inu_cgs'),('native_LTE','LTE_Inu_cgs')]:
            ax.plot(wave[mask]-centre,(native[key]/native['continuum_Inu_cgs'])[mask],label='Reference '+name)
            ax.plot(wave[mask]-centre,cases[name]['normalized'][mask],ls='--',label='LW/HOP '+name)
        ax.set_title(str(centre)+' A');ax.set_xlabel('AIR offset / A');ax.set_ylabel('I/Ic');ax.legend(fontsize=7)
    fig.suptitle('Fixed complete native O state / HOP background: remaining transfer/profile conventions')
    fig.tight_layout();fig.savefig(a.output/'comparison.png',dpi=150);plt.close(fig)
    for path,receipt in inputs.items():
        if fingerprint(path)!=receipt:raise ValueError('Input changed during transfer analysis')
    doc=dict(generator=fingerprint(__file__),scientific_admission='HOLD',inputs=inputs,transfer_seconds=r['seconds'],windows=rows,limitations=['Fixed native O populations, not independent coupled SE or collision qualification.','Finite +/-0.5A areas, not complete EWs or abundance corrections; no observed fits.','Native and LW profile/atomic constants, lower ghost and formal-solver differences remain explicit; no result-selected adjustment.','6363 transition absent; control cannot validate that diagnostic.'])
    (a.output/'report.json').write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n');print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
