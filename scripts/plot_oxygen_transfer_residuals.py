"""Plot own-state LW-minus-native residuals so close profile overlays remain reviewable."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={};curves=[]
    cases=[('QLG-bot, table matched','rya1371-oxygen-linear-no-stark-analysis-v1','rya1371-hop-table-combined-v1'),('Native fp64 kernel','rya1371-oxygen-accurate-profile-analysis-v1','rya1371-hop-accurate-profile-transfer-v1'),('VDW-off diagnostic','rya1371-oxygen-no-vdw-analysis-v1','rya1371-hop-no-vdw-transfer-v1'),('Native vacuum-centre packet','rya1371-oxygen-no-vdw-analysis-v1','rya1371-hop-native-vacuum-centres-v1')]
    for label,native,transfer in cases:
        ds=[]
        for path in [a.raw_root/native/'spectra.npz',a.raw_root/transfer/'native_NLTE.npz']:
            inputs[str(path)]=fingerprint(path)
            with np.load(path)as d:ds.append({k:d[k].copy()for k in d.files})
        n,t=ds;np.testing.assert_array_equal(n['wave_vac_A'],t['wave_vac_A']);np.testing.assert_array_equal(n['wave_air_A'],t['wave_air_A']);res=t['normalized']-n['Inu_cgs']/n['continuum_Inu_cgs']
        if not np.isfinite(res).all():raise ValueError('Nonfinite profile residual')
        curves.append((label,n['wave_air_A'],res))
    fig,axes=plt.subplots(3,2,figsize=(12,10),constrained_layout=True)
    for ax,centre in zip(axes.ravel(),[6158.17,6300.304,7771.944,7774.166,7775.388,8446.4]):
        for label,wave,res in curves:
            mask=abs(wave-centre)<=.5;ax.plot(wave[mask]-centre,res[mask]*1e4,label=label,lw=1.2)
        ax.axhline(0,color='0.6',lw=.6);ax.set_title(f'{centre:.3f} AIR A');ax.set_xlabel('AIR offset / A');ax.set_ylabel('LW/HOP minus native, delta(I/Ic) x 10000')
    axes.ravel()[0].legend(fontsize=8);fig.suptitle('Each residual uses its own complete native population state; intrinsic continuum, no fitted shift/scale\nVDW-off is a diagnostic, not an adopted physical model',fontsize=11)
    a.output.mkdir(parents=True,exist_ok=False);fig.savefig(a.output/'residuals.png',dpi=170);plt.close(fig)
    for path,f in inputs.items():
        if fingerprint(path)!=f:raise ValueError('Changed plot input')
    (a.output/'report.json').write_text(json.dumps(dict(generator=fingerprint(__file__),inputs=inputs,plot=fingerprint(a.output/'residuals.png'),scientific_admission='HOLD',scope='Own-state six-window NLTE normalized residuals; no observational fit, cross-state population attribution or scientific admission.'),indent=2)+'\n')
if __name__=='__main__':main()
