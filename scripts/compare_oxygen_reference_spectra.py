"""Matched vertical intensity comparison, with explicit sampled-window refusals."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.bride_observables import window_integral
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['native','bride','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();inputs={}
    def read(path):
        inputs[str(path)]=fingerprint(path)
        with np.load(path)as d:return {k:d[k].copy()for k in d.files}
    native=read(a.native/'spectra.npz');cases={name:read(a.bride/(name+'.npz'))for name in ['local','held_multi_b','LTE']}
    wave=cases['local']['wave_air_A'];vac=cases['local']['wave_vac_A']
    selected=np.searchsorted(native['wave_vac_A'],vac)
    if np.any(selected>=len(native['wave_vac_A'])):raise ValueError('Missing native wavelength')
    native={k:v[selected]for k,v in native.items()}
    # The native array contains exactly the requested segments. Require exact
    # wavelength identity before point-by-point engine comparisons.
    for d in [native,*cases.values()]:
        np.testing.assert_array_equal(d['wave_vac_A'],vac)
        np.testing.assert_array_equal(d['wave_air_A'],wave)
    profiles={name:d['normalized']for name,d in cases.items()}
    profiles['native_NLTE']=native['Inu_cgs']/native['continuum_Inu_cgs']
    profiles['native_LTE']=native['LTE_Inu_cgs']/native['continuum_Inu_cgs']
    if any(x.shape!=wave.shape or not np.isfinite(x).all()for x in profiles.values()):raise ValueError('Invalid profile array')
    rows=[]
    for centre in [6158.17,6300.304,6363.776,7771.944,7774.166,7775.388,8446.4]:
        lo,hi=centre-.5,centre+.5;mask=(wave>=lo)&(wave<=hi)
        row=dict(centre_air_A=centre,lo_air_A=lo,hi_air_A=hi,coverage='absent atomic transition control'if centre==6363.776 else 'theoretical O-only window')
        for name,profile in profiles.items():row[name+'_area_mA']=1000*window_integral(wave,profile,lo,hi,max_gap=.0051)
        row['max_abs_normalized_native_minus_bride']=float(np.max(abs(profiles['native_NLTE'][mask]-profiles['local'][mask])))
        row['max_abs_normalized_multi_shape_minus_bride']=float(np.max(abs(profiles['held_multi_b'][mask]-profiles['local'][mask])))
        row['max_relative_native_continuum_minus_bride']=float(np.max(abs(native['continuum_Inu_cgs'][mask]*1e-3/cases['local']['continuum_Inu_SI'][mask]-1)))
        rows.append(row)
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'windows.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,figsize=(12,10))
    for ax,centre in zip(axes.flat,[6158.17,6300.304,7771.944,7774.166,7775.388,8446.4]):
        mask=abs(wave-centre)<=.5
        for name in ['local','held_multi_b','LTE','native_NLTE','native_LTE']:ax.plot(wave[mask]-centre,profiles[name][mask],label=name)
        ax.set_title(str(centre)+' A / vertical intensity');ax.set_xlabel('AIR wavelength offset / A');ax.set_ylabel('I / Ic');ax.legend(fontsize=7)
    fig.suptitle('5750 K / logg4.5 / A(O)8.70: same source atom; native physics differences remain')
    fig.tight_layout();fig.savefig(a.output/'comparison.png',dpi=150);plt.close(fig)
    for path,receipt in inputs.items():
        if fingerprint(path)!=receipt:raise ValueError('Input changed during comparison')
    report=dict(schema='codex.oxygen_reference_spectra/1',scientific_admission='HOLD',generator=fingerprint(__file__),inputs=inputs,windows=rows,limitations=['Finite +/-0.5A areas are not full equivalent widths or abundance corrections.','Vertical intensity only; no flux/intensity substitution or observation fit.','Same source atmosphere/atom does not imply identical consumed mesh/EOS/opacity/collisions/boundary.','6363A transition absent from held atom; continuum/tail signal is not line validation.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
