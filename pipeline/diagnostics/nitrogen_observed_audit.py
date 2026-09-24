"""RYA-1220: observed feature absorption, never a deblended N I EW."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--engine-root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(a.engine_root));sys.path.insert(0,str(a.engine_root/'scripts'))
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from measure_band_ew import load_window_ex
    from pipeline.wavelength_util import vac_to_air
    lines=[7442.29,8216.33,8629.23,8683.40,10108.90]
    holdings=[('iag_fts_solar_atlas','solar_iag'),
              ('kpno_solar_atlas','solar_kpno_kurucz2005_corrected'),
              ('kpno_solar_atlas','solar_kpno_molecfit_corrected')]
    rows=[];fig,axes=plt.subplots(5,1,figsize=(9,13))
    for i,wave in enumerate(lines):
        for inst,holding in holdings:
            try:
                win=load_window_ex(inst,wave,2.,holding=holding)
                w,f=np.asarray(win.wave),np.asarray(win.flux)
                finite=np.isfinite(w)&np.isfinite(f)
                w,f=w[finite],f[finite]
                if len(w)<5: raise ValueError('insufficient observed pixels')
                pd.DataFrame({'wavelength_air_A':w,'flux':f}).to_csv(a.out/f'{holding}_{wave:.2f}.csv',index=False)
                axes[i].plot(w-wave,f,label=holding,lw=.8)
                for half in [.15,.25,.5,1.1]:
                    m=abs(w-wave)<=half
                    if m.sum()<3: continue
                    absorbed=float(np.trapezoid(1-f[m],w[m])*1000)
                    rows.append(dict(line_air_A=wave,holding=holding,half_width_A=half,
                        n_pixels=int(m.sum()),feature_absorption_mA=absorbed,
                        minimum_flux=float(f[m].min()),geometry='disk_integrated_flux',
                        disposition='BLENDED_FEATURE_NOT_ATOMIC_EW',error=''))
            except Exception as e:
                rows.append(dict(line_air_A=wave,holding=holding,half_width_A=None,
                    n_pixels=0,feature_absorption_mA=None,minimum_flux=None,
                    geometry='disk_integrated_flux',disposition='NOT_MEASURED',error=str(e)))
        axes[i].set(xlabel=f'Wavelength offset from {wave:.2f} A',ylabel='Normalized flux',xlim=(-1.1,1.1))
        axes[i].axvline(0,color='black',ls=':',lw=.5)
    axes[0].legend(fontsize=7);fig.tight_layout()
    fig.savefig(a.out/'ni_observed_profiles.pdf');fig.savefig(a.out/'ni_observed_profiles.png',dpi=150)
    pd.DataFrame(rows).to_csv(a.out/'ni_observed_feature_audit.csv',index=False)
    jpath=a.engine_root/'data/results/rya1214_crires_jk/solar_crires_plus_j_rya1219_rest.csv'
    j=pd.read_csv(jpath)
    # Fetch the actual edge pixels, rather than assuming a broad instrument span.
    iw=load_window_ex('iag_fts_solar_atlas',11080.,3.,holding='solar_iag')
    edge=float(np.max(iw.wave))
    mol=pd.read_csv(a.engine_root/'data/reference/amarsi2021_cno/derived/amarsi2021_cno_molecular_lines.csv')
    cn=mol[(mol.species=='CN') & (mol.element_parameter=='logepsN')].copy()
    cn['air_A']=vac_to_air(cn.wavelength_vac_nm.to_numpy()*10)
    overlap=cn[(cn.air_A<=edge)&(cn.air_A>=j.wavelength_air_A.min())]
    summary={'iag_verified_pixels_through_A':edge,'crires_j_start_A':float(j.wavelength_air_A.min()),
             'matched_cn_transitions':len(overlap),'comparison_status':'NO_SHARED_WAVELENGTH_COVERAGE',
             'j_source_sha256':hashlib.sha256(jpath.read_bytes()).hexdigest(),
             'note':'IAG upper edge probe is conservative; declared holding end 11083.46 A also lies below J.'}
    (a.out/'cn_overlap_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary));print('Read original spectra; wrote diagnostic CSV and plots only.')


if __name__=='__main__':main()
