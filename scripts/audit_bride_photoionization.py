"""Independent UV photoionization quadrature on verified fixed-population replays.

Uses linear interpolation of raw cross-sections, not backend rate integration.
Differences from backend WENO/grid weights are reported, never calibrated away.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rate-report',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.multi_atom import read_multi_atom
    import lightweaver.constants as c
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report=json.loads(args.rate_report.read_text())
    if report['status']!='DIAGNOSTIC_COMPLETE':raise ValueError('Incomplete rate replay')
    source_path=require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(source_path)!=report['atom']:raise ValueError('Changed source atom')
    source=read_multi_atom(source_path)
    pairs=[(x.j,x.i)for x in source.continua]
    if len(set(pairs))!=len(pairs):raise ValueError('Ambiguous continuum endpoints')
    args.output.mkdir(parents=True,exist_ok=False)
    rows=[];kernels={};fig,axes=plt.subplots(1,2,figsize=(12,4))
    result=dict(schema='codex.bride_photoionization_audit/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),source_report=fingerprint(args.rate_report),
        experimental_collision_convention=report.get('experimental_collision_convention'),
        lte_uv_background_lines=report.get('lte_uv_background_lines',False),
        finite_column_top_cap=report.get('finite_column_top_cap',False),
        experimental_fe_background=report.get('experimental_fe_background'),
        atom=report['atom'],arrays=[],
        method='4*pi/h integral sigma(lambda)*Jnu(lambda)/lambda d(lambda), raw linear sigma; backend WENO/grid-weight differences retained',
        interpretation='J/B indicates nonlocal radiation, not proof of opacity error. No rate, cross-section, boundary or abundance change.')
    for entry in report['modes']:
        path=args.rate_report.parent/entry['array']['name']
        if fingerprint(path)!=entry['array']:raise ValueError('Changed replay array')
        result['arrays'].append(fingerprint(path))
        data=np.load(path);wave=data['wavelength_vac_nm'];J=data['Jnu_W_m2_sr_Hz'];B=data['Bnu_W_m2_sr_Hz']
        for continuum in source.continua:
            raw_wave=continuum.wavelength[::-1]*.1
            raw_sigma=continuum.alpha[::-1]*1e-4
            mask=(wave>=raw_wave[0])&(wave<=raw_wave[-1])
            w=wave[mask]
            if len(w)<2:raise ValueError('Unresolved continuum support')
            sigma=np.interp(w,raw_wave,raw_sigma)
            kernel=4*np.pi/c.HPlanck*sigma[:,None]*J[mask]/w[:,None]
            thermal=4*np.pi/c.HPlanck*sigma[:,None]*B[mask]/w[:,None]
            photo=np.trapezoid(kernel,w,axis=0);planck=np.trapezoid(thermal,w,axis=0)
            backend=data['radiative_rates_s_inv'][continuum.j-1,continuum.i-1]
            for target in (-4.9,-4.,-3.,-2.,-1.2621,0.):
                depth=int(np.argmin(abs(data['logtau500']-target)))
                rows.append(dict(mode=entry['name'],lower_1based=continuum.i,upper_1based=continuum.j,
                    logtau500=data['logtau500'][depth],temperature_K=data['temperature_K'][depth],
                    threshold_vac_nm=raw_wave[-1],support_min_vac_nm=raw_wave[0],
                    independent_photo_s_inv=photo[depth],backend_photo_s_inv=backend[depth],
                    independent_over_backend=photo[depth]/backend[depth],
                    independent_planck_photo_s_inv=planck[depth],J_weighted_over_B=photo[depth]/planck[depth]))
            if continuum.i==6 and continuum.j==40:
                # Diagnostic top point and previously declared photospheric sample.
                depths=[0,int(np.argmin(abs(data['logtau500']+1.2621)))]
                if entry['name']=='local_nlte':
                    for ax,depth in zip(axes,depths):
                        ax.plot(w,kernel[:,depth],label='Local Jnu')
                        ax.plot(w,thermal[:,depth],label='Local Planck Bnu',ls='--')
                        ax.set(xlabel='Vacuum wavelength (nm)',ylabel='Photoionization kernel (s^-1 nm^-1)',
                               title=f'Lower triplet continuum, logtau={data["logtau500"][depth]:.4f}')
                        ax.set_yscale('log');ax.legend()
                kernels[entry['name']]=dict(wavelength_vac_nm=w.tolist(),
                    selected_depth_indices=depths,kernel_s_inv_nm=kernel[:,depths].T.tolist(),
                    planck_kernel_s_inv_nm=thermal[:,depths].T.tolist())
    with (args.output/'photoionization.csv').open('x')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    result['triplet_lower_selected']=[r for r in rows if r['lower_1based']==6 and r['upper_1based']==40]
    result['comparison']=fingerprint(args.output/'photoionization.csv')
    result['selected_kernels']=kernels
    fig.tight_layout();fig.savefig(args.output/'photoionization_kernel.png',dpi=150);plt.close(fig)
    result['figure']=fingerprint(args.output/'photoionization_kernel.png')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['triplet_lower_selected'][:6],indent=2))


if __name__=='__main__':main()
