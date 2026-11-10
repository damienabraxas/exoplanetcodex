"""Unfitted shipped-atlas overlay of the native LTE molecular-line control."""
import argparse,json,os,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    os.environ['CODEX_CONT_REF']='iag'
    from pipeline import local_continuum as lc
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    analysis=json.loads((a.analysis/'report.json').read_text());run=json.loads((a.source/'report.json').read_text())
    if fingerprint(a.source/'report.json')!=analysis['source']:raise ValueError('Changed native source')
    comparison_path=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json';comparison=json.loads(comparison_path.read_text());b=comparison['broadening']
    if fingerprint(comparison_path)!=analysis['comparison']:raise ValueError('Changed reference broadening')
    observation=require_codex_path('data.kurucz2005_residual')
    if fingerprint(observation)!=analysis['observation']:raise ValueError('Changed Kitt source')
    kw,kf=read_kurucz2005_residual(observation,7765,7782)
    iw,ifl=lc.solar_atlas_reference(7765,7782,b['R'])
    if fingerprint(require_codex_path('data.solar_iag_baker2020'))!=analysis['iag']:raise ValueError('Changed IAG source')
    curves={}
    for name in ('off','all'):
        path=a.source/name/'flux.csv'
        if fingerprint(path)!=next(c['flux']for c in run['cases']if c['name']==name):raise ValueError('Changed native flux')
        raw=np.genfromtxt(path,delimiter=',',names=True);w=raw['wavelength_air_A'];segment=np.array([(w[0]/10,w[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        curves[name]=(w,ispec.apply_post_fundamental_effects(w/10,raw['normalized_flux'].copy(),segment,
            macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(SOLAR_GRAVITATIONAL_REDSHIFT,)))
    a.output.mkdir(parents=True,exist_ok=False);fig,axes=plt.subplots(2,2,figsize=(13,8))
    for axis,(left,right,label)in zip(axes.ravel(),[(7771.2,7772.7,'O7771'),(7773.4,7774.9,'O7774'),(7774.65,7776.15,'O7775'),(7776.4,7777.8,'FeII neighborhood')]):
        for wave,flux,name,color in [(kw,kf,'Kitt shipped','black'),(iw,ifl,'IAG shipped helper (Gaussian R)','grey')]:
            keep=(wave>=left)&(wave<=right);axis.plot(wave[keep],flux[keep],label=name,color=color,lw=1.)
        for name,color in [('off','tab:blue'),('all','tab:orange')]:
            wave,flux=curves[name];keep=(wave>=left)&(wave<=right);axis.plot(wave[keep],flux[keep],label='Native LTE: molecular lines '+name,color=color,lw=1.3)
        axis.set(title=label,xlabel='Observer AIR A',ylabel='Normalized flux');axis.grid(alpha=.2);axis.legend(fontsize=7)
    fig.suptitle('Native5750/4.5 LTE/O8.70 molecular-line control, HOLD\nFixed broadening + gravity; shipped atlas levels, no observed refit or abundance/width/RV fit')
    fig.tight_layout();figure=a.output/'profiles.png';fig.savefig(figure,dpi=150);plt.close(fig)
    (a.output/'report.json').write_text(json.dumps(dict(schema='codex.bride_native_molecules_plot/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        source=fingerprint(a.source/'report.json'),analysis=fingerprint(a.analysis/'report.json'),reference_helper=fingerprint(lc.__file__),
        figure=fingerprint(figure),scope='Illustrative native LTE control, not matched Bride NLTE, isolated blend validation or calibrated FTS LSF; both RV analyses retained separately.'),indent=2)+'\n');print(a.output)


if __name__=='__main__':main()
