"""Plot finite crossed-continuum accounting without adopting a continuum."""
import argparse,csv,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    import numpy as np
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    from pipeline.bride_reference import fingerprint
    report=json.loads((a.source/'report.json').read_text())
    if report['schema']!='codex.bride_crossed_continuum/1':raise ValueError('Wrong crossed source')
    records=[r for r in report['records'] if r['reference']=='synthesis' and r['centre_air_A']==7775.388]
    a.output.mkdir(parents=True,exist_ok=False)
    fig,axes=plt.subplots(2,1,figsize=(12,8))
    for axis,case in zip(axes,('canonical_refined','lower_gravity_refined')):
        selected=[r for r in records if r['case']==case];x=np.arange(len(selected));width=.25
        for offset,key,label in [(-width,'support_path_average','Support path average'),(0,'denominator_path_average','Denominator path average'),(width,'total','Endpoint total')]:
            axis.bar(x+offset,[r['accounting']['observed_area_mA'][key] for r in selected],width,label=label)
        labels=[]
        for r in selected:
            design={'selected_with_clipping':'Selected + clipping','selected_no_clipping':'Selected, no clipping','coefficient_support_no_clipping':'Coefficient support, fixed'}[r['design']]
            labels.append(r['stark_mode']+'\n'+design)
        axis.set_xticks(x,labels,fontsize=8);axis.axhline(0,color='black',lw=.6)
        axis.set(title=case,ylabel='Observed 7775 window area change (mA)');axis.grid(axis='y',alpha=.2);axis.legend(fontsize=8)
    fig.suptitle('Fixed-O FeII width controls: crossed historical continuum, HOLD\nFinite path accounting; correlated diagnostics, not uncertainty terms or adopted fits')
    fig.tight_layout();fig.savefig(a.output/'accounting.png',dpi=150);plt.close(fig)
    (a.output/'report.json').write_text(json.dumps(dict(schema='codex.bride_crossed_continuum_plot/1',generator=fingerprint(__file__),
        source=fingerprint(a.source/'report.json'),figure=fingerprint(a.output/'accounting.png'),scientific_admission='HOLD'),indent=2)+'\n')
    print(a.output)


if __name__=='__main__':main()
