"""Plot the full-column coefficient screen without inferring emergent flux."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from pipeline.bride_reference import fingerprint
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    source=json.loads((args.source/'report.json').read_text())
    with (args.source/'species_depths.csv').open() as f:rows=list(csv.DictReader(f))
    result=dict(schema='codex.bride_lte_screen_analysis/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        source=fingerprint(args.source/'report.json'),depths=fingerprint(args.source/'species_depths.csv'),
        scope='Absolute-weight footprint averages of local natural-only opacity, no flux/EW or formation-depth inference.',summary=[])
    fig,axes=plt.subplots(3,2,figsize=(12,10))
    for ci,case in enumerate(source['cases']):
        for li,centre in enumerate((7771.944,7774.166,7775.388)):
            fp=f'{centre}_synthesis';rr=[r for r in rows if r['case']==case['name'] and r['footprint']==fp]
            totals={key:np.zeros(56) for key in ('other_candidates','modeled_FeI','Fe_II')};tau=np.full(56,np.nan)
            for r in rr:
                d=int(r['depth_index']);tau[d]=float(r['logtau500']);v=float(r['mean_ratio'])
                totals['modeled_FeI' if r['group']=='modeled_FeI' else 'other_candidates'][d]+=v
                if r['group']=='Fe II':totals['Fe_II'][d]+=v
            if not np.isfinite(tau).all():raise ValueError('Incomplete depth screen')
            ax=axes[li,ci]
            for name,value in totals.items():ax.semilogy(tau,value,label=name)
            ax.axvline(case['reference_logtau500_actual'],color='grey',ls=':',label='nearest logtau500=0')
            ax.set(title=f"{case['name']}: {centre:.3f} AIR A",xlabel='Native logtau500 (all56 depths)',ylabel='Weighted local line / background absorption');ax.grid(alpha=.2);ax.legend(fontsize=8)
            d=case['reference_depth_index'];result['summary'].append(dict(case=case['name'],footprint=fp,
                reference_depth_index=d,logtau500=float(tau[d]),ratios={name:float(value[d]) for name,value in totals.items()}))
    fig.suptitle('Natural-only LTE opacity screen at historical synthesis anchors; HOLD, no predicted flux');fig.tight_layout()
    args.output.mkdir(parents=True,exist_ok=False);fig.savefig(args.output/'depth_ratios.png',dpi=140);plt.close(fig)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
