"""Plot receipted conditional continuum controls without choosing an OPA frame."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from pipeline.bride_reference import fingerprint
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report=json.loads(args.audit.read_text());data={}
    for case in report['cases']:
        path=args.audit.parent/case['array']['name']
        if fingerprint(path)!=case['array']:raise ValueError('Changed coefficient arrays')
        data[case['name']]=dict(np.load(path))
    a=data['source_is_vacuum_matched_MARCS'];b=data['source_is_air_matched_MARCS']
    o=data['source_is_vacuum_benchmark_O870'];wave=a['wavelength_source_A']
    fig,axes=plt.subplots(2,2,figsize=(12,7),sharex=True)
    for row,target in enumerate((-4.9,-1.)):
        depth=int(np.argmin(abs(a['logtau500']-target)))
        for column,kind in enumerate(('absorption','scattering')):
            ax=axes[row,column]
            for d,label,style in ((a,'OPA coordinate treated as vacuum','-'),(b,'OPA coordinate treated as air','--')):
                ax.plot(wave*.1,d['local_'+kind+'_m2_kg'][:,depth]/d['marcs_'+kind+'_m2_kg'][:,depth],style,label=label,lw=1)
            ax.axhline(1,color='black',lw=.7);ax.set_yscale('log');ax.set_xscale('log')
            ax.set(ylabel='Local / MARCS '+kind,title=f'logtau500={a["logtau500"][depth]:.4f}')
            ax.grid(alpha=.2);ax.legend(fontsize=7)
    for ax in axes[-1]:ax.set_xlabel('OPA source coordinate / 10 (medium unverified)')
    fig.suptitle('LTE continuum diagnostics: matched MARCS composition; no spectral admission')
    fig.tight_layout();args.output.mkdir(parents=True,exist_ok=False)
    figure=args.output/'conditional_continuum_ratios.png';fig.savefig(figure,dpi=160);plt.close(fig)
    i=int(np.argmin(abs(wave-2772.477)))
    result=dict(schema='codex.bride_marcs_continuum_analysis/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),audit=fingerprint(args.audit),figure=fingerprint(figure),
        wavelength_source_A_277_sample=float(wave[i]),
        top_277_frame_fractional_effect=float(b['local_absorption_m2_kg'][i,0]/a['local_absorption_m2_kg'][i,0]-1),
        top_277_O866_to_870_fractional_effect=float(o['local_absorption_m2_kg'][i,0]/a['local_absorption_m2_kg'][i,0]-1),
        caveat='Both frame hypotheses retained. FUV below2000A is vacuum in both project-convention hypotheses. Sampling ratios do not certify matched EOS/opacity partition or native producer inputs.')
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
