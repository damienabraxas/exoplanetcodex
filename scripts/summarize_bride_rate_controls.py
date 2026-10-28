"""Archive outcomes of explicit collision-convention controls, including failures."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',nargs='+',type=Path,required=True)
    p.add_argument('--rate-audit',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--gencol-source',type=Path,required=True)
    args=p.parse_args()
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_observables import window_integral
    args.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema='codex.bride_convention_controls/1',scientific_admission='HOLD',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),gencol=fingerprint(args.gencol_source),
        gencol_url='https://github.com/sjbruce3301/FORWARD-coronal-modeling/blob/15145239ffdbb62d19dca0749b11c5920f9fdf13/idl/FORCOMP/NOTUSEDBYFORWARD/EXTRAS/gencol.f',
        exact_grid_producer_build='UNKNOWN; controls are not admitted native MULTI replication',runs=[])
    rate=json.loads((args.rate_audit/'report.json').read_text())
    if rate['status']!='DIAGNOSTIC_COMPLETE':raise ValueError('Rate audit incomplete')
    report['rate_audit_receipt']=fingerprint(args.rate_audit/'report.json')
    report['rate_replay_physical_imbalance']={m['name']:m['physical_imbalance_max'] for m in rate['modes']}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,4))
    prior=ROOT/'data/results/rya1281/population_audit_v1'
    for mode in ('local_lte','local_nlte','external_b_local_lte'):
        path=prior/f'native5750_{mode}_triplet.csv';d=np.genfromtxt(path,delimiter=',',names=True)
        ax.plot(d['air_A'],d['normalized_flux'],label=mode,lw=1)
    ref=np.array(json.loads((ROOT/'data/results/rya1281/population_audit_v1/report.json').read_text())['oxygen_transitions'])
    for directory in args.runs:
        summary=json.loads((directory/'summary.json').read_text())
        if summary['status']!='EXPERIMENTS_FINISHED':raise ValueError('Run set not finished')
        run=dict(name=directory.name,summary=fingerprint(directory/'summary.json'),cases=[])
        for entry in summary['cases']:
            path=directory/entry['name'];m=json.loads((path/'manifest.json').read_text())
            if fingerprint(path/'manifest.json')!=entry['manifest']:raise ValueError('Changed source manifest')
            if m['status'] not in ('FAILED','CONVERGED_EXPERIMENT'):raise ValueError('Nonterminal source case')
            result=dict(name=entry['name'],manifest=entry['manifest'],status=m['status'],
                error=m.get('error'),elapsed_seconds=m['elapsed_seconds'],controls=m['controls'])
            history=[json.loads(line) for line in (path/'history.jsonl').read_text().splitlines()]
            result['history']=fingerprint(path/'history.jsonl');result['recorded_iterations']=len(history)
            result['rejected_invalid_ng_count']=sum(x.get('rejected_invalid_ng',False) for x in history)
            result['last_recorded_diagnostics']=history[-1] if history else None
            if m['status']=='CONVERGED_EXPERIMENT':
                if fingerprint(path/'flux.csv')!=m['flux']:raise ValueError('Changed convention spectrum')
                d=np.genfromtxt(path/'flux.csv',delimiter=',',names=True)
                result['flux']=m['flux'];result['populations']=fingerprint(path/'populations.npz')
                result['window_area_mA']=[1000*window_integral(d['wavelength_air_A'],d['normalized_flux'],c-.5,c+.5) for c in (7771.944,7774.166,7775.388)]
                stem=directory.name+'_'+entry['name']
                np.savetxt(args.output/f'{stem}_profiles.csv',np.column_stack((d['wavelength_air_A'],d['normalized_flux'])),
                    delimiter=',',header='air_A,normalized_flux',comments='')
                ax.plot(d['wavelength_air_A'],d['normalized_flux'],label=entry['name'],ls='--')
            run['cases'].append(result)
        report['runs'].append(run)
    ax.set(xlabel='Air wavelength (Angstrom)',ylabel='Intrinsic normalized flux',
        title='Native5750/4.5/O8.70: explicit non-adopted collision conventions',xlim=(7770.8,7776.2))
    ax.legend(fontsize=7);fig.tight_layout();fig.savefig(args.output/'convention_profiles.png',dpi=150);plt.close(fig)
    report['status']='DIAGNOSTIC_COMPLETE'
    (args.output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report['runs'],indent=2))


if __name__=='__main__':main()
