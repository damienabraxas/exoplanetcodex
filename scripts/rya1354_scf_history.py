"""Native SCF history and timestamped log snapshots; never state qualification."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import shutil


def history(text):
    """Keep unfinished outer cycles distinct; ignore inner `Iteration number:`."""
    matches=list(re.finditer(r'^\s*Iteration number\s+(\d+)\s*$',text,re.M))
    records=[]
    for index,match in enumerate(matches):
        end=matches[index+1].start() if index+1<len(matches) else len(text)
        body=text[match.end():end];rows={}
        for line in body.splitlines():
            fields=line.split()
            if len(fields)<11 or not re.fullmatch(r'\d+[spdfgh]-?',fields[0]) or fields[2] not in ('1','2','3','4'):
                continue
            value=float(fields[4].replace('D','E'))
            if not math.isfinite(value) or value<0:
                raise ValueError('Invalid native orbital residual')
            rows[fields[0]]=dict(residual=value,method=int(fields[2]),nodes=int(fields[-1]))
        energy=re.findall(r'Weighted average energy of these levels\s*=\s*([\d.DEd+\-]+)',body)
        complete=bool(rows and energy and 'Generalised occupation numbers:' in body)
        records.append(dict(iteration=int(match.group(1)),complete_printed_cycle=complete,
            orbitals=rows,max_reported_residual=max((r['residual'] for r in rows.values()),default=None),
            weighted_energy_hartree=float(energy[-1].replace('D','E')) if energy else None))
    if not records or any(a['iteration']>=b['iteration'] for a,b in zip(records,records[1:])):
        raise ValueError('Missing or nonmonotonic outer SCF cycle inventory')
    return records


def snapshot(run,out):
    """Capture a read-only prefix of live logs; timing is asynchronous, not atomic."""
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Capture native live logs on Sirius')
    if (run/'receipt.json').exists():
        raise ValueError('Terminal receipt exists; collect the full native campaign instead')
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'driver.py')
    for name in ('rmcdhf.stdin','rmcdhf.stdout','rmcdhf.stderr'):
        (out/name).write_bytes((run/'ground'/name).read_bytes())
    rows=history((out/'rmcdhf.stdout').read_text())
    capture=dict(schema='rya1354.running-scf-log-snapshot.v8',captured_utc=datetime.now(timezone.utc).isoformat(),
        source=str(run),scope='POSTHOC_ASYNCHRONOUS_LOG_PREFIX_NOT_TERMINAL_RECEIPT_NOT_ACCEPTED_STATE',
        scientific_admission='HOLD',adoption_allowed=False,accepted_lifetime=None,
        latest_printed_outer_iteration=rows[-1]['iteration'],history=rows,
        files={p.name:dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size)
               for p in sorted(out.iterdir()) if p.is_file()})
    (out/'capture.json').write_text(json.dumps(capture,indent=2)+'\n')


def plot(bundle,out):
    """Plot captured native residuals; the criterion is not an accuracy band."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from rya1354_native_parts import load_parts
    files,_=load_parts(bundle)
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for name,label in (('nodebudget','Original selection; NSOLV24'),('pcomponents','Explicit both2p; NSOLV24')):
        prefix=f'evidence/{name}_log_snapshot_v8/'
        capture=json.loads(files[prefix+'capture.json'])
        rows=history(files[prefix+'rmcdhf.stdout'].decode())
        if rows!=capture['history']:
            raise ValueError('Capture history differs from native log prefix')
        complete=[r for r in rows if r['complete_printed_cycle']]
        axes[0].semilogy([r['iteration'] for r in complete],[r['max_reported_residual'] for r in complete],label=label)
        first=complete[:7]
        axes[1].semilogy([r['iteration'] for r in first],[r['max_reported_residual'] for r in first],'-o',label=label)
    for ax in axes:
        ax.axhline(.0625**6,color='black',linestyle=':',label='Original orbital criterion')
        ax.set(xlabel='Completed printed SCF cycle',ylabel='Largest last-reported orbital residual')
        ax.grid(alpha=.2)
    axes[0].set_title('Captured histories; neither state qualified')
    axes[1].set_title('First seven completed cycles')
    axes[0].legend(fontsize=8)
    fig.suptitle('O I n4 ground optimization: diagnostic snapshots, not accuracy estimates')
    fig.savefig(out,dpi=150,metadata={'Description':'Asynchronous captured native log prefixes; no terminal state, lifetime, error bar or accuracy claim.'})
    plt.close(fig)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run',type=Path)
    mode.add_argument('--plot-bundle',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.run:snapshot(args.run,args.out)
    else:plot(args.plot_bundle,args.out)
