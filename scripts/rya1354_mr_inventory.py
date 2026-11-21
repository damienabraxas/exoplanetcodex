"""Generate and audit a bounded oxygen MR-SD inventory, without diagonalization.

Reference families are the n<=4 subset of Li2023 Table1 MCDHF references.
This is a new reduced model inventory, not reproduction of their calculation.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import re
import shutil

from rya1354_solver_smoke import run, sha


def references(parity, core, multireference):
    prefix='1s(2,i)2s(2,'+('*' if core else 'i')+')'
    if parity=='even':
        tails=['2p(4,*)']
        if multireference:
            tails += ['2p(3,*)'+orb+'(1,*)' for orb in ('3p','4p','4f')]
    elif parity=='odd':
        tails=['2p(3,*)3s(1,*)']
        if multireference:
            tails += ['2p(3,*)'+orb+'(1,*)' for orb in ('4s','3d','4d')]
    else:
        raise ValueError('Unknown parity')
    return [prefix+tail for tail in tails]


def inventory(text, parity, core):
    """Check native three-record CSFs, eight electrons, frozen1s and parity."""
    lines=text.splitlines()
    if lines[0]!='Core subshells:' or lines[1].strip() or lines[4]!='CSF(s):':
        raise ValueError('Unexpected native CSF header/core representation')
    body=[line for line in lines[5:] if not line.startswith(' *') and not line.startswith('*')]
    if not body or len(body)%3:
        raise ValueError('Incomplete native CSF records')
    blocks=Counter(); seen=set()
    for i in range(0,len(body),3):
        record=tuple(body[i:i+3])
        if record in seen:
            raise ValueError('Duplicate native CSF')
        seen.add(record)
        fields=re.findall(r'(\d+[spdfgh]-?)\s*\(\s*(\d+)\)',record[0])
        occ={orb:int(number) for orb,number in fields}
        if (len(occ)!=len(fields) or sum(occ.values())!=8 or occ.get('1s')!=2
                or (not core and occ.get('2s')!=2)):
            raise ValueError('Electron/core occupation differs from declared model')
        odd=sum('spdfgh'.index(re.search('[spdfgh]',orb)[0])*number for orb,number in occ.items())%2
        label=record[2].split()[-1]
        expected='+' if parity=='even' else '-'
        if odd!=(parity=='odd') or not label.endswith(expected):
            raise ValueError('Wrong native occupation/block parity')
        blocks[label]+=1
    return dict(CSFs=sum(blocks.values()),blocks=dict(blocks),duplicate_CSFs=0,
                electron_count=8,closed_1s=True,closed_2s=not core)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic native generators belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.mr-inventory.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],cases={},
                 source='Li2023 arXiv2304.02310 Table1 MCDHF reference families restricted to n<=4; no fitted strengths or observed-energy inputs',
                 scope='CSF_GENERATION_ONLY_NO_WAVEFUNCTION_ENERGY_RATE_OR_GRADE')
    try:
        for parity in ('even','odd'):
            for core in (False,True):
                for mr in (False,True):
                    name=f'{parity}_'+('CV' if core else 'VV')+'_'+('MR' if mr else 'SR')
                    out=args.out/name;out.mkdir()
                    configs=references(parity,core,mr)
                    # Match all ground J=0,1,2 blocks and odd target J=1,2 blocks.
                    js='0,4' if parity=='even' else '2,4'
                    stdin='*\n0\n'+'\n'.join(configs)+'\n\n4s,4p,4d,4f\n'+js+'\n2\nn\n'
                    receipt['steps'].append(run(args.grasp_bin/'rcsfgenerate',stdin,out,1800))
                    if receipt['steps'][-1]['semantic_errors']:
                        raise ValueError('Native generator refusal')
                    receipt['cases'][name]=dict(references=configs,active='4s,4p,4d,4f',twoJ=js,
                                                **inventory((out/'rcsf.out').read_text(),parity,core))
                    (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='INVENTORY_COMPLETED_NO_ATOMIC_PREDICTION'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
