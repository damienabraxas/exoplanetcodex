"""Matched DC/static-Breit CI on qualified n5 orbitals; no failed-wave input."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil

from rya1354_joint_outer import qualified_source
from rya1354_lifetime_ladder import STATES,identity
from rya1354_solver_smoke import run,sha
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order

PHYSICS={'DC':'n\nn\nn\nn\nn\nn\n','Breit':'y\ny\n1.d-6\nn\nn\nn\nn\n'}


def controls(state,case):
    roots=STATES[state]['roots'].removesuffix('5\n') if state=='ground' else STATES[state]['roots']
    return 'y\nstate\n'+PHYSICS[case]+roots


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('source','n4-source','parent','n4-core','core','grasp-bin','out'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for name,value in vars(args).items():
        setattr(args,name,value.resolve())
    audit=qualified_source(args)
    args.out.mkdir(parents=True,exist_ok=False)
    drivers=args.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):
        shutil.copyfile(p,drivers/p.name)
    receipt=dict(schema='rya1354.qualified-n5-ci.v16',status='RUNNING',n=5,
        started_utc=datetime.now(timezone.utc).isoformat(),source_receipt_sha256=sha(args.source/'receipt.json'),
        source_audit=audit,cases=list(PHYSICS),steps=[],scientific_admission='HOLD',
        exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,
        model='Same qualified sequential n5 full2s-SD/DC source orbitals, CSFs and nucleus for all cases. Native RCI DC control and same-space static-limit Breit. No orbital optimization or failed joint/core-release waves. Native internal RCI grid retained; no observed-energy substitution, gauge selection or accuracy claim. Historical n3 component limitation retained.')
    try:
        for case in PHYSICS:
            root=args.out/case;root.mkdir();selected={}
            for state in STATES:
                seed=args.source/('outer_n5_'+state+'_verify');out=root/state;out.mkdir()
                for name in ('state.c','state.w','isodata'):
                    shutil.copyfile(seed/name,out/name)
                shutil.copyfile(seed/'state.lsj.lbl',out/'source.lsj.lbl')
                identity((out/'source.lsj.lbl').read_text(),state)
                shutil.copyfile(out/'state.c',out/'rcsf.inp')
                for exe,stdin in [('rangular','y\n'),('rci',controls(state,case)),('jj2lsj','state\ny\nn\ny\n')]:
                    receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,43200))
                    if receipt['steps'][-1]['semantic_errors']:
                        raise ValueError('Native CI semantic refusal')
                for name in ('state.c','state.w','isodata'):
                    if (out/name).read_bytes()!=(seed/name).read_bytes():
                        raise ValueError('Matched CI source CSF/orbital/nucleus changed')
                if not (out/'state.cm').is_file():
                    raise ValueError('Native CI mixing missing')
                identity((out/'state.lsj.lbl').read_text(),state);selected[state]=out
            require_triplet_quintet_order((selected['triplet']/'state.lsj.lbl').read_text(),(selected['quintet']/'state.lsj.lbl').read_text())
            for name,lower,operators,same in [('to_ground','ground','E1,M2,E3',False),('to_5S','quintet','M1,E2,M3',True)]:
                receipt['steps']+=transitions(args.grasp_bin,selected[lower],selected['triplet'],root/name,operators,same,ci=True)
            receipt['last_completed_case']=case
            (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
