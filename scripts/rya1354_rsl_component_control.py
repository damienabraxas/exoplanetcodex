"""Exercise pinned native GETRSL; compile/run only on Sirius, no atomic solve."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess

from rya1354_solver_smoke import sha

FIXTURE='''program component_control
use orb_C, only: NW, NP, NH
implicit none
integer :: indices(4), count, i
NW=4
NP(1:4)=[2,2,2,3]
NH(1:4)=[character(len=2) :: 's ', 'p-', 'p ', 's ']
do i=1,3
  call GETRSL(indices,count)
  write(*,'(*(I0,1X))') count,indices(1:count)
end do
end program component_control
'''
INPUT='2s,2p\n2s,2p-,2p\n2p*\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Compile native controls on Sirius')
    root=args.grasp_root.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'driver.py')
    source=out/'component_control.f90';source.write_text(FIXTURE)
    compiler=Path(shutil.which('gfortran')).resolve()
    libraries=[root/'build/lib'/name for name in ('lib9290.a','libmod.a')]
    receipt=dict(schema='rya1354.native-rsl-component-control.v8',status='RUNNING',
        captured_utc=datetime.now(timezone.utc).isoformat(),scientific_admission='HOLD',
        compiler_sha256=sha(compiler),compiler_version=subprocess.check_output([str(compiler),'--version'],text=True),
        source_sha256=sha(root/'src/lib/lib9290/getrsl.f90'),
        libraries={p.name:sha(p) for p in libraries},atomic_prediction=False)
    try:
        cmd=[str(compiler),'-I'+str(root/'build/modules/mod'),str(source),*(str(p) for p in libraries),'-o',str(out/'component_control')]
        receipt['compile_command']=cmd
        built=subprocess.run(cmd,cwd=out,capture_output=True,text=True)
        (out/'compile.stdout').write_text(built.stdout);(out/'compile.stderr').write_text(built.stderr)
        built.check_returncode()
        (out/'component.stdin').write_text(INPUT)
        tested=subprocess.run([str(out/'component_control')],cwd=out,input=INPUT,capture_output=True,text=True,timeout=30)
        (out/'component.stdout').write_text(tested.stdout);(out/'component.stderr').write_text(tested.stderr)
        tested.check_returncode()
        actual=[[int(v) for v in line.split()] for line in tested.stdout.splitlines()]
        if actual!=[[2,1,3],[3,1,2,3],[2,2,3]]:
            raise ValueError('Native GETRSL component selection differs')
        receipt.update(status='NATIVE_COMPONENT_SELECTION_VERIFIED_NO_ATOMIC_PREDICTION',selections=actual,
            interpretation='2p excludes2p-; explicit2p-,2p or2p* includes both. GETOLD applies METHOD1 only to selected spectroscopic subshells; correlation default METHOD3 is not node-checked in the same way.')
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        receipt['files']={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
