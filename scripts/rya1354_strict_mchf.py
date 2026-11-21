"""Isolated native MCHF orbital-stop diagnostic; never modify upstream objects."""
import argparse
from datetime import datetime,timezone
import difflib
import json
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess

from rya1354_solver_smoke import sha


def strict_source(text):
    energy='      CONV = ECONV .or. DP1 .LE. Z2'
    stop='      IF ( LAST ) GO TO 17'
    if text.count(energy)!=1 or text.count(stop)!=1:
        raise ValueError('Pinned MCHF stopping statements differ')
    return text.replace(energy,'      CONV = DP1 .LE. Z2').replace(stop,
        "      IF (DP1 .LE. Z2) THEN\n"
        "         WRITE(OUT,*) 'RYA1354_ORBITAL_STOP',DP1,Z2\n"
        "         GO TO 17\n"
        "      ENDIF\n"
        "      LAST = .FALSE.")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--atsp',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Native solver builds belong on Sirius')
    args.atsp=args.atsp.resolve()
    args.out=args.out.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'build_driver.py')
    src=args.atsp/'src/appl/mchf'
    original=(src/'scf.f').read_text()
    changed=strict_source(original)
    (args.out/'scf.original.f').write_text(original)
    (args.out/'scf.f').write_text(changed)
    (args.out/'strict.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),fromfile='scf.original.f',tofile='scf.f')))
    make=(src/'Makefile').read_text().replace('\\\n',' ')
    objects=[]
    for variable in ('USPOBJ','SYSOBJ'):
        body=re.search(r'^'+variable+r'\s*=([^\n]*)',make,re.M)[1]
        names=shlex.split(body.replace('\\\n',' ').replace('${FC_MALLOC}','LINUX'))
        objects.extend(str(args.out/'scf.o') if name=='scf.o' else str(src/name) for name in names)
    compile=['gfortran','-O2','-fno-automatic','-fcray-pointer','-std=legacy','-fallow-argument-mismatch','-c',str(args.out/'scf.f'),'-o',str(args.out/'scf.o')]
    link=['gfortran','-o',str(args.out/'mchf'),'-static-libgfortran','-L'+str(args.atsp/'lib'),*objects,
          '-ldvdson','-lrad','-lcom','-lnet','/usr/lib/x86_64-linux-gnu/liblapack.so.3','/usr/lib/x86_64-linux-gnu/libblas.so.3']
    with (args.out/'build.log').open('x') as log:
        for command in (compile,link):
            subprocess.run(command,cwd=args.out,stdout=log,stderr=log,check=True)
    receipt=dict(schema='rya1354.strict-mchf-build.v1',started_utc=datetime.now(timezone.utc).isoformat(),
                 compile=compile,link=link,original_executable_sha256=sha(args.atsp/'bin/mchf'),
                 upstream_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.atsp,text=True).strip(),
                 upstream_dirty=subprocess.check_output(['git','status','--porcelain'],cwd=args.atsp,text=True),
                 change='Require native DP1<=Z2 at the actual function-stop branch; disable energy-only acceptance; print DP1/Z2 at stop. Orbital/CI update equations unchanged.',
                 files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in args.out.iterdir() if p.is_file()})
    (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
