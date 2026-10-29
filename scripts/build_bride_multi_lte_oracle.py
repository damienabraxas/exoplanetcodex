"""Build the actual held MULTI LTEPOP routine, without RT/EOS/opacity inputs.

Native DP include conversion follows make_mul23_dbl.csh. Only workspace capacity
is enlarged. Molecular QPART is an unreachable guard in this O-only wrapper.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
PIN='4d6aa3c39d634b38da772e469a5692435eba1354'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Authorized Mac oracle only')
    from pipeline.bride_reference import fingerprint
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.source_root,text=True).strip()!=PIN:
        raise ValueError('Unpinned source checkout')
    relative=Path('pymulti/files/multi/mul23/source');source=args.source_root/relative
    receipts=[]
    def read(name):
        path=source/name
        raw=subprocess.check_output(['git','show',f'{PIN}:{(relative/name).as_posix()}'],cwd=args.source_root)
        if raw!=path.read_bytes():raise ValueError('Modified native input: '+name)
        receipts.append(fingerprint(path));return raw.decode('latin1')
    text=read('mul23_sub.f');blocks=[];routines=[]
    for name in ('LTEPOP','GETWRD'):
        start=re.search(rf'^      SUBROUTINE {name}\b',text,re.M)
        end=re.search(r'^      END\s*$',text[start.start():],re.M)
        block=text[start.start():start.start()+end.end()]+'\n'
        routines.append(dict(name=name,sha256=hashlib.sha256(block.encode('latin1')).hexdigest()))
        blocks.append(block)
    args.output.mkdir(parents=True,exist_ok=False)
    includes=['CATOM','CATMOS','CCONST','CLU','CTRAN','CSLINE','CBCNST','C3D']
    for name in includes:
        # Same explicit REAL -> REAL(8) include conversion as the native DP script.
        converted=re.sub(r'\bREAL ', 'REAL(8) ',read(name))
        (args.output/name).write_text(converted)
    (args.output/'PREC').write_text(read('PREC.dp'))
    param=read('PARAM')
    if 'MK=30'not in param:raise ValueError('Unexpected native workspace dimensions')
    (args.output/'PARAM').write_text(param.replace('MK=30','MK=64'))
    blocks.append(read('mul23_bconst.f'))
    wrapper="""      SUBROUTINE BRIDE_LTE(NLEV,ND,T_IN,NE_IN,E_IN,G_IN,
     *                     STAGE_IN,FRAC,CONST_OUT)
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CCONST'
      INCLUDE 'CLU'
      INTEGER NLEV,ND,STAGE_IN(NLEV)
      REAL(8) T_IN(ND),NE_IN(ND),E_IN(NLEV),G_IN(NLEV)
      REAL(8) FRAC(NLEV,ND),CONST_OUT(8)
      NK=NLEV
      NDEP=ND
      ATOMID='O'
      LJOBLO=6
      HCE=HH*CC/EE
      EK=EE/BK
      DO I=1,NK
        EV(I)=E_IN(I)*HCE
        G(I)=G_IN(I)
        ION(I)=STAGE_IN(I)
      ENDDO
      DO K=1,NDEP
        TEMP(K)=T_IN(K)
        NE(K)=NE_IN(K)
        TOTN(K)=1.0D0
      ENDDO
      CALL LTEPOP
      DO K=1,NDEP
        DO I=1,NK
          FRAC(I,K)=NSTAR(I,K)
        ENDDO
      ENDDO
      CONST_OUT(1)=EE
      CONST_OUT(2)=HH
      CONST_OUT(3)=CC
      CONST_OUT(4)=BK
      CONST_OUT(5)=EM
      CONST_OUT(6)=UU
      CONST_OUT(7)=PI
      CONST_OUT(8)=EK
      END

      REAL(8) FUNCTION QPART(ELEMID,K)
      CHARACTER*4 ELEMID
      INTEGER K
      STOP 'Molecular branch forbidden in O-only LTE oracle'
      END
"""
    generated=args.output/'multi_lte_dp.f';generated.write_text('\n'.join(blocks)+wrapper,encoding='latin1')
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Existing gfortran required; no install fallback')
    library=args.output/'libmulti_lte.dylib'
    command=[compiler,'-std=legacy','-O2','-fcheck=bounds','-fPIC','-dynamiclib',str(generated),'-o',str(library)]
    run=subprocess.run(command,cwd=args.output,capture_output=True,text=True)
    (args.output/'compiler.log').write_text(run.stdout+run.stderr)
    if run.returncode:raise RuntimeError('Native LTE build failed; inspect compiler.log')
    report=dict(schema='codex.bride_multi_lte_oracle/1',source_commit=PIN,
        generator=fingerprint(__file__),sources=receipts,routines=routines,
        compiler=compiler,compiler_version=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],
        command=command,generated=fingerprint(generated),library=fingerprint(library),
        interpretation='Actual atomic LTEPOP plus GETWRD and BCONST; native DP includes; MK30->64 capacity only; O-only unreachable molecular guard. No full MULTI RT/EOS or exact grid producer identity.')
    (args.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Built native LTEPOP oracle',library)


if __name__=='__main__':main()
