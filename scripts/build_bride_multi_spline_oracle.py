"""Build a local diagnostic oracle from pinned, already-held MULTI source.

No downloads/installations or full MULTI build. Original files are read-only.
Only TAUTSP, PPVALU and INTERV are extracted; DP conversion follows the source's
PREC mechanism. This does not identify the producer of a departure grid.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Authorized Mac oracle only')
    from pipeline.bride_reference import fingerprint
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.source_root,text=True).strip()
    if commit!='4d6aa3c39d634b38da772e469a5692435eba1354':raise ValueError('Unpinned MULTI source')
    relative='pymulti/files/multi/mul23/source/mul23_opacu.f'
    source=args.source_root/relative
    pinned=subprocess.check_output(['git','show',f'{commit}:{relative}'],cwd=args.source_root)
    if pinned!=source.read_bytes():raise ValueError('Modified MULTI source')
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Held gfortran required; no installation fallback')
    text=pinned.decode();blocks=[];receipts=[]
    for name in ('TAUTSP','PPVALU','INTERV'):
        start=re.search(rf'^      (?:SUBROUTINE|FUNCTION) {name}\b',text,re.M)
        if start is None:raise ValueError(f'Missing native routine {name}')
        end=re.search(r'^      END\s*$',text[start.start():],re.M)
        if end is None:raise ValueError('Missing routine terminator')
        block=text[start.start():start.start()+end.end()]+'\n'
        receipts.append(dict(name=name,source_sha256=hashlib.sha256(block.encode()).hexdigest()))
        block=block.replace("      INCLUDE 'PREC'","      IMPLICIT REAL(8) (A-H,O-Z)")
        blocks.append(block)
    args.output.mkdir(parents=True,exist_ok=False)
    generated=args.output/'multi_spline_dp.f';generated.write_text('\n'.join(blocks))
    library=args.output/'libmulti_spline.dylib'
    command=[compiler,'-std=legacy','-O2','-fPIC','-dynamiclib',str(generated),'-o',str(library)]
    result=subprocess.run(command,capture_output=True,text=True)
    (args.output/'compiler.log').write_text(result.stdout+result.stderr)
    if result.returncode:raise RuntimeError('Native oracle compilation failed; see compiler.log')
    report=dict(schema='codex.bride_multi_spline_oracle/1',source_commit=commit,
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source=fingerprint(source),source_relative=relative,routines=receipts,
        generator=fingerprint(__file__),compiler=compiler,
        compiler_version=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],
        command=command,precision='Unmodified routines with PREC expanded as IMPLICIT REAL(8), matching DP mechanism',
        generated=fingerprint(generated),library=fingerprint(library),
        producer_identity='This held source is not proven identical to the departure-grid producer')
    (args.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
