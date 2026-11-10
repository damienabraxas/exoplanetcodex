"""Build pinned Mac TS baseline and output-only component-reporting derivative."""
import argparse
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
PIN='21269948750e2177fe32c5df38e1d9bb81363374'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac-only diagnostic build')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    ispec=require_codex_path('engines.ispec')
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=ispec,text=True).strip()!=PIN:
        raise ValueError('Changed iSpec source checkout')
    base=Path('synthesizer/turbospectrum');source=ispec/base/'source'
    makefile=ispec/base/'exec-gf/Makefile'
    text=makefile.read_text()
    block=text.split('babsma_luobj =',1)[1].split('eqwidt_luobj =',1)[0]
    filenames=re.findall(r'\$\(SOURCE\)/([\w.]+)',block)
    if len(filenames)<30:raise ValueError('Unexpected native build recipe')
    receipts=[]
    for path in [makefile]+[p for p in source.iterdir()if p.suffix in ('.f','.f90','.inc','.par')]:
        raw=subprocess.check_output(['git','show',f'{PIN}:{path.relative_to(ispec)}'],cwd=ispec)
        if raw!=path.read_bytes():raise ValueError('Modified original source: '+path.name)
        receipts.append(fingerprint(path))
    args.output.mkdir(parents=True,exist_ok=False)
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Existing gfortran required')
    report=dict(schema='codex.bride_ts_budget_build/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),source_commit=PIN,source_files=receipts,
        compiler=compiler,compiler_version=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],
        interpretation='NTO/NTPO print selection and diagnostic field formatting only in isolated derivative; source files protected. No rate or opacity formula edit.',builds=[])
    for name in ('baseline','reporting'):
        directory=args.output/name;shutil.copytree(source,directory)
        if name=='reporting':
            path=directory/'babsma.f';text=path.read_text()
            if text.count('      NTO=0')!=1:raise ValueError('Unexpected print-selection site')
            path.write_text(text.replace('      NTO=0','      NTO=1\n      NTPO(1)=1'))
            path=directory/'detabs.f';text=path.read_text()
            if text.count('36A6')!=1 or text.count('36F6.4')!=1:
                raise ValueError('Unexpected component print formats')
            # Native print path overflows its 36 fields with today's >36
            # components and reverts to the wavelength's REAL descriptor.
            path.write_text(text.replace('36A6','100A10').replace('36F6.4','100E16.8').replace('2E11.4,','2E24.16,'))
            text=path.read_text();anchor='      fakt(17)=presneutral(ntp,26)/(rhokt*part(15,1))'
            if text.count(anchor)!=1:raise ValueError('Unexpected Fe population reporting site')
            text=text.replace(anchor,anchor+"\n      write(iwrit,*) 'BRIDE_FE_STATE',ntp,t(ntp),pe(ntp),ro,\n     & part(15,1),presneutral(ntp,26),presion(ntp,26)")
            path.write_text(text)
            path=directory/'absko.f';text=path.read_text();anchor='! NH cont included !'
            if text.count(anchor)!=1:raise ValueError('Unexpected post-DETABS site')
            text=text.replace(anchor,anchor+"\n      write(iwrit,*) 'BRIDE_POST_STATE',ntp,xla(jp),\n     & cont/rosav(ntp),opnh*1.d-16*partryck(ntp,13)/\n     & (t(ntp)*1.38066d-16)/rosav(ntp),sumabs,sums ca".replace('sums ca','sumsca'))
            path.write_text(text)
        command=[compiler,'-fconvert=big-endian','-O3','-fno-automatic','-fbacktrace',
                 '-std=legacy','-mcmodel=small','-fallow-argument-mismatch',
                 '-o','babsma_lu',*filenames]
        run=subprocess.run(command,cwd=directory,capture_output=True,text=True,timeout=90)
        (directory/'compiler.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError('Build failed: '+name)
        report['builds'].append(dict(name=name,command=command,program=fingerprint(directory/'babsma_lu'),
            babsma_source=fingerprint(directory/'babsma.f'),detabs_source=fingerprint(directory/'detabs.f'),
            absko_source=fingerprint(directory/'absko.f'),
            compiler_log=fingerprint(directory/'compiler.log')))
    (args.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Built baseline and reporting TS binaries:',args.output)


if __name__=='__main__':main()
