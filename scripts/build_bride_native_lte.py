"""Isolated unmodified native TS opacity and LTE-spectrum build for Mac diagnostics."""
import argparse,json,platform,re,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac-only build')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    plan=json.loads((ROOT/'config/experiments/rya1281_native_molecules.json').read_text());pin=plan['source_commit']
    engine=require_codex_path('engines.ispec');base=engine/'synthesizer/turbospectrum';source=base/'source';recipe=base/'exec-gf/Makefile'
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=engine,text=True).strip()!=pin:raise ValueError('Wrong native source pin')
    receipts=[]
    for path in [recipe]+[p for p in source.iterdir()if p.suffix in ('.f','.f90','.inc','.par')]:
        if path.read_bytes()!=subprocess.check_output(['git','show',pin+':'+str(path.relative_to(engine))],cwd=engine):raise ValueError('Changed original source '+path.name)
        receipts.append(dict(relative=str(path.relative_to(engine)),receipt=fingerprint(path)))
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Existing compiler required; no install')
    a.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema='codex.bride_native_lte_build/1',scientific_admission='HOLD',generator=fingerprint(__file__),source_commit=pin,
        sources=receipts,compiler=compiler,compiler_version=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],builds=[])
    text=recipe.read_text()
    for target,following in [('babsma_lu','eqwidt_lu'),('bsyn_lu',None)]:
        block=text.split(target+'obj =',1)[1]
        block=block.split(following+'obj =',1)[0]if following else block.split('#clean',1)[0]
        names=re.findall(r'\$\(SOURCE\)/([\w.]+)',block)
        if len(names)<30 or len(set(names))!=len(names):raise ValueError('Unexpected native build recipe')
        # Make builds modules before users; retain exactly the same source files.
        names=[n for n in names if n.endswith('.f90')]+[n for n in names if not n.endswith('.f90')]
        directory=a.output/target;shutil.copytree(source,directory)
        command=[compiler,'-fconvert=big-endian','-O3','-fno-automatic','-fbacktrace','-std=legacy','-mcmodel=small','-fallow-argument-mismatch','-o',target,*names]
        r=subprocess.run(command,cwd=directory,capture_output=True,text=True,timeout=plan['build_seconds'])
        (directory/'compiler.log').write_text(r.stdout+r.stderr)
        if r.returncode:raise RuntimeError('Native build failed: '+target)
        report['builds'].append(dict(target=target,command=command,program=fingerprint(directory/target),compiler_log=fingerprint(directory/'compiler.log')))
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(a.output)


if __name__=='__main__':main()
