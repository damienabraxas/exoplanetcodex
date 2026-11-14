"""Both-gauge minimal FAC rates with unpruned diagnostic queries; Sirius only."""
import argparse
import difflib
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess

from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_coherent import input_text as numerical_input
from rya1354_fac_controls import diagnostics
from rya1354_fac_oxygen_levels import inspect_inventory
from rya1354_fac_mapping import verify_transport
from rya1354_fac_recoupling_audit import review_case

OPS=(('E1',-1,'ground'),('M2',2,'ground'),('E3',-3,'ground'),
     ('M1',1,'quintet'),('E2',-2,'quintet'),('M3',3,'quintet'))


def read_files(root):
    return {p.name:p.read_bytes() for p in root.iterdir() if p.is_file()}


def verify_receipt(root):
    receipt=json.loads((root/'receipt.json').read_text())
    for name,r in receipt['files'].items():
        p=root/name
        if sha(p)!=r['sha256'] or p.stat().st_size!=r['bytes']:
            raise ValueError('FAC source artifact changed: '+name)
    return receipt


def rate_input(multiplier,scale,gauge,mapping,levels):
    upper,=mapping['targets']['3S'];quintet,=mapping['targets']['5S']
    by_id={r['index']:r for r in levels}
    if set(by_id)!=set(range(15)):raise ValueError('FAC rate level IDs differ')
    even=sorted(int(i) for i,r in mapping['levels'].items() if r['parity']=='+')
    lower=[r['index'] for r in levels if r['energy_eV']<by_id[upper]['energy_eV']]
    if set(lower)!=set(even+[quintet]) or len(even)!=5:
        raise ValueError('Minimal FAC lower-state inventory differs')
    text=numerical_input(multiplier,scale)+f'SetTransitionOptions({gauge},0)\nSetTransitionCut(0,0)\n'
    requests=[]
    for op,m,dest in OPS:
        lows=even if dest=='ground' else [quintet]
        # Pinned FAC requires a single database filename per DB type/session.
        text+=f"TransitionTable('rates.tr.b',{lows},[{upper}],{m})\n"
        for low in lows:
            jl=by_id[low]['two_J'];ju=by_id[upper]['two_J'];rank=2*abs(m)
            if abs(jl-ju)<=rank<=jl+ju:
                text+=f"RawTransition('raw.tsv',{low},{upper},{m})\n"
                requests.append(dict(operator=op,multipole=m,lower=low,upper=upper))
    if len(requests)!=13:raise ValueError('FAC declared channel count differs')
    text+="PrintTable('rates.tr.b','rates.tr',1)\n"
    return text,requests


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('source','bootstrap','controls','mapping','probe','out'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic builds/calculations belong on Sirius')
    for name,value in vars(args).items():setattr(args,name,value.resolve())
    build=verify_receipt(args.bootstrap);control=verify_receipt(args.controls);mapped=verify_receipt(args.mapping)
    if (subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.source,text=True).strip()!=PIN
            or subprocess.check_output(['git','diff','--name-only'],cwd=args.source,text=True).strip()
            or sha(args.source/'sfac/sfac')!=build['binary']['sha256']
            or control['status']!='FAC_COHERENT_CONTROLS_COMPLETED_REVIEW_REQUIRED'
            or mapped['status']!='FAC_ANGULAR_TRANSFORM_COMPLETED_REVIEW_REQUIRED'
            or mapped['source_receipt_sha256']!=sha(args.controls/'receipt.json')):
        raise ValueError('FAC rate source qualification differs')
    prepared={}
    for case,multiplier,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
        source=read_files(args.controls/case);converted=read_files(args.mapping/case)
        if (source['input.sf'].decode()!=numerical_input(multiplier,scale)
                or not diagnostics(source['sfac.stdout'].decode(),multiplier,256)['printed_rule_verified']
                or converted['state.c']!=source['o_export_000.c'] or converted['state.cm']!=source['o_export_000.cm']):
            raise ValueError('FAC numerical/mapping input differs')
        verify_transport(source['o.basis'].decode(),source['o_export_000.c'].decode(),source['o_export_000.cm'])
        prepared[case]=(source,review_case(source,converted))
    args.out.mkdir(parents=True,exist_ok=False);args.probe.mkdir(parents=True,exist_ok=False)
    drivers=args.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*'):
        if p.is_file():shutil.copyfile(p,drivers/p.name)
    receipt=dict(schema='rya1354.fac-model-rates.v21',status='RUNNING',pin=PIN,cases={},steps=[],
        started_utc=datetime.now(timezone.utc).isoformat(),source_receipt_sha256=sha(args.controls/'receipt.json'),
        mapping_receipt_sha256=sha(args.mapping/'receipt.json'),scientific_admission='HOLD',
        exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Minimal FAC model both-gauge relativistic one-photon channel experiment. Isolated SFAC command queries unchanged library before writer pruning. Nonzero native status is never silently zero-filled. No accuracy/calibrated grade.')
    def execute(label,command,cwd,timeout):
        with (cwd/(label+'.stdout')).open('xb') as out,(cwd/(label+'.stderr')).open('xb') as err:
            p=subprocess.run(command,cwd=cwd,stdout=out,stderr=err,timeout=timeout,check=False)
        receipt['steps'].append(dict(label=label,command=command,cwd=str(cwd),returncode=p.returncode))
        if p.returncode:raise RuntimeError(label+' exit '+str(p.returncode))
    try:
        original=(args.source/'sfac/sfac.c').read_text()
        snippet=Path(__file__).with_name('rya1354_fac_raw_transition.c').read_text()
        marker='#include "global.h"\n';registration='  {"TransitionTable", PTransitionTable, METH_VARARGS},'
        if original.count(marker)!=1 or original.count(registration)!=1:
            raise ValueError('Pinned SFAC insertion points differ')
        patched=original.replace(marker,marker+'\n'+snippet+'\n',1).replace(registration,
            '  {"RawTransition", PRawTransition, METH_VARARGS},\n'+registration,1)
        (args.probe/'sfac.c').write_text(patched);(args.probe/'upstream.sfac.c').write_text(original)
        (args.probe/'observability.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),fromfile='pinned/sfac.c',tofile='diagnostic/sfac.c')))
        for name in ('stoken.c','stoken.h'):shutil.copyfile(args.source/'sfac'/name,args.probe/name)
        shutil.copyfile(args.source/'LICENSE',args.probe/'LICENSE')
        shutil.copyfile(args.source/'sfac/Makefile',args.probe/'source.Makefile')
        library=args.source/'libfac.a';before=sha(library)
        execute('make',['make','-f',str(args.source/'sfac/Makefile'),'sfac'],args.probe,180)
        binary=args.probe/'sfac'
        if sha(library)!=before:raise ValueError('FAC mathematical library changed')
        receipt['probe']=dict(path=str(binary),sha256=sha(binary),bytes=binary.stat().st_size,libfac_sha256=before,
                              upstream_binary_sha256=build['binary']['sha256'])
        proof=args.out/'frontend_evidence';proof.mkdir()
        for p in args.probe.iterdir():
            if p.is_file() and p.name!='sfac':shutil.copyfile(p,proof/p.name)
        for name in ('faclib/transition.c','faclib/dbase.c','faclib/consts.h','doc/manual.tex'):
            dest=args.out/'source_documents'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.source/name,dest)
        for case,multiplier,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
            source,mapping=prepared[case];levels=inspect_inventory(source['o.lev'].decode())
            for gauge_name,gauge in [('B',2),('C',1)]:
                name=case+'_'+gauge_name;out=args.out/name;out.mkdir()
                text,requests=rate_input(multiplier,scale,gauge,mapping,levels)
                (out/'input.sf').write_text(text)
                execute('sfac',[str(binary),'input.sf'],out,300)
                d=diagnostics((out/'sfac.stdout').read_text(),multiplier,256)
                if not d['printed_rule_verified']:raise ValueError('FAC rate numerical guard failed')
                for file in ('o_export_000.c','o_export_000.cm'):
                    if (out/file).read_bytes()!=source[file]:raise ValueError('FAC regenerated state differs from mapping')
                rows=[l for l in (out/'raw.tsv').read_text().splitlines() if not l.startswith('#') and l.strip()]
                if len(rows)!=13:raise ValueError('FAC raw channel queries incomplete')
                receipt['cases'][name]=dict(gauge=gauge,mode=0,requests=requests,diagnostics=d,
                                            source_case=case,raw_records=len(rows))
                (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_MODEL_RATES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
