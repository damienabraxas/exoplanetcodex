"""v44: signed native E1 component records for the sealed v41 CI/rate models.

Re-runs the EXACT v41 rate input per case/gauge with a read-only frontend that
appends RawE1Components calls immediately after the E1 RawTransition calls, i.e.
in the actual native photon-grid context. All v41 physical outputs, raw rates and
the independent rate audit must reproduce byte-for-byte/exactly. No new SCF, CI,
Hamiltonian, rates, lifetime, uncertainty or grade.
"""
import argparse
from datetime import datetime,timezone
import difflib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_rates import verify_receipt
from rya1354_fac_radial_angular_remaining_metadata import hold
from rya1354_fac_radial_recouple import CASES
from rya1354_fac_radial_rates import rate_input,run_native
from rya1354_fac_radial_ci_v36 import review as ci_review
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_radial_operator import prepare as operator_prepare
from rya1354_fac_e1_components_audit import parse,verify,grid_context

FRONTEND='toolchains/fac-config-probe-v27'
STATE_FILES=('o_export_000.c','o_export_000.cm','source.pot','ci.pot','raw.tsv')
# Measured null: v36 CI vs its v41 regeneration differ ONLY in input.sf/native.time/o.lev/o.lev.b
# (session stamp/timing). rates.tr(.b) carry the same stamp. Everything else must be byte-identical.
UNSTABLE={'input.sf','native.time','o.lev','o.lev.b','rates.tr','rates.tr.b','components.tsv'}
FORBIDDEN=('AngZ','angz','mix_cut','MixCut','PrepAngular','max_nbs','SetTransitionMaxE','SetTransitionMaxM')


def input_text(scale,gauge,mapping):
    """Exact v41 input plus component calls after the last E1 RawTransition."""
    text,requests,lower=rate_input(scale,gauge,mapping)
    if any(x in text for x in FORBIDDEN):raise ValueError('Native angular/mixing cut controls differ from defaults')
    upper,=mapping['targets']['3S']
    e1=[r['lower'] for r in requests if r['multipole']==-1]
    if not e1 or len(set(e1))!=len(e1):raise ValueError('E1 channel inventory differs')
    last=f"RawTransition('raw.tsv',{e1[-1]},{upper},-1)\n"
    if text.count(last)!=1:raise ValueError('E1 insertion point differs')
    calls=''.join(f"RawE1Components('components.tsv',{low},{upper},{gauge})\n" for low in e1)
    head,tail=text.split(last)
    if not tail.startswith('TransitionTable('):raise ValueError('E1 photon-grid context would not be the native one')
    return head+last+calls+tail,requests,e1


def prepare(root):
    base=operator_prepare(root)
    parent=root/'runs/fac_radial_rates_v41';receipt=verify_receipt(parent);hold(receipt)
    index=json.loads((parent/'index.json').read_text())
    probe=index['probe'];frontend=root/FRONTEND
    if Path(probe['path']).resolve()!=(frontend/'sfac').resolve() or sha(frontend/'sfac')!=probe['sha256']:
        raise ValueError('v41 rate frontend identity differs')
    v43=root/'runs/fac_radial_operator_v43';r43=verify_receipt(v43);hold(r43)
    if r43['status']!='FAC_RADIAL_OPERATOR_COMPLETED_REVIEW_REQUIRED':raise ValueError('v43 radial checkpoint differs')
    result=dict(schema='rya1354.e1-components-input.v44',parent_receipt_sha256=base['parent_receipt_sha256'],
        parent_index_sha256=base['parent_index_sha256'],v43_receipt_sha256=sha(v43/'receipt.json'),fac_pin=index['fac_pin'],
        rate_frontend=dict(path=FRONTEND,sfac_c_sha256=sha(frontend/'sfac.c'),sfac_sha256=probe['sha256'],
                           libfac_sha256=probe['libfac_sha256'],upstream_binary_sha256=probe['upstream_binary_sha256']),
        snippet_sha256=sha(Path(__file__).with_suffix('.c')),potential_layout=index['potential_layout'],
        constants=base['constants'],cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    for case in CASES:
        entry=index['cases'][case];mapping=entry['mapping'];scale=100 if case.endswith('100') else 10
        upper,=mapping['targets']['3S'];gauges={}
        for suffix,gauge in [('B',2),('C',1)]:
            src=parent/(case+'_'+suffix);text,requests,e1=input_text(scale,gauge,mapping)
            if (src/'input.sf').read_text()!=rate_input(scale,gauge,mapping)[0]:raise ValueError('Parent input differs')
            gauges[suffix]=dict(parent=str(src.relative_to(root)),gauge=gauge,e1_lower=e1,input_sha256=sha_text(text),
                                parent_raw_sha256=sha(src/'raw.tsv'),parent_audit=receipt['cases'][case+'_'+suffix]['audit'])
        result['cases'][case]=dict(scale=scale,upper=upper,reference_lower=base['cases'][case]['lower'],
            qualified_CI=entry['qualified_CI'],gauges=gauges)
    return result


def sha_text(text):
    import hashlib
    return hashlib.sha256(text.encode()).hexdigest()


def build(root,source,probe,evidence,index):
    probe.mkdir(parents=True,exist_ok=False);evidence.mkdir(parents=True,exist_ok=False)
    frontend=root/FRONTEND;fe=index['rate_frontend']
    if sha(frontend/'sfac.c')!=fe['sfac_c_sha256'] or sha(Path(__file__).with_suffix('.c'))!=index['snippet_sha256']:
        raise ValueError('Frontend source or snippet differs from frozen index')
    original=(frontend/'sfac.c').read_text();snippet=Path(__file__).with_suffix('.c').read_text()
    marker='#include "global.h"\n';registration='  {"TransitionTable", PTransitionTable, METH_VARARGS},'
    if original.count(marker)!=1 or original.count(registration)!=1:raise ValueError('Probe insertion points differ')
    patched=original.replace(marker,marker+'\n'+snippet+'\n',1).replace(registration,
        '  {"RawE1Components", PRawE1Components, METH_VARARGS},\n'+registration,1)
    (probe/'sfac.c').write_text(patched);(probe/'rate_frontend.sfac.c').write_text(original)
    (probe/'observability.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),
        fromfile='v27-rate-frontend/sfac.c',tofile='e1-components/sfac.c')))
    for n in ('stoken.c','stoken.h'):shutil.copyfile(frontend/n,probe/n)
    shutil.copyfile(source/'sfac/Makefile',probe/'source.Makefile');shutil.copyfile(source/'LICENSE',probe/'LICENSE')
    before=sha(source/'libfac.a')
    with (probe/'make.stdout').open('xb') as out,(probe/'make.stderr').open('xb') as err:
        p=subprocess.run(['make','-f',str(source/'sfac/Makefile'),'sfac'],cwd=probe,stdout=out,stderr=err,timeout=180,check=False)
    if p.returncode or sha(source/'libfac.a')!=before:raise ValueError('Component frontend build failed/changed library')
    for path in probe.iterdir():
        if path.is_file() and path.name!='sfac':shutil.copyfile(path,evidence/path.name)
    return dict(path=str(probe/'sfac'),sha256=sha(probe/'sfac'),libfac_sha256=before,rate_frontend_sha256=fe['sfac_sha256'])


def identical(out,parent):
    """Every stable parent file must exist and match byte-for-byte."""
    names=sorted(p.name for p in parent.iterdir() if p.is_file() and p.name not in UNSTABLE)
    bad=[n for n in names if not (out/n).is_file() or (out/n).read_bytes()!=(parent/n).read_bytes()]
    if bad:raise ValueError('Native output differs from sealed v41: '+', '.join(bad[:10]))
    if any(n not in names for n in STATE_FILES):raise ValueError('Required state file not compared')
    return names


def review(out,root,entry,g,index):
    parent=root/g['parent'];scale=entry['scale']
    text=(out/'input.sf').read_text()
    if sha_text(text)!=g['input_sha256']:raise ValueError('Declared component input differs')
    if ci_review(out,index['potential_layout'],scale,declared_input=text)!=entry['qualified_CI']:raise ValueError('CI replay differs')
    return identical(out,parent)


def component_review(out,scale,g,case_mapping,c):
    raw_text=(out/'raw.tsv').read_text()
    requests=rate_input(scale,g['gauge'],case_mapping)[1]
    audit=rate_review(raw_text,(out/'rates.tr').read_text(),requests,g['gauge'],case_mapping,c)
    if audit!=g['parent_audit']:raise ValueError('Independent rate audit differs from sealed v41')
    raw={}
    for line in raw_text.splitlines():
        if line.startswith('#') or not line.strip():continue
        f=line.split()
        if int(f[5])==-1:raw[int(f[1])]=float(f[12])
    header,contexts=parse((out/'components.tsv').read_text())
    if [x['lower'] for x in contexts]!=g['e1_lower'] or any(x['gauge']!=g['gauge'] for x in contexts):
        raise ValueError('Component context inventory differs')
    cut=float(header['ANGZCUT'])
    return {str(x['lower']):dict(verify(x,raw[x['lower']],cut),grid=grid_context(x)) for x in contexts}


def main():
    ap=argparse.ArgumentParser()
    for key in ('root','source','index','probe','out'):ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--index-sha256',required=True);a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Native component work belongs on Sirius')
    for k,v in vars(a).items():
        if isinstance(v,Path):setattr(a,k,v.resolve())
    if sha(a.index)!=a.index_sha256:raise ValueError('Declared component index differs')
    index=json.loads(a.index.read_text());hold(index)
    if index!=prepare(a.root):raise ValueError('Frozen component inputs differ')
    fe=index['rate_frontend']
    if (index['fac_pin']!=PIN or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or sha(a.source/'libfac.a')!=fe['libfac_sha256'] or sha(a.source/'sfac/sfac')!=fe['upstream_binary_sha256']
        or any(os.environ.get(k)!='1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'))):
        raise ValueError('Component source/runtime differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for path in Path(__file__).parent.glob('rya1354_*'):
        if path.is_file():shutil.copyfile(path,drivers/path.name)
    shutil.copyfile(a.index,a.out/'index.json')
    docs=a.out/'source_documents';docs.mkdir()
    for n in ('faclib/structure.c','faclib/structure.h','faclib/transition.c','faclib/radial.c','faclib/consts.h','faclib/config.h'):
        dest=docs/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.source/n,dest)
    mappings=json.loads((a.root/'runs/fac_radial_rates_v41/index.json').read_text())['cases']
    receipt=dict(schema='rya1354.e1-components.v44',status='RUNNING',fac_pin=PIN,started_utc=datetime.now(timezone.utc).isoformat(),
        index_sha256=a.index_sha256,cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    try:
        receipt['probe']=build(a.root,a.source,a.probe,a.out/'frontend_evidence',index);binary=Path(receipt['probe']['path'])
        entry=index['cases']['ground_coherent10'];g=entry['gauges']['B']
        out=a.out/'canary_unchanged_v41_input';out.mkdir();parent=a.root/g['parent']
        shutil.copyfile(parent/'source.pot',out/'source.pot');shutil.copyfile(parent/'input.sf',out/'input.sf')
        run_native(binary,out,1800)
        receipt['canary']=dict(case='ground_coherent10_B',byte_identical=identical(out,parent))
        (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        for case in CASES:
            entry=index['cases'][case]
            for suffix,g in entry['gauges'].items():
                name=case+'_'+suffix;out=a.out/name;out.mkdir();parent=a.root/g['parent']
                shutil.copyfile(parent/'source.pot',out/'source.pot')
                text=input_text(entry['scale'],g['gauge'],mappings[case]['mapping'])[0];(out/'input.sf').write_text(text)
                run_native(binary,out,1800)
                compared=review(out,a.root,entry,g,index)
                receipt['cases'][name]=dict(compared_byte_identical=compared,
                    components=component_review(out,entry['scale'],g,mappings[case]['mapping'],index['constants']))
                (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_E1_COMPONENTS_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
