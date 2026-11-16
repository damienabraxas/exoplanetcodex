"""Cold radial E1 frequency/gauge diagnostics of sealed v41 orbital models."""
import argparse
from datetime import datetime,timezone
import difflib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess

from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_rates import verify_receipt
from rya1354_fac_radial_preflight import OLD,NEW,name,wave,input_text as radial_input,review as radial_review
from rya1354_fac_potential_fields import compare
from rya1354_fac_radial_recouple import CASES
from rya1354_fac_radial_angular_remaining_metadata import hold
from rya1354_fac_radial_rates import run_native
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants
from rya1354_fac_radial_rates import rate_input

ACTIVE=tuple(x for x in OLD+NEW if x[0]>=2 and x not in ((2,-1),))


def pairs():
    def ell(k):return k if k>0 else -k-1
    return [(a,b) for a in ACTIVE for b in ACTIVE
            if (ell(a[1])+ell(b[1]))%2==1
            and abs((2*abs(a[1])-1)-(2*abs(b[1])-1))<=2<=2*abs(a[1])+2*abs(b[1])-2]


def input_text(scale,gauge,energy,photon_floor=None):
    if scale not in (10,100) or gauge not in (1,2) or not math.isfinite(energy) or energy<0:
        raise ValueError('Undeclared radial operator controls')
    text=radial_input(scale)
    if photon_floor is not None:
        if not math.isfinite(photon_floor) or photon_floor<=0:raise ValueError('Invalid photon-grid floor')
        text+=f"SetOption('radial:awmin',{photon_floor:.17g})\n"
    text+=f'SetTransitionOptions({gauge},0)\n'
    for (n,k),(n2,k2) in pairs():
        text+=f"RawRadialE1('radial.tsv',{n},{k},{n2},{k2},{energy:.17g},{gauge})\n"
    for n,k in OLD+NEW:text+=f"WaveFuncTable('probe_after_{name(n,k)}.wav',{n},{k})\n"
    return text+"SavePotential('probe.pot')\n"


def raw_rows(text,gauge,request,c):
    rows={}
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):continue
        f=line.split()
        if len(f)!=16:raise ValueError('Radial operator column count differs')
        n,k,n2,k2,g,mode,gap,ng=map(int,f[:8]);key=((n,k),(n2,k2));v=list(map(float,f[8:]))
        if key not in pairs() or key in rows or (g,mode,gap,ng)!=(gauge,0,int(request==0),1):
            raise ValueError('Radial operator pair/gauge/context differs')
        if not all(math.isfinite(x) for x in v):raise ValueError('Nonfinite radial operator')
        req,e,ep1,ep2,aw,native_aw,r,nr=v
        if req!=request or e<=0 or max(ep1,ep2)>=0 or (request>0 and e!=request) or (request==0 and e!=abs(ep1-ep2)):
            raise ValueError('Radial operator frequency/orbital energy differs')
        if aw!=c['FINE_STRUCTURE_CONST']*e or native_aw!=aw:raise ValueError('Native photon grid differs')
        rows[key]=dict(requested_au=req,used_au=e,eps1_au=ep1,eps2_au=ep2,radial_FR=r,radial_NR_length=nr)
    if set(rows)!=set(pairs()):raise ValueError('Incomplete radial pair coverage')
    return rows


def prepare(root,*,photon_floor=None):
    parent=root/'runs/fac_radial_rates_v41';receipt=verify_receipt(parent);hold(receipt)
    index=json.loads((parent/'index.json').read_text());hold(index)
    if receipt['status']!='FAC_RADIAL_RATES_COMPLETED_REVIEW_REQUIRED' or sha(parent/'index.json')!=receipt['index_sha256']:
        raise ValueError('Original rate checkpoint differs')
    c=constants((parent/'source_documents/consts.h').read_text())
    result=dict(schema='rya1354.radial-operator-input.v42',parent_receipt_sha256=sha(parent/'receipt.json'),
        parent_index_sha256=sha(parent/'index.json'),fac_pin=index['fac_pin'],probe=index['probe'],
        potential_layout=index['potential_layout'],constants=c,cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    for case in CASES:
        entry=index['cases'][case];mapping=entry['mapping'];scale=100 if case.endswith('100') else 10
        upper,=mapping['targets']['3S'];low=min(mapping['levels'],key=lambda i:mapping['levels'][i]['energy_au'])
        level=mapping['levels'][low]
        if (level['dominant'],level['two_J'],level['parity'])!=('1s(2).2s(2).2p(4)_3P',4,'+'):
            raise ValueError('Reference ground-state identity differs')
        raw={}
        for suffix,gauge in [('B',2),('C',1)]:
            path=parent/(case+'_'+suffix);text,requests,lower=rate_input(scale,gauge,mapping)
            if (path/'input.sf').read_text()!=text:raise ValueError('Parent input differs')
            audit=rate_review((path/'raw.tsv').read_text(),(path/'rates.tr').read_text(),requests,gauge,mapping,c)
            if audit!=receipt['cases'][case+'_'+suffix]['audit']:raise ValueError('Parent independent raw audit differs')
            raw[suffix]=(path/'raw.tsv').read_text()
        reference=[line.split() for line in raw['B'].splitlines() if line.strip() and not line.startswith('#')]
        row,=[f for f in reference if tuple(map(int,(f[0],f[1],f[5])))==(upper,int(low),-1)]
        result['cases'][case]=dict(reference=str((parent/(case+'_B')).relative_to(root)),scale=scale,upper=upper,lower=int(low),
            reference_level=level,upper_level=mapping['levels'][str(upper)],energy_au=float(row[11]),
            raw=raw,levels={i:mapping['levels'][i] for i in map(str,lower)},
            radial_preflight=entry['qualified_CI']['radial_preflight'])
    for case,value in result['cases'].items():
        partner=('equal' if case.startswith('ground') else 'ground')+'_'+case.split('_',1)[1]
        value['partner']=partner;value['paired_energy_au']=result['cases'][partner]['energy_au']
    if photon_floor is not None:
        if not math.isfinite(photon_floor) or photon_floor<=0:raise ValueError('Invalid photon-grid floor')
        previous=root/'runs/fac_radial_operator_v42';prior=verify_receipt(previous);hold(prior)
        if prior['status']!='FAILED_RETAINED' or set(prior['cases'])!={'ground_coherent10_B_own','ground_coherent10_B_paired'}:
            raise ValueError('Retained first diagnostic failure differs')
        result.update(schema='rya1354.radial-operator-input.v43',photon_grid_min=photon_floor,
            retained_failure=str(previous.relative_to(root)),retained_failure_receipt_sha256=sha(previous/'receipt.json'))
    return result


def build(source,probe,evidence):
    probe.mkdir(parents=True,exist_ok=False);evidence.mkdir(parents=True,exist_ok=False)
    original=(source/'sfac/sfac.c').read_text();snippet=Path(__file__).with_suffix('.c').read_text()
    marker='#include "global.h"\n';registration='  {"TransitionTable", PTransitionTable, METH_VARARGS},'
    if original.count(marker)!=1 or original.count(registration)!=1:raise ValueError('Probe insertion points differ')
    patched=original.replace(marker,marker+'\n'+snippet+'\n',1).replace(registration,
        '  {"RawRadialE1", PRawRadialE1, METH_VARARGS},\n'+registration,1)
    (probe/'sfac.c').write_text(patched);(probe/'upstream.sfac.c').write_text(original)
    (probe/'observability.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),fromfile='pinned/sfac.c',tofile='read-only-radial/sfac.c')))
    for n in ('stoken.c','stoken.h'):shutil.copyfile(source/'sfac'/n,probe/n)
    shutil.copyfile(source/'sfac/Makefile',probe/'source.Makefile');shutil.copyfile(source/'LICENSE',probe/'LICENSE')
    before=sha(source/'libfac.a')
    with (probe/'make.stdout').open('xb') as out,(probe/'make.stderr').open('xb') as err:
        p=subprocess.run(['make','-f',str(source/'sfac/Makefile'),'sfac'],cwd=probe,stdout=out,stderr=err,timeout=180,check=False)
    if p.returncode or sha(source/'libfac.a')!=before:raise ValueError('Read-only radial build failed/changed library')
    for path in probe.iterdir():
        if path.is_file() and path.name!='sfac':shutil.copyfile(path,evidence/path.name)
    return dict(path=str(probe/'sfac'),sha256=sha(probe/'sfac'),libfac_sha256=before,upstream_binary_sha256=sha(source/'sfac/sfac'))


def review_files(files,entry,index,gauge,energy):
    text=input_text(entry['scale'],gauge,energy,index.get('photon_grid_min'))
    if files['input.sf'].decode()!=text:raise ValueError('Declared radial operator input differs')
    files=dict(files)
    files['input.sf']=radial_input(entry['scale']).encode()
    preflight=radial_review(files,entry['scale'],index['potential_layout'])
    if preflight!=entry['radial_preflight']:raise ValueError('Original radial preflight differs')
    for n,k in OLD+NEW:
        prefix='before_' if (n,k) in OLD else 'new_'
        if wave(files[prefix+name(n,k)+'.wav'].decode(),n,k)!=wave(files['probe_after_'+name(n,k)+'.wav'].decode(),n,k):
            raise ValueError('Probe changed physical orbital')
    potential=compare(files['source.pot'],files['probe.pot'],index['potential_layout'])
    rows=raw_rows(files['radial.tsv'].decode(),gauge,energy,index['constants'])
    return dict(radial_preflight=preflight,post_probe_orbitals_unchanged=True,potential=potential,pairs=len(rows),
        rows=[dict(n1=a[0],kappa1=a[1],n2=b[0],kappa2=b[1],**v) for (a,b),v in rows.items()])


def review(out,entry,index,gauge,energy):
    return review_files({p.name:p.read_bytes() for p in out.iterdir() if p.is_file()},entry,index,gauge,energy)


def main():
    ap=argparse.ArgumentParser()
    for key in ('root','source','index','probe','out'):ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--index-sha256',required=True);a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Radial diagnostic work belongs on Sirius')
    for k,v in vars(a).items():
        if isinstance(v,Path):setattr(a,k,v.resolve())
    if sha(a.index)!=a.index_sha256:raise ValueError('Declared radial diagnostic index differs')
    index=json.loads(a.index.read_text());hold(index)
    if index!=prepare(a.root,photon_floor=index.get('photon_grid_min')):raise ValueError('Original frozen rate inputs differ')
    old=index['probe']
    if (index['fac_pin']!=PIN or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or sha(a.source/'libfac.a')!=old['libfac_sha256'] or sha(a.source/'sfac/sfac')!=old['upstream_binary_sha256']
        or any(os.environ.get(k)!='1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'))):
        raise ValueError('Diagnostic source/runtime differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for path in Path(__file__).parent.glob('rya1354_*'):
        if path.is_file():shutil.copyfile(path,drivers/path.name)
    shutil.copyfile(a.index,a.out/'index.json')
    docs=a.out/'source_documents';docs.mkdir()
    for n in ('faclib/radial.c','faclib/transition.c','faclib/orbital.c','faclib/consts.h','faclib/radial.h','sfac/sfac.c','doc/manual.tex'):
        dest=docs/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.source/n,dest)
    version='v43' if 'photon_grid_min' in index else 'v42'
    receipt=dict(schema='rya1354.radial-operator.'+version,status='RUNNING',fac_pin=PIN,started_utc=datetime.now(timezone.utc).isoformat(),
        index_sha256=a.index_sha256,cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    try:
        receipt['probe']=build(a.source,a.probe,a.out/'frontend_evidence');binary=Path(receipt['probe']['path'])
        canaries=[]
        entry=index['cases']['ground_coherent10']
        for label,exe in [('original',a.source/'sfac/sfac'),('diagnostic',binary)]:
            out=a.out/('canary_'+label);out.mkdir();shutil.copyfile(a.root/entry['reference']/'source.pot',out/'source.pot')
            (out/'input.sf').write_text(radial_input(10));run_native(exe,out,300);canaries.append(out)
        relevant=[p.name for p in canaries[0].iterdir() if p.suffix in ('.wav','.pot','.overlap')]
        if any((canaries[0]/n).read_bytes()!=(canaries[1]/n).read_bytes() for n in relevant):
            raise ValueError('Read-only frontend canary changed radial output')
        receipt['canary']=dict(all_physical_files_byte_identical=True,files=sorted(relevant))
        if 'retained_failure' in index:
            previous=a.root/index['retained_failure'];prior=verify_receipt(previous)
            receipt['retained_failure_receipt_sha256']=sha(previous/'receipt.json');receipt['previous_successful_controls_exact']={}
        for case,entry in index['cases'].items():
            for suffix,gauge in [('B',2),('C',1)]:
                for kind,energy in [('own',entry['energy_au']),('paired',entry['paired_energy_au']),('orbital_gap',0.0)]:
                    label=case+'_'+suffix+'_'+kind;out=a.out/label;out.mkdir()
                    shutil.copyfile(a.root/entry['reference']/'source.pot',out/'source.pot')
                    (out/'input.sf').write_text(input_text(entry['scale'],gauge,energy,index.get('photon_grid_min')));run_native(binary,out,300)
                    receipt['cases'][label]=review(out,entry,index,gauge,energy)
                    if 'retained_failure' in index and label in prior['cases']:
                        if receipt['cases'][label]!=prior['cases'][label] or (out/'radial.tsv').read_bytes()!=(previous/label/'radial.tsv').read_bytes():
                            raise ValueError('New photon floor changed a prior successful control')
                        receipt['previous_successful_controls_exact'][label]=True
                    (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_RADIAL_OPERATOR_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
