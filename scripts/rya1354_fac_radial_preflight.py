"""Bound4s/4p radial preflight in existing saved potentials; no CI/rates."""
import argparse
from collections import defaultdict
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import platform
import re
import shutil
import struct
import subprocess
from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_rates import verify_receipt

OLD=((1,-1),(2,-1),(2,1),(2,-2),(3,-1),(3,1),(3,-2),(3,2),(3,-3))
NEW=((4,-1),(4,1),(4,-2))
KAPPAS=(-1,1,-2,2,-3)


def name(n,k):return f'n{n}_k{k}'


def input_text(scale):
    text="SetAtom('O')\nSetBreit(0)\nSetSE(0)\nSetVP(0)\nSetMS(0,0)\nSetConfigEnergyMode(-1)\nSetOption('radial:config_energy',-1)\n"
    for key,value in [('enerelerr',1e-4),('enerelerr1',1e-8),('eneabserr',1e-5)]:text+=f"SetOption('orbital:{key}',{value/scale:.17g})\n"
    text+="RestorePotential('source.pot')\nSavePotential('restored.pot')\n"
    for n,k in OLD:text+=f"WaveFuncTable('before_{name(n,k)}.wav',{n},{k})\n"
    for k in KAPPAS:text+=f"RadialOverlaps('before_k{k}.overlap',{k})\n"
    text+="SavePotential('old_solved.pot')\n"
    for n,k in NEW:text+=f"WaveFuncTable('new_{name(n,k)}.wav',{n},{k})\n"
    for n,k in OLD:text+=f"WaveFuncTable('after_{name(n,k)}.wav',{n},{k})\n"
    for k in KAPPAS:text+=f"RadialOverlaps('after_k{k}.overlap',{k})\n"
    return text+"SavePotential('extended.pot')\n"


def wave(text,n,k):
    headers={}
    for line in text.splitlines():
        m=re.match(r'^#\s*(\w+)\s*=\s*([-+\d.Ee]+)\s*$',line)
        if m:headers[m[1]]=float(m[2])
    if headers.get('n')!=n or headers.get('kappa')!=k or not math.isfinite(headers.get('energy',math.nan)) or headers['energy']>=0:raise ValueError('Requested bound orbital/header/energy differs')
    rows=[]
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):continue
        f=line.split()
        if len(f)!=12:raise ValueError('Pinned bound wavefunction column count differs')
        vals=list(map(float,f[1:]))
        if not f[0].isdigit() or int(f[0])!=len(rows) or not all(math.isfinite(x) for x in vals):raise ValueError('Invalid bound radial sample/index')
        if vals[0]<=0 or (rows and vals[0]<=rows[-1][1]):raise ValueError('Radial grid is not positive/increasing')
        rows.append((int(f[0]),*vals))
    if not rows or len(rows)!=int(headers.get('ilast',-2))+1 or not any(r[4]!=0 for r in rows):raise ValueError('Incomplete/zero bound radial function')
    # Columns are pinned source: index,r,rVT0,rVeff(cached),P,Q,... .
    # Effective-potential scratch column is not a per-orbital observable on cached re-export.
    stable=[(r[0],r[1],r[2],r[4],r[5]) for r in rows]
    return dict(n=n,kappa=k,energy_eV=headers['energy'],samples=len(rows),qrnorm=headers.get('qrnorm'),ilast=int(headers['ilast'])),stable


def overlaps(text,k,orbitals):
    required={(n1,n2) for n1,k1 in orbitals for n2,k2 in orbitals if k1==k2==k};seen={};maxoff=0.;maxh=0.
    for line in text.splitlines():
        if not line.strip():continue
        f=line.split()
        if len(f)!=11:raise ValueError('Native overlap layout differs')
        n1,k1=int(f[0]),int(f[1]);n2,k2=int(f[4]),int(f[5]);key=(n1,n2)
        if key not in required or key in seen or k1!=k or k2!=k:raise ValueError('Unexpected/duplicate overlap identity')
        values=[float(f[i]) for i in (2,3,6,7,8,9,10)]
        if not all(math.isfinite(x) for x in values) or values[0]>=0 or values[2]>=0:raise ValueError('Unbound/nonfinite overlap record')
        if n1==n2 and f[8]!='1.00000E+00':raise ValueError('Native normalization not unity at native print precision')
        if n1!=n2:maxoff=max(maxoff,abs(values[4]))
        maxh=max(maxh,abs(values[5]-values[6]));seen[key]=values
    if set(seen)!=required:raise ValueError('Missing native overlap pair')
    return dict(pairs=len(seen),diagonal_printed_unity=True,max_abs_offdiagonal_overlap=maxoff,max_abs_printed_H_asymmetry_au=maxh,
                scope='Native printed normalization/overlap diagnostic; no physical accuracy or grid-convergence tolerance inferred.')


def review(files,scale,spec):
    if files['input.sf'].decode()!=input_text(scale):raise ValueError('Radial probe input differs')
    from rya1354_fac_potential_fields import compare
    comparisons={fn:compare(files['source.pot'],files[fn],spec) for fn in ('restored.pot','old_solved.pot','extended.pot')}
    if not comparisons['restored.pot']['serialized_byte_identical']:raise ValueError('Potential restore roundtrip changed bytes')
    mode,flag,core,nmax,ngrid=struct.unpack('<5i',files['source.pot'][:20])
    if ngrid!=1500 or nmax<4:raise ValueError('Paired radial grid/bound-state envelope differs')
    result=dict(potential_comparisons=comparisons,grid_points=ngrid,grid_supported_nmax=nmax,old_orbitals={},new_orbitals={},overlaps={})
    for n,k in OLD:
        a,ar=wave(files['before_'+name(n,k)+'.wav'].decode(),n,k);b,br=wave(files['after_'+name(n,k)+'.wav'].decode(),n,k)
        if a!=b or ar!=br:raise ValueError('Old orbital physical radial samples changed')
        result['old_orbitals'][name(n,k)]=a
    for n,k in NEW:result['new_orbitals'][name(n,k)]=wave(files['new_'+name(n,k)+'.wav'].decode(),n,k)[0]
    for k in KAPPAS:
        result['overlaps'][str(k)]=dict(before=overlaps(files['before_k'+str(k)+'.overlap'].decode(),k,OLD),after=overlaps(files['after_k'+str(k)+'.overlap'].decode(),k,OLD+NEW))
    return result


def main():
    ap=argparse.ArgumentParser()
    for n in ('source','baseline','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Radial bound-state calculations belong on Sirius')
    for n,v in vars(a).items():setattr(a,n,v.resolve())
    parent=verify_receipt(a.baseline);binary=a.source/'sfac/sfac'
    if (parent['status']!='FAC_FULL_VALENCE_LEVELS_COMPLETED_REVIEW_REQUIRED' or parent['fac_pin']!=PIN
        or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or sha(binary)!=parent['probe']['upstream_binary_sha256']):raise ValueError('Radial baseline/source/binary differs')
    a.out.mkdir(parents=True,exist_ok=False)
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,a.out/p.name)
    for n in ('doc/manual.tex','faclib/radial.c','faclib/orbital.c','sfac/sfac.c','faclib/consts.h'):
        p=a.out/'source_documents'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.source/n,p)
    from rya1354_fac_potential_fields import layout
    spec=layout((a.source/'faclib/radial.c').read_text(),(a.source/'faclib/orbital.c').read_text(),(a.source/'faclib/consts.h').read_text())
    receipt=dict(potential_layout=spec,schema='rya1354.fac-radial-preflight.v32',status='RUNNING',fac_pin=PIN,binary_sha256=sha(binary),baseline_receipt_sha256=sha(a.baseline/'receipt.json'),started_utc=datetime.now(timezone.utc).isoformat(),
        cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Restored qualified potential;bound4s/4p preflight only;no freshSCF,CI,LSor rates. Native overlaps and unchanged old physical P/Q/r/rVT0 samples;all potential/grid/header bytes unchanged exceptsource-derivedW/dW/dW2 caches;not radial grid convergence or accuracy.')
    try:
        for policy in ('ground','equal'):
            for label,scale in [('coherent10',10),('coherent100',100)]:
                name=policy+'_'+label;out=a.out/name;out.mkdir();shutil.copyfile(a.baseline/name/'o.pot',out/'source.pot');(out/'input.sf').write_text(input_text(scale))
                with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
                    p=subprocess.run([str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=300,check=False)
                if p.returncode:raise ValueError(f'Radial native exit{p.returncode}: '+name)
                files={p.name:p.read_bytes() for p in out.iterdir() if p.is_file()};receipt['cases'][name]=review(files,scale,spec)
                (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_RADIAL_PREFLIGHT_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
