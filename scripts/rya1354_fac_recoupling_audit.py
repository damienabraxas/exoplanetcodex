"""Replay exact FAC->LS angular transport, with limited analytic term controls."""
import argparse
from collections import Counter
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re
import struct

from rya1354_native_parts import load_parts
from rya1354_fac_mapping import verify_transport
from rya1354_fac_control_audit import mixing
from rya1354_fac_recouple import CONTROLS,GRASP_PIN
from rya1354_fac_controls import diagnostics
from rya1354_fac_coherent import input_text as numerical_input

ROOT=Path(__file__).resolve().parents[1]
RUN='runs/fac_recoupling_v19/'
SOURCE='runs/fac_coherent_controls_v18/'


def records(data):
    result=[];offset=0
    while offset<len(data):
        if offset+8>len(data):raise ValueError('Truncated angular record')
        n=struct.unpack_from('<i',data,offset)[0]
        if n<0 or offset+n+8>len(data) or struct.unpack_from('<i',data,offset+n+4)[0]!=n:
            raise ValueError('Angular record marker differs')
        result.append(data[offset+4:offset+n+4]);offset+=n+8
    return result


def ls_basis(text):
    blocks=[];current=[];lines=[l for l in text.splitlines() if l.strip()]
    i=0
    while i<len(lines):
        if lines[i].strip()=='*':
            if current:blocks.append(current);current=[]
            i+=1;continue
        if i+1>=len(lines):raise ValueError('Incomplete LS basis')
        occupations=[(n,l,int(q)) for n,l,q in re.findall(r'(\d+)([spdfg])\(\s*(\d+)\)',lines[i])]
        terms=lines[i+1].split()
        if sum(q for n,l,q in occupations)!=8 or len(terms)!=2*len(occupations)-1:
            raise ValueError('LS basis electron/coupling inventory differs')
        if occupations==[('1','s',2),('2','s',2),('2','p',4)]:
            name='1s(2).2s(2).2p(4)_'+terms[-1]
        elif occupations==[('1','s',2),('2','s',2),('2','p',3),('3','s',1)]:
            name='1s(2).2s(2).2p(3)'+terms[-2]+'.3s_'+terms[-1]
        else:raise ValueError('Unexpected LS configuration')
        current.append(name);i+=2
    if current:blocks.append(current)
    return blocks


def labels(text):
    rows={};current=None
    for line in text.splitlines():
        m=re.match(r'^\s*(\d+)\s+(\d+)\s+([+-])\s+(-[\d.]+)\s+([\d.]+)%',line)
        if m:
            ident,j,p,e,norm=m.groups();ident=int(ident)
            if ident in rows:raise ValueError('Duplicate converted FAC level ID')
            current=dict(index=ident,two_J=2*int(j),parity=p,energy=e,components=[])
            rows[ident]=current
        else:
            m=re.match(r'^\s*(-?[\d.]+)\s+([\d.]+)\s+(1s\(.+)$',line)
            if m:
                if current is None:raise ValueError('Unheaded LS coefficient')
                current['components'].append(dict(coefficient=m[1],weight=m[2],name=m[3]))
    if set(rows)!=set(range(15)):raise ValueError('Converted FAC level IDs incomplete')
    return rows


def review_case(source,converted):
    """Replay matrices and printed coefficients; no threshold fit to energies."""
    meta=mixing(source['o_export_000.cm']);cm=records(source['o_export_000.cm'])
    ls=ls_basis(converted['state.lsj.c'].decode());lab=labels(converted['state.lsj.lbl'].decode())
    t=records(converted['state.lsj.T'])
    if t[0]!=b'jj2lsj' or struct.unpack('<4i',t[1])!=(8,15,5,7) or len(ls)!=7:
        raise ValueError('Transformation source inventory differs')
    pos=9;checked=0;max_orth=0.0;result={};jlines=converted['state.lsj.j'].decode().splitlines();jpos=0
    for block_index,block in enumerate(meta['blocks']):
        nb=block['basis'];nv=block['levels'];names=ls[block_index]
        if len(names)!=nb or struct.unpack('<4i',t[2+block_index])!=(block_index+1,nb,nv,block['two_J']+1):
            raise ValueError('Transformation block dimensions differ')
        marker=t[pos];pos+=1
        if marker[:18]!=b' *   Block Number=' or struct.unpack('<i',marker[18:])[0]!=block_index+1:
            raise ValueError('Transformation block marker differs')
        matrix=[[None]*nb for _ in range(nb)]
        for _ in range(nb*nb):
            jj,l,coef=struct.unpack('<iid',t[pos]);pos+=1
            if not 1<=jj<=nb or not 1<=l<=nb or matrix[jj-1][l-1] is not None or not math.isfinite(coef):
                raise ValueError('Transformation coefficient identity differs')
            matrix[jj-1][l-1]=coef
        for i in range(nb):
            for j in range(nb):
                dot=sum(matrix[k][i]*matrix[k][j] for k in range(nb))
                max_orth=max(max_orth,abs(dot-(1 if i==j else 0)))
        values=struct.unpack('<'+'d'*(nb*nv),cm[2+4*block_index+3])
        energies=struct.unpack('<'+'d'*(nv+1),cm[2+4*block_index+2])
        while jpos<len(jlines) and '2*J =' not in jlines[jpos]:jpos+=1
        if jpos==len(jlines):raise ValueError('Full LS coefficient block absent')
        header=re.search(r'2\*J =\s*(\d+)\s+NUMBER =\s*(\d+)',jlines[jpos]);jpos+=1
        if not header or tuple(map(int,header.groups()))!=(block['two_J'],nv):
            raise ValueError('Full LS coefficient header differs')
        for i,ident in enumerate(block['indices']):
            while jpos<len(jlines) and 'Ssms=' not in jlines[jpos]:jpos+=1
            if jpos+2>=len(jlines):raise ValueError('Full LS vector missing')
            jpos+=1;row=re.match(r'^\s*\d+\s+(-[\d.]+)\s+(1s\(.+)$',jlines[jpos]);jpos+=1
            if row is None:raise ValueError('Full LS vector identity missing')
            printed=re.findall(r'[-+]?\d+\.\d{8}',jlines[jpos]);jpos+=1
            if len(printed)!=nb:raise ValueError('Full LS coefficient inventory differs')
            expected=[sum(values[i*nb+k]*matrix[k][l] for k in range(nb)) for l in range(nb)]
            if any(Decimal(format(x,'.8f'))!=Decimal(y) for x,y in zip(expected,printed)):
                raise ValueError('LS coefficient/sign/order reconstruction differs')
            checked+=nb;entry=lab[ident]
            p='+' if block['parity']==1 else '-'
            energy=format(energies[0]+energies[i+1],'.9f')
            if entry['two_J']!=block['two_J'] or entry['parity']!=p or entry['energy']!=energy or row[1]!=energy:
                raise ValueError('FAC ID/parity/J/native energy transport differs')
            weights={name:float(value)**2 for name,value in zip(names,printed)}
            dominant=max(weights,key=weights.get)
            if row[2]!=dominant or entry['components'][0]['name']!=dominant:
                raise ValueError('Native dominant LS identity differs')
            for component in entry['components']:
                if component['name'] not in names:raise ValueError('Unknown printed LS component')
                l=names.index(component['name'])
                if Decimal(component['coefficient'])!=Decimal(printed[l]):
                    raise ValueError('Printed LS-label coefficient differs')
                if Decimal(component['weight'])!=Decimal(format(expected[l]**2,'.8f')):
                    raise ValueError('Printed LS-label weight differs')
            result[ident]=dict(two_J=block['two_J'],parity=p,dominant=dominant,
                dominant_printed_weight=weights[dominant],components=entry['components'])
    if pos!=len(t):raise ValueError('Unconsumed angular transformation records')
    expected=Counter([('+',0,'3P'),('+',0,'1S'),('+',2,'3P'),('+',4,'3P'),('+',4,'1D'),
        ('-',0,'3P'),('-',2,'3S'),('-',2,'3D'),('-',2,'3P'),('-',2,'1P'),
        ('-',4,'5S'),('-',4,'3D'),('-',4,'1D'),('-',4,'3P'),('-',6,'3D')])
    if Counter((r['parity'],r['two_J'],r['dominant'].rsplit('_',1)[1]) for r in result.values())!=expected:
        raise ValueError('Minimal LS term inventory differs')
    for p,j,term in [('+',2,'3P'),('-',0,'3P'),('-',6,'3D')]:
        one,=[r for r in result.values() if (r['parity'],r['two_J'])==(p,j)]
        if one['dominant'].rsplit('_',1)[1]!=term or one['dominant_printed_weight']!=1.0:
            raise ValueError('Analytic unique-term control failed')
    targets={term:[i for i,r in result.items() if r['dominant']=='1s(2).2s(2).2p(3)4S.3s_'+term] for term in ('3S','5S')}
    if any(len(v)!=1 for v in targets.values()):raise ValueError('Target parent/term mapping ambiguous')
    return dict(levels={str(k):v for k,v in result.items()},targets=targets,full_LS_coefficients_reconstructed=checked,
                max_transformation_orthonormality_residual=max_orth,unique_term_controls=3)


def summarize(base):
    source,_=load_parts(base/'fac_coherent_complete_v18');files,manifest=load_parts(base/'fac_recoupling_complete_v19')
    receipt=json.loads(files[RUN+'receipt.json'])
    if receipt['status']!='FAC_ANGULAR_TRANSFORM_COMPLETED_REVIEW_REQUIRED' or receipt['grasp_pin']!=GRASP_PIN:
        raise ValueError('Recoupling terminal scope differs')
    if receipt['source_receipt_sha256']!=hashlib.sha256(source[SOURCE+'receipt.json']).hexdigest():
        raise ValueError('Recoupling numerical source differs')
    cases={}
    for case,multiplier,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
        original={k.removeprefix(SOURCE+case+'/'):v for k,v in source.items() if k.startswith(SOURCE+case+'/')}
        converted={k.removeprefix(RUN+case+'/'):v for k,v in files.items() if k.startswith(RUN+case+'/')}
        if (original['input.sf'].decode()!=numerical_input(multiplier,scale)
                or not diagnostics(original['sfac.stdout'].decode(),multiplier,256)['printed_rule_verified']):
            raise ValueError('FAC numerical source requalification failed')
        if (converted['state.c']!=original['o_export_000.c'] or converted['state.cm']!=original['o_export_000.cm']
                or converted['jj2lsj.stdin'].decode()!=CONTROLS):
            raise ValueError('Recoupling exact input/controls differ')
        transport=verify_transport(original['o.basis'].decode(),original['o_export_000.c'].decode(),original['o_export_000.cm'])
        if transport!=receipt['transport'][case]:raise ValueError('Native/export transport receipt differs')
        cases[case]=dict(transport=transport,recoupling=review_case(original,converted))
    return dict(schema='rya1354.fac-recoupling-audit.v19',cases=cases,parts=manifest['parts'],
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        rates_admitted=False,independent_angular_implementation=False,
        scope='Exact native/export/LS representation consistency with limited analytic term controls for this minimal configuration only. FAC Hamiltonian independent; angular converter shared with GRASP. No physical accuracy/calibrated benchmark.',
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')
