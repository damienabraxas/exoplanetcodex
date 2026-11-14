"""Full expanded FAC LS matrix replay; heavy replay runs on Sirius only."""
from collections import Counter
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import struct
from rya1354_fac_recoupling_audit import records
from rya1354_fac_space import configurations


def packed(occ,terms):
    """Pinned PACKLS rules restricted to this n3 s/p/d (d occupancy<=2) space."""
    result='';holes=0;ns=len(occ)
    for i,(n,l,q) in enumerate(occ):
        ell='spd'.index(l);full=4*ell+2;last=i==ns-1
        if ell==2 and q>2:raise ValueError('Unsupported d-shell label')
        prevfull=4*'spd'.index(occ[i-1][1])+2 if i else 0;prevq=occ[i-1][2] if i else 0
        if i and prevq!=prevfull:holes+=prevfull-prevq
        result+=str(n)+l+(f'({q})' if q!=1 else '')
        omit=(i and ((holes<2 and prevq<=1 and q==2)
            or (holes==0 and prevq==prevfull and last and (q<3 or ell<2 or (ell==2 and q>7)))))
        if q!=1 and not omit and q<full-1 and ns!=1:result+=terms[i][:2]
        skip=(i and prevq==prevfull and not last and ((q in (1,full-1) and holes<2) or holes==0))
        if not skip and ((i and q!=full) or last):result+='_'+terms[ns+i-1][:2]
        if not last:result+='.'
    return result


def ls_blocks(text,rank):
    lines=[l for l in text.splitlines() if l.strip()];blocks=[];block=[];i=0
    while i<len(lines):
        if lines[i].strip()=='*':
            if block:blocks.append(block);block=[]
            i+=1;continue
        occ=[(int(n),l,int(q)) for n,l,q in re.findall(r'(\d+)([spdfg])\(\s*(\d+)\)',lines[i])]
        if i+1>=len(lines):raise ValueError('Missing LS coupling line')
        terms=lines[i+1].split();i+=2
        if sum(q for n,l,q in occ)!=8 or len(terms)!=2*len(occ)-1:raise ValueError('LS occupation/coupling count differs')
        q={(n,l):v for n,l,v in occ}
        if q.pop((1,'s'),0)!=2 or q.pop((2,'s'),0)!=2:raise ValueError('LS fixed occupations differ')
        val=tuple(q.pop(key,0) for key in ((2,'p'),(3,'s'),(3,'p'),(3,'d')))
        if q or val not in configurations(rank):raise ValueError('Undeclared LS occupation')
        block.append(dict(name=packed(occ,terms),occupations=occ,terms=terms))
    if block:blocks.append(block)
    return blocks


def label_rows(text,nlev):
    rows={};current=None
    for line in text.splitlines():
        m=re.match(r'^\s*(\d+)\s+(\d+)\s+([+-])\s+(-[\d.]+)\s+([\d.]+)%',line)
        if m:
            ident,j,p,e,norm=m.groups();ident=int(ident)
            if ident in rows:raise ValueError('Duplicate original LS ID')
            current=dict(two_J=2*int(j),parity=p,energy=e,components=[]);rows[ident]=current
        else:
            m=re.match(r'^\s*(-?[\d.]+)\s+([\d.]+)\s+(1s\(.+)$',line)
            if m:
                if current is None:raise ValueError('Unheaded component')
                current['components'].append(m.groups())
    if set(rows)!=set(range(nlev)):raise ValueError('Incomplete LS IDs')
    return rows


def review(source,converted,rank):
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Expanded angular replay belongs on Sirius')
    cm=records(source['o_export_000.cm']);nele,ncsf,norb,nlev,nvalues,nblocks=struct.unpack('<6i',cm[1])
    ls=ls_blocks(converted['state.lsj.c'].decode(),rank);lab=label_rows(converted['state.lsj.lbl'].decode(),nlev)
    t=records(converted['state.lsj.T'])
    if t[0]!=b'jj2lsj' or struct.unpack('<4i',t[1])!=(nele,ncsf,norb,nblocks) or len(ls)!=nblocks:
        raise ValueError('Expanded angular header differs')
    pos=2+nblocks;checked=0;maxorth=0.;result={};jlines=converted['state.lsj.j'].decode().splitlines();jpos=0
    for b in range(nblocks):
        _,nb,nv,jp,p=struct.unpack('<5i',cm[2+4*b]);names=[r['name'] for r in ls[b]]
        if len(names)!=nb or len(set(names))!=nb:raise ValueError('LS labels ambiguous or incomplete')
        if struct.unpack('<4i',t[2+b])!=(b+1,nb,nv,jp):raise ValueError('Angular block differs')
        marker=t[pos];pos+=1
        if marker[:18]!=b' *   Block Number=' or struct.unpack('<i',marker[18:])[0]!=b+1:raise ValueError('Angular marker differs')
        matrix=[[None]*nb for _ in range(nb)]
        for _ in range(nb*nb):
            jj,l,coef=struct.unpack('<iid',t[pos]);pos+=1
            if not 1<=jj<=nb or not 1<=l<=nb or matrix[jj-1][l-1] is not None or not math.isfinite(coef):raise ValueError('Angular coefficient identity differs')
            matrix[jj-1][l-1]=coef
        for i in range(nb):
            for j in range(nb):maxorth=max(maxorth,abs(sum(matrix[k][i]*matrix[k][j] for k in range(nb))-(i==j)))
        ids=struct.unpack('<'+'i'*nv,cm[3+4*b]);energies=struct.unpack('<'+'d'*(nv+1),cm[4+4*b]);values=struct.unpack('<'+'d'*(nb*nv),cm[5+4*b])
        while jpos<len(jlines) and '2*J =' not in jlines[jpos]:jpos+=1
        if jpos==len(jlines):raise ValueError('Missing LS vector block')
        h=re.search(r'2\*J =\s*(\d+)\s+NUMBER =\s*(\d+)',jlines[jpos]);jpos+=1
        if not h or tuple(map(int,h.groups()))!=(jp-1,nv):raise ValueError('LS vector header differs')
        for i,ident in enumerate(ids):
            while jpos<len(jlines) and 'Ssms=' not in jlines[jpos]:jpos+=1
            if jpos+2>=len(jlines):raise ValueError('Missing wrapped LS vector')
            jpos+=1;row=re.match(r'^\s*\d+\s+(-[\d.]+)\s+(1s\(.+)$',jlines[jpos]);jpos+=1
            if not row:raise ValueError('Missing LS vector identity')
            printed=[]
            while jpos<len(jlines) and 'Ssms=' not in jlines[jpos] and '2*J =' not in jlines[jpos]:
                printed.extend(re.findall(r'[-+]?\d+\.\d{8}',jlines[jpos]));jpos+=1
            if len(printed)!=nb:raise ValueError('Wrapped LS coefficient count differs')
            expected=[sum(values[i*nb+k]*matrix[k][l] for k in range(nb)) for l in range(nb)]
            if any(Decimal(format(x,'.8f'))!=Decimal(y) for x,y in zip(expected,printed)):raise ValueError('Expanded LS coefficient/sign replay differs')
            checked+=nb;entry=lab[ident];parity='+' if p==1 else '-';energy=format(energies[0]+energies[i+1],'.9f')
            if (entry['two_J'],entry['parity'],entry['energy'],row[1])!=(jp-1,parity,energy,energy):raise ValueError('LS original ID/J/parity/energy differs')
            dominant=max(range(nb),key=lambda l:float(printed[l])**2)
            if row[2]!=names[dominant] or entry['components'][0][2]!=names[dominant]:raise ValueError('Dominant packed LS name differs')
            for coef,weight,name in entry['components']:
                if name not in names:raise ValueError('Unknown LS component label')
                l=names.index(name)
                if Decimal(coef)!=Decimal(printed[l]) or Decimal(weight)!=Decimal(format(expected[l]**2,'.8f')):raise ValueError('LS label coefficient/weight differs')
            result[str(ident)]=dict(two_J=jp-1,parity=parity,energy_au=energies[0]+energies[i+1],dominant=names[dominant],dominant_printed_weight=float(printed[dominant])**2,components=[dict(coefficient=c,weight=w,name=n) for c,w,n in entry['components']])
    if pos!=len(t) or len(result)!=nlev:raise ValueError('Unconsumed angular records/levels')
    targets={term:[int(i) for i,r in result.items() if r['dominant']=='1s(2).2s(2).2p(3)4S.3s_'+term] for term in ('3S','5S')}
    if any(len(v)!=1 for v in targets.values()):raise ValueError('Expanded target parent/term mapping ambiguous')
    return dict(levels=result,targets=targets,full_LS_coefficients_reconstructed=checked,max_transformation_orthonormality_residual=maxorth)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--converted',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    from rya1354_fac_rates import verify_receipt,read_files
    from rya1354_fac_space_audit import transport
    from rya1354_fac_space import inspect,input_text
    from rya1354_fac_controls import diagnostics
    from rya1354_fac_recouple import CONTROLS,GRASP_PIN
    parent=verify_receipt(a.source);converted=verify_receipt(a.converted)
    if (parent['status']!='FAC_EXPANDED_LEVELS_COMPLETED_REVIEW_REQUIRED' or converted['status']!='FAC_ANGULAR_TRANSFORM_COMPLETED_REVIEW_REQUIRED' or converted['grasp_pin']!=GRASP_PIN or converted['source_receipt_sha256']!=hashlib.sha256((a.source/'receipt.json').read_bytes()).hexdigest()):raise ValueError('Expanded source qualification differs')
    report=dict(schema='rya1354.fac-expanded-ls-audit.v24',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,independent_angular_implementation=False,cases={})
    for rank in (1,2):
        for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
            name=f'n3_r{rank}_{label}';s=read_files(a.source/name);c=read_files(a.converted/name)
            if (s['input.sf'].decode()!=input_text(rank,mult,scale) or not diagnostics(s['sfac.stdout'].decode(),mult,256)['printed_rule_verified'] or c['state.c']!=s['o_export_000.c'] or c['state.cm']!=s['o_export_000.cm'] or c['jj2lsj.stdin'].decode()!=CONTROLS):raise ValueError('Expanded replay input differs')
            inspect(s['o.lev'].decode(),rank);transport(s['o.basis'].decode(),s['o_export_000.c'].decode(),s['o_export_000.cm'],rank)
            report['cases'][name]=review(s,c,rank)
    a.out.write_text(json.dumps(report,indent=2)+'\n')
