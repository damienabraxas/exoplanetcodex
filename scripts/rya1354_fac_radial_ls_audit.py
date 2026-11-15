"""BLAS-assisted full signed LS replay with exact legacy-sum rounding fallback."""
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
from rya1354_fac_stream_transport_v35 import configurations,N3,EXTENDED


def nk_d_term(term):
    # getchLS writes rawseniority; getchLS_lbl NK convention maps it to indexw.
    parent=term[:2];v=int(term[2:])
    maps={'1S':{0:1,4:2},'1D':{2:1,4:2},'1G':{2:1,4:2},
          '2D':{1:1,3:2,5:3},'2F':{3:1,5:2},'2G':{3:1,5:2},
          '3P':{2:1,4:2},'3F':{2:1,4:2}}
    if parent in maps:
        if v not in maps[parent]:raise ValueError('Unsupported NK d-shell seniority')
        w=maps[parent][v]
    else:w=0
    return parent+str(w)


def packed(occ,terms):
    """Pinned PACKLS rules restricted to this n3 s/p/d (d occupancy<=4) space."""
    terms=list(terms)
    for i,(n,l,q) in enumerate(occ):
        if l=='d':terms[i]=nk_d_term(terms[i])
    result='';holes=0;ns=len(occ)
    for i,(n,l,q) in enumerate(occ):
        ell='spd'.index(l);full=4*ell+2;last=i==ns-1
        if ell==2 and q>4:raise ValueError('Unsupported d-shell label')
        prevfull=4*'spd'.index(occ[i-1][1])+2 if i else 0;prevq=occ[i-1][2] if i else 0
        if i and prevq!=prevfull:holes+=prevfull-prevq
        result+=str(n)+l+(f'({q})' if q!=1 else '')
        omit=(i and ((holes<2 and prevq<=1 and q==2)
            or (holes==0 and prevq==prevfull and last and (q<3 or ell<2 or (ell==2 and q>7)))))
        if q!=1 and not omit and q<full-1 and ns!=1:
            use3=ell==2 and ((q==3 and terms[i][:2]=='2D') or (q==4 and terms[i][:2] in ('3F','3P','1G','1D','1S')))
            result+=terms[i][:3] if use3 else terms[i][:2]
        skip=(i and prevq==prevfull and not last and ((q in (1,full-1) and holes<2) or holes==0))
        if not skip and ((i and q!=full) or last):result+='_'+terms[ns+i-1][:2]
        if not last:result+='.'
    return result


def ls_blocks(text,shells=N3):
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
        val=tuple(q.pop(key,0) for key in tuple((int(name[:-1]),name[-1]) for name,_,_ in shells))
        if q or val not in configurations(shells):raise ValueError('Undeclared LS occupation')
        block.append(dict(name=packed(occ,terms),occupations=occ,terms=terms))
    if block:blocks.append(block)
    return blocks


def stdout_ids(text):
    result=[];waiting=False
    for line in text.splitlines():
        if 'Weights of major contributors to ASF in LS-coupling:' in line:
            waiting=True;continue
        if waiting:
            m=re.match(r'^\s*(\d+)\s+(\d+)\s+([+-])\s+.*\bof\b',line)
            if m:
                ident,j,p=m.groups();result.append((int(ident),2*int(j),p));waiting=False
    if waiting:raise ValueError('Missing explicit stdout LS identity')
    return result


def label_rows(text,nlev,expected=None,summary=None):
    if expected is not None and (len(expected)!=nlev or summary!=expected):raise ValueError('Explicit stdout/binary ID order differs')
    rows={};current=None;position=0
    for line in text.splitlines():
        m=re.match(r'^\s*(\d+|\*{3})\s+(\d+)\s+([+-])\s+(-[\d.]+)\s+([\d.]+)%',line)
        if m:
            token,j,p,e,norm=m.groups();j=2*int(j)
            if expected is None:
                if token=='***':raise ValueError('Unqualified overflowed display ID')
                ident=int(token)
            else:
                if position>=nlev:raise ValueError('Extra LS label header')
                ident,ej,ep=expected[position]
                if (j,p)!=(ej,ep):raise ValueError('LS label/binary J/parity order differs')
                if (token=='***' and ident<=999) or (token!='***' and int(token)!=ident):raise ValueError('LS display ID inconsistent with verified native ID')
            if ident in rows:raise ValueError('Duplicate original LS ID')
            position+=1
            current=dict(two_J=j,parity=p,energy=e,display_overflow=token=='***',components=[]);rows[ident]=current
        else:
            m=re.match(r'^\s*(-?[\d.]+)\s+([\d.]+)\s+(1s\(.+)$',line)
            if m:
                if current is None:raise ValueError('Unheaded component')
                current['components'].append(m.groups())
    if set(rows)!=set(range(nlev)):raise ValueError('Incomplete LS IDs')
    return rows


def product(values,matrix):
    return values@matrix


def review(source,converted,shells=N3):
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Expanded angular replay belongs on Sirius')
    import numpy as np
    signed_fallbacks=0;weight_fallbacks=0
    cm=records(source['o_export_000.cm']);nele,ncsf,norb,nlev,nvalues,nblocks=struct.unpack('<6i',cm[1])
    ls=ls_blocks(converted['state.lsj.c'].decode(),shells)
    expected_ids=[]
    for b in range(nblocks):
        _,nb,nv,jp,p=struct.unpack('<5i',cm[2+4*b])
        expected_ids.extend((i,jp-1,'+' if p==1 else '-') for i in struct.unpack('<'+'i'*nv,cm[3+4*b]))
    explicit_ids=stdout_ids(converted['jj2lsj.stdout'].decode())
    lab=label_rows(converted['state.lsj.lbl'].decode(),nlev,expected_ids,explicit_ids)
    t=records(converted['state.lsj.T'])
    if t[0]!=b'jj2lsj' or struct.unpack('<4i',t[1])!=(nele,ncsf,norb,nblocks) or len(ls)!=nblocks:
        raise ValueError('Expanded angular header differs')
    pos=2+nblocks;checked=0;maxorth=0.;result={};jlines=converted['state.lsj.j'].decode().splitlines();jpos=0
    for b in range(nblocks):
        _,nb,nv,jp,p=struct.unpack('<5i',cm[2+4*b]);names=[r['name'] for r in ls[b]]
        if len(names)!=nb or len(set(names))!=nb:raise ValueError('LS labels ambiguous or incomplete')
        columns={name:i for i,name in enumerate(names)}
        if struct.unpack('<4i',t[2+b])!=(b+1,nb,nv,jp):raise ValueError('Angular block differs')
        marker=t[pos];pos+=1
        if marker[:18]!=b' *   Block Number=' or struct.unpack('<i',marker[18:])[0]!=b+1:raise ValueError('Angular marker differs')
        matrix=np.full((nb,nb),np.nan)
        for _ in range(nb*nb):
            jj,l,coef=struct.unpack('<iid',t[pos]);pos+=1
            if not 1<=jj<=nb or not 1<=l<=nb or not np.isnan(matrix[jj-1,l-1]) or not math.isfinite(coef):raise ValueError('Angular coefficient identity differs')
            matrix[jj-1,l-1]=coef
        maxorth=max(maxorth,float(np.max(np.abs(matrix.T@matrix-np.eye(nb)))))
        ids=struct.unpack('<'+'i'*nv,cm[3+4*b]);energies=struct.unpack('<'+'d'*(nv+1),cm[4+4*b]);values=np.frombuffer(cm[5+4*b],dtype='<f8').reshape(nv,nb)
        products=product(values,matrix)
        def legacy_dot(row,column):
            return sum(float(values[row,k])*float(matrix[k,column]) for k in range(nb))
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
            expected=products[i].tolist()
            for l,y in enumerate(printed):
                if Decimal(format(expected[l],'.8f'))!=Decimal(y):
                    expected[l]=legacy_dot(i,l);signed_fallbacks+=1
                    if Decimal(format(expected[l],'.8f'))!=Decimal(y):raise ValueError('Expanded LS coefficient/sign replay differs')
            checked+=nb;entry=lab[ident];parity='+' if p==1 else '-';energy=format(energies[0]+energies[i+1],'.9f')
            if (entry['two_J'],entry['parity'],entry['energy'],row[1])!=(jp-1,parity,energy,energy):raise ValueError('LS original ID/J/parity/energy differs')
            dominant=max(range(nb),key=lambda l:float(printed[l])**2)
            if row[2]!=names[dominant] or entry['components'][0][2]!=names[dominant]:raise ValueError('Dominant packed LS name differs')
            for coef,weight,name in entry['components']:
                if name not in columns:raise ValueError('Unknown LS component label')
                l=columns[name]
                if Decimal(coef)!=Decimal(printed[l]):raise ValueError('LS label coefficient/weight differs')
                if Decimal(weight)!=Decimal(format(expected[l]**2,'.8f')):
                    expected[l]=legacy_dot(i,l);weight_fallbacks+=1
                    if Decimal(weight)!=Decimal(format(expected[l]**2,'.8f')):raise ValueError('LS label coefficient/weight differs')
            result[str(ident)]=dict(two_J=jp-1,parity=parity,energy_au=energies[0]+energies[i+1],dominant=names[dominant],dominant_printed_weight=float(printed[dominant])**2,native_component_count=len(entry['components']),label_display_id_overflow=entry['display_overflow'])
    if pos!=len(t) or len(result)!=nlev:raise ValueError('Unconsumed angular records/levels')
    targets={term:[int(i) for i,r in result.items() if r['dominant']=='1s(2).2s(2).2p(3)4S.3s_'+term] for term in ('3S','5S')}
    if any(len(v)!=1 for v in targets.values()):raise ValueError('Expanded target parent/term mapping ambiguous')
    return dict(levels=result,targets=targets,full_LS_coefficients_reconstructed=checked,max_transformation_orthonormality_residual=maxorth,label_ID_overflows=sum(r['label_display_id_overflow'] for r in result.values()),original_ID_authority='Native binary CM IDs, independently matched to explicit I4 LS stdout IDs/order/J/parity; I3 label overflow preserved',signed_rounding_fallbacks=signed_fallbacks,weight_rounding_fallbacks=weight_fallbacks,orthonormality_arithmetic='NumPy BLAS Gram diagnostic only; no acceptance threshold. Signed coefficients/weights retain exact old Python-sum fallback on any native-precision discrepancy.')
