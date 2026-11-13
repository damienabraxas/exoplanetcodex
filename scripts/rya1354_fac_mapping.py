"""Exact FAC native-to-export coefficient transport before any LS recoupling."""
from pathlib import Path
import re
import struct

from rya1354_fac_control_audit import mixing


def fac_name(name):
    result=[]
    for item in name.split('.'):
        m=re.fullmatch(r'(\d+)([spdfg])([+-])(\d+)\((\d+)\)(\d+)',item)
        if m is None:
            raise ValueError('Unsupported FAC coupled-name layout')
        n,l,sign,q,sj,tj=m.groups()
        result.append((int(n),l,sign,int(q),int(sj),int(tj)))
    return tuple(result)


def twice_j(text):
    text=text.strip()
    if '/' in text:
        a,b=map(int,text.split('/'))
        if 2*a%b:
            raise ValueError('Non half-integer CSF angular momentum')
        return 2*a//b
    return 2*int(text)


def csfs(text):
    lines=[x for x in text.splitlines()[5:] if not x.startswith(' *') and not x.startswith('*')]
    if len(lines)%3:
        raise ValueError('Incomplete exported CSF triples')
    result=[]
    for offset in range(0,len(lines),3):
        a,b,c=lines[offset:offset+3];items=[];closed_last=None
        if len(a)%9:
            raise ValueError('Exported FAC CSF orbital width differs')
        for index in range(len(a)//9):
            m=re.fullmatch(r'\s*(\d+)([spdfg])([ -])\(\s*(\d+)\)',a[index*9:(index+1)*9])
            if m is None:
                raise ValueError('Exported FAC CSF orbital field differs')
            n,l,minus,q=m.groups();n=int(n);q=int(q);ell='spdfg'.index(l)
            sjmax=2*ell+(-1 if minus=='-' else 1)
            if q==sjmax+1:
                closed_last=(n,l,'-' if minus=='-' else '+',q,0,0)
                continue
            sj=twice_j(b[index*9:(index+1)*9])
            total=c[1+index*9:1+(index+1)*9].strip()
            tj=(twice_j(total) if total else sj if not items else items[-1][-1] if sj==0 else None)
            if tj is None:
                raise ValueError('Exported open-shell coupling absent')
            items.append((n,l,'-' if minus=='-' else '+',q,sj,tj))
        # ConstructLevelName retains the last full shell when no open shell exists.
        result.append(tuple(items) if items else (closed_last,))
    return result


def verify_transport(basis_text,c_text,cm_bytes):
    meta=mixing(cm_bytes);basis={};coefficients={};mode=''
    for line in basis_text.splitlines():
        if 'Basis Table' in line:mode='basis';continue
        if 'Mixing Coefficients' in line:mode='mixing';continue
        if not line.strip() or line.startswith('#'):continue
        f=line.split()
        if mode=='basis':
            if len(f)!=10:
                raise ValueError('Native FAC basis field count differs')
            symmetry,parity,twoj,k=map(int,f[:4]);key=(symmetry,fac_name(f[-1]))
            if key in basis:
                raise ValueError('Duplicate native FAC basis identity')
            basis[key]=k
        elif mode=='mixing':
            if len(f)!=11:
                raise ValueError('Native FAC mixing field count differs')
            level,symmetry,parity,twoj,rank,k=map(int,f[:6])
            key=(level,symmetry,k)
            if key in coefficients:raise ValueError('Duplicate native FAC mixing entry')
            coefficients[key]=f[9]
    records=[];offset=0
    while offset<len(cm_bytes):
        size=struct.unpack_from('<i',cm_bytes,offset)[0]
        records.append(cm_bytes[offset+4:offset+size+4]);offset+=size+8
    exported=csfs(c_text);start=0;checked=0;used=set()
    for index,block in enumerate(meta['blocks']):
        nb=block['basis'];nv=block['levels'];p=0 if block['parity']==1 else 1
        symmetry=2*block['two_J']+p
        order=[]
        for descriptor in exported[start:start+nb]:
            key=(symmetry,descriptor)
            if key not in basis:raise ValueError('Exported CSF differs from native FAC basis')
            order.append(basis[key])
        if len(set(order))!=nb:raise ValueError('Exported FAC CSF basis order is not bijective')
        values=struct.unpack('<'+'d'*(nb*nv),records[2+4*index+3])
        for i,level in enumerate(block['indices']):
            for j,k in enumerate(order):
                value=values[i*nb+j];key=(level,symmetry,k)
                if key not in coefficients:
                    if value!=0.0:raise ValueError('Nonzero exported FAC coefficient lacks native source')
                    continue
                if format(value,'.8E')!=coefficients[key]:
                    raise ValueError('FAC export coefficient/sign/order differs from native')
                checked+=1;used.add(key)
        start+=nb
    if start!=len(exported) or used!=set(coefficients):
        raise ValueError('FAC basis/mixing transport incomplete')
    return dict(native_basis_count=len(basis),exported_CSFs=len(exported),
        signed_coefficients_checked=checked,level_indices_preserved=15,
        precision='Exact FAC %.8E native-print formatting equality against exported binary doubles',
        scope='Coupled-basis identity/order/sign transport only; no LS or physical accuracy admission')
