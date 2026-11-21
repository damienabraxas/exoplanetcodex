"""Expanded-space native/export transport and inventory replay; no LS/rates."""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import struct
from rya1354_native_parts import load_parts
from rya1354_fac_recoupling_audit import records
from rya1354_fac_mapping import fac_name,csfs
from rya1354_fac_space import configurations,expected_inventory,input_text,inspect
from rya1354_fac_controls import diagnostics


def transport(basis_text,c_text,cm,rank):
    rr=records(cm)
    if rr[0]!=b'G92MIX':raise ValueError('Mixing signature differs')
    nele,ncsf,norb,nlev,nvalues,nblocks=struct.unpack('<6i',rr[1])
    expected=expected_inventory(configurations(rank))
    if nele!=8 or nlev!=sum(expected.values()) or ncsf!=nlev or len(rr)!=2+4*nblocks:
        raise ValueError('Expanded mixing dimensions differ')
    basis={};coefficients={};mode=''
    for line in basis_text.splitlines():
        if 'Basis Table' in line:mode='basis';continue
        if 'Mixing Coefficients' in line:mode='mixing';continue
        if not line.strip() or line.startswith('#'):continue
        f=line.split()
        if mode=='basis':
            if len(f)!=10:raise ValueError('Native basis field count differs')
            symmetry,parity,j,k=map(int,f[:4]);key=(symmetry,f[-2],fac_name(f[-1]))
            if key in basis:raise ValueError('Duplicate coupled native basis')
            basis[key]=k
        elif mode=='mixing':
            if len(f)!=11:raise ValueError('Native coefficient field count differs')
            level,symmetry,parity,j,order,k=map(int,f[:6]);key=(level,symmetry,k)
            if key in coefficients:raise ValueError('Duplicate native coefficient')
            coefficients[key]=f[9]
    exported=csfs(c_text);lines=[x for x in c_text.splitlines()[5:] if not x.startswith(' *') and not x.startswith('*')]
    if len(exported)!=ncsf or len(lines)!=3*ncsf:raise ValueError('Expanded CSF count differs')
    allowed=set(configurations(rank));occupation_counts=Counter()
    for i in range(0,len(lines),3):
        q=Counter()
        for n,l,num in re.findall(r'(\d+)([spdfg])[ -]\(\s*(\d+)\)',lines[i]):q[(int(n),l)]+=int(num)
        if q.pop((1,'s'),0)!=2 or q.pop((2,'s'),0)!=2:raise ValueError('Fixed core occupations differ')
        val=tuple(q.pop(key,0) for key in ((2,'p'),(3,'s'),(3,'p'),(3,'d')))
        if q or val not in allowed:raise ValueError('Undeclared expanded occupation')
        occupation_counts[val]+=1
    # ConstructLevelName suppresses closed NR shells in sname and closed jj shells
    # in name. Their pair distinguishes closed-shell occupations lost by name alone.
    nrnames=[]
    for i in range(0,len(lines),3):
        occ=Counter()
        for n,l,num in re.findall(r'(\d+)([spdfg])[ -]\(\s*(\d+)\)',lines[i]):occ[(int(n),l)]+=int(num)
        nrnames.append('.'.join(str(n)+l+str(q) for (n,l),q in sorted(occ.items(),key=lambda x:(x[0][0],'spdfg'.index(x[0][1])))
                                if q!=2*(2*'spdfg'.index(l)+1)))
    offset=0;ids=[];used=set();checked=0;actual=Counter();totalvalues=0
    for block in range(nblocks):
        pos=2+4*block;number,nb,nv,jp,parity=struct.unpack('<5i',rr[pos])
        if number!=block+1 or nb!=nv or nb<=0 or jp<1 or parity not in (-1,1):raise ValueError('Expanded block header differs')
        if [len(rr[pos+i]) for i in (1,2,3)]!=[4*nv,8*(nv+1),8*nb*nv]:raise ValueError('Expanded block array sizes differ')
        indices=struct.unpack('<'+'i'*nv,rr[pos+1]);ids.extend(indices)
        energies=struct.unpack('<'+'d'*(nv+1),rr[pos+2]);values=struct.unpack('<'+'d'*(nb*nv),rr[pos+3])
        if not all(math.isfinite(x) for x in (*energies,*values)):raise ValueError('Nonfinite expanded eigenvector')
        p=0 if parity==1 else 1;j=jp-1;symmetry=2*j+p;actual[(p,j)]+=nv
        order=[]
        for local,descriptor in enumerate(exported[offset:offset+nb]):
            key=(symmetry,nrnames[offset+local],descriptor)
            if key not in basis:raise ValueError('Expanded CSF identity transport differs')
            order.append(basis[key])
        if len(set(order))!=nb:raise ValueError('Expanded basis permutation not bijective')
        for i,ident in enumerate(indices):
            for column,k in enumerate(order):
                key=(ident,symmetry,k);value=values[i*nb+column]
                if key not in coefficients:
                    if value!=0:raise ValueError('Unprinted nonzero expanded coefficient')
                else:
                    if format(value,'.8E')!=coefficients[key]:raise ValueError('Expanded signed coefficient transport differs')
                    used.add(key);checked+=1
        offset+=nb;totalvalues+=nb*nv
    if (offset!=ncsf or totalvalues!=nvalues or actual!=expected or len(ids)!=nlev
            or set(ids)!=set(range(nlev)) or used!=set(coefficients)):
        raise ValueError('Expanded aggregate/per-symmetry inventory differs')
    return dict(levels=nlev,CSFs=ncsf,orbitals=norb,blocks=nblocks,signed_printed_coefficients_verified=checked,
        declared_configurations_present=len(occupation_counts),LS_mapping_admitted=False)


def summarize(base):
    files,manifest=load_parts(base/'fac_space_complete_v23');root='runs/fac_space_v23/'
    receipt=json.loads(files[root+'receipt.json'])
    if receipt['status']!='FAC_EXPANDED_LEVELS_COMPLETED_REVIEW_REQUIRED':raise ValueError('Expanded run incomplete')
    control,_=load_parts(base/'fac_coherent_complete_v18')
    source=json.loads(control['runs/fac_coherent_controls_v18/receipt.json'])
    if receipt['pin']!=source['pin'] or receipt['binary_sha256']!=source['binary_sha256']:
        raise ValueError('Expanded FAC source/executable provenance differs')
    result=dict(schema='rya1354.fac-space-audit.v23',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,cases={})
    for rank in (1,2):
        for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
            name=f'n3_r{rank}_{label}';p=root+name+'/'
            if files[p+'input.sf'].decode()!=input_text(rank,mult,scale):raise ValueError('Expanded input differs')
            d=diagnostics(files[p+'sfac.stdout'].decode(),mult,256)
            if not d['printed_rule_verified']:raise ValueError('Expanded stopping guard failed')
            inspect(files[p+'o.lev'].decode(),rank)
            result['cases'][name]=dict(diagnostics=d,transport=transport(files[p+'o.basis'].decode(),
                files[p+'o_export_000.c'].decode(),files[p+'o_export_000.cm'],rank))
    result['receipt_sha256']=hashlib.sha256(files[root+'receipt.json']).hexdigest()
    return result

if __name__=='__main__':
    print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
