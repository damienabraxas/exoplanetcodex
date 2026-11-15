"""Source-mapped principal-complex identity plus full streamed coefficient transport."""
from collections import Counter
from itertools import combinations, product
import math
import mmap
from pathlib import Path
import re
import struct
from rya1354_fac_mapping import csfs, fac_name

N3=(('2p',1,6),('3s',0,2),('3p',1,6),('3d',2,10))
EXTENDED=N3+(('4s',0,2),('4p',1,6))


def configurations(shells):
    return [q for q in product(*(range(min(cap,4)+1) for _,_,cap in shells)) if sum(q)==4]


def inventory(shells,declared):
    counts=Counter()
    for q in declared:
        states=Counter({0:1});parity=0
        for (_,ell,_),n in zip(shells,q):
            parity+=ell*n
            shell=Counter(sum(c) for c in combinations([2*m+s for m in range(-ell,ell+1) for s in (-1,1)],n))
            out=Counter()
            for x,a in states.items():
                for y,b in shell.items():out[x+y]+=a*b
            states=out
        for j in range(0,max(states)+1,2):
            if states[j]!=states[j+2]:counts[(parity%2,j)]+=states[j]-states[j+2]
    return counts


def record_ranges(data):
    pos=0;result=[]
    while pos<len(data):
        if pos+8>len(data):raise ValueError('Truncated record marker')
        size=struct.unpack_from('<i',data,pos)[0]
        if size<0 or pos+size+8>len(data) or struct.unpack_from('<i',data,pos+size+4)[0]!=size:raise ValueError('Invalid record marker')
        result.append((pos+4,size));pos+=size+8
    return result


def transport(basis_path,c_path,cm_path,shells=N3):
    declared=configurations(shells);expected=inventory(shells,declared)
    c_text=Path(c_path).read_text();exported=csfs(c_text)
    lines=[x for x in c_text.splitlines()[5:] if not x.startswith(' *') and not x.startswith('*')]
    nrnames=[];complexes=[];occupations=Counter()
    for pos in range(0,len(lines),3):
        # Relativistic +/- components must be added, not overwritten.
        occ=Counter()
        for n,l,q in re.findall(r'(\d+)([spdfg])[ -]\(\s*(\d+)\)',lines[pos]):occ[int(n),l]+=int(q)
        principal=Counter()
        for (n,l),q in occ.items():principal[n]+=q
        complexes.append('.'.join(str(n)+'*'+str(q) for n,q in sorted(principal.items()) if q!=2*n*n))
        nrnames.append('.'.join(str(n)+l+str(q) for (n,l),q in sorted(occ.items(),key=lambda x:(x[0][0],'spdfg'.index(x[0][1]))) if q!=4*'spdfg'.index(l)+2))
        if not nrnames[-1]:
            (n,l),q=max(occ.items(),key=lambda x:(x[0][0],'spdfg'.index(x[0][1])))
            nrnames[-1]=str(n)+l+str(q)
        if occ.pop((1,'s'),0)!=2 or occ.pop((2,'s'),0)!=2:raise ValueError('Fixed core occupation differs')
        val=tuple(occ.pop((int(name[:-1]),name[-1]),0) for name,_,_ in shells)
        if occ or val not in declared:raise ValueError('Undeclared occupation')
        occupations[val]+=1
    basis={};native_slots=set();seen_header=False
    with Path(basis_path).open() as src:
        for line in src:
            if 'Basis Table' in line:seen_header=True;continue
            if 'Mixing Coefficients' in line:break
            if not line.strip() or line.startswith('#'):continue
            f=line.split()
            if not seen_header or len(f)!=10:raise ValueError('Native basis field count differs')
            sym,p,j,k=map(int,f[:4]);key=(sym,f[7],f[-2],fac_name(f[-1]))
            if sym!=2*j+p or key in basis or (sym,k) in native_slots:raise ValueError('Native coupled basis duplicate/identity differs')
            basis[key]=k;native_slots.add((sym,k))
    with Path(cm_path).open('rb') as handle,mmap.mmap(handle.fileno(),0,access=mmap.ACCESS_READ) as data:
        rr=record_ranges(data)
        if len(rr)<2 or data[rr[0][0]:rr[0][0]+rr[0][1]]!=b'G92MIX' or rr[1][1]!=24:raise ValueError('Mixing signature/header differs')
        nele,ncsf,norb,nlev,nvalues,nblocks=struct.unpack_from('<6i',data,rr[1][0])
        if nele!=8 or ncsf!=nlev or nlev!=sum(expected.values()) or len(rr)!=2+4*nblocks or len(exported)!=ncsf or len(lines)!=3*ncsf:raise ValueError('Mixing aggregate dimensions differ')
        level_map={};blocks=[];actual=Counter();offset=0;total=0;used_native=set()
        for b in range(nblocks):
            h=2+4*b
            if rr[h][1]!=20:raise ValueError('Block header size differs')
            number,nb,nv,jp,parity=struct.unpack_from('<5i',data,rr[h][0])
            if number!=b+1 or nb!=nv or nb<=0 or jp<1 or parity not in (-1,1) or [rr[h+i][1] for i in (1,2,3)]!=[4*nv,8*(nv+1),8*nb*nv]:raise ValueError('Block dimensions differ')
            p=0 if parity==1 else 1;j=jp-1;sym=2*j+p;actual[p,j]+=nv
            columns={}
            for column,descriptor in enumerate(exported[offset:offset+nb]):
                key=(sym,complexes[offset+column],nrnames[offset+column],descriptor)
                if key not in basis or basis[key] in columns:raise ValueError('CSF identity permutation differs')
                columns[basis[key]]=column;used_native.add(key)
            if set(columns)!=set(range(nb)):raise ValueError('Native basis permutation not bijective')
            for x in struct.iter_unpack('<d',data[rr[h+2][0]:rr[h+2][0]+rr[h+2][1]]):
                if not math.isfinite(x[0]):raise ValueError('Nonfinite energy')
            seen=bytearray(nb*nv);start=rr[h+3][0]
            blocks.append((start,nb,nv,seen))
            ids=struct.unpack_from('<'+'i'*nv,data,rr[h+1][0])
            for row,ident in enumerate(ids):
                if ident in level_map:raise ValueError('Duplicate exported level identity')
                level_map[ident]=(b,row,sym,p,j,columns)
            offset+=nb;total+=nb*nv
        if offset!=ncsf or total!=nvalues or actual!=expected or set(level_map)!=set(range(nlev)) or used_native!=set(basis) or set(occupations)!=set(declared):raise ValueError('Aggregate/per-symmetry identity inventory differs')
        checked=0;mixing=False
        with Path(basis_path).open() as src:
            for line in src:
                if 'Mixing Coefficients' in line:mixing=True;continue
                if not mixing or not line.strip() or line.startswith('#'):continue
                f=line.split()
                if len(f)!=11:raise ValueError('Native coefficient field count differs')
                ident,sym,p,j,rank,k=map(int,f[:6])
                if ident not in level_map:raise ValueError('Unknown native level identity')
                b,row,es,ep,ej,columns=level_map[ident];start,nb,nv,seen=blocks[b]
                if (sym,p,j)!=(es,ep,ej) or k not in columns or not 0<=rank<nb:raise ValueError('Native coefficient identity differs')
                index=row*nb+columns[k]
                if seen[index]:raise ValueError('Duplicate native coefficient')
                value=struct.unpack_from('<d',data,start+8*index)[0]
                if not math.isfinite(value) or format(value,'.8E')!=f[9]:raise ValueError('Signed coefficient transport differs')
                seen[index]=1;checked+=1
        for start,nb,nv,seen in blocks:
            for index in range(nb*nv):
                value=struct.unpack_from('<d',data,start+8*index)[0]
                if not math.isfinite(value) or (not seen[index] and value!=0):raise ValueError('Unprinted nonzero/nonfinite coefficient')
    return dict(levels=nlev,CSFs=ncsf,orbitals=norb,blocks=nblocks,signed_printed_coefficients_verified=checked,declared_configurations_present=len(occupations),LS_mapping_admitted=False)
