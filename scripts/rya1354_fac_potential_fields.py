"""Pinned SavePotential/SetPotDP layout: protect all except named W scratch arrays."""
import hashlib
import math
import re
import struct


def layout(radial,orbital,consts):
    writer=radial.split('int SavePotential(char *fn) {',1)[1].split('\nint ModifyPotential(',1)[0]
    fields=re.findall(r'fwrite\(&p->(\w+), sizeof\((int|double)\), 1, f\)',writer)
    if fields[:5]!=[(k,'int') for k in ('mode','flag','r_core','nmax','maxrp')]:raise ValueError('Pinned potential scalar layout differs')
    header=sum(4 if kind=='int' else 8 for _,kind in fields)
    if header!=348 or 'fwrite(p->dws, sizeof(double), p->nws, f)' not in writer:raise ValueError('Pinned potential serialized header differs')
    dp=radial.split('void SetPotDP(POTENTIAL *p) {',1)[1].split('\nvoid AllocDWS',1)[0]
    arrays=re.findall(r'p->(\w+) = x;',dp)
    if len(arrays)!=36 or arrays[22:25]!=['W','dW','dW2']:raise ValueError('Pinned potential work array layout differs')
    nk=int(re.search(r'^#define NKSEP\s+(\d+)',consts,re.M)[1]);nk1=int(re.search(r'^#define NKSEP1\s+(\d+)',consts,re.M)[1])
    if (nk,nk1)!=(5,6):raise ValueError('Pinned potential array count differs')
    func=orbital.split('int SetPotentialW(POTENTIAL *pot, double e, int kappa, int k) {',1)[1].split('\nint SetPotentialWZMC',1)[0]
    if any(s not in func for s in ('xi = e - w*pot->VT[k][i];','- 2.0*kappa*x/pot->rad[i]','Differential(pot->W, pot->dW','Differential(pot->dW, pot->dW2')):raise ValueError('Energy/kappa-dependent W derivation differs')
    return dict(header=header,arrays=36+3*nk+3*nk1,work_arrays={k:arrays.index(k) for k in ('W','dW','dW2')})


def compare(source,other,spec):
    if len(source)<20 or len(other)!=len(source):raise ValueError('Potential lengths differ')
    n=struct.unpack_from('<i',source,16)[0];h=spec['header'];count=spec['arrays']
    if n<=0 or len(source)!=h+8*n*count:raise ValueError('Potential payload length/grid differs')
    if source[:h]!=other[:h]:raise ValueError('Protected potential header changed')
    changed=[];protected=bytearray(source[:h])
    allowed={v:k for k,v in spec['work_arrays'].items()}
    for i in range(count):
        a=source[h+i*n*8:h+(i+1)*n*8];b=other[h+i*n*8:h+(i+1)*n*8]
        if i not in allowed:
            if a!=b:raise ValueError('Protected potential/grid array changed: '+str(i))
            protected.extend(a)
        else:
            if not all(math.isfinite(x[0]) for x in struct.iter_unpack('<d',a+b)):raise ValueError('Nonfinite derived W scratch')
            if a!=b:changed.append(allowed[i])
    return dict(serialized_byte_identical=source==other,protected_header_grid_potential_byte_identical=True,
                derived_energy_kappa_work_arrays_changed=changed,protected_sha256=hashlib.sha256(protected).hexdigest(),
                scope='Only source-proven energy/kappa-dependent W,dW,dW2 may vary; every other byte protected. Raw saved files unchanged/retained.')
