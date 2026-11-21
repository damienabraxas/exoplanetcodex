"""Diagnostic snapshots of complete native ORBOUT wavefunctions; never a stop waiver.

Pinned orbout.f90: G92RWF header, (n,kappa,E,MF), (PZ,P,Q), radius
records. Native little-endian gfortran with four-byte sequential record markers.
Compare saved post-normalization radial functions, not a claimed physical error.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import re
import struct


def decode(data):
    records=[];offset=0
    while offset<len(data):
        if offset+8>len(data):
            raise ValueError('Truncated native record')
        size=struct.unpack_from('<i',data,offset)[0]
        if size<0 or offset+8+size>len(data) or struct.unpack_from('<i',data,offset+4+size)[0]!=size:
            raise ValueError('Incomplete or incompatible native sequential record')
        records.append(data[offset+4:offset+4+size]);offset+=size+8
    if not records or records[0]!=b'G92RWF' or (len(records)-1)%3:
        raise ValueError('Wrong native orbital-file identity')
    orbitals={}
    for i in range(1,len(records),3):
        if len(records[i])!=20:
            raise ValueError('Unexpected native orbital header')
        n,k,e,mf=struct.unpack('<iidi',records[i])
        if n<1 or k==0 or abs(k)>n or mf<=0 or len(records[i+1])!=8*(1+2*mf) or len(records[i+2])!=8*mf:
            raise ValueError('Wrong native orbital-array shape')
        values=struct.unpack('<'+'d'*(1+2*mf),records[i+1])
        radius=struct.unpack('<'+'d'*mf,records[i+2])
        if not all(math.isfinite(v) for v in (e,*values,*radius)) or (n,k) in orbitals:
            raise ValueError('Nonfinite or duplicate native orbital')
        orbitals[(n,k)]=dict(energy_parameter=e,pz=values[0],p=values[1:mf+1],q=values[mf+1:],radius=radius)
    return orbitals


def difference(a,b):
    if set(a)!=set(b):
        raise ValueError('Orbital inventories differ')
    rows=[]
    for key in a:
        x,y=a[key],b[key];common=min(len(x['radius']),len(y['radius']))
        if x['radius'][:common]!=y['radius'][:common]:
            raise ValueError('Radial grid changed across snapshots')
        m=max(len(x['p']),len(y['p']))
        xp=(*x['p'],*((0.,)*(m-len(x['p']))));xq=(*x['q'],*((0.,)*(m-len(x['q']))))
        yp=(*y['p'],*((0.,)*(m-len(y['p']))));yq=(*y['q'],*((0.,)*(m-len(y['q']))))
        delta=max(abs(v-w) for v,w in zip((*xp,*xq),(*yp,*yq)))
        rows.append(dict(n=key[0],kappa=key[1],max_saved_radial_change=delta,pz_change=y['pz']-x['pz'],
                         energy_parameter_change=y['energy_parameter']-x['energy_parameter']))
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    raw=None
    for attempt in range(10):
        try:
            data=(args.source/'rwfn.out').read_bytes()
            decode(data)
            if data==(args.source/'rwfn.out').read_bytes():
                raw=data;break
        except ValueError:
            continue
    if raw is None:
        raise ValueError('No complete stable native wavefunction read; no snapshot written')
    index=len(list(args.out.glob('snapshot*.w')))+1
    name=f'snapshot{index:03}.w'
    with (args.out/name).open('xb') as file:file.write(raw)
    iterations=re.findall(r'Iteration number\s+(\d+)',(args.source/'rmcdhf.stdout').read_text())
    record=dict(schema='rya1354.native-wave-snapshot.v1',source=str(args.source),captured_utc=datetime.now(timezone.utc).isoformat(),
                sha256=hashlib.sha256(raw).hexdigest(),printed_iteration_at_capture=int(iterations[-1]) if iterations else None,
                timing='Asynchronous complete/stable ORBOUT file read; exact cycle-to-wave timing not asserted',
                scope='POST_NORMALIZATION_WAVE_CHANGE_DIAGNOSTIC_ONLY_NOT_AN_ACCURACY_OR_STOPPING_WAIVER')
    if index>1:
        previous=args.out/f'snapshot{index-1:03}.w'
        record['versus_previous']=difference(decode(previous.read_bytes()),decode(raw))
    (args.out/f'snapshot{index:03}.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record))


if __name__=='__main__':
    main()
