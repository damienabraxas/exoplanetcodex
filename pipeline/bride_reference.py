"""Read exact held comparison inputs; no RT/SE or atmosphere interpolation."""
from pathlib import Path
import hashlib
import numpy as np


def fingerprint(path):
    path=Path(path)
    return dict(name=path.name,bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def read_marcs_structure(path):
    """Read the two standard MARCS structure tables in their printed precision.

    Source Depth is inward cm; output height is outward m. Pe is dyn/cm2;
    density is g/cm3. Hydrogen total is rho*X/mH, not neutral H or NLTE pops.
    Only data/source units are parsed here; no EOS or physics is substituted.
    """
    text=Path(path).read_text();lines=text.splitlines()
    model=lines[0].strip()
    count_line=next(i for i,l in enumerate(lines)if 'Number of depth points' in l)
    n=int(lines[count_line].split()[0])
    first=next(i for i,l in enumerate(lines)if 'lgTauR' in l and 'lgTau5' in l)
    second=next(i for i,l in enumerate(lines)if 'KappaRoss' in l and 'Density' in l)
    a=np.asarray([[float(x)for x in l.split()]for l in lines[first+1:first+1+n]])
    b=np.asarray([[float(x)for x in l.split()]for l in lines[second+1:second+1+n]])
    if a.shape!=(n,9)or b.shape!=(n,8)or not np.isfinite(a).all()or not np.isfinite(b).all():
        raise ValueError('Malformed/nonfinite MARCS structure tables')
    if not np.array_equal(a[:,0],np.arange(1,n+1))or not np.array_equal(a[:,0],b[:,0]):
        raise ValueError('MARCS structure row identities differ')
    if not np.array_equal(a[:,1],b[:,1]):
        raise ValueError('MARCS optical-depth identities differ across tables')
    abundance_start=next(i for i,l in enumerate(lines)if 'Logarithmic chemical number abundances' in l)
    abundances=[float(v)for line in lines[abundance_start+1:count_line]for v in line.split()]
    if len(abundances)!=92 or not np.isfinite(abundances).all()or abundances[0]!=12.:
        raise ValueError('Expected 92 finite MARCS elemental abundances with H=12')
    X=float(next(l for l in lines if 'are X, Y' in l).split()[0])
    gravity=float(next(l for l in lines if 'Surface gravity' in l).split()[0])
    vturb=float(next(l for l in lines if 'Microturbulence parameter' in l).split()[0])*1000.
    if not (0<X<1 and gravity>0 and vturb>=0 and np.all(a[:,4:7]>0) and np.all(b[:,3]>0)):
        raise ValueError('Invalid MARCS physical fields')
    if not np.all(np.diff(a[:,3])>0):
        raise ValueError('MARCS depth must increase inward')
    return dict(model_id=model,ndep=n,height_m=-a[:,3]*.01,
        temperature_K=a[:,4].copy(),electron_pressure_Pa=a[:,5]*.1,
        rho_kg_m3=b[:,3]*1000.,hydrogen_mass_fraction=X,
        abundances_dex=abundances,logg_cgs=float(np.log10(gravity)),vturb_m_s=vturb,
        logtau_ross=a[:,1].copy(),logtau_500=a[:,2].copy())


def read_exact_departure_record(binary,aux_text,*,model_id,abundance,nlevels,ndep):
    """Read one native MULTI grid record with exact identity, never nearest-node.

    The verified Gerber layout is 500 ASCII bytes, two little-endian i4,
    ndep linear tau doubles, then (level,depth) little-endian doubles.
    Return raw-record hash; do not copy/hash the multi-GB grid onto the Mac.
    """
    hits=[]
    for line in aux_text.splitlines():
        fields=line.split()
        if len(fields)>=9 and not line.lstrip().startswith('#'):
            if fields[0].strip("'")==model_id and float(fields[7])==abundance:
                hits.append(fields)
    if len(hits)!=1:
        raise ValueError(f'Expected one exact model/abundance record, found {len(hits)}')
    row=hits[0];offset=int(row[8])-1
    if offset<1000:raise ValueError('Invalid one-based record pointer')
    size=500+8+8*ndep*(1+nlevels)
    with Path(binary).open('rb')as f:f.seek(offset);raw=f.read(size)
    if len(raw)!=size:raise ValueError('Truncated native departure record')
    identity=raw[:500].decode('ascii').strip()
    dimensions=np.frombuffer(raw[500:508],dtype='<i4')
    if identity!=model_id or tuple(dimensions)!=(ndep,nlevels):
        raise ValueError('Native record identity or dimensions do not match requested input')
    tau=np.frombuffer(raw[508:508+8*ndep],dtype='<f8').copy()
    populations=np.frombuffer(raw[508+8*ndep:],dtype='<f8').reshape(nlevels,ndep).copy()
    if not np.isfinite(tau).all()or np.any(tau<=0)or not np.all(np.diff(tau)>0):
        raise ValueError('Invalid native linear optical-depth grid')
    if not np.isfinite(populations).all()or np.any(populations<0):
        raise ValueError('Invalid native departure coefficients')
    return dict(model_id=identity,abundance=abundance,pointer_1based=offset+1,
        record_sha256=hashlib.sha256(raw).hexdigest(),tau=tau,departures=populations)
