"""Read native TS opacity output (no appended MARCS abundance vector)."""
from pathlib import Path

import numpy as np


def native_real32_print(values):
    """Input default REAL, then native 1P E14.7 output; no fitted tolerance."""
    values=np.asarray(values,dtype=np.float32)
    return np.asarray([float(f'{x:.7E}') for x in values])


def read_ts_opacity(path):
    fields=Path(path).read_text().split()
    if len(fields)<4 or fields[0].strip("'")!='MRXF':raise ValueError('Not TS MRXF')
    n=int(fields[1]);standard=float(fields[2]);nw=int(fields[3])
    if not (1<=n<=400 and 1<=nw<=10000 and standard>0):raise ValueError('Bad TS dimensions')
    a=np.asarray([float(x.replace('D','E')) for x in fields[4:]])
    if len(a)!=nw+n*(8+2*nw) or not np.isfinite(a).all():raise ValueError('Malformed TS opacity')
    wave=a[:nw];rows=a[nw:].reshape(n,8+2*nw);s=rows[:,:8]
    # Plane-parallel output stores a signed geometric depth in column one;
    # it is not the positive normalized radius of the public MARCS OPA file.
    if (np.any(np.diff(wave)<=0) or np.any(wave<=0) or np.any(s[:,1:6]<=0)
            or np.any(s[:,6]<0) or np.any(s[:,7]<=0) or np.any(rows[:,8:]<0)):
        raise ValueError('Invalid TS physical output')
    absolute=rows[:,8:].reshape(n,nw,2)*s[:,7,None,None]*.1
    return dict(ndep=n,nwave=nw,wavelength_air_A=wave,standard_wavelength_air_A=standard,
        tau_standard=s[:,1],temperature_K=s[:,2],electron_pressure_Pa=s[:,3]*.1,
        gas_pressure_Pa=s[:,4]*.1,rho_kg_m3=s[:,5]*1000,vturb_m_s=s[:,6]*1000,
        absorption_m2_kg=absolute[:,:,0].T,scattering_m2_kg=absolute[:,:,1].T)


def read_ts_components(path,ndep,nwave):
    """Read isolated reporting-build output; require complete ordered blocks."""
    states=[];records=[];post=[];labels=None;pending=None
    for line in Path(path).read_text().splitlines():
        if 'BRIDE_POST_STATE' in line:
            values=[float(x.replace('D','E')) for x in line.split()[1:]]
            if len(values)!=6 or values[0]!=1:raise ValueError('Malformed post-DETABS state')
            post.append(values[1:])
        elif 'BRIDE_FE_STATE' in line:
            fields=line.split();values=[float(x.replace('D','E')) for x in fields[1:]]
            if len(values)!=7 or values[0]!=1:raise ValueError('Malformed Fe state')
            states.append(values[1:])
        elif line.startswith('0WAVEL.='):
            prefix,names=line.split('SCAT  ',1)
            names=[names[i:i+10].strip() for i in range(0,len(names),10)]
            names=[n for n in names if n]
            if labels is None:labels=names
            if names!=labels or len(set(labels))!=len(labels):raise ValueError('Ambiguous component labels')
            pending=float(prefix.split('=',1)[1].split()[0])
        elif line.startswith(' T='):
            values=[float(x.replace('D','E')) for x in line.split()[1:]]
            if pending is None or labels is None or len(values)!=3+len(labels):
                raise ValueError('Incomplete component record')
            records.append([pending,*values]);pending=None
    if len(states)!=ndep or len(records)!=ndep*(nwave+1) or pending is not None:
        raise ValueError('Incomplete component depth/wavelength blocks')
    states=np.asarray(states);records=np.asarray(records).reshape(ndep,nwave+1,-1)
    if not np.isfinite(states).all() or not np.isfinite(records).all() or np.any(states<=0):
        raise ValueError('Invalid component state')
    if np.any(records[:,:,2:]<0):raise ValueError('Negative component coefficient')
    fractions=records[:,:,4:]
    # Native PROV fractions are default REAL; bounded summation error, not an
    # empirically chosen normalization or correction to the retained values.
    eps=np.finfo(np.float32).eps;bound=2*len(labels)*eps/(1-len(labels)*eps)
    for subset in (fractions[:,:,:-4],fractions[:,:,-4:]):
        np.testing.assert_allclose(subset.sum(axis=2),1.,rtol=0,atol=bound)
    if post:
        post=np.asarray(post).reshape(ndep,nwave+1,5)
        if (not np.isfinite(post).all() or np.any(post[:,:,0]<=0)
                or np.any(post[:,:,1:]<0) or np.any(post[:,:,3]<=0)):
            raise ValueError('Invalid post-DETABS state')
    else:post=None
    return dict(labels=labels,states=states,records=records,fractions=fractions,
        post=post,
        state_fields=['temperature_K','electron_pressure_dyn_cm2','rho_g_cm3',
                      'FeI_partition','FeI_pressure_dyn_cm2','FeII_pressure_dyn_cm2'])
