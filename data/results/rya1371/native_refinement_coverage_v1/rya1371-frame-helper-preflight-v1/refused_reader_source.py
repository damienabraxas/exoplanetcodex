"""Independent strict reader for the pinned OpenMULTI3D1D diagnostic formats.

Layout owner: atom_mod.f90%write_rates and patch_t%save, at the engine pin in
config/reference_engines.yaml. No inferred endian, axis, element or LTE fallback.
"""
from pathlib import Path
import numpy as np


def take(stream,dtype,count):
    data=np.fromfile(stream,dtype=dtype,count=int(count))
    if data.size!=count:raise ValueError('Truncated native reference array')
    return data


def read_rate_pair(radiative,collision):
    with Path(radiative).open('rb')as f:
        nk,nrad,nx,ny,nz,iteration=map(int,take(f,'<i4',6))
        if not(1<nk<=1000 and 0<nrad<=100000 and nx==ny==1 and 0<nz<=10000 and iteration>0):raise ValueError('Unexpected1D rate header')
        lower=take(f,'<i4',nrad)-1;upper=take(f,'<i4',nrad)-1
        invalid=(lower<0)|(lower>=nk)|(upper<0)|(upper>=nk)|(lower==upper)
        sentinel=(upper==-1)&(lower>=0)&(lower<nk)
        if np.any(invalid&~sentinel):raise ValueError('Corrupt native transition identity')
        n=take(f,'<f8',nk*nz).reshape(nz,nk,order='F').T;s=take(f,'<f8',nk*nz).reshape(nz,nk,order='F').T
        up=take(f,'<f4',nrad*nz).reshape(nz,nrad,order='F').T.astype(float);down=take(f,'<f4',nrad*nz).reshape(nz,nrad,order='F').T.astype(float)
        if f.read():raise ValueError('Unexpected extra radiative bytes')
    with Path(collision).open('rb')as f:
        if list(take(f,'<i4',4))!=[nk,nx,ny,nz]:raise ValueError('Radiative/collision identity mismatch')
        cup=take(f,'<f4',nz*nk*(nk-1)//2).reshape(nz,nk*(nk-1)//2,order='F').T.astype(float)
        if f.read():raise ValueError('Unexpected extra collision bytes')
    if any(not np.isfinite(x).all()for x in [n,s,up,down,cup])or np.any(n<=0)or np.any(s<=0)or any(np.any(x<0)for x in [up,down,cup]):raise ValueError('Invalid reference populations/rates')
    # atom_mod creates one continuum slot per level. A top-ion level with
    # no next ion has j=0 and is explicitly deactivated; only zero rates qualify.
    if np.any(up[invalid]!=0)or np.any(down[invalid]!=0):raise ValueError('Active rate has an invalid/sentinel level identity')
    c=np.zeros((nk,nk,nz));r=np.zeros_like(c);k=0
    for j in range(1,nk):
        for i in range(j):c[j,i]=cup[k];c[i,j]=cup[k]*s[i]/s[j];k+=1
    for k,(i,j)in enumerate(zip(lower,upper)):
        if not invalid[k]:r[j,i]+=up[k];r[i,j]+=down[k]
    return dict(populations_cm3=n,nstar_cm3=s,collision_rates_s_inv=c,radiative_rates_s_inv=r,lower=lower,upper=upper,inactive_sentinel_slots=np.flatnonzero(invalid),iteration=iteration)


def read_patch_set(folder,kind):
    folder=Path(folder);meta=folder/(kind+'_meta.txt');result={}
    for line in meta.read_text().splitlines():
        if not line.strip()or line.lstrip().startswith('*'):continue
        name,bits,*shape=line.split();shape=tuple(map(int,shape));bits=int(bits)
        if bits not in [32,64]or not shape or any(v<=0 for v in shape)or Path(name).name!=name:raise ValueError('Invalid patch metadata')
        path=folder/name
        if path.stat().st_size!=np.prod(shape)*(bits//8):raise ValueError('Native patch byte-count mismatch')
        a=np.fromfile(path,dtype='<f'+str(bits//8)).reshape(shape,order='F')
        if not np.isfinite(a).all():raise ValueError('Nonfinite native patch')
        key=tuple(int(name.rsplit('_',1)[1].split('.')[0][i:i+3])for i in range(0,6,3))
        if key in result and not np.array_equal(result[key],a):raise ValueError('Contradictory duplicate patch identity')
        result[key]=a
    if not result:raise ValueError('Empty native patch set')
    return result


def interpolate_positive(source_height,values,target_height):
    """Explicit log-positive height interpolation, with only float32 endpoint rounding."""
    x=np.asarray(source_height,float);y=np.asarray(values,float);t=np.asarray(target_height,float)
    if x.ndim!=1 or y.shape[-1]!=len(x)or not np.all(np.diff(x)>0)or any(not np.isfinite(v).all()for v in [x,y,t])or np.any(y<=0):raise ValueError('Invalid positive interpolation layout')
    rounding=float(len(x)*np.finfo(np.float32).eps*max(1.,float(np.ptp(x))))
    if t.min()<x[0]-rounding or t.max()>x[-1]+rounding:raise ValueError('Physical extrapolation of independent reference prohibited')
    clipped=np.clip(t,x[0],x[-1]);out=np.stack([np.exp(np.interp(clipped,x,np.log(v)))for v in y.reshape(-1,len(x))]).reshape(y.shape[:-1]+t.shape)
    return out,dict(method='log-positive linear height interpolation; not an independent interpolation-accuracy bound',max_endpoint_rounding_clamp_m=float(np.max(abs(clipped-t))),endpoint_roundoff_bound_m=rounding)
def per_angstrom_to_per_hz(wavelength_vac_A, intensity):
    """Native erg/s/cm2/A/sr to erg/s/cm2/Hz/sr (vacuum Jacobian)."""
    import numpy as np
    wave=np.asarray(wavelength_vac_A,dtype=float);values=np.asarray(intensity)
    if wave.shape!=values.shape or np.any(wave<=0) or not np.isfinite(wave).all() or not np.isfinite(values).all():
        raise ValueError('Invalid wavelength/intensity conversion arrays')
    return values*wave**2/2.99792458e18


def supported_air_axis(wavelength_vac_A):
    """AIR coordinates only above200nm; keep vacuum primary for UV/IR output.

    NaN denotes the explicitly unsupported AIR convention, never missing intensity.
    """
    from pipeline.wavelength_util import vac_to_air
    wave=np.asarray(wavelength_vac_A,dtype=float)
    if wave.ndim!=1 or not np.isfinite(wave).all()or np.any(wave<=0)or np.any(np.diff(wave)<=0):raise ValueError('Invalid explicit vacuum wavelength axis')
    supported=wave>2000.;air=np.full(wave.shape,np.nan)
    if np.any(supported):air[supported]=vac_to_air(wave[supported])
    return air,supported
