"""Strict reader for the isolated, stage-labelled native physical-rate sidecars.

Rate sums retain native float32 integration. BB stores thermal/deviation
integrals; BF downward sums still require the source LTE population ratio.
These files are NOT the original preconditioned rates_O exports.
"""
from pathlib import Path
import numpy as np
from pipeline.openmulti3d_output import take


def _array(f,dtype,nz,n):
    value=take(f,dtype,nz*n).reshape(nz,n,order='F').T
    if not np.isfinite(value).all():raise ValueError('Nonfinite observer array')
    return value


def read_state(path):
    with Path(path).open('rb')as f:
        magic,version,iteration,zs,nx,ny,nz,nk,nr,nline=map(int,take(f,'<i4',10))
        if (magic,version,nx,ny)!=(1371,1,1,1)or not(0<iteration<=1500 and 0<=zs<10000 and 0<nz<=10000 and 1<nk<=1000 and 0<nline<=nr<=100000):raise ValueError('Invalid observer state header')
        lower=take(f,'<i4',nr)-1;upper=take(f,'<i4',nr)-1
        n=_array(f,'<f8',nz,nk);star=_array(f,'<f8',nz,nk);norm=_array(f,'<f4',nz,nr)
        if f.read():raise ValueError('Trailing observer state bytes')
    if np.any(n<=0)or np.any(star<=0)or np.any(norm<=0):raise ValueError('Invalid observer physical state')
    return dict(iteration=iteration,zs=zs,nz=nz,nk=nk,nr=nr,nline=nline,lower=lower,upper=upper,populations_cm3=n,nstar_cm3=star,norm=norm)


def read_integrals(path,state,*,radiation_iteration):
    up=np.zeros((state['nr'],state['nz']),dtype=np.float64);down=up.copy();count=0;identity=[]
    with Path(path).open('rb')as f:
        while f.peek(1):
            magic,version,kind,it,zs,nx,ny,nz,nr,nk,fs,direction=map(int,take(f,'<i4',12))
            if (magic,version,nx,ny)!=(1371,1,1,1)or kind not in [1,2]or direction not in [-1,1]or not(0<nr<=state['nr']and fs>0):raise ValueError('Invalid observer integral header')
            if (it,zs,nz,nk)!=(radiation_iteration,state['zs'],state['nz'],state['nk']):raise ValueError('Observer stage mismatch')
            indices=take(f,'<i4',nr)-1
            if len(set(indices))!=nr or np.any(indices<0)or np.any(indices>=state['nr']):raise ValueError('Invalid observer transition map')
            a=_array(f,'<f4',nz,nr);b=_array(f,'<f4',nz,nr)
            n=_array(f,'<f8',nz,nk);star=_array(f,'<f8',nz,nk)
            if not np.array_equal(n,state['populations_cm3'])or not np.array_equal(star,state['nstar_cm3']):raise ValueError('Radiation sweep used a different population stage')
            # Keep double accumulation to avoid making cancellation in the diagnostic
            # depend on the observer call order; source inputs/integrals remain fp32.
            up[indices]+=a;down[indices]+=b;count+=1;identity.append((kind,fs,direction))
    if not count:raise ValueError('Empty observer integrals')
    return up,down,dict(records=count,kind_window_direction=identity,accumulation='float64 sum of source float32 contributions; not source summation-order replay')


def physical_matrix(state,integrals,weights):
    th,dev,audit=integrals;nk=state['nk'];nz=state['nz'];nr=state['nr'];w=np.asarray(weights,float)
    if w.shape!=(nk,)or not np.isfinite(w).all()or np.any(w<=0):raise ValueError('Invalid statistical weights')
    R=np.zeros((nk,nk,nz))
    for k,(i,j)in enumerate(zip(state['lower'],state['upper'])):
        if i<0 or i>=nk or j<0 or j>=nk or i==j:
            if np.any(th[k]!=0)or np.any(dev[k]!=0):raise ValueError('Nonzero inactive observer rate')
            continue
        if k<state['nline']:
            a=th[k]*state['norm'][k];b=dev[k]*state['norm'][k];g=w[i]/w[j];beta=state['nstar_cm3'][j]*g/state['nstar_cm3'][i]
            if np.any(beta<=0)or np.any(beta>=1):raise ValueError('Invalid BB LTE ratio')
            u=a*beta/(1-beta)+b;d=g*(a/(1-beta)+b)
        else:u=th[k];d=dev[k]*state['nstar_cm3'][i]/state['nstar_cm3'][j]
        if np.any(u<0)or np.any(d<0)or not np.isfinite(u+d).all():raise ValueError('Invalid gross physical rate')
        R[j,i]+=u;R[i,j]+=d
    return R
