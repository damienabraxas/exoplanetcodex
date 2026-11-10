"""Small-memory inverse secant updates in number-conserving log coordinates."""
import numpy as np


def log_coordinates(populations):
    n=np.asarray(populations)
    if n.ndim!=2 or not np.isfinite(n).all()or np.any(n<=0):raise ValueError('Positive finite populations required')
    return np.log(n[1:]/n[0])


def normalized_log_populations(coordinates):
    x=np.asarray(coordinates)
    if x.ndim!=2 or not np.isfinite(x).all():raise ValueError('Finite log coordinates required')
    logs=np.vstack([np.zeros((1,x.shape[1])),x])
    return logs-np.logaddexp.reduce(logs,axis=0)


def populations_from_coordinates(coordinates,total):
    logs=normalized_log_populations(coordinates);n=np.exp(logs)*np.asarray(total)
    if not np.isfinite(n).all()or np.any(n<=0):raise ValueError('Invalid decoded populations; no clipping or floor')
    return n


class InverseSecant:
    """H starts at I; H_new=H+(s-Hy)(s^T H)/(s^T H y).

    Store low-rank outer products, restart at the declared memory limit or an
    unresolved secant denominator. No approximate Jacobian alters physical rates.
    """
    def __init__(self, *, memory=20, denominator_rtol=1e-12):
        if memory<1 or not 0<denominator_rtol<1:raise ValueError('Invalid secant configuration')
        self.memory=memory;self.denominator_rtol=denominator_rtol;self.u=[];self.v=[]

    def apply(self,vector,*,transpose=False):
        x=np.asarray(vector);out=x.copy()
        for u,v in zip(self.u,self.v):out+=(v*np.dot(u,x)if transpose else u*np.dot(v,x))
        return out

    def update(self,step,residual_change):
        s=np.asarray(step);y=np.asarray(residual_change)
        if s.shape!=y.shape or not np.isfinite(s).all()or not np.isfinite(y).all():raise ValueError('Invalid secant vectors')
        restart=len(self.u)>=self.memory
        if restart:self.u=[];self.v=[]
        hy=self.apply(y);v=self.apply(s,transpose=True);denominator=float(np.dot(v,y));scale=float(np.linalg.norm(v)*np.linalg.norm(y))
        if not np.isfinite(denominator)or scale==0 or abs(denominator)<=self.denominator_rtol*scale:
            self.u=[];self.v=[];return dict(updated=False,restarted=True,rank=0)
        u=(s-hy)/denominator
        if not np.isfinite(u).all():raise ValueError('Nonfinite inverse secant update')
        self.u.append(u);self.v.append(v)
        return dict(updated=True,restarted=restart,rank=len(self.u))
