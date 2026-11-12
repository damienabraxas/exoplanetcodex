"""Analytic line-centre/per-frequency LTE identities; no RT/SE or adoption."""
import numpy as np


def occupation(x):
    x=np.asarray(x,dtype=float)
    if not np.isfinite(x).all()or np.any(x<=0)or np.any(x>700):raise ValueError('Invalid bounded thermal exponent')
    return 1/np.expm1(x)


def balanced_defect(up,down):
    up=np.asarray(up,dtype=float);down=np.asarray(down,dtype=float)
    if up.shape!=down.shape or not np.isfinite(up).all()or not np.isfinite(down).all()or np.any(up<=0)or np.any(down<=0):raise ValueError('Invalid thermal balance pair')
    return (up-down)/(up+down)


def bb_thermal_defect(beta_star,photon_occupation):
    """beta_star=(n*_j*g_i)/(n*_i*g_j); A and population factors cancel."""
    beta=np.asarray(beta_star,dtype=float);n=np.asarray(photon_occupation,dtype=float)
    if beta.shape!=n.shape or np.any(beta<=0)or np.any(beta>=1):raise ValueError('Invalid BB LTE ratio')
    return balanced_defect(n,beta*(1+n))


def bf_thermal_defect(intensity,einstein_prefactor,beta):
    """At n=n*, alpha/n* factors cancel; frequency-level identity only."""
    J=np.asarray(intensity,dtype=float);bare=np.asarray(einstein_prefactor,dtype=float);b=np.asarray(beta,dtype=float)
    if J.shape!=bare.shape or b.shape!=J.shape or not np.isfinite(b).all()or np.any(b<=0)or np.any(b>=1):raise ValueError('Invalid BF thermal fields')
    return balanced_defect(J,b*(bare+J))
