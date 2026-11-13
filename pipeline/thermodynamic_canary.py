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


def collision_balance(matrix, populations):
    """Audit C[destination,source,depth] without altering coefficients or states.

    Return signed pair defects only where a rate exists. Disconnected pairs are
    absent, not successful zero-defect measurements. One-way edges are reported
    at +/-1 rather than hidden by a logarithm or positive-only mask.
    """
    C=np.asarray(matrix,dtype=float);n=np.asarray(populations,dtype=float)
    if n.ndim!=2 or C.shape!=(n.shape[0],n.shape[0],n.shape[1]):
        raise ValueError('Collision matrix/state dimensions differ')
    if not np.isfinite(C).all()or np.any(C<0)or not np.isfinite(n).all()or np.any(n<=0):
        raise ValueError('Invalid collision matrix/state')
    if np.any(C[np.arange(len(n)),np.arange(len(n))]!=0):
        raise ValueError('Expected off-diagonal collision coefficients')
    flow=C*n[None,:,:];reverse=flow.transpose(1,0,2);gross=flow+reverse
    present=gross>0
    defect=np.full_like(gross,np.nan)
    np.divide(flow-reverse,gross,out=defect,where=present)
    return defect,present


def integrated_defect(up,down,weights):
    """Integrate positive directional terms on [wavelength,depth] grids.

    Common transition/depth factors cancel in the normalized defect. Weights
    must already carry the intended frequency measure and profile/cross section.
    """
    up=np.asarray(up,dtype=float);down=np.asarray(down,dtype=float);w=np.asarray(weights,dtype=float)
    if up.ndim!=2 or down.shape!=up.shape or w.shape!=up.shape:
        raise ValueError('Integrated rate shapes differ')
    if any(not np.isfinite(x).all()for x in [up,down,w])or np.any(up<=0)or np.any(down<=0)or np.any(w<0)or np.any(w.sum(axis=0)<=0):
        raise ValueError('Invalid integrated rate terms/measure')
    return balanced_defect(np.sum(w*up,axis=0),np.sum(w*down,axis=0))


def native_split_planck_defect(beta,thermal_integral):
    """Native line_rate_coefficients identity at J=B (dev=0).

    The common g_i/g_j factor cancels against the LTE population ratio.
    This checks a source-derived formula; it is not a native solver run.
    """
    beta=np.asarray(beta,dtype=float);th=np.asarray(thermal_integral,dtype=float)
    if beta.shape!=th.shape or np.any(beta<=0)or np.any(beta>=1):
        raise ValueError('Invalid native thermal split ratio')
    nref=beta/(1-beta)
    return balanced_defect(nref*th,beta*((nref+1)*th))
