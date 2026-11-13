"""Independent stationary, isotropic-profile rate integrals in SI units.

J is the saved angle-weighted intensity (not renormalized). Native ray weights
need not sum to exactly one; preserve that sum in spontaneous terms.
"""
import numpy as np


def trapezoid_measure(wavelength_nm):
    wave=np.asarray(wavelength_nm,dtype=float)
    if wave.ndim!=1 or len(wave)<2 or not np.isfinite(wave).all()or np.any(wave<=0)or np.any(np.diff(wave)<=0):
        raise ValueError('Invalid increasing vacuum wavelength grid')
    out=np.empty_like(wave);out[0]=(wave[1]-wave[0])*.5;out[-1]=(wave[-1]-wave[-2])*.5;out[1:-1]=(wave[2:]-wave[:-2])*.5
    return out


def bb_rates(wave,profile,J,*,lambda0_nm,Aji,Bij,Bji,angular_sum):
    measure=trapezoid_measure(wave)[:,None];wave=np.asarray(wave,dtype=float)
    phi=np.asarray(profile,dtype=float);J=np.asarray(J,dtype=float)
    if phi.ndim!=2 or J.shape!=phi.shape or len(wave)!=len(phi):raise ValueError('Invalid BB field/profile shapes')
    if any(not np.isfinite(x).all()or np.any(x<0)for x in [phi,J])or any(not np.isfinite(x)or x<=0 for x in [lambda0_nm,Aji,Bij,Bji,angular_sum]):raise ValueError('Invalid BB integral inputs')
    normalization=np.sum(measure*phi,axis=0)*angular_sum
    if np.any(normalization<=0):raise ValueError('Empty BB profile integral')
    weight=measure*phi/normalization[None,:]*(lambda0_nm/wave[:,None])
    field=np.sum(weight*J,axis=0)
    return Bij*field,Bji*field+Aji*np.sum(weight,axis=0)*angular_sum


def bf_rates(wave,alpha,J,*,nstar_ratio,T_K,h,k,c,angular_sum):
    measure=trapezoid_measure(wave)[:,None];wave=np.asarray(wave,dtype=float)
    alpha=np.asarray(alpha,dtype=float);J=np.asarray(J,dtype=float);T=np.asarray(T_K,dtype=float);ratio=np.asarray(nstar_ratio,dtype=float)
    if J.shape!=(len(wave),len(T))or alpha.shape!=wave.shape or ratio.shape!=T.shape:raise ValueError('Invalid BF field/state shapes')
    if any(not np.isfinite(x).all()or np.any(x<0)for x in [J,alpha])or any(not np.isfinite(x).all()or np.any(x<=0)for x in [T,ratio])or any(not np.isfinite(x)or x<=0 for x in [h,k,c,angular_sum]):raise ValueError('Invalid BF integral inputs')
    nu=c/(wave[:,None]*1e-9);bare=2*h*nu**3/c**2;beta=np.exp(-h*nu/(k*T[None,:]))
    weight=4*np.pi/h*measure/wave[:,None]*alpha[:,None]
    return np.sum(weight*J,axis=0),ratio*np.sum(weight*beta*(bare*angular_sum+J),axis=0)
