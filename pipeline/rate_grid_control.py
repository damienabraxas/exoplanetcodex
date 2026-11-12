"""Explicit numerical midpoint refinement; no line endpoints or physics changes."""
import numpy as np
from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS


def midpoint_grid(wavelength_vac_nm,passes,*,interval_vac_nm=None):
    wave=np.asarray(wavelength_vac_nm,dtype=float)
    if isinstance(passes,bool)or not isinstance(passes,int)or passes not in [1,2]:raise ValueError('Require one or two predeclared midpoint passes')
    if wave.ndim!=1 or len(wave)<2 or not np.isfinite(wave).all()or np.any(wave<=0)or np.any(np.diff(wave)<=0):raise ValueError('Require strictly increasing finite vacuum nm grid')
    interval=None if interval_vac_nm is None else np.asarray(interval_vac_nm,dtype=float)
    if interval is not None and(interval.shape!=(2,)or not np.isfinite(interval).all()or interval[0]>=interval[1]or interval[0]<wave[0]or interval[1]>wave[-1]):raise ValueError('Invalid fixed refinement interval')
    for _ in range(passes):
        gaps=np.ones(len(wave)-1,dtype=bool)if interval is None else(wave[:-1]>=interval[0])&(wave[1:]<=interval[1])
        if len(wave)+np.count_nonzero(gaps)>MAX_PACKET_WAVELENGTHS:raise ValueError('Midpoint diagnostic exceeds original 20000-point packet scope')
        wave=np.sort(np.concatenate([wave,wave[:-1][gaps]+np.diff(wave)[gaps]/2]))
    if np.any(np.diff(wave)<=0):raise ValueError('Midpoint is not representable')
    return wave
