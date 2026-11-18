"""Catalog wavelength-envelope accounting, never instrument/science admission."""
import numpy as np


def catalog_interval_vac_nm(lo_nm,hi_nm,frame):
    bounds=np.asarray([lo_nm,hi_nm],dtype=float)
    if not np.isfinite(bounds).all()or not 0<bounds[0]<bounds[1]:raise ValueError('Invalid catalog wavelength-nm bounds')
    if frame=='air':
        if bounds[0]<=200:raise ValueError('Catalog AIR interval crosses unsupported UV conversion domain')
        from pipeline.wavelength_util import air_to_vac
        return air_to_vac(bounds*10.)/10.
    if frame in ['vacuum','vacuum_wavenumber']:
        # Catalog columns are explicitly wavelength_nm, including IAG. The latter
        # describes raw PRODUCT coordinates, not the units of these catalog bounds.
        return bounds.copy()
    raise ValueError('Unresolved catalog wavelength frame')


def interval_intersections(bounds,windows):
    lo,hi=map(float,bounds)
    if not np.isfinite([lo,hi]).all()or not lo<hi:raise ValueError('Invalid target interval')
    result=[]
    for a,b in windows:
        if not np.isfinite([a,b]).all()or not a<b:raise ValueError('Invalid declared diagnostic interval')
        if max(lo,a)<min(hi,b):result.append([max(lo,float(a)),min(hi,float(b))])
    return result
