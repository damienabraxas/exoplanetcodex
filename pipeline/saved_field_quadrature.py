"""Bounded quadrature of declared saved-field interpolants, not new radiation."""
import numpy as np
from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS


def gaussian_panel_batches(panels,order,*,batch_panels=256):
    x=np.asarray(panels,float)
    if x.ndim!=1 or len(x)<2 or not np.isfinite(x).all()or np.any(x<=0)or np.any(np.diff(x)<=0):raise ValueError('Invalid Gaussian panel endpoints')
    if isinstance(order,bool)or order not in [4,8,16]or isinstance(batch_panels,bool)or not isinstance(batch_panels,int)or batch_panels<=0 or batch_panels*order>MAX_PACKET_WAVELENGTHS:raise ValueError('Gaussian query exceeds original order/packet contract')
    nodes,weights=np.polynomial.legendre.leggauss(order)
    for first in range(0,len(x)-1,batch_panels):
        lower=x[first:first+batch_panels];upper=x[first+1:first+batch_panels+1];lower=lower[:len(upper)];span=upper-lower
        query=(lower[:,None]+span[:,None]*(nodes[None,:]+1)*.5).ravel()
        measure=(span[:,None]*weights[None,:]*.5).ravel()
        if np.any(query.reshape(-1,order)<=lower[:,None])or np.any(query.reshape(-1,order)>=upper[:,None])or np.any(measure<=0):raise ValueError('Unrepresentable positive interior Gaussian panel')
        yield query,measure


def interpolate_saved_field(wave,J,query):
    """Physical J_nu linearly interpolated in vacuum wavelength, no extrapolation."""
    wave=np.asarray(wave,float);J=np.asarray(J,float);q=np.asarray(query,float)
    if wave.ndim!=1 or len(wave)<2 or J.ndim!=2 or J.shape[0]!=len(wave)or q.ndim!=1 or not len(q)or len(q)>MAX_PACKET_WAVELENGTHS or any(not np.isfinite(v).all()for v in [wave,J,q])or np.any(wave<=0)or np.any(np.diff(wave)<=0)or np.any(J<0)or np.any(q<wave[0])or np.any(q>wave[-1]):raise ValueError('Invalid saved-field interpolation or extrapolation')
    return np.stack([np.interp(q,wave,J[:,depth])for depth in range(J.shape[1])],axis=1)
