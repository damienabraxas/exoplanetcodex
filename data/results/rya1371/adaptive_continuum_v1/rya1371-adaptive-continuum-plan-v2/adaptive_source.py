"""Bounded nested trapezoid error indicators; never remove original nodes."""
import numpy as np
from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS


def nested_indices(base,fine):
    base=np.asarray(base,float);fine=np.asarray(fine,float)
    if any(x.ndim!=1 or len(x)<2 or not np.isfinite(x).all()or np.any(x<=0)or np.any(np.diff(x)<=0)for x in [base,fine]):raise ValueError('Invalid nested wavelength grids')
    idx=np.searchsorted(fine,base)
    if np.any(idx>=len(fine))or not np.array_equal(fine[idx],base)or idx[0]!=0 or idx[-1]!=len(fine)-1:raise ValueError('Inherited nodes/endpoints missing')
    cost=np.diff(idx)
    if np.any(cost<2)or np.any(cost&(cost-1)):raise ValueError('Require at least one nested dyadic refinement per base panel')
    return idx,cost


def panel_indicators(base,fine,integrand):
    """Sum absolute finest/coarsened panel differences, normalized per depth.

    An error indicator for point allocation, not a rigorous quadrature bound.
    """
    idx,_=nested_indices(base,fine);F=np.asarray(integrand,float);fine=np.asarray(fine,float)
    if F.ndim!=2 or F.shape[0]!=len(fine)or not np.isfinite(F).all()or np.any(F<0):raise ValueError('Invalid positive continuum integrand')
    total=np.sum((F[1:]+F[:-1])*.5*np.diff(fine)[:,None],axis=0)
    if np.any(total<=0):raise ValueError('Empty continuum rate support')
    scores=[]
    for start,end in zip(idx[:-1],idx[1:]):
        x=fine[start:end+1];f=F[start:end+1];width=np.diff(x)
        small=(f[1:]+f[:-1])*.5*width[:,None]
        fine_pairs=small[::2]+small[1::2]
        coarse_pairs=(f[::2][1:]+f[::2][:-1])*.5*(x[2::2]-x[:-2:2])[:,None]
        scores.append(float(np.max(np.sum(abs(fine_pairs-coarse_pairs),axis=0)/total)))
    return np.asarray(scores)


def adaptive_grid(base,fine,scores,additional_points):
    idx,cost=nested_indices(base,fine);fine=np.asarray(fine,float);scores=np.asarray(scores,float)
    if scores.shape!=cost.shape or not np.isfinite(scores).all()or np.any(scores<0):raise ValueError('Invalid continuum panel scores')
    if isinstance(additional_points,bool)or not isinstance(additional_points,int)or additional_points<=0 or len(fine)+additional_points>MAX_PACKET_WAVELENGTHS:raise ValueError('Adaptive budget exceeds original 20000-point scope')
    remaining=additional_points;chosen=[];nodes=[]
    # Largest estimated relative error per added node first. Stable ties preserve
    # source panel order, independent of desired reference agreement.
    for panel in np.argsort(-scores/cost,kind='stable'):
        if scores[panel]<=0 or cost[panel]>remaining:continue
        local=fine[idx[panel]:idx[panel+1]+1];nodes.extend(local[:-1]+np.diff(local)*.5);chosen.append(int(panel));remaining-=int(cost[panel])
    if not nodes:raise ValueError('No positive, affordable adaptive panel')
    out=np.sort(np.concatenate([fine,np.asarray(nodes)]))
    if np.any(np.diff(out)<=0):raise ValueError('Adaptive midpoint not representable')
    np.testing.assert_array_equal(out[np.searchsorted(out,fine)],fine)
    return out,np.asarray(chosen,dtype=int)
