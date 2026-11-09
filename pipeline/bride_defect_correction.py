"""Number-conserving defect correction; approximate operators do not set roots."""
import numpy as np

def corrected_populations(operator, physical, populations, total):
    """Solve P delta=-R n, sum(delta)=N-sum(n), with row equilibration.

    All matrices use [destination,source,depth]. The largest-population row is
    replaced by conservation, matching the inspected LW statistical-equilibrium
    kernel. No inverse, floor, clipping, damping or altered physical rates.
    """
    p=np.asarray(operator);r=np.asarray(physical);n=np.asarray(populations)
    if p.shape!=r.shape or p.shape!=(n.shape[0],n.shape[0],n.shape[1]):
        raise ValueError('Inconsistent operator/population shapes')
    if not all(np.isfinite(x).all()for x in [p,r,n,total]) or np.any(n<=0):
        raise ValueError('Invalid defect-correction inputs')
    matrix=np.moveaxis(p.copy(),2,0)
    rhs=-np.einsum('ijd,jd->di',r,n)
    scale=np.max(abs(matrix),axis=2)
    if np.any(scale==0):raise ValueError('Disconnected preconditioner equation')
    matrix/=scale[:,:,None];rhs/=scale
    depth=np.arange(n.shape[1]);row=np.argmax(n,axis=0)
    matrix[depth,row,:]=1.;rhs[depth,row]=np.asarray(total)-n.sum(axis=0)
    delta=np.linalg.solve(matrix,rhs[:,:,None])[:,:,0].T
    return n+delta
