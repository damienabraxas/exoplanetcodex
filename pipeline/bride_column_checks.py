"""Pure diagnostic column checks; does not admit interpolated atmosphere inputs."""
import numpy as np


def integrated_column(height_m,rho_kg_m3,top_column_kg_m2):
    z=np.asarray(height_m,dtype=float);rho=np.asarray(rho_kg_m3,dtype=float)
    if (z.ndim!=1 or len(z)<2 or rho.shape!=z.shape or not np.isfinite(z).all()
            or not np.isfinite(rho).all() or np.any(rho<=0) or np.any(np.diff(z)>=0)
            or not np.isfinite(top_column_kg_m2) or top_column_kg_m2<=0):
        raise ValueError('Column check requires ordered outward height, positive density/top column in SI')
    return top_column_kg_m2+np.r_[0.,np.cumsum(.5*(rho[:-1]+rho[1:])*(-np.diff(z)))]


def candidate_column(columns,weights,*,scale):
    """Convex native-row candidate only, using the existing simplex's weights."""
    values=np.asarray(columns,dtype=float);w=np.asarray(weights,dtype=float)
    if (values.ndim!=2 or w.shape!=(values.shape[0],) or not np.isfinite(values).all()
            or not np.isfinite(w).all() or np.any(w<0) or not np.isclose(w.sum(),1.,rtol=0,atol=1e-12)
            or np.any(values<=0) or np.any(np.diff(values,axis=1)<=0)):
        raise ValueError('Invalid native columns or nonconvex interpolation weights')
    if scale=='linear':return w@values
    if scale=='log':return np.exp(w@np.log(values))
    raise ValueError('Declare linear or log column interpolation')
