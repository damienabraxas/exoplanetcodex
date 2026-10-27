"""Pure input controls for the explicitly authorized Mac validation experiment."""
import copy
import numpy as np


def diagnostic_collision_scaling(atom, scales):
    """Copy and scale named collision classes for causal controls, never adoption."""
    allowed = {c.kind for c in atom.collisions}
    if set(scales)-allowed:
        raise ValueError('Unknown collision class in diagnostic scaling')
    if any(not np.isfinite(v) or v <= 0 for v in scales.values()):
        raise ValueError('Collision scales must be positive and finite')
    result = copy.deepcopy(atom)
    for collision in result.collisions:
        collision.rates *= scales.get(collision.kind, 1.)
    return result


def marcs_filename(teff,logg,feh=0.,alpha=0.):
    return f'p{teff}_g+{logg:.1f}_m0.0_t01_st_z{feh:+.2f}_a{alpha:+.2f}_c+0.00_n+0.00_o{alpha:+.2f}_r+0.00_s+0.00.mod'


def interpolate_marcs_structure(nodes,target,*,thermodynamic_scale):
    """Experimental linear-simplex interpolation on native common tauRoss layers.

    Uses the SciPy primitive used by iSpec, retaining density as well as T/Pe/
    geometry required by Bride. No EOS, atmosphere-equilibrium solve, temperature
    offset, nearest-node fallback or extrapolation is invented. Sampling density
    is NOT physical interpolation accuracy; test held-out native models.
    """
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import Delaunay
    points=np.asarray([p for p,s in nodes],dtype=float)
    reference=nodes[0][1]
    names=('height_m','temperature_K','electron_pressure_Pa','rho_kg_m3','logtau_ross','logtau_500')
    for p,s in nodes:
        if s['ndep']!=reference['ndep'] or not np.array_equal(s['logtau_ross'],reference['logtau_ross']):
            raise ValueError('Atmospheres do not share the native tauRoss layer grid')
        if s['abundances_dex']!=reference['abundances_dex']or s['hydrogen_mass_fraction']!=reference['hydrogen_mass_fraction']:
            raise ValueError('Cannot interpolate across unmodeled composition changes')
    target=np.asarray(target,dtype=float);tri=Delaunay(points)
    simplex=int(tri.find_simplex(target))
    if simplex<0:raise ValueError('Target outside atmosphere interpolation hull; no nearest fallback')
    if thermodynamic_scale not in ('linear','log_positive','log_positive_layers'):raise ValueError('Unknown thermodynamic interpolation convention')
    values=np.stack([np.column_stack([s[k]for k in names])for p,s in nodes])
    if thermodynamic_scale.startswith('log_positive'):
        if np.any(values[:,:,[2,3]]<=0):raise ValueError('Nonpositive pressure/density input')
        values[:,:,[2,3]]=np.log(values[:,:,[2,3]])
    actual=LinearNDInterpolator(tri,values)([target])[0]
    if thermodynamic_scale.startswith('log_positive'):actual[:,[2,3]]=np.exp(actual[:,[2,3]])
    if thermodynamic_scale=='log_positive_layers':
        thickness=np.stack([-np.diff(s['height_m'])for p,s in nodes])
        if np.any(thickness<=0):raise ValueError('Nonpositive layer thickness')
        dz=np.exp(LinearNDInterpolator(tri,np.log(thickness))([target])[0])
        anchor=int(np.argmin(abs(reference['logtau_ross'])))
        height=np.empty(reference['ndep']);height[anchor]=actual[anchor,0]
        for i in range(anchor,len(height)-1):height[i+1]=height[i]-dz[i]
        for i in range(anchor-1,-1,-1):height[i]=height[i+1]+dz[i]
        actual[:,0]=height
    if not np.isfinite(actual).all():raise ValueError('Nonfinite interpolated structure')
    result=copy.deepcopy(reference)
    for i,name in enumerate(names):result[name]=actual[:,i].copy()
    if any(np.any(result[k]<=0)for k in ('temperature_K','electron_pressure_Pa','rho_kg_m3')):
        raise ValueError('Unphysical interpolated field')
    if not np.all(np.diff(result['height_m'])<0):raise ValueError('Invalid interpolated geometry')
    result['logg_cgs']=float(target[1]);result['model_id']=f'EXPERIMENTAL_INTERPOLATED_MARCS_T{target[0]}_g{target[1]}'
    partial=tri.transform[simplex,:2].dot(target-tri.transform[simplex,2])
    weights=np.r_[partial,1-partial.sum()]
    return result,dict(method='SciPy LinearNDInterpolator, common native logtauRoss rows; density retained; no EOS re-solve',
                       thermodynamic_scale=thermodynamic_scale,
                       target=target.tolist(),simplex_source_indices=tri.simplices[simplex].tolist(),weights=weights.tolist())


def controlled_structure(original,*,xi,drop_bottom):
    if not isinstance(drop_bottom,int) or drop_bottom<0 or drop_bottom>original['ndep']-4:
        raise ValueError('Invalid bottom-depth selection')
    if not np.isfinite(xi) or xi<0:raise ValueError('Invalid microturbulence')
    result=copy.deepcopy(original);n=original['ndep']-drop_bottom
    for name in ('height_m','temperature_K','electron_pressure_Pa','rho_kg_m3','logtau_ross','logtau_500'):
        result[name]=result[name][:n].copy()
    result['ndep']=n;result['vturb_m_s']=xi*1000
    return result


def assert_collision_coverage(atom,temperature):
    """No implicit WENO extrapolation/clamping beyond any supplied rate table."""
    lo=float(np.min(temperature));hi=float(np.max(temperature))
    bad=[]
    for c in atom.collisions:
        grid=atom.temperature if c.temperature is None else c.temperature
        if lo<grid[0] or hi>grid[-1]:bad.append((c.kind,c.i,c.j,float(grid[0]),float(grid[-1])))
    if bad:raise ValueError(f'Temperature range {lo}..{hi} exceeds {len(bad)} collision tables; first={bad[0]}')
