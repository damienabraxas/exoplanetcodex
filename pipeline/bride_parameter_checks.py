"""Pure input controls for the explicitly authorized Mac validation experiment."""
import copy
import numpy as np


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
