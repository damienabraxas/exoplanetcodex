"""Explicit experimental native-O collision recipe; no implicit H levels or rate adoption."""
import numpy as np

def collision_matrix(records,coefficients_cm3,temperature_K,ne_cm3,H_ground_raw_cm3,H_I_total_cm3,nstar_cm3,energy_eV,g,e_k,*,control='native_contract'):
    """Return C[destination,source,depth] in s^-1 for the declared fixed state.

    Coefficients must come from the source-qualified native log-Hermite oracle.
    Raw Hydrogen level density deliberately excludes occupation weights, matching
    native atom_mod. Alternative controls isolate semantics; none is a default.
    Native CH0 inverse uses the O LTE ratio, not an independently sourced MN rate.
    """
    if control not in ['native_contract','total_H_substitution','directed_CH0_only','LW_CH_direction']:raise ValueError('Unknown explicit collision control')
    T=np.asarray(temperature_K);ne=np.asarray(ne_cm3);H=np.asarray(H_ground_raw_cm3);HI=np.asarray(H_I_total_cm3);star=np.asarray(nstar_cm3);q=np.asarray(coefficients_cm3);energy=np.asarray(energy_eV);weight=np.asarray(g)
    if T.ndim!=1 or star.shape!=(41,len(T))or q.shape!=(len(records),len(T))or any(x.shape!=T.shape for x in [ne,H,HI])or energy.shape!=(41,)or weight.shape!=(41,):raise ValueError('Collision state dimensions differ')
    if any(not np.isfinite(x).all()or np.any(x<=0)for x in [T,ne,H,HI,star,q,weight])or not np.isfinite(energy).all()or np.any(np.diff(energy)<0)or not np.isfinite(e_k)or e_k<=0:raise ValueError('Invalid finite positive collision state')
    matrix=np.zeros((41,41,len(T)))
    for index,row in enumerate(records):
        idx,kind,lo,up=row[:4];lo=int(lo);up=int(up)
        if int(idx)!=index or not 0<=lo<up<41 or kind not in ['CE','CI','CH','CH0']:raise ValueError('Unsupported native collision record')
        ratio=star[up]/star[lo];coefficient=q[index]
        if kind=='CE':down=coefficient*ne*np.sqrt(T)*weight[lo]/weight[up];upward=down*ratio
        elif kind=='CI':upward=coefficient*ne*np.sqrt(T)*np.exp(-(energy[up]-energy[lo])*e_k/T);down=upward/ratio
        elif kind=='CH':
            down=coefficient*(HI if control=='total_H_substitution'else H);upward=down*ratio
            if control=='LW_CH_direction':upward=coefficient*H;down=upward/ratio
        else:
            upward=coefficient*(HI if control=='total_H_substitution'else H);down=upward/ratio
            if control=='directed_CH0_only':down=np.zeros_like(upward)
        matrix[up,lo]+=upward;matrix[lo,up]+=down
    if not np.isfinite(matrix).all()or np.any(matrix<0):raise ValueError('Nonfinite/negative collision matrix')
    return matrix
