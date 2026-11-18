"""Read-only finite-atom Saha/Boltzmann canaries; never solver populations."""
import numpy as np


def normalized_weights(log_weights,total):
    log=np.asarray(log_weights,float);total=np.asarray(total,float)
    if log.ndim!=2 or total.shape!=(log.shape[1],)or not np.isfinite(log).all()or not np.isfinite(total).all()or np.any(total<=0):raise ValueError('Invalid LTE normalization state')
    weights=np.exp(log-np.max(log,axis=0))
    result=weights/weights.sum(axis=0)*total
    if not np.isfinite(result).all()or np.any(result<=0):raise ValueError('Nonpositive LTE population; no clipping')
    return result


def native_lte_weights(energy_eV,g,stage,T,ne_cm3,constants,*,round_native=True,shield=True):
    """Source atom_mod REAL4 e_kT/lsaha/dxi; two ion stages only.

    Energies and statistical weights are actual consumed native-table values.
    Returns unnormalized log weights and source composition auxiliaries. Absolute
    normalization/partition/EOS and total element density are not established.
    """
    E=np.asarray(energy_eV,float);g=np.asarray(g,float);stage=np.asarray(stage)
    T=np.asarray(T,float);ne=np.asarray(ne_cm3,float)
    if E.ndim!=1 or g.shape!=E.shape or stage.shape!=E.shape or T.ndim!=1 or ne.shape!=T.shape or not np.array_equal(np.unique(stage),[1,2])or np.any(np.diff(E)<0)or any(not np.isfinite(x).all()for x in [E,g,T,ne])or any(np.any(x<=0)for x in [g,T,ne]):raise ValueError('Invalid two-stage native LTE state')
    if not isinstance(round_native,bool)or not isinstance(shield,bool):raise ValueError('Explicit LTE precision/shield booleans required')
    if set(constants)!={'hh','bk','ee','em','ec','pi'}or any(not np.isfinite(v)or v<=0 for v in constants.values()):raise ValueError('Invalid source LTE constants')
    h,k,e,m,q,pi=[constants[x]for x in ['hh','bk','ee','em','ec','pi']]
    saha=(2*pi*m*k/h**2)**1.5
    lsaha=np.log(2*saha*T*np.sqrt(T))-np.log(ne)
    theta=np.log10(np.exp(1.))*e/k
    dxi_const=q**3*np.sqrt(8*pi)/(theta*k*e)
    dxi=dxi_const*(theta/T)*np.sqrt(T*k*ne)if shield else np.zeros_like(T)
    ek=e/k/T
    if round_native:ek,lsaha,dxi=[v.astype(np.float32).astype(float)for v in [ek,lsaha,dxi]]
    charge=stage[:,None]-1
    log=np.log(g[:,None])-ek*(E[:,None]-dxi*charge*stage[:,None]*.5)+charge*lsaha
    return log,dict(e_kT=ek,lsaha=lsaha,dxi_eV=dxi)


def held_lw_weights(energy_J,g,stage,T,ne_m3,constants,*,shield=True):
    """Independent log form of held-LW finite-level LTE, including Debye term.

    Stage numbering here is one-based, unlike the zero-based LW AtomicLevel.
    """
    E=np.asarray(energy_J,float);g=np.asarray(g,float);stage=np.asarray(stage);T=np.asarray(T,float);ne=np.asarray(ne_m3,float)
    if E.ndim!=1 or g.shape!=E.shape or stage.shape!=E.shape or T.ndim!=1 or ne.shape!=T.shape or not np.array_equal(np.unique(stage),[1,2])or any(not np.isfinite(x).all()for x in [E,g,T,ne])or any(np.any(x<=0)for x in [g,T,ne]):raise ValueError('Invalid held-LW two-stage state')
    if set(constants)!={'h','k','m','q','epsilon'}or any(not np.isfinite(v)or v<=0 for v in constants.values()):raise ValueError('Invalid held-LW constants')
    h,k,m,q,epsilon=[constants[x]for x in ['h','k','m','q','epsilon']]
    saha=np.log(2*(2*np.pi*m*k*T/h**2)**1.5/ne)
    lowering=np.sqrt(8*np.pi/k)*(q*q/(4*np.pi*epsilon))**1.5*np.sqrt(ne/T)if shield else np.zeros_like(T)
    charge=stage[:,None]-stage[0]
    return np.log(g[:,None]/g[0])-(E[:,None]-E[0]-charge*lowering)/(k*T)+charge*saha
