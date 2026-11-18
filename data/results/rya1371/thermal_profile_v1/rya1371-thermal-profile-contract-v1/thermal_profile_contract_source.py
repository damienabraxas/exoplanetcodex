"""Read-only component equilibrium checks; absent channels are not passes."""
import numpy as np


def background_balance(chi,eta,scattering,blackbody):
    arrays=[np.asarray(x,dtype=float)for x in [chi,eta,scattering,blackbody]]
    chi,eta,sca,B=arrays
    if chi.ndim!=2 or any(x.shape!=chi.shape or not np.isfinite(x).all()for x in arrays):raise ValueError('Invalid background grid')
    if any(np.any(x<0)for x in [chi,eta,sca])or np.any(B<=0)or np.any(sca>chi):raise ValueError('Invalid background coefficients')
    absorption=chi-sca;expected=absorption*B
    # eta is thermal emission only: coherent scattering enters exactly once.
    emission=eta+sca*B;loss=chi*B
    def summarize(up,down):
        gross=up+down;mask=gross>0
        if not np.any(mask):return dict(present_cells=0,max_absolute_signed_defect=None,worst_wave_depth=None)
        defect=np.zeros_like(gross);np.divide(up-down,gross,out=defect,where=mask)
        index=np.unravel_index(np.argmax(abs(defect)),defect.shape)
        return dict(present_cells=int(mask.sum()),max_absolute_signed_defect=float(abs(defect[index])),worst_wave_depth=list(map(int,index)),signed_defect_at_worst=float(defect[index]))
    return dict(thermal_absorption=summarize(eta,expected),total_at_J_equals_B=summarize(emission,loss),
        pure_scattering_cells=int(np.count_nonzero((absorption==0)&(sca>0))),zero_extinction_cells=int(np.count_nonzero(chi==0)),
        scattering_fraction_max=float(np.max(np.divide(sca,chi,out=np.zeros_like(chi),where=chi>0))))


def collision_summary(C,nstar):
    from pipeline.thermodynamic_canary import collision_balance
    defect,present=collision_balance(C,nstar)
    # Report each unordered pair once; do not double the inventory.
    select=np.triu(np.ones(present.shape[:2],dtype=bool),1)[:,:,None]&present
    if not np.any(select):return dict(pair_depths=0,max_absolute_signed_defect=None)
    return dict(pair_depths=int(select.sum()),one_way_pair_depths=int(np.count_nonzero(abs(defect[select])==1)),max_absolute_signed_defect=float(np.max(abs(defect[select]))))
