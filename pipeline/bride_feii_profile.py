"""Fixed-window shape diagnostics and explicit stage-II native Stark branch."""
import numpy as np


def feii_stark_gamma(temperature,ne_m3,row,mode):
    from pipeline.bride_fe_stark import stark_gamma
    if mode=='provided':return stark_gamma(temperature,ne_m3,row,'provided')
    if mode!='cowley':raise ValueError('Require explicit FeII Stark mode')
    T,ne=np.broadcast_arrays(temperature,ne_m3)
    if not np.isfinite(T).all() or not np.isfinite(ne).all() or np.any(T<=0) or np.any(ne<0):raise ValueError('Invalid Stark fields')
    gap=row['ionization_eV']-(row['excitation_eV']+3.40*3647./row['wavelength_air_A'])
    if not np.isfinite(gap) or gap<=0:raise ValueError('Invalid FeII upper gap')
    return 5.42e-7*13.598*(2./gap)**2*(ne*1e-6)/3.**2


def feature_metrics(wave,flux,centre,*,core_half_A,shoulder_inner_A,shoulder_outer_A):
    from scripts.audit_bride_kitt_residuals import integral
    w,f=np.asarray(wave,float),np.asarray(flux,float)
    if w.ndim!=1 or w.shape!=f.shape or not np.isfinite(w).all() or not np.isfinite(f).all() or np.any(np.diff(w)<=0):raise ValueError('Invalid feature arrays')
    if not 0<core_half_A<shoulder_inner_A<shoulder_outer_A:raise ValueError('Invalid fixed windows')
    if w[0]>centre-shoulder_outer_A or w[-1]<centre+shoulder_outer_A:raise ValueError('Truncated feature')
    width=shoulder_outer_A-shoulder_inner_A
    left=integral(w,f,centre-shoulder_outer_A,centre-shoulder_inner_A)/width
    right=integral(w,f,centre+shoulder_inner_A,centre+shoulder_outer_A)/width
    level=(left+right)/2
    if level<=0:raise ValueError('Nonpositive shoulder statistic')
    core=integral(w,f,centre-core_half_A,centre+core_half_A)/(2*core_half_A)
    return dict(shoulder_level=level,left_shoulder=left,right_shoulder=right,core_mean=core,
        contrast=1-core/level,relative_core_area_mA=1000*(2*core_half_A-integral(w,f/level,centre-core_half_A,centre+core_half_A)),
        scope='Fixed shoulder contrast, not an isolated EW, fitted/adopted continuum, formation diagnostic or uncertainty estimate')
