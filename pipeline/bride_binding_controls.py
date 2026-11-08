"""Explicit immutable state controls; no EOS substitution or provider admission."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class BindingControl:
    case:str
    height_m:tuple
    temperature_K:tuple
    electron_pressure_Pa:tuple
    ne_m3:tuple
    nH_m3:tuple
    vturb_m_s:tuple
    ne_owner:str
    nH_owner:str


def control_state(structure,records,case,*,kBoltzmann_SI):
    choices={'model_mass_H':'MARCS_LW_nH_m3','native_network_H':'nH_network_m3'}
    if case not in choices:raise ValueError('Unknown explicit H-normalization control')
    n=structure['ndep'];names=records.dtype.names or ()
    required=['depth_index','temperature_K','electron_pressure_dyn_cm2','xi_cm_s',*choices.values()]
    if any(k not in names for k in required)or len(records)!=n or not np.array_equal(records['depth_index'],np.arange(1,n+1)):raise ValueError('Invalid binding depth identity')
    fields=[np.asarray(records[k],float)for k in required[1:]];height=np.asarray(structure['height_m'],float)
    if not np.isfinite(kBoltzmann_SI)or kBoltzmann_SI<=0 or height.shape!=(n,)or not np.isfinite(height).all()or np.any(np.diff(height)>=0):raise ValueError('Invalid unit/geometry binding')
    if any(not np.isfinite(x).all()for x in fields)or any(np.any(records[k]<=0)for k in ['temperature_K','electron_pressure_dyn_cm2',*choices.values()])or np.any(records['xi_cm_s']<0):raise ValueError('Invalid exported thermodynamic state')
    T=records['temperature_K'].astype(float).copy();Pe=records['electron_pressure_dyn_cm2'].astype(float)*.1;ne=Pe/(kBoltzmann_SI*T)
    return BindingControl(case,tuple(height),tuple(T),tuple(Pe),tuple(ne),tuple(records[choices[case]].astype(float)),tuple(records['xi_cm_s'].astype(float)*.01),
        'exported native Pe/(LW k * exported native T)',choices[case])
