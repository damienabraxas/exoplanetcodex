"""One-variable test of the turnover. If raising A(Fe) suppresses its own line depths
because Fe I BOUND-FREE opacity scales with A(Fe), then removing that component must
remove the turnover. Same ladder, same window, only the opacity table changes."""
import os, sys, numpy as np
sys.path.insert(0,'/Users/ryanschmitt/codex/ispec'); sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195')
exec(open('synth_pair.py').read().split("from pipeline.abundances_derive")[0])
from pipeline.abundances_derive import _ISPEC_SOLAR_ABUND_FILE, _SYNTH_ISOTOPE_FILE
TEFF,LOGG,MH,VT=5772.0,4.44,0.0,1.0; LO,HI=335.0,337.0
pack=ispec.load_modeled_layers_pack('/Users/ryanschmitt/codex/ispec/input/atmospheres/MARCS.GES/')
atm=ispec.interpolate_atmosphere_layers(pack,{'teff':TEFF,'logg':LOGG,'MH':MH,'alpha':0.0},code='turbospectrum')
ab=ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE); iso=ispec.read_isotope_data(_SYNTH_ISOTOPE_FILE)
ll=ispec.read_atomic_linelist(os.environ['NEARUV_LL']); ll=ll[(ll['wave_nm']>=LO-1)&(ll['wave_nm']<=HI+1)]
wave=np.arange(LO,HI,0.0002)
def run(table,dFe=0.0):
    TABLE['name']=table; fixed=[]
    if dFe:
        fixed=ab[ab['code']==26].copy(); fixed['Abund']=fixed['Abund']+dFe
    return np.asarray(ispec.generate_spectrum(wave,atm,TEFF,LOGG,MH,0.0,ll,iso,ab,fixed,
        microturbulence_vel=VT,macroturbulence=0.0,vsini=0.0,R=0,verbose=0,
        code='turbospectrum',tmp_dir=SP))
S,Z=f'{SP}/jonabs_STOCK.dat',f'{SP}/jonabs_FeI_ZEROED.dat'
b_s,b_z=run(S),run(Z); m=(b_s>0)&(b_z>0); d=lambda x:1.0-x[m].mean()
print(f"{'dA(Fe)':>8}  {'STOCK (Fe I bf on)':>20}  {'ZEROED (Fe I bf off)':>22}")
for x in (0.0,0.2,0.4,0.6,0.8,1.0):
    ds=d(run(S,x))-d(b_s) if x else 0.0
    dz=d(run(Z,x))-d(b_z) if x else 0.0
    print(f"{x:+8.2f}  {ds:+20.5f}  {dz:+22.5f}")
