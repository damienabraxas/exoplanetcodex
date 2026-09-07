"""Bracket, never extrapolate. The abundance response of a blanketed near-UV window
saturates, so a linear fit over 0-0.2 dex cannot be run out to the Fe I bf effect."""
import os, sys, numpy as np
sys.path.insert(0,'/Users/ryanschmitt/codex/ispec'); sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195')
exec(open('synth_pair.py').read().split("from pipeline.abundances_derive")[0])
from pipeline.abundances_derive import _ISPEC_SOLAR_ABUND_FILE, _SYNTH_ISOTOPE_FILE
TEFF, LOGG, MH, VT = 5772.0, 4.44, 0.0, 1.0
LO_nm, HI_nm = 335.0, 337.0
pack = ispec.load_modeled_layers_pack('/Users/ryanschmitt/codex/ispec/input/atmospheres/MARCS.GES/')
atm  = ispec.interpolate_atmosphere_layers(pack,{'teff':TEFF,'logg':LOGG,'MH':MH,'alpha':0.0},code='turbospectrum')
ab   = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
iso  = ispec.read_isotope_data(_SYNTH_ISOTOPE_FILE)
ll   = ispec.read_atomic_linelist(os.environ['NEARUV_LL'])
ll   = ll[(ll['wave_nm']>=LO_nm-1)&(ll['wave_nm']<=HI_nm+1)]
wave = np.arange(LO_nm, HI_nm, 0.0002)
def run(table, dFe=0.0):
    TABLE['name']=table
    fixed=[]
    if dFe:
        fixed=ab[ab['code']==26].copy(); fixed['Abund']=fixed['Abund']+dFe
    return np.asarray(ispec.generate_spectrum(wave, atm, TEFF, LOGG, MH, 0.0, ll, iso, ab,
        fixed, microturbulence_vel=VT, macroturbulence=0.0, vsini=0.0, R=0, verbose=0,
        code='turbospectrum', tmp_dir=SP))
S,Z = f'{SP}/jonabs_STOCK.dat', f'{SP}/jonabs_FeI_ZEROED.dat'
base=run(S); fz=run(Z); m=(base>0)&(fz>0)
d=lambda x: 1.0-x[m].mean()
target=d(fz)-d(base)
print(f"target depth change from removing Fe I bf: {target:+.5f}\n")
lad=[]
for x in (0.1,0.2,0.4,0.6,0.8,1.0,1.5):
    dd=d(run(S,x))-d(base); lad.append((x,dd))
    print(f"  A(Fe) {x:+.2f} dex -> depth change {dd:+.5f}{'   <-- brackets target' if dd>=target else ''}")
xs=np.array([a for a,_ in lad]); ys=np.array([b for _,b in lad])
if ys.max()>=target:
    i=int(np.argmax(ys>=target)); lo=(0.0,0.0) if i==0 else lad[i-1]
    eq=lo[0]+(target-lo[1])*(lad[i][0]-lo[0])/(lad[i][1]-lo[1])
    print(f"\n=> Fe I bf is worth {eq:.3f} dex in this window (BRACKETED, interpolated)")
else:
    print(f"\n=> NOT BRACKETED even at {xs.max()} dex: the window's abundance response "
          f"saturates below the Fe I bf effect. An opacity change and an abundance change "
          f"are not interchangeable here; report the opacity number, not a dex equivalent.")
