"""RYA-1204 Lever B — add near-UV molecular opacity and re-measure.

The molecules directory is swapped by monkeypatching iSpec's own glob root, so model,
abundances, atomic line list, opacity table and wave grid are the production call
unchanged and the ONLY difference is the molecular list.
"""
import os, sys, glob, re, tempfile, shutil
import numpy as np
sys.path.insert(0,'/Users/ryanschmitt/codex/ispec'); sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195')
import ispec, ispec.synth.turbospectrum as tsmod
from pipeline.abundances_derive import _ISPEC_SOLAR_ABUND_FILE, _SYNTH_ISOTOPE_FILE
SP=os.path.dirname(os.path.abspath(__file__))
MOL=os.path.join(SP,'molecules')

# iSpec symlinks `molecules_dir` into the run dir then globs "molecules/*.bsyn".
# Point that symlink at ours by patching os.symlink for that one call.
_orig_symlink=os.symlink
def symlink(src,dst,*a,**k):
    if str(dst).endswith('/molecules'): src=MOL
    return _orig_symlink(src,dst,*a,**k)
os.symlink=symlink

TEFF,LOGG,MH,VT=5772.0,4.44,0.0,1.0
LO_nm,HI_nm=335.0,337.0                       # the same window Lever A was measured in
pack=ispec.load_modeled_layers_pack('/Users/ryanschmitt/codex/ispec/input/atmospheres/MARCS.GES/')
atm=ispec.interpolate_atmosphere_layers(pack,{'teff':TEFF,'logg':LOGG,'MH':MH,'alpha':0.0},code='turbospectrum')
ab=ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE); iso=ispec.read_isotope_data(_SYNTH_ISOTOPE_FILE)
ll=ispec.read_atomic_linelist('/Users/ryanschmitt/codex/rya1195/data/linelists/ispec_nearuv_3000_3780/atomic_lines.tsv')
ll=ll[(ll['wave_nm']>=LO_nm-1)&(ll['wave_nm']<=HI_nm+1)]
wave=np.arange(LO_nm,HI_nm,0.0002)

def run(use_mol,dFe=0.0):
    fixed=[]
    if dFe:
        fixed=ab[ab['code']==26].copy(); fixed['Abund']=fixed['Abund']+dFe
    return np.asarray(ispec.generate_spectrum(wave,atm,TEFF,LOGG,MH,0.0,ll,iso,ab,fixed,
        microturbulence_vel=VT,macroturbulence=0.0,vsini=0.0,R=0,verbose=0,
        code='turbospectrum',use_molecules=use_mol,tmp_dir=SP))

no_mol=run(False); with_mol=run(True)
m=(no_mol>0)&(with_mol>0)
d=lambda x:1.0-x[m].mean()
print(f'valid pixels {m.sum()}/{m.size}')
print(f'  molecules OFF (production today): depth {d(no_mol):.5f}')
print(f'  molecules ON  (Lever B)         : depth {d(with_mol):.5f}')
print(f'  depth change from molecules     : {d(with_mol)-d(no_mol):+.5f}')

# Calibrate into dex on the SAME footing, and BRACKET rather than extrapolate (RYA-1202).
base=no_mol; tgt=d(with_mol)-d(base); lad=[]
print('\nabundance ladder at molecules OFF (for the dex conversion):')
for x in (0.05,0.10,0.20,0.40):
    dd=d(run(False,x))-d(base); lad.append((x,dd)); print(f'   A(Fe) {x:+.2f} -> {dd:+.5f}')
xs=np.array([a for a,_ in lad]); ys=np.array([b for _,b in lad])
print(f'\ntarget (molecular depth change) {tgt:+.5f}')
if ys.max()>=tgt>=0:
    i=int(np.argmax(ys>=tgt)); lo=(0.0,0.0) if i==0 else lad[i-1]
    eq=lo[0]+(tgt-lo[1])*(lad[i][0]-lo[0])/(lad[i][1]-lo[1])
    print(f'=> LEVER B is worth {eq:.3f} dex (BRACKETED). Adding opacity LOWERS the derived')
    print(f'   abundance, so the near-UV A(Fe) moves DOWN by ~{eq:.3f} dex.')
else:
    print('=> not bracketed by this ladder; report the depth change, not a dex equivalent.')
np.save('leverB.npy',np.array([d(no_mol),d(with_mol),tgt]))
