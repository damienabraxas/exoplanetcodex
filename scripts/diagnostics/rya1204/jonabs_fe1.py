"""Extract jonabs component 17 (FeI bf) as a 38-wavelength x 18-temperature array.
inabs.f reads, per temperature K: READ(IREAT,104)(XKAP(JJ,K),JJ=1,NLATB) -- so the
values are temperature-major, 18 blocks of 38."""
import numpy as np, re
src='/Users/ryanschmitt/codex/ispec/synthesizer/turbospectrum/DATA/jonabs_vac_v19.2.dat'
L=open(src,errors='replace').read().split('\n')
h=467                                   # 1-based line 468: the FeI header
assert L[h].startswith('FeI'), L[h]
nlatb=int(re.search(r'NLATB=\s*(\d+)',L[h+1]).group(1))
# wavelength grid
tok=[]; i=h+2
while len(tok)<nlatb: tok+=L[i].split(); i+=1
lam=np.array(tok[:nlatb],float)
tline=next(k for k in range(i,i+6) if 'ILOGT' in L[k])
ntetb=int(re.search(r'NTETB=\s*(\d+)',L[tline]).group(1))
ilogt=int(re.search(r'ILOGT=\s*(\d+)',L[tline]).group(1))
kvadl=int(re.search(r'KVADL=\s*(\d+)',L[h+1]).group(1))
tok=[]; i=tline+1
while len(tok)<ntetb: tok+=L[i].split(); i+=1
T=np.array(tok[:ntetb],float)
tok=[]; 
while len(tok)<nlatb*ntetb and i<len(L):
    s=L[i].split()
    if s and re.match(r'^[\d.]+E[+-]\d+$',s[0]): tok+=s
    elif tok: break
    i+=1
X=np.array(tok[:nlatb*ntetb],float).reshape(ntetb,nlatb)   # [temperature, wavelength]
print(f'component 17: NLATB={nlatb} NTETB={ntetb} ILOGT={ilogt} KVADL={kvadl} '
      f'(KVADL=0 => LINEAR interpolation in wavelength)')
print(f'wavelength nodes {lam.min():.0f}-{lam.max():.0f} A; temperatures {T.min():.0f}-{T.max():.0f} K')
print(f'values {X.min():.3e} .. {X.max():.3e}  (cm^2 per Fe I atom)')
np.savez('jonabs_fe1.npz', lam=lam, T=T, X=X)
# the band's nodes, at the temperatures nearest the continuum-forming layer
j=[k for k,w in enumerate(lam) if 2800<=w<=4200]
it=int(np.argmin(abs(T-5723)))
print(f'\nnearest tabulated temperature to 5723 K: {T[it]:.0f} K')
print('  lambda   sigma(cm^2)')
for k in j: print(f'  {lam[k]:7.0f}  {X[it,k]:.4e}')
