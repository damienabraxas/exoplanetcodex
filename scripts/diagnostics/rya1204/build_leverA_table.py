"""RYA-1204 Lever A — the PHYSICAL re-tabulation, built and measured rather than argued.

Replace jonabs component 17 across 3000-3780 A with the resonance-averaged Bautista/NORAD
cross-section on a FINE grid (20 A spacing), keeping every node outside the band untouched.
That is both halves of the ticket's Lever A: the correct values AND the sampling the coarse
grid cannot represent.

🔴 NOTHING IS SCALED TO A TARGET. The values are Bautista's own, resonance-AVERAGED onto
each grid interval (Bautista & Pradhan's own prescription), converted into jonabs's
convention -- sum_i g_i exp(-E_i/kT) sigma_i, no partition function -- which was
established by measurement in step0_final.py, not assumed.
"""
import numpy as np, pickle, re
RY_EV,HC_EVA,K_EV=13.605693,12398.419,8.617333262e-5
states=pickle.load(open('fe1_states.pkl','rb')); BE_g=max(s['BE'] for s in states)
SRC='/Users/ryanschmitt/codex/ispec/synthesizer/turbospectrum/DATA/jonabs_vac_v19.2.dat'
L=open(SRC,errors='replace').read().split('\n')
h=467; assert L[h].startswith('FeI')
nlatb=int(re.search(r'NLATB=\s*(\d+)',L[h+1]).group(1))
i=h+2; tok=[]
while len(tok)<nlatb: tok+=L[i].split(); i+=1
lam_old=np.array(tok[:nlatb],float)
tline=next(k for k in range(i,i+6) if 'ILOGT' in L[k])
ntetb=int(re.search(r'NTETB=\s*(\d+)',L[tline]).group(1))
j=tline+1; tok=[]
while len(tok)<ntetb: tok+=L[j].split(); j+=1
T=np.array(tok[:ntetb],float)
end=j
tok=[]
while len(tok)<nlatb*ntetb and end<len(L):
    s=L[end].split()
    if s and re.match(r'^[\d.]+E[+-]\d+$',s[0]): tok+=s
    elif tok: break
    end+=1
X_old=np.array(tok[:nlatb*ntetb],float).reshape(ntetb,nlatb)

def sigsum(lam,t):
    E=HC_EVA/np.asarray(lam,float)/RY_EV; kT=K_EV*t/RY_EV
    tot=np.zeros_like(E)
    for s in states:
        w=s['g']*np.exp(-(BE_g-s['BE'])/kT)
        sg=np.interp(E,s['E'],s['sig'],left=0.,right=s['sig'][-1])
        tot+=w*np.where(E>=s['BE'],sg,0.)
    return tot*1e-18

# fine grid inside the band; keep every node outside it exactly as shipped
NEW=np.arange(3000.,3781.,20.)
keep=lam_old[(lam_old<3000-1)|(lam_old>3780+1)]
lam_new=np.sort(np.concatenate([keep,NEW]))
X_new=np.zeros((ntetb,len(lam_new)))
for k,t in enumerate(T):
    # outside the band: the shipped values, interpolated onto any node we did not move
    out=np.interp(lam_new,lam_old,X_old[k])
    # inside the band: Bautista, resonance-AVERAGED over each 20 A interval rather than
    # sampled at a point -- a point sample of a resonance forest is a lottery.
    b=np.empty(len(lam_new))
    for m,w in enumerate(lam_new):
        if 3000-1<=w<=3780+1:
            fine=np.arange(max(3000.,w-10.),min(3780.,w+10.)+0.25,0.25)
            b[m]=sigsum(fine,t).mean()
        else:
            b[m]=out[m]
    X_new[k]=np.where((lam_new>=3000-1)&(lam_new<=3780+1),b,out)
print(f'nodes {nlatb} -> {len(lam_new)}  ({len(NEW)} inside 3000-3780 A, was 3)')
i0=(lam_new>=3000)&(lam_new<=3780); it=int(np.argmin(abs(T-5900)))
old_band=np.interp(np.arange(3000.,3781.,1.),lam_old,X_old[it])
new_band=np.interp(np.arange(3000.,3781.,1.),lam_new,X_new[it])
print(f'band mean at {T[it]:.0f} K: shipped {old_band.mean():.4e} -> re-tabulated '
      f'{new_band.mean():.4e}  ({new_band.mean()/old_band.mean():.3f}x)')

def fmt(a,per=6,w=10,p=3):
    out=[]
    for s in range(0,len(a),per):
        out.append(''.join(f'{v:{w}.{p}f}' if v>=1 else f'{v:{w}.{p}f}' for v in a[s:s+per]))
    return out
def fmte(a,per=6):
    return [''.join(f' {v:9.3E}' for v in a[s:s+per]) for s in range(0,len(a),per)]
block=[L[h]]
block.append(re.sub(r'NLATB=\s*\d+',f'NLATB={len(lam_new)}',L[h+1]))
block+= [''.join(f'{v:10.0f}.' if False else f'{v:10.0f}' for v in lam_new[s:s+6]) for s in range(0,len(lam_new),6)]
block.append(L[tline])
block+= [''.join(f'{v:10.1f}' for v in T[s:s+6]) for s in range(0,len(T),6)]
for k in range(ntetb): block+=fmte(X_new[k])
open('jonabs_LEVER_A.dat','w').write('\n'.join(L[:h]+block+L[end:]))
print('wrote jonabs_LEVER_A.dat')
