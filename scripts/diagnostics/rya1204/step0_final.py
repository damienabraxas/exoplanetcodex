"""RYA-1204 Step 0 — the decisive plot and verdict.

🔴 NORMALISATION, ESTABLISHED BY MEASUREMENT NOT ASSUMPTION. jonabs component 17 stores
   sum_i g_i exp(-E_i/kT) sigma_i, i.e. WITHOUT dividing by the partition function -- the
   ratio jonabs/(Bautista/U) tracks U across the whole temperature grid (25.06 -> 20.4 at
   1000 K, 47.14 -> ~52 at 8000 K). Comparing before fixing this would have shown a
   spurious factor ~40 and "proved" a lever that does not exist.
"""
import numpy as np, pickle
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
RY_EV,HC_EVA,K_EV=13.605693,12398.419,8.617333262e-5
states=pickle.load(open('fe1_states.pkl','rb')); BE_g=max(s['BE'] for s in states)
J=np.load('jonabs_fe1.npz'); jlam,jT,X=J['lam'],J['T'],J['X']

def bautista(l,t):
    """sum_i g_i exp(-Eexc/kT) sigma_i, in cm^2 -- jonabs's own convention."""
    E=HC_EVA/np.asarray(l,float)/RY_EV; kT=K_EV*t/RY_EV
    tot=np.zeros_like(E)
    for s in states:
        w=s['g']*np.exp(-(BE_g-s['BE'])/kT)
        sg=np.interp(E,s['E'],s['sig'],left=0.,right=s['sig'][-1])
        tot+=w*np.where(E>=s['BE'],sg,0.)
    return tot*1e-18

T=5900.0; it=int(np.argmin(abs(jT-T)))
lam=np.arange(3000.,3781.,0.5)
full=bautista(lam,T)
band=(jlam>=2800)&(jlam<=4200); nodes=jlam[band]
jn=X[it][band]                                   # what jonabs actually tabulates
bn=bautista(nodes,T)                             # Bautista at the same nodes
coarse_j=np.interp(lam,nodes,jn)                 # jonabs as the code uses it (KVADL=0)
coarse_b=np.interp(lam,nodes,bn)                 # Bautista sampled the same way

print(f"T = {jT[it]:.0f} K (jonabs grid), band 3000-3780 A\n")
print(f"{'lambda':>8} {'jonabs':>12} {'Bautista':>12} {'jonabs/Baut':>12}")
for k,w in enumerate(nodes): print(f"{w:8.0f} {jn[k]:12.4e} {bn[k]:12.4e} {jn[k]/bn[k]:12.3f}")
r=jn/bn
print(f"\n(b) ABSOLUTE SCALE at the nodes: median {np.median(r):.3f}x  range {r.min():.3f}-{r.max():.3f}")
print(f"(a) UNDER-SAMPLING (Bautista full vs Bautista at these nodes, linear):")
print(f"    band mean ratio {full.mean()/coarse_b.mean():.4f}x")
print(f"(a+b) what the CODE uses vs full-resolution Bautista:")
print(f"    band mean ratio {coarse_j.mean()/full.mean():.4f}x  "
      f"(jonabs band mean {coarse_j.mean():.4e}, Bautista full {full.mean():.4e})")

fig,ax=plt.subplots(2,1,figsize=(11,8),sharex=True,gridspec_kw={'height_ratios':[3,2]})
ax[0].semilogy(lam,full,lw=.6,color='#4477dd',label='Bautista / NORAD, full resolution')
ax[0].semilogy(lam,coarse_j,lw=2,color='#dd6644',label='jonabs component 17, as the code uses it (4 nodes, linear)')
ax[0].plot(nodes,jn,'o',color='#dd6644',ms=8,zorder=5,label='jonabs tabulated nodes')
ax[0].plot(nodes,bn,'s',color='#4477dd',ms=7,mfc='none',zorder=5,label='Bautista at the same nodes')
ax[0].set_ylabel(r'$\sum_i g_i e^{-E_i/kT}\sigma_i$   (cm$^2$)')
ax[0].set_title(f'RYA-1204 Step 0 — Fe I bound-free, near-UV, T = {jT[it]:.0f} K\n'
                f'band-mean agreement {coarse_j.mean()/full.mean():.3f}x — the sampling is NOT the lever',fontsize=11)
ax[0].legend(fontsize=8.5); ax[0].grid(alpha=.25)
ax[1].semilogy(lam,coarse_j/full,lw=.6,color='#666')
ax[1].axhline(1,color='k',lw=1); ax[1].axhline(2,color='#cc3333',ls='--',lw=1,label='Bell 2001 needed 2x')
ax[1].set_ylabel('jonabs / Bautista'); ax[1].set_xlabel(r'wavelength ($\AA$)')
ax[1].legend(fontsize=8.5); ax[1].grid(alpha=.25); ax[1].set_ylim(0.02,60)
plt.tight_layout(); plt.savefig('rya1204_step0_fe1bf.png',dpi=140)
print('\nplot -> rya1204_step0_fe1bf.png')
np.savez('step0_final.npz',lam=lam,full=full,coarse_j=coarse_j,nodes=nodes,jn=jn,bn=bn,T=jT[it])

# ── the verdict is not a property of one temperature ─────────────────────────────────
# Swept over the line-forming range, so a single-T artifact cannot carry the conclusion:
#   T(K)   undersampling   code-vs-Bautista
#   3800      0.935            1.874
#   4500      0.956            1.663
#   5200      0.977            1.539
#   5900      1.000            1.441
#   6600      1.025            1.364
#   7300      1.052            1.332
#   8000      1.078            1.299
# Under-sampling stays within +/-7% of unity everywhere -- the 4-node grid is not the
# defect -- and the shipped table sits ABOVE Bautista at EVERY line-forming temperature,
# never the 2x BELOW that Bell 2001 needed.
