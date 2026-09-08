"""RYA-1204 Step 0 — parse NORAD/Bautista Fe I level-resolved photoionization.

Per LS bound state: (2S+1, L, parity, ns) / (nr, ntot) / (BE_Ry, ac) / ntot x (E_Ry, sig_Mb).
Walked as a flat token stream — blank lines between records make line-based parsing fragile.
"""
import numpy as np, pickle, re

txt = open('fe1.px.txt', errors='replace').read()
m = re.search(r'^\s*26\s+25\s+P\s*$', txt, re.M)
tok = txt[m.end():].split()
p, states = 0, []
while p + 8 <= len(tok):
    S2, L, ip, ns = (int(tok[p+k]) for k in range(4))
    nr, ntot = int(tok[p+4]), int(tok[p+5])
    BE = float(tok[p+6])
    q = p + 8
    arr = np.array(tok[q:q+2*ntot], float).reshape(ntot, 2)
    states.append({'S2':S2,'L':L,'ip':ip,'ns':ns,'g':S2*(2*L+1),'BE':BE,
                   'E':arr[:,0].copy(),'sig':arr[:,1].copy()})
    p = q + 2*ntot
print(f'parsed {len(states)} LS bound states')
be = np.array([s['BE'] for s in states])
print(f'binding energies {be.min():.5f}-{be.max():.5f} Ry; ground (max BE) = {be.max():.6f} Ry '
      f'= {be.max()*13.605693:.3f} eV  [Fe I IP is 7.902 eV]')
print(f'total (E,sigma) points: {sum(len(s["E"]) for s in states):,}')

RY_EV, HC_EVA = 13.605693, 12398.419
e_lo, e_hi = HC_EVA/3780.0/RY_EV, HC_EVA/3000.0/RY_EV
print(f'\nnear-UV 3000-3780 A photon energy = {e_lo:.4f}-{e_hi:.4f} Ry '
      f'({e_lo*RY_EV:.3f}-{e_hi*RY_EV:.3f} eV)')
reach = [s for s in states if s['BE'] <= e_hi]
print(f'states ionisable somewhere in the band (BE <= {e_hi:.4f} Ry): {len(reach)} of {len(states)}')
print(f'  -> their excitation energies: '
      f'{(be.max()-np.array([s["BE"] for s in reach])).min()*RY_EV:.2f}-'
      f'{(be.max()-np.array([s["BE"] for s in reach])).max()*RY_EV:.2f} eV above ground')
print('\nGROUND STATE CANNOT CONTRIBUTE: its BE is '
      f'{be.max()*RY_EV:.3f} eV, threshold {HC_EVA/(be.max()*RY_EV):.0f} A -- far blueward of the band.')
pickle.dump(states, open('fe1_states.pkl','wb'))
