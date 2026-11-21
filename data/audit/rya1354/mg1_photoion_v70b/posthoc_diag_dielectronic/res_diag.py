"""POST-HOC diagnostic (v70b): 3p2 1S resonance width vs the dielectronic polarisation scale. Not a declared test."""
import sys, json, numpy as np
sys.path.insert(0, '/Users/ryanschmitt/codex/rya1354/scripts')
import rya1354_photoion as P, rya1354_twovalence_model as T
scale, rbox = float(sys.argv[1]), int(sys.argv[2])
T.ALPHA2_SCALE = scale
P.syms = lambda: [(0, 0, 0), (1, 0, 1)]
sp, B = P.setup(); P.add_continuum(B, rbox)
st, thr = P.solve_all(sp, B)
d1 = st[(1, 0, 1)]; j = next(j for j in d1['match'] if d1['terms'][j][1] == P.INIT_1P); i1 = d1['match'][j]
f, fv, de = P.f_and_fv(B, d1, (1, 0, 1), st[(0, 0, 0)], (0, 0, 0), i1)
out = '/private/tmp/claude-501/-Users-ryanschmitt/1e68bc2c-fab3-4f2d-8d42-1e6f87e7f11a/scratchpad/diag'
np.savez(f'{out}/res_{scale}_{rbox}.npz', e=st[(0, 0, 0)]['e'] - thr, f=f, fv=fv)
print('done', scale, rbox)
