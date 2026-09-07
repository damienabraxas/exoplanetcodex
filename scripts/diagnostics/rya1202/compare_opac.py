"""RYA-1202 — what does Fe I bound-free photoionization actually carry in our near-UV
continuous opacity? Measured by running babsma twice through the real code path with only
the Fe I bf cross-sections zeroed.

MODELOPAC layout: "'MRXF' ndep Teff" / nlam / nlam wavelengths / then per depth an
8-number header (index 2 is T) followed by nlam interleaved (absorption, scattering).
"""
import numpy as np

def read(path):
    t = open(path).read().split()
    assert t[0] == "'MRXF'"
    ndep, nlam = int(t[1]), int(t[2 + 1])
    i = 4
    lam = np.array(t[i:i+nlam], float); i += nlam
    T = np.empty(ndep); A = np.empty((ndep, nlam)); S = np.empty((ndep, nlam))
    for d in range(ndep):
        hdr = np.array(t[i:i+8], float); i += 8
        T[d] = hdr[2]
        blk = np.array(t[i:i+2*nlam], float).reshape(nlam, 2); i += 2*nlam
        A[d], S[d] = blk[:,0], blk[:,1]
    assert i == len(t), (i, len(t))
    return lam, T, A, S

lam, T, a_s, _ = read('run_stock/opac.out')
l2,  T2, a_z, _ = read('run_fezero/opac.out')
assert np.array_equal(lam, l2) and np.allclose(T, T2)
band = (lam >= 3000) & (lam <= 3800)
print(f"grid {lam.size} pts {lam.min():.0f}-{lam.max():.0f} A ({band.sum()} in the near-UV "
      f"policy window), {T.size} depths, T {T.min():.0f}-{T.max():.0f} K\n")

frac = (a_s - a_z) / a_s                       # Fe I bf share of continuous absorption
# The near-UV continuum forms around T ~ Teff; report there and across the profile.
d_cont = int(np.argmin(abs(T - 5772.0)))
print(f"CONTINUUM-FORMING LAYER (T={T[d_cont]:.0f} K, nearest Teff):")
f = frac[d_cont][band]
print(f"  Fe I bf share of kappa_cont: median {np.median(f)*100:.2f}%  "
      f"min {f.min()*100:.2f}%  max {f.max()*100:.2f}%")

print("\nby depth (median over 3000-3800 A):")
for d in range(0, T.size, 7):
    print(f"  depth {d:2d}  T={T[d]:6.0f} K   {np.median(frac[d][band])*100:6.2f}%")

print("\nby wavelength (continuum-forming layer):")
for w in (3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800):
    j = int(np.argmin(abs(lam - w)))
    print(f"  {lam[j]:7.0f} A   {frac[d_cont][j]*100:6.2f}%")

# First-order lever: for a weak line, depth ~ kappa_line/kappa_cont, so a change in
# continuous opacity maps to d log A of the same size and opposite sign.
dlog = np.log10(a_s / a_z)
sel = dlog[:, band]
print(f"\nd log10(kappa_cont) attributable to Fe I bf, 3000-3800 A:")
print(f"  continuum-forming layer: {np.median(dlog[d_cont][band]):+.4f} dex")
print(f"  over all depths: median {np.median(sel):+.4f}  max {sel.max():+.4f} dex")
