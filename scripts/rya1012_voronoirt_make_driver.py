"""
scripts/rya1012_voronoirt_make_driver.py
========================================
RYA-1012 -- build a runnable driver from VoronoiRT's compare_searchlight.jl.

WHY: the shipped script cannot run its own test. All three test calls at the
bottom of compare_searchlight.jl are commented out and the live entry point is
do_timing(DATA, QUADRATURE), which needs a Bifrost cube we do not have. The
npzwrite calls that produced the bundled reference arrays are commented out too,
so the script emits nothing even if it did run.

This patcher touches ONLY plumbing -- plot backend, quadrature input file, the
author's own commented-out npzwrite lines (redirected to a NEW output dir so the
bundled references are never overwritten), and the entry point. No algorithm,
solver or tolerance line is modified: nothing is tuned toward the target
(RYA-161 spirit). Run from the VoronoiRT repo root; VoronoiRT itself is NOT
committed to this repo (RYA-1012 G3 -- it is external MIT code we evaluated).

Result it reproduces: I_160_45_regular.npy bitwise identical; I_20_15_regular.npy
to 9 ULP (max rel 1.16e-15) -- at phi=195 deg, not the 15 deg its filename claims.
"""

import sys, io

src = "src/compare_searchlight.jl"
dst = "src/rya1012_driver.jl"
quad = sys.argv[1]
mode = sys.argv[2]  # regular | irregular | both

L = io.open(src, encoding="utf-8").read().split("\n")


def setl(n, text):  # 1-indexed
    L[n - 1] = text


NPZ = "../data/rya1012_out/"

# (1) plot backend: pyplot -> gr (avoids PyPlot/matplotlib; no numerical effect)
setl(8, "gr()")

# (2) quadrature INPUT file (angles come from data, not from a code edit)
setl(104, '    weights, θ_array, ϕ_array, n_angles = VoronoiRT.read_quadrature("%s")' % quad)
setl(192, '    _, θ_array, ϕ_array, _ = VoronoiRT.read_quadrature("%s")' % quad)

# (3) restore the author's OWN commented-out npzwrite calls, redirected to a NEW
#     output dir so the bundled reference arrays are never overwritten.
#     Integer angle naming matches the bundled filenames (I_20_15_*.npy).
setl(143, '            npzwrite("%sI_$(floor(Int,θ))_$(floor(Int,ϕ))_voronoi.npy", ustrip.(bottom_I))' % NPZ)
# symmetric write for the theta>90 branch (author only wrote the theta<90 one),
# needed to produce the I_160_45_voronoi counterpart.
setl(126, '            npzwrite("%sI_$(floor(Int,θ))_$(floor(Int,ϕ))_voronoi.npy", ustrip.(top_I))\n' % NPZ + L[125])
setl(148, '        npzwrite("%sx_voronoi.npy", x)' % NPZ)
setl(149, '        npzwrite("%sy_voronoi.npy", y)' % NPZ)
setl(161, '    npzwrite("%sx_regular.npy", ustrip.(x))' % NPZ)
setl(162, '    npzwrite("%sy_regular.npy", ustrip.(y))' % NPZ)
setl(207, '            npzwrite("%sI_$(floor(Int,θ))_$(floor(Int,ϕ))_regular.npy", ustrip.(I))' % NPZ)
setl(215, '            npzwrite("%sI_$(floor(Int,θ))_$(floor(Int,ϕ))_regular.npy", ustrip.(I))' % NPZ)

# (4) entry point: drop the shipped `do_timing(DATA, QUADRATURE)` tail (which
#     requires the Bifrost cube) and call the searchlight tests instead.
calls = {
    "regular": ["searchlight_regular()"],
    "irregular": ["searchlight_irregular()"],
    "both": ["searchlight_regular()", "searchlight_irregular()"],
}[mode]
L = L[:491] + calls

io.open(dst, "w", encoding="utf-8").write("\n".join(L) + "\n")
print("wrote %s  quad=%s  mode=%s" % (dst, quad, mode))
