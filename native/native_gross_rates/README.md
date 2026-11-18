# Native physical-rate observer (experimental only)

Derivative of OpenMULTI3D source32085d26c46e9919a986397b5514a9728abe2074,
with the separately held accurate profile patch. Upstream GPL notices remain in
all derivative source files. No installed source/executable/default changes.

`gross_rates.patch` repeats BB/BF rate integration against the SAME existing
formal radiation field, explicitly disabling preconditioning only in that local
rate-integral call. Solver accumulators are saved/restored exactly; no population,
radiation, opacity, global preconditioning or convergence control changes.
Serial one-thread use only. The source run must reproduce ALL original native
binary outputs byte-for-byte before accepting the observer.

Sidecar v1 `audit_gross` records:12 little-endian int32 fields
(magic1371,version1,kind BB1/BF2,radiation iteration,zs,nx,ny,nz,nr,nk,fs,dir),
then int32 global transition indices, float32 upward/thermal and downward/deviation
integrals, float64 actual n and nStar. Arrays are Fortran contiguous.
`audit_state` header:10 int32 fields(magic,version,population iteration,zs,
nx,ny,nz,nk,nrad,nline), lower/upper int32 indices, n/nStar float64,
normalized profile factors float32. All rates use s^-1; populations cm^-3.
Native radiation iteration12 precedes population update13. The reader additionally
requires exact n/nStar identity across all contributions; no nearest-stage fallback.

BB thermal/deviation sums are normalized by the saved source profile factors and
converted with source statistical weights/nStar in FP64. BF downward sums require
the native LTE ratio. The diagnostic sums source FP32 contributions in FP64, not
an exact replay of native accumulator rounding order. Native effective SE weights
in the ORIGINAL rates_O files remain semantically distinct.

This exports physical rate integrals at PRE-update radiation populations. It is
not a new native formal solve at POST-update populations, nor independently
validated solar/3D/all-band physics. No source or observational fitting/adoption.
