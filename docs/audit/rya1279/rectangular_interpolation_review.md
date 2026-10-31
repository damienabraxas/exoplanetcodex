# Simplex versus rectangular interpolation

RYA-1281 / parent RYA-1278, 2026-10-31. Mac input diagnostics only.
Baseline branch commit ae6601fbe2f6d8fe9bcf279fc564e275c128e82c;
fetched main f8c723a85f905425b1f7ae189f15a8763998ffc0 unchanged.

## Question and design

Does using all four corners reduce held-out atmosphere errors without sacrificing
joint thermodynamic/geometry/column consistency? Compare the existing SciPy
simplex implementation with a separate diagnostic bilinear operator. Production
interpolation is unchanged. Neither operator solves EOS or atmosphere equilibrium.

Fifteen native MARCS files span 5250/5500/5750/6000/6250 K and logg 4/4.5/5.
Three overlapping 500 K by 1 dex rectangles use outer corners and hold out the
midpoints of four edges and the center. These are 15 cell-target comparisons,
**11 distinct held-out model identities**, not 15 independent models/cells.
For each, retain both methods, all three existing field conventions and all three
column constructions: **270 cases**. Same composition and common native tauRoss
rows; no source edits, nearest fallback, extrapolation, fitted adjustment or
temperature offset. The diagonal/source weights are retained in every case.

Source hashes, code hashes and native closure checks are in report.json;
summary.csv records every combination; profiles.csv retains candidate/native
T, Pe, rho and outward height at each depth for all 90 field interpolants.
Electron-density error is calculated from Pe/T (kB cancels). This is an input
comparison, not an EOS validation or a coupled NLTE electron-density calculation.

## Results and tradeoff

Central targets, log Pe/rho and log layer thickness with log-column interpolation:

| Cell Teff range | Held-out Teff/logg | Method | Max T-profile error K | Max Pe error | Max rho error | Max column error | Max mass closure discrepancy |
|---|---|---|---:|---:|---:|---:|---:|
|5250-5750|5500/4.5|Simplex|22.45|5.29%|2.30%|2.17%|0.500%|
|5250-5750|5500/4.5|Bilinear|9.98|5.51%|5.13%|5.16%|1.509%|
|5500-6000|5750/4.5|Simplex|19.80|5.00%|0.87%|0.72%|0.545%|
|5500-6000|5750/4.5|Bilinear|7.88|4.67%|3.15%|3.12%|1.513%|
|5750-6250|6000/4.5|Simplex|20.35|4.76%|0.75%|0.94%|0.581%|
|5750-6250|6000/4.5|Bilinear|13.35|3.80%|2.59%|2.66%|1.414%|

Bilinear reduces temperature error at all three centers but worsens density,
column and mass-closure errors. Pe improves at two centers and worsens at one.
The methods agree on the edges (both reduce to the same 1D interpolation;
reported metric differences are zero). This is a useful implementation control,
not an independent atmosphere check.

Across all cell-targets with log-layer geometry, maximum T-profile error ranges
4.15-22.45 K for simplex and 4.15-14.30 K for bilinear. This is **not a stellar
Teff uncertainty**. Native rho*dz closure baselines span 0.0518-0.3073% for these
15 models, before interpolation. No acceptance threshold is inferred from this.

Density-integrated candidate columns close by construction, yet their held-out
native column errors span 0.668-5.176% for both methods across these targets;
sharing endpoints of this range does not imply identical interior results.
Exact self-consistency of one constructed integral is not physical certification.
No operator is promoted on the basis of a single metric or selected target.

## Verification and disposition

Sixteen focused tests pass, including an analytic bilinear cross-term field,
input-order permutation, exact corner weights, incomplete-cell/extrapolation
refusal and a known log-bilinear density/pressure/layer geometry. Existing
mass-integration and simplex input tests remain passing. No new verified software
defect was discovered; the prior characterized model-domain ledger entry is
extended with this evidence. The production interpolator still withholds columns;
the diagnostic bilinear result also withholds them. Native-only Fe opacity packet
checks, all previous O solves, Kitt Peak AIR observations, gold and adopted values
are unchanged. **Scientific HOLD** remains; no RT/SE or production output here.

Reproduce with the Bride Python environment and Mac host profile:
`python scripts/audit_bride_rectangular_interpolation.py --output <new-immutable-directory>`.
Original: experiments.bride_mac_grid/rya1281-rectangular-interpolation-v1.
Three byte-identical review copies: data/results/rya1281/rectangular_interpolation_v1.

## Next bounded task

Compare the actual solar/10 K targets using their nearest native rectangular
cells, rather than these deliberately coarser held-out cells, and retain method
spread as sensitivity only, not independent accuracy. Compare field errors by
depth before interpreting oxygen sensitivity; a declared photospheric subset
is not a computed line contribution function. Qualify EOS/pressure/geometry and
column support jointly before any interpolated cap is admitted. Then regenerate
background packets for the exact admitted structures; the existing native packet
must not be reused on changed T/ne/nH fields. Continue visible residual/blend and
producer-contract work independently. O then C/N, Fe, another star; full 3D and
all-band claims remain separate gates.
