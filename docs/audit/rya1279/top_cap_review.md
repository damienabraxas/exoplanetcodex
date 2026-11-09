# Finite-column top boundary: the triplet surface collapse is resolved in controls

RYA-1281 / RYA-1278, continuing from63e1dc6e on stacked draft PR #598. All work
below is serial Mac diagnostics. Scientific admission remains **HOLD**: the
triplet excitation result improves substantially, but O II and producer/background
contracts remain unresolved. No default solver, adopted abundance or source atom
was changed.

## Result at unchanged rates and atmosphere

A native-column, constant-source cap removes most of the sharp first-point
triplet-population discrepancy for both tested collision contracts. The published
reference departures were never supplied to the solve or used to choose a cap.

| Control | Top lower b | Top upper b | Top O II b | Upper b maximum difference, all depths |
|---|---:|---:|---:|---:|
| Independent departure reference |0.988478|1.281391|8.168057|--|
| Sourced MN pair, zero incoming boundary (prior phase) |0.025152|0.031916|0.393150|1.249476|
| Same sourced pair, finite-column cap |0.973814|1.279690|16.252022|0.026197|
| Raw LTE-restoring pair, finite-column cap |0.966449|1.269186|18.546761|0.036755|

For the sourced pair the all-depth upper-triplet difference falls about48-fold.
Its upper difference over the predeclared logtau[-3,0] subset remains0.005560;
the raw-pair cap result is0.005396. This subset is not a contribution function.
[Full-depth departure plot](../../../data/output/rya1281/top_cap_v1/departures.png) and
[fixed triplet profiles](../../../data/output/rya1281/top_cap_v1/flux.png) retain both controls.
The O II discrepancy is material: all-level/depth maximum differences are8.08396
and10.37870. These are **not** full-population-parity results.

Intrinsic fixed +/-0.5 A triplet areas barely move:

| Control | Areas (mA) |
|---|---|
| Sourced pair, prior zero boundary |68.777573 /59.230425 /46.280367|
| Sourced pair, cap outer face |68.777633 /59.230434 /46.280340|
| Raw pair, cap outer face |68.777897 /59.229484 /46.278823|

Thus a nearly unchanged visible spectrum can accompany a dramatically different
surface population. No abundance fit, Kitt Peak renormalization, second telluric
correction or wavelength-frame change was made. These new observables cover the
777 nm VIS triplet; no all-band or new-star validation is claimed.

## Native input and the boundary approximation

The [official MARCS documentation](https://marcs.astro.uu.se/docs.html) and
[2018 format note](https://marcs.astro.uu.se/) establish RHOX as column mass,
expressed in g/cm^2. The existing printed native column is now preserved as
kg/m^2 (multiply by10); it is not reconstructed from pressure. Native density
converts g/cm^3 to kg/m^3 (multiply by1000). At the benchmark top:

    m_top =0.09156734 kg/m^2
    rho_top =1.0854e-6 kg/m^3
    L_effective =m_top/rho_top =84362.76027 m

This is an effective opacity-per-mass cap length, not an invented geometric mesh
or an independently calculated atmosphere above the grid. At each iteration,
actual top extinction and emissivity define a constant-source transfer solution:

    tau_cap(lambda,mu) =chi_top(lambda)*L_effective/mu
    I_in =S_top(lambda)*(1-exp(-tau_cap)),  S_top=eta_total/chi_top.

The source includes thermal/atomic emission and background scattering using the
same lagged J as the backend formal solution. It is not a fitted Planck
irradiation field. The expression matches the finite-depth prescription in held
pymulti4d6aa3c TRANSP/TRANF for a mass-column interpretation, but the actual
modified grid producer's input-depth convention is still unidentified.

The sourced case has cap vertical optical depths3.7655 at130.204 nm,4.4328 at
135.560 nm,9.3411e-5 at277.248 nm, and2.1846e-5 at500 nm. Small visible optical
depth does not justify zero incoming radiation for those UV wavelengths.
The printed MARCS tau500 at the same point is1.20865e-5; it is not equal to this
constant-opacity-per-mass cap's tau500. No scale factor was introduced to force
them equal. A native-optical-depth boundary and a uniform mass-column cap must
not be silently called identical; exact producer/background input remains needed.

The cap is opt-in. It requires native column metadata, nonnegative extinction
and emission, and angle-independent static coefficients. Unsupported conditions
fail. The existing atmosphere interpolation does not validate column-mass
interpolation: that field is withheld on interpolated structures, and the runner
refuses applying this cap to such a case. No borrowed nearest-node column.

## Numerical controls and failed attempts

The first two predeclared runs refreshed the cap once per SE population step.
Both reached the1500-iteration limit without joint convergence:

| Original attempt | Seconds | Disposition |
|---|---:|---|
| Sourced MN + cap |110.6826|FAILED; no admitted population/flux output|
| Raw LTE pair + cap |121.3685|FAILED; no admitted population/flux output|

Their complete manifests/histories remain immutable. The analysis explicitly
records failed controls; it does not silently select successful curves.

A separate predeclared numerical plan solves radiation and the cap at **fixed
populations** before each SE step. The equations, source coefficients, native
column, abundance and outer tolerances are unchanged. Inner J and cap changes
must each reach1e-8; the outer population/J/cap/rate tolerances remain1e-6,
conservation1e-10. The600-second/case and1500-outer-iteration budgets were retained.

| Nested control | Outer steps | Inner formal calls | Seconds | Final outer J change | Final cap change |
|---|---:|---:|---:|---:|---:|
| Sourced MN + cap |291|2485|199.8076|9.3512e-7|2.4812e-7|
| Raw LTE pair + cap |836|5491|487.5533|4.0423e-7|2.8011e-7|

Terminal per-step physical rate residuals are2.61e-15/2.64e-15 and particle errors
6.66e-16/8.88e-16. No invalid Ng update was rejected in the successful runs.
No tolerance, population floor, abundance or physical rate was adjusted to force
convergence. Both are full O41, native MARCS5750/4.5, xi1, A(O)8.70, rays3/qscale1.

## Independent checks of the implementation and final states

Tests cover zero-column, thin/thick, pure-emission and angular slab limits;
actual backend inward/outward intensities against the analytic isothermal slab;
scattering inclusion; deep-copy isolation; and real native-column ray sampling.
The MARCS parser's column conversion is checked independently on synthetic input.

The backend's simple formal_sol path intentionally does not record DepthData.
Coefficient capture therefore uses the full formal-solution path without any SE
update. Diagnostic ray tracing also leaves J-quadrature weights unset; explicit
normalized flux-quadrature weights are supplied for that capture. Nonfinite
updates and missing coefficient captures fail rather than producing an empty cap.
These details were found by failing backend tests and corrected before the
nested science runs. The original failed science cases never reached flux export.

Outgoing rays are transferred through the same constant-source cap to its outer
face, including emission and attenuation. Separate grid-point fluxes are saved;
sharp diagnostic-line coefficients are evaluated on the actual diagnostic grid,
not interpolated from coarse SE line sampling. Line and continuum rays use their
respective solved cap sources.

A fresh fixed-population replay reconstructs J/cap and R+C without SE updates:

| Collision contract | Own NLTE maximum imbalance | Local LTE | Independent b on local LTE |
|---|---:|---:|---:|
| Sourced MN + cap |4.645995e-7|0.877765|0.329475|
| Raw pair + cap |2.008346e-7|0.893265|0.389598|

The final populations remain equilibria within the declared outer tolerance when
radiation is recomputed. Foreign departures remain non-equilibria of the tested
local physics; they are not native MULTI rate dumps.

Independent quadrature of the sourced-cap lower-triplet photoionization rate
still agrees with the backend to0.0274% at the top:16248.9998 versus16244.5438
s^-1. Its J-weighted rate is64.6221 times the local-Planck rate. The boundary
change therefore does not remove the remaining nonlocal UV/background question.

## Limits and next work

The cap resolves a demonstrated surface-transfer artifact in these controlled
comparisons. It does not validate the complete atom or identify the exact
producer. Native LTE populations, electron density, background opacities and
transition-rate/J outputs are the next important reference evidence. Departure
coefficients alone do not certify identical LTE denominators across engines.
The inspected held LTEPOP includes the usual0.5 electron phase-space factor;
no missing-factor-of-two claim follows from the O II discrepancy.

Next: check the producer's actual mass/optical-depth boundary convention and LTE
normalization before assigning the O II difference uniquely to photoionization.
Use the preserved rate budgets to separate background/EOS/atomic-input effects.
Do not fit an opacity multiplier or cap length to the reference. Matched solar
parameter work with a cap also needs a validated column interpolation contract.
Then finish solar O diagnostics, proceed to C/N, then Fe, then another star.

## Evidence and validation

Original stores under CODEX_WORK_ROOT/bride_mac_checks:
rya1281-top-cap-v1 (two failed controls), rya1281-top-cap-nested-v1 (two converged),
rya1281-top-cap-analysis-v2, rya1281-top-cap-author-rates-v1,
rya1281-top-cap-raw-rates-v1, rya1281-top-cap-photoionization-v1.
Earlier analysis versions remain. Final generator/inputs/population/cap snapshots
are hash-receipted. Full rate arrays remain external.

Review package:28 registered artifacts in data/results/rya1281/top_cap_v1;3 inspected figures under
data/output/rya1281/top_cap_v1. Reproduction uses the two committed top-cap plans,
recorded code commits5ddb7571 (single-refresh) and3aca3092 (nested), and the CLI
commands in GENERATORS.yaml with **new immutable output directories**.

171 focused science tests plus3 installed iSpec checks pass. Fixed-physics grid
admission refuses these experimental boundary cases. Detailed guard/main/library
status is in progress.md and the issue handoff. No full-suite/CI-green claim,
production emission, scientific signoff, merge or publication.
