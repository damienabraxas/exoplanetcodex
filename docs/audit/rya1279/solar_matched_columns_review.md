# RYA-1281 matched native and canonical-solar O columns

Scientific HOLD. This is the next physical-input control after the completed
[O777shape milestone](O_shape_milestone_review.md). It tests native5750/logg4.5,
temperature-only5772/4.5, gravity-only5750/4.438and combined5772/4.438, with
the SAME current O atom/collision/background/xi conventions. Canonical solar
labels come from current `get_star_params('solar')`; they are not optimized.

## Source ownership, interpolation and cap

All four held native nodes in the nearest5750–6000/g4–4.5cell are fingerprinted.
They share56native tauRoss rows,92element composition and hydrogen mass fraction.
The baseline uses native MARCS source precision. Other columns use the previously
audited SciPy simplex/log-positive Pe/rho/log-layer-thickness interpolation,
with all56depths, no extrapolation/nearest fallback or new equilibrium solve.
This is an experimental input construction, not independent atmosphere/EOS or
physical interpolation accuracy. It changes whole structures, not a constant
kinetic-temperature offset at every depth.

The source-owned ne=Pe/(LW k*T) and nH=rho*X/(LW Amu*H mass) boundary is used
for ALL four columns. xi remains1km/s. No exported5750native-network H/CN/Pe
packet is transplanted onto an interpolated atmosphere. Original MARCS files,
native inventory and prior spectra remain unchanged. No molecular provider is
inserted; the native background molecule table stays explicitly empty, while
existing H/proton background continuum remains active.

The interpolator deliberately removes column mass. No interpolated column mass
is invented for the boundary. ALL four cases retain the original geometric cap
length84362.76027271051m, explicitly experimental and not a validated new solar
column mass. This freezes one boundary intervention while structures change.
Thermal/charge/chemistry equilibrium is not re-solved.

The native-node interpolation identity screen passes. Maximum source-node
differences:height1.16415e-10m,Pe2.27374e-13Pa,rho2.71051e-19kg/m3,
T/tauRoss/tau500exact. The common-field/source/geometry guards are unchanged.

## Preflight input responses

Four dry-run Contexts pass, full6796point source O SE layout and56depths retained.
Maximum all-depth changes relative to the fresh source-precision baseline:

|Column|Max kinetic T change(K)|Max relative Pe change|Max relative nH change|Max height change(m)|
|---|---:|---:|---:|---:|
|T-only|28.2832|0.0406300|0.0187219|900.2517|
|g-only|8.1592|0.0684314|0.0935210|92058.9093|
|Solar combined|33.4268|0.0563149|0.1104920|93104.1439|

These are input sensitivities, not model errors or oxygen contribution-weighted
values. The solar simplex weights are0.088/0.788/0.124on6000g4.5/5750g4.5/
5750g4.0. No higher-dimensional or validated3Dphysical interpolation claim.

The baseline's source-precision versus previous native-export differences are
Tmax0.000390625K,Pe1.953125e-5Pa,nerelative7.82788e-8; nH,height andxi
are exact. These differences are recorded before assigning T/g responses.
The fresh baseline replays the earlier normalized flux within1.05286e-8absolute
and true continuum1.62180e-9relative, under prospective1e-6numerical screens.
Those screens/identities are not physical accuracy estimates.

## Solves and independent admission

Attempt1 uses the SAME saved completed v2 model-H fractional populations,
separately renormalized to each O total. No sequential matrix warm starts,
imported J/rates, abundance/RV/continuum tuning or source correction.
Source O8.70/41levels/256lines/22continua, bare native background/passive atoms,
source TAUTSP/directional collisions and author MN LTE-pair hypotheses remain.
No added Fe/CN line provider or independent physical-source admission.

Inherited controls:5rays/q2, joint population/rate/J/cap1e-6, nested radiation
1e-8,max100inner/1500outer,600seconds per case, Ng with invalid-population
rejection,7emergent rays/4001AIR points7763–7783A. All cases serial on Mac,
OMP/OPENBLAS1. Source/helper/backend/constant/atom receipts and starting layouts
are checked. Original physical numerical kernels are not changed.

Each completed case requires independently rebuilt per-class collision sum
against active.C, transition radiative rates, net/gross flow and O conservation,
plus a fresh standalone continuum. Final populations must remain exact in both
emergent-ray contexts. The whole-case600second deadline is also checked after
flux extraction. A timeout is preserved, never promoted to convergence.

The downstream analyzer requires ALL four completed manifests, verifies their
common source/physics/owners/cap, replays fresh saved C+R rates independently,
checks conservation and the declared common initialization invariant, then measures the same
native-grid fixed-shoulder O metrics against the frozen Kitt/IAG observations.
Both RVs/two shoulder schemes/nominal solar broadening stay fixed; no new
continuum fit. T/g finite path averages and interaction are correlated controlled
differences, not an inverse fit, preferred atmosphere or independent error budget.

## Reproduction

Load `config/hosts/mac.env`; held LW virtualenv with OMP/OPENBLAS1, cases serial:

```
python scripts/run_bride_solar_matched_columns.py --case baseline --output <fresh-matrix>/baseline
python scripts/run_bride_solar_matched_columns.py --case temperature --output <fresh-matrix>/temperature
python scripts/run_bride_solar_matched_columns.py --case gravity --output <fresh-matrix>/gravity
python scripts/run_bride_solar_matched_columns.py --case solar --output <fresh-matrix>/solar
```

Use the project analysis virtualenv for the postprocessing:

```
python scripts/analyze_bride_solar_matched_columns.py --source <completed-matrix> --output <fresh-matched-analysis>
```

Raw logs are `/private/tmp/rya1281-solar-matched-<case>-v1.log`; input-only logs
use `rya1281-solar-matched-dry-<case>.log`. No overwrite of previous attempts.
Dry Contexts completed in0.63–0.78s. Baseline real solve passed3outer iterations,
16.85937s,peak696254464bytes/one thread, fresh independent rate residual
3.95623e-7 and conservation6.66134e-16.

Tests so far:7held-LW matched-input/H/cap tests pass1.84s;22interpolator/column
tests pass0.37s;6source-rate/spline/charge tests pass1.90s;4analysis-shape/
finite-accounting tests pass in the project environment. Source/code receipts
remain fixed during the serial matrix. No install/download/benchmark/Sirius/CI.

## Preserved attempt1 and predeclared convergence controls

The first matrix finishes with baseline and gravity admitted. Gravity passes
176outer iterations/265.300727s,peak685604864bytes, independent rate3.45056e-7.
Temperature hits600seconds with last complete iteration303:J/cap3.25505e-6,
population3.25276e-6. Combined solar hits600seconds at iteration312:J/cap
1.08537e-6,population1.08462e-6. These near-threshold histories are not final
audits: no final population/rate/continuum/spectrum/report is emitted for either.
The all-four admission condition is not met; no partial factorial is promoted.

Attempt2 changes only initialization, prospectively committed before its runs.
All four cases use the same saved departure-pattern shape b=source_n/source_nstar,
multiplied by each case's LOCAL nstar and normalized per depth to its O total.
Thus common fractions are no longer asserted; common b-pattern SHAPE up to the
per-depth conservation factor is checked instead. No J/rates or physical
departure solution is imported/admitted. Source rates, atom, columns, controls,
Ng, cap and600second per-case budgets remain unchanged. Attempt1 is retained.
Source-domain tests for the new positive seed and invariance pass;8matched-input/
H/cap tests now pass1.85s. No active old run was edited or interrupted.

A third, separately labeled MALI population-operator control is predeclared in
`config/experiments/rya1281_solar_matched_mali.json`, with its own runner/helper.
It leaves the original numerical helper byte-for-byte intact and retains ALL
physical coefficients/contexts, seed policy, tolerances and budgets. After the
fixed-population J/cap field converges, an explicit `lambdaIterate=False` formal
refresh constructs LW's MALI Gamma; J/cap must remain within1e-8. Exact outer
J/cap changes are measured across the whole step. That Gamma is used ONLY for
the iterative population update. Physical admission still rebuilds off-diagonal
C+transition-R independently, with fresh radiation, positivity and conservation.
An accelerated operator is never mislabeled as the physical rate matrix.

This is an experimental preconditioner comparison, not a changed rate equation,
physical-model admission or numerical default. Completed legacy baseline/gravity
states provide same-physics endpoint controls. MALI dry Context passes. Existing
92numerical/iteration tests pass2.60s, including actual MALI versus physical-matrix
distinctions and analytic fixtures. Existing library warnings remain visible.
Full O-source admission remains required before any MALI result is compared to
the solar observations. Owned solves stay serial; no process was killed.

The installed0.16.1method signature was inspected directly:
`formal_sol_gamma_matrices(...lambdaIterate=False...)`. Its docstring identifies
`True` as setting the approximate Lambda term to zero, and the installed
`iterate_ctx.py` calls the default formal update before `stat_equil`.
This establishes the API/operator meaning, not a source rebuild attestation or
independent physical validation. Original compiled backend receipts are retained.
No numerical/default helper is silently replaced.


## Completed solver-admission diagnostic milestone

The four original common-fraction cases completed their bounded attempts. Baseline
passes in 3 outer iterations (16.8594 s), fresh physical residual
3.9562315360044704e-7; gravity passes in 176 (265.3007 s), residual
3.450558193900417e-7. Temperature and combined solar time out at 600 s; no final
spectrum or rate admission is claimed for them. The fresh baseline agrees with
the prior model-H spectrum within 1.05286e-8 normalized flux and 1.62180e-9
relative continuum, below the prospective 1e-6 numerical parity screens.

The common-departure-pattern restart does not solve the convergence problem:
baseline passes, all other cases time out under the original budget. These are
initialization controls, not independent physical results or a faster validated
solver. Every dry/run attempt remains preserved. The subsequent homogeneous
MALI control also times out in all four cases. Its baseline physical residual
stays near 0.005738 even when population updates become small. No MALI spectrum
is admitted. The synthetic numerical tests passing does not establish validity
for this source atom and custom collision/cap setup.

A bounded audit freezes the already qualified common-fraction baseline populations,
converges the cap/radiation, then samples direct/MALI/direct operators without any
population update. Physical residual remains 3.956231546503086e-7 in all three.
For direct iteration, (Gamma n - R n)/gross differs by at most
5.991892822873286e-15; for MALI it differs by 0.005703133226556716. The return to
direct iteration reproduces its original values. The maximum occurs at zero-based
level 40, depth 54; neighboring deep rows show the same pattern. Thus this
homogeneous approximate operator does not preserve the physical fixed point in
this held configuration. This is an observed interface/algorithm failure, not
proof of an upstream Lightweaver defect or a source-to-binary rebuild attestation.
Source inspection: LwMiddleLayer.pyx formal_sol_gamma_matrices initializes Gamma
with collisions; SimdFullIterationTemplates.hpp accumulates physical transition
rates separately from the approximate operator; UpdatePopulations.cpp stat_eq_impl
solves homogeneous Gamma populations with the largest-population row replaced by
number conservation. Row equilibration leaves those homogeneous equations intact.

A separate prospective correction uses P delta = -R n with sum(delta)=N-sum(n),
where P is MALI and R is independently reconstructed physical C+transition-R.
The root is therefore defined by the physical rates, not homogeneous P. No inverse,
clipping, floors, damping or relaxed tolerance is introduced. Ng remains enabled;
invalid accelerated populations restore the positive unaccelerated candidate.
The two targeted tests discriminate a deliberately wrong approximate stationary
state from the physical root, and recover/conserve the independent two-level
solution at two depths. They pass. The actual full-oxygen baseline is a distinct
required qualification; mathematical correctness of the correction alone does
not establish convergence or scientific admission.

91 focused source-column/defect/numerical tests pass (666 existing visible library
warnings). Earlier 92 numerical/iteration tests are complementary, not an
independent full-oxygen qualification. No adopted abundance, grid, continuum,
reference, source collision convention or production default changes. Wavelength
comparisons remain AIR with the existing single canonical conversion; no new
observation shift or telluric operation is applied.

Next: qualify the physical-residual update on baseline, compare to the archived
baseline and gravity endpoints, then complete all four cases before running the
predeclared factorial shape analysis. If the correction fails its original budget,
characterize its conditioning and acceleration before attempting more physical
columns. Only a complete admitted matrix can answer whether canonical solar T/g
reduces the observed O777 strength discrepancy. Broader xi/T/g and CNO/Fe, full3D
NLTE, source/EOS and independent reference validation remain open.

The full-oxygen defect-correction baseline also exhausts its original 600 s
budget (403 completed outer iterations), with last population change 1.72508e-4,
physical residual 9.00823e-5, J change 2.04845e-5. No final physical admission or
spectrum. This establishes that the mathematically root-preserving correction
is not yet a qualified acceleration for this configuration. Its failure seed is
initialization only. Next work is conditioning/acceleration diagnosis on a fixed
baseline, not additional solar parameter solves with this control. The completed
milestone is solver-admission diagnosis and failure preservation; the physical
canonical-solar matrix remains incomplete.

Handoff checks: generated target mirror and v238 register structure pass; staged
whitespace check passes. Generator guard retains exactly the same 2436 diagnostic
strings as the prior O-shape baseline (zero added/removed); its global exit remains
nonzero. Earlier global library/ledger/tracker debts remain, not a green global
science verdict. Raw archived outputs were byte-verified against their external
attempts. No active RT job remains.
