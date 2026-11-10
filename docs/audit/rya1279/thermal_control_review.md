# RYA-1281: source-matched temperature/gravity slice

Scientific HOLD. This advances the broader parameter program with six declared
physical cases: Teff5572/5672/5772/5822 at g4.438, and g4.338/4.538 at5772.
The nominal point is shared. Xi1/composition/O8.70 stay fixed. One separate
nominal numerical refinement changes angular, wavelength and flux-ray sampling.
These interpolated static1D O777 VIS diagnostics do not establish an abundance,
a fresh atmosphere-equilibrium solution, full3D or all-band admission.

## Source and warm-start contract

`rya1281_mac_thermal_slice.json` declares each axis and nearest native brackets.
Changing an undeclared axis, composition, O abundance, depth cut, collision
prescription or convergence budget refuses. Both previously refused hotter
candidates remain in the broader domain ledger; no source extrapolation or
row truncation introduced. Native-only provider/cap admission remains unchanged.

Fresh TS point evaluations: native identity plus six physical cases, baseline
and reporting builds. Seven exact identities/110input/schema/neighbor refusals
pass. Full56 rows/1601packet wavelengths, requested T/Pe/xi and serialized default
REAL output checked. Each complete certificate binds local ne/nH/composition,
thermal/geometry/velocity/column identity. Producer density is separate: maximum
producer/local rho discrepancy ranges0.414-1.383% over the physical cases.
Local density is retained, with the existing per-neutral Fe conversion; no raw
producer per-mass density substitution or fitted opacity factor.

Same full O atom/source MN-TAUTSP collisions,41 source FeI pressure lines and
Cowley Stark; atomic EOS/continuum/O-collision baseline retained, H2 omitted.
New contexts recompute local LTE populations, densities and pressure widths.
Warm starts preserve the checked nominal departure coefficients at the new LTE
state, then normalize to the new O total. The exact nominal shortcut preserves
its previous populations. Warm starts are subsequently solved, never accepted
as converged outputs. Positive/conserved projection and wrong-input tests pass.

Only the derived native CH negative-spline clipping count may vary with local T;
source oracle, directions and prescriptions must match. Author MN provenance
must match exactly. No collision multiplier or clipping rule fitted to Kitt.
Original coupled solver/helper remains unchanged. Source preparationd3121d07
and runner45af6341 precede the respective producer/O work.

## Results

| Case | Iterations | Seconds | Fresh independent rate residual | 7771 /7774 /7775 area mA |
|---|---:|---:|---:|---|
| thermal_T5572 | 273 | 277.038 | 4.460633e-07 | 57.53155/49.41209/38.18722 |
| thermal_T5672 | 246 | 356.511 | 4.312502e-07 | 64.70097/55.77166/43.41635 |
| thermal_nominal | 16 | 24.815 | 4.137014e-07 | 72.12826/62.41175/48.93311 |
| thermal_T5822 | 283 | 383.483 | 4.846015e-07 | 75.88629/65.79010/51.76184 |
| thermal_g4338 | 118 | 172.434 | 4.392779e-07 | 74.57871/64.69266/50.93255 |
| thermal_g4538 | 178 | 240.258 | 4.609377e-07 | 69.82493/60.26935/47.05611 |

The nominal raw flux/continuum arrays reproduce the previous xi nominal exactly.
Each new physical case uses original1500iterations/600s, joint1e-6, conservation
1e-10 and inner1e-8 controls. Fresh independent reconstruction from saved R+C
is required after final RT. Initial/final populations, nstar, rates/J and complete
histories are retained; in-loop roundoff is not the fresh residual.

Nominal refinement changes rays3->5, line quadrature scale1->2 and flux rays5->7,
with the same1500iteration/600s, joint1e-6/inner1e-8 controls. Converges
153iterations/484.133s, fresh independently reconstructed rate residual
3.5728843e-07. Aperture changes are
-0.00039747/+0.00197945/+0.00288037mA; largest absolute0.002881mA.
Max normalized profile changes inside the three O apertures are
6.800e-5/2.852e-5/1.845e-5. This combined nominal-only refinement cannot apportion
individual angular/wavelength/flux errors or bound all parameter extremes. It is
far smaller in area than the several-mA nominal observational gap.

All seven independent saved-R+C reconstructions pass. Total replay1938.67s
(32.31min), including refinement484.13s; process peak1007468544bytes for refinement,
including imports/setup. One context thread; host loads in reports. Not full-grid/
3D throughput or a physical error budget.

Finite three-aperture response vectors over Teff5672->5822, g4.338->4.538 and
xi0.9->1.1 have signed direction cosines T/g=-0.999386, T/xi=+0.997084 and
g/xi=-0.993800. Prior xi source/kernel/control/thermal compatibility is checked;
these are unweighted finite-direction comparisons, not Pearson correlations,
derivatives, inverse covariance or proof of exact degeneracy. Similar strengths
make parameter discrimination fragile; full profiles/independent constraints
and additional diagnostics remain necessary. No parameter selected to fit Kitt.

[Response plot](../../../data/output/rya1281/thermal_slice_v1/analysis/thermal_response.png),
[spectral comparison](../../../data/output/rya1281/thermal_slice_v1/analysis/thermal_control.png),
[numerical differences](../../../data/output/rya1281/thermal_slice_v1/analysis/numerical_refinement.png).


## Interpretation and observational conditioning

Temperature/gravity changes are causal model-response controls, not revised solar
parameters. Nearest-cell simplex interpolation can change source cells across
the sweep; responses include this declared grid/operator choice. Numerical ray/
wavelength refinement does not bound atmosphere interpolation, EOS or atomic
physics uncertainty. No covariance/abundance error budget or preferred parameter
is inferred from these points.

Kitt remains AIR, internal RT vacuum, with the original observation and both
canonical broadening/gravitational-shift and zero-RV variants retained. Fixed
solar post-effects isolate atmosphere response, not a self-consistent stellar
mass/radius/velocity variation. ORCH synthesis/IAG continuum controls remain
separate from shipped baseline, with full anchor/reference spread preserved;
no favorable continuum, second telluric correction, RV/LSF or abundance fit.

All six physical cases plus refinement have synthesis/IAG continuum controls.
Guided observed areas span74.57242-74.59816 /65.84765-65.98026 /
46.28040-46.92637mA. IAG81.73581/70.81716/57.70338 and p95-envelope
74.10889/63.73224/49.94709 unchanged. This reference spread and one-sided-anchor
qualification remain much larger than the nominal numerical refinement; no
favorable correction is chosen or claimed physically qualified.


## Validation and resume

45tests pass: five thermal-axis/bracket/source-domain/warm-total/refinement/
metadata controls plus40earlier source/backend/curve/provider/cap/observable/
reference/continuum regressions. Deliberate altered fields/constructed column
refuse even with a recomputed digest. Fresh parameter-domain checks still refuse
5872/5972 at976tables; no zero/extrapolation/depth-cut fallback. An exploratory
system-Python integration command refused unsupported NumPy; rerun in the held
project environment, no output/data changed. No new verified accepted science
or production bug; no fabricated bug-ledger row for intended refusals.

88result copies plus one external-archive receipt (89artifacts total), and17
figure copies are preserved under thermal_slice_v1 result/output roots. Actual generator/helper writers and commands are registered in
GENERATORS.yaml. Full producer arrays remain at experiments.bride_thermal_coefficients;
metadata copies are not standalone packets. The refined rates.npz is6382021bytes
(6.1MiB), above the repository5MiB ceiling. Its complete unchanged R/C/J archive
stays at experiments.bride_thermal_refined_rates; rates_external.json verifies
SHA256 944ca6614c00e616ffc82c85ab795a87bcbd08c4619f31157d1acdb509576316.
The local commit hook refused the archive copy, resolved by this registered
external receipt; numerical content/report hashes are unchanged. Registered Mac-store roots:
rya1281-thermal-source-v1, thermal-coefficients-v1, thermal-coupled-v1,
thermal-refined-v1, thermal-continuum-v1 and thermal-review-v2. Reviewv1 retained;
v2 adds aperture profile metrics and qualifies the direction interpretation.

Science/tests use held Bride Python/LW0.16.1; analysis uses existing project
Python/iSpec/pandas, no installs. An optional nice request was denied by sandbox;
analysis commands still completed at their existing priority. No priority-control
success claimed and no permission escalation needed for the successful analysis.
All owned jobs finished. No source spectrum/atom/grid/gold/adopted product changed.


Next: representative refinement at parameter extremes, additional gravity/
composition points and predeclared interactions. Resolve hotter-source coverage
explicitly. Independent matched1D O/atom lineage, full EOS/chemistry/opacity,
continuum/LSF and dimensionality gates remain open. C/N before Fe abundance or
another star. Mac only, no Sirius science/CI load, no adoption/production/merge/
publication/full3D/all-band claim.
