# RYA-1281: wider xi sweep and broader parameter-domain preflight

Scientific HOLD. Five fresh coupled O cases extend the same solar experiment to
xi0.5/0.75/1/1.25/1.5km/s. Together with the previous slice there are seven
unique xi points. All five new solves converge, and the repeated nominal flux
is array-identical. This is a wider1D O777 response test, not completion of the
requested multi-parameter/full-spectrum/3D program or a measured abundance.

## Inputs, convergence and source controls

`rya1281_mac_xi_wide.json` fixes solar5772/g4.438/O8.70, full56-depth thermal
fields/composition, experimental column/cap, bound FeI continuum,41 source FeI
pressure lines, Cowley Stark and original O collision conventions. H2 stays
omitted to retain the atomic EOS/continuum/O-collision baseline. New native
producer evaluations and exact certificates pass6identities/94refusals; actual
serialized xi checked, no changed-xi certificate relabel or flipped admission.
No physics kernel changed for this batch. Plan authored before producer calls;
656921de precedes all five O runs. Later analysis-only commits do not change
solver/kernel bytes. Main f8c723a8 refreshed unchanged; stacked draft #598/#597.

Each case starts from the same checked earlier nominal populations and performs
new coupled SE/radiation. Original1500iterations/600s, joint1e-6, conservation1e-10
and inner1e-8 controls retained. Separate saved radiative and collisional rates
independently reconstruct the final physical matrix/residual after fresh RT.

| xi km/s | Iterations | Seconds | Independent fresh rate residual |
|---|---:|---:|---:|
| 0.5 | 73 | 101.485 | 2.466418e-7 |
| 0.75 | 93 | 123.238 | 3.647904e-7 |
| 1.0 | 16 | 24.273 | 4.137014e-7 |
| 1.25 | 61 | 84.834 | 3.463197e-7 |
| 1.5 | 138 | 159.411 | 3.223193e-7 |

Total493.24s (8.22min), Mac serial1thread; maximum process peak750927872bytes
includes imports/setup. Host loads captured. No full3D/grid-scaling measurement.
Some intermediate Ng steps have elevated physical residuals; histories retained,
joint convergence prevents accepting those steps. No tolerance relaxed.

## Wider oxygen response

Intrinsic +/-0.5AIR A diagnostic aperture areas in mA, not isolated oxygen EWs:

| xi km/s | 7771.944 | 7774.166 | 7775.388 |
|---|---:|---:|---:|
| 0.5 | 70.38364 | 60.96249 | 47.86601 |
| 0.75 | 71.12707 | 61.58032 | 48.32156 |
| 1.0 | 72.12826 | 62.41175 | 48.93311 |
| 1.25 | 73.35446 | 63.42896 | 49.67894 |
| 1.5 | 74.77040 | 64.60199 | 50.53589 |

xi1.5-minus-nominal adds2.64214/2.19024/1.60278mA; xi0.5 subtracts
1.74462/1.44927/1.06709mA. Response is monotonic and clearly nonlinear across
this interval. Do not extrapolate a0.1km/s derivative to the full range or infer
a preferred solar xi from these fixedO8.70 comparisons. Thermal structure is
fixed; self-consistent atmosphere response remains a separate question.

Kitt remains AIR, internal RT vacuum. Both canonical broadening/gravitational
shift and zero-RV variants retained. Original shipped observation, telluric
state, abundances and all earlier spectra remain unchanged. Synthesis/IAG ORCH
continuum legs rerun for every case. Across xi0.5-1.5, guided observed areas span
approximately74.58716-74.58933 /65.93192-65.93911 /46.73492-46.79456mA;
IAG81.73581/70.81716/57.70338 and p95 envelope74.10889/63.73224/49.94709
unchanged. Reference/anchor qualification still dominates interpretation, no
favorable continuum or xi selected, no second telluric correction or LSF/RV fit.

[Seven-point response](../../../data/output/rya1281/xi_wide_v1/combined_v3/combined_xi_response.png)
and [wider spectral comparison](../../../data/output/rya1281/xi_wide_v1/analysis/xi_control.png).

## Adaptive-grid reconstruction diagnostic

Compatible-input checks require matching actual thermal fields/composition,
source kernels/line packets, collision receipts and solve controls before joining
the two slices. Duplicate nominal spectra must be array-identical. Negative
controls reject changed pressure source, one-ULP T, convergence tolerance or
fixed-population mode. Seven distinct points retained; no best-case filtering.

Using the coarse0.5/0.75/1/1.25/1.5 curve to reconstruct the known0.9/1.1 points:

- **Interpolating aperture areas:** max absolute discrepancy0.028993mA with
  linear interpolation and0.001157mA with shape-preserving cubic (PCHIP).
- **Separately interpolating emergent Fnu and continuum, then dividing:**
  max profile discrepancy across7763-7783A is0.0012244 linear and0.00010752
  PCHIP. Within the three O apertures, worst discrepancy is0.00019787 linear
  and0.00002272 PCHIP. These are distinct from interpolating integrated areas.

Thus good integrated-area agreement alone does not establish profile accuracy.
This is a known-point reconstruction check, not blind validation, a physical
error bound, a certified interpolated SE solution or adoption of either operator.
The raw-flux-ratio test distinguishes this from interpolating normalized spectra.
Representative numerical refinements at sweep extremes remain owed.

## Temperature/gravity preflight for the next axes

Read-only candidates declared in `rya1281_parameter_domain.json`. Actual held
collision tables span1000-10000K or2000-200000K. All56 candidate atmosphere
layers are checked against every table; source models/interpolator/atom hashed.
This is necessary temperature-domain coverage only, not EOS/opacity/physical
atmosphere/provider or RT admission.

| Candidate Teff K | Layer T min-max K | Disposition |
|---|---|---|
| 5572 | 3904.2384-9736.7484 | domain covered only |
| 5672 | 3983.2384-9833.1084 | domain covered only |
| 5772 | 4062.4452-9930.0128 | domain covered only |
| 5822 | 4100.8652-9980.4728 | domain covered only |
| 5872 | 4139.2852-10030.9328 | REFUSED:976tables uncovered |
| 5972 | 4216.1708-10131.8648 | REFUSED:976tables uncovered |

At both refused temperatures, uncovered source tables are611CE/39CI/302CH/24CH0.
No out-of-domain coefficients extrapolated, top/bottom rows cut or rates zeroed.
This is a held-source coverage gap, not a physical stellar-temperature ceiling.
At5772K, gravity candidates4.338/4.388/4.438/4.488/4.538 all pass this domain
check. In total9/11 candidate checks covered,2refused; duplicate nominal checks
belong to different axes, not different stars. No temperature/gravity O solve
has been run by this preflight. Full candidate fields/certificates still owed.

## Validation and resume

40tests pass: new compatible-physics and raw-flux-ratio controls plus existing
sweep/binding/pressure/Stark/trace/cap/observable/reference/continuum regressions.
Additional6producer identities/94runtime refusals and5independent rate checks
are separate evidence. One monitoring read briefly encountered an empty history
at case startup; it did not alter a job/result and subsequent reads report the
preparing state. No new verified accepted science/production bug; no bug-ledger
row fabricated for expected source-domain refusals.

71artifacts/14figures are byte-identical copies in xi_wide_v1 result/output roots.
Full producer coefficient arrays remain in registered external holdings; copied
certificate metadata is not a standalone packet. Actual generator/helper writers
and commands registered. External Mac roots:
`rya1281-xi-wide-source-v1`, `rya1281-xi-wide-coefficients-v1`,
`rya1281-xi-wide-coupled-v1`, `rya1281-xi-wide-review-v1`,
`rya1281-xi-wide-continuum-v1`, `rya1281-xi-combined-v3`,
`rya1281-parameter-domain-v1`. Combinedv1/v2 retained. v3 makes CLI declarations explicit for the static generator checker; all comparison metrics, curve CSV and figure remain identical to v2. Initial checker added two false flag diagnostics for the dynamic declarations; corrected before push.

Next: source-matched Teff/gravity inputs and representative numerical refinement
at extremes. Preserve both refused hotter candidates until collision-source
coverage is resolved; no silent contraction of the requested wider sweep.
Then composition and declared interactions alongside independent matched1D O,
full EOS/opacity/continuum/model-dimensionality gates. C/N before Fe abundance
or another star. Mac only, no Sirius science/CI load. Scientific HOLD; no
adopted abundance/production/merge/publication/full3D/all-band claim.
