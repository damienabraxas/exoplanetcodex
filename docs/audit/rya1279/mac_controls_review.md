# Mac controlled checks: input sensitivity,10K sampling and adaptive-grid prototype

> 2026-10-27 correction: these are RH-convention experimental controls. The later [rate-convention audit](rate_conventions_review.md) identifies unresolved CH/CH0 source semantics; native MULTI physics parity is not established. Original results remain preserved.

Ryan authorized moving this validation work to the Mac, exploring temperature,
microturbulence and other inputs, checking the air-wavelength convention, and
keeping the broader adaptive/full3D/all-band goal explicit. **All36 new RT/SE
calculations ran on the Mac**, serially with one thread/nice15. No new Sirius
science job, dependency installation, raw-data change or adopted abundance was
needed. Combined case wall time was1218.39s (20.3minutes); interactive development
and analysis time is separate. All36 met the joint population/J/physical-rate
criterion1e-6 and conservation gate1e-10. Largest terminal physical-rate residual
was1.88e-14. This is experimental1D work, not a full3D or production result.

## What was varied

| Input | Completed controls | Interpretation |
|---|---|---|
| Teff | Native5500/5750/6000K;5752/5762/5772/5782/5792K at logg4.438 |10K fine sampling uses full thermodynamic-field interpolation, then a new NLTE solve at every node; not a temperature offset or interpolation of output spectra |
| Microturbulence |0.5,0.9,1.0,1.1,1.5km/s |Line/SE xi at fixed native MARCS t01 structure;0.9/1.1 reproduces the orchestrator paired perturbation size |
| Gravity |Native logg4.0/4.5/5.0; solar4.438 interpolation |Wide diagnostic drift, not the Sun's parameter uncertainty |
| Trace oxygen |8.6/8.7/8.8 |Response probes, not fits or changed adopted values |
| Metal mixture |[M/H]-0.25/+0.25 native MARCS structures, trace O fixed8.70 |Low-metal source includes its native alpha+0.10 pattern; a mixture change, not an isolated Fe derivative |
| Bottom boundary |56/55/54 native5750K depths |Tests the layer selection required by the10000K collision-table ceiling |
| Macroturbulence |Canonical solar3.8km/s and±1km/s probes |Same installed iSpec post-effects kernel as synthesis; independently specified, not fitted |
| Rotation |Canonical vsini1.8km/s and±1km/s probes |Same kernel and its explicit limb-darkening0.6 convention |
| Instrumental width |R300000/500000 catalog interval |Gaussian approximation in iSpec; not a measured local FTS sinc LSF |
| Alignment |Zero,±0.5km/s, known solar gravitational redshift+0.633km/s |Model-side diagnostic shifts, not changes to frozen observed pixels |
| Grid representation |Linear fields, log Pe/rho, log Pe/rho plus positive layer thickness |Independently held-out native temperature and gravity nodes test interpolation |

Atomic gf/rates/model completeness, a complete weak-blend/molecular spectrum,
magnetic fields, granulation/horizontal transfer, a fresh atmosphere/EOS solution,
and a full observational covariance budget are **not** covered by this table.
This is a response grid, not the RYA-587 abundance-uncertainty budget.

## Main parameter results

All areas below are the same predeclared±0.5Angstrom windows used in the first
Kitt Peak check; they contain blends/pseudo-continuum absorption and are not
isolated oxygen-only total EWs. Intrinsic native-model areas in milliangstrom:

| Control |7771.944 |7774.166 |7775.388 |
|---|---:|---:|---:|
| Kitt Peak2005 observation |77.56|67.20|53.79|
|5750/4.5, xi1, O8.70 baseline |61.31|53.08|41.84|
| xi0.9 |60.96|52.78|41.62|
| xi1.1 |61.70|53.40|42.08|
| xi0.5 |59.88|51.87|40.95|
| xi1.5 |63.50|54.91|43.18|
|5500K (common bottom cut) |46.71|39.94|30.88|
|6000K (common bottom cut) |76.34|66.86|53.68|
| logg4.0 |69.55|61.03|49.19|
| logg5.0 |54.61|46.53|35.71|
| O8.60 |55.52|47.76|37.24|
| O8.80 |67.47|58.75|46.78|
| Lower-metal native mixture |61.64|53.20|41.70|
| Higher-metal native mixture |60.18|52.26|41.45|

The ordinary±0.1km/s xi change affects areas by less than1%; it cannot supply
the original21–22% missing absorption. Temperature has a much larger effect:
6000K nearly matches the areas, but canonical solar Teff is5772K. Choosing6000K
because it fits would replace a known input with a compensating error. Gravity,
abundance and mixture probes likewise describe sensitivity, not new best values.

Canonical broadening changes the first-line minimum from0.7146 to0.7593 but
changes its window area only from61.315 to61.284milliangstrom. The known redshift
control improves alignment/RMS while leaving areas effectively unchanged. A
visually close core or improved RMS is not proof of correct absorption physics.

The full-depth and one/two-layer-cut native areas differ by at most
0.0000365milliangstrom. That supports these particular lower-boundary selections;
it is not a license to truncate arbitrary atmospheres. The6000K full structure
reaches10161.4K, outside supplied collision tables; the retained55-depth structure
stays within domain. No rate clamping or extrapolation was used.

##10K resolution: an independently tested grid-representation improvement

The existing iSpec interpolation primitive was inspected. Its packaged atmosphere
output does not retain density needed by this Bride interface, and its general
interpolator can choose a nearest model outside a hull. The experimental adapter
instead retains the native MARCS density/pressure/geometry, requires common
tauRoss layers and composition, and **rejects out-of-hull requests**. It uses
SciPy LinearNDInterpolator without creating an EOS or a new atmosphere-equilibrium
solution. Source files, interpolation simplices/weights and transforms are logged.

Holding out the5750/4.5 atmosphere from5500/6000 models tests temperature
interpolation; holding it out from logg4.0/5.0 models tests gravity interpolation.
Neither test uses Kitt Peak as an optimization target.

| Interpolation coordinates | Temperature holdout area error | Gravity holdout area error |
|---|---:|---:|
| Linear thermodynamic fields/heights |7.38–8.67% low|48.26–52.73% low|
| Log Pe/rho, linear heights |0.78–0.93% low|31.01–34.09% low|
| Log Pe/rho and positive geometric layer thickness |0.28–0.35% low|0.18–0.19% high|

The gravity holdout's temperature error is only1.25K in the diagnostic native
logtau500[-3,0] window, yet direct-height interpolation misstates layer thickness
by as much as77.1%. The positive-layer transform reduces that maximum to1.84%.
Thus a close temperature profile alone cannot validate an interpolated atmosphere.
The prototype preserves positive layers and monotone geometry instead of making
the spectrum agree by altering abundance. Log-pressure interpolation is motivated
by the held native MARCS interpolator; density and layer-coordinate extensions
are experimental and explicitly tested here. No Fortran source was vendored.

The improved10K sequence, with canonical solar logg4.438, xi1 and fixed O8.70,
gives these broadened window areas:

| Teff K |7771.944 |7774.166 |7775.388 |
|---|---:|---:|---:|
|5752|62.344|54.081|42.770|
|5762|62.943|54.626|43.233|
|5772|63.542|55.172|43.698|
|5782|64.142|55.719|44.164|
|5792|64.743|56.268|44.632|

Each10K step changes areas by about0.60/0.55/0.47milliangstrom. The tested native
holdout errors are smaller than that step signal. This supports resolving these
responses, **not universal10K physical accuracy** or validation of every pressure,
EOS, composition,3D or spectral domain. More native holdouts/refinement are needed.

The5772/4.438 experimental spectrum still falls18.07/17.90/18.77% below Kitt Peak
in these windows. Fixing the known solar input mismatch and the interpolation
representation does not remove the original model/observation discrepancy.

## Air wavelengths, Earth observations and canonical holding handling

The selected file is the same frozen **solar_kpno_kurucz2005_corrected** product
as the earlier check; its hash is verified before analysis. Its header explicitly
declares **vacuum nm including gravitational redshift**. The shared reader converts
once to **air Angstrom**, which is the comparison/plot coordinate. RT is evaluated
at vacuum nm using the canonical conversion in the other direction. At the7774.166
air anchor, the air-to-vacuum displacement is2.13941Angstrom; a missing or double
conversion is independently tested. The solar gravitational shift is only
0.01641Angstrom there and is a separate model-side control.

The raw model Fnu is stellar **surface** flux, not absolute irradiance at Earth.
This comparison uses continuum-normalized flux, so distance dilution cancels.
An absolute Earth-irradiance comparison would additionally need the geometric
dilution and spectral-density/observation-response contract. It was not performed.
Air here means the published standard-air wavelength convention, not a simulation
of local observatory weather or a reason to reintroduce corrected tellurics.

The same owners as ORCH are called: telluric_policy.gate_holding/holding_basis/
applied_state, band_policy.check_intake, and the normalization declaration check.
The holding passes, has corrected/applied state, and the red-optical policy says
to trust the supplied atlas continuum. No second correction, continuum placement,
pixel masking based on bad model residuals, or observed renormalization occurs.
The diagnostic scan reports normalized flux, median local p95=0.9963525 and
slope=-0.001016 across its22Angstrom padded span. Those numbers are not permission
to shift the continuum; weak blends/pseudo-continuum and a missing residual-error
budget remain. The canonical gate's legacy instrument-level reason string is
preserved verbatim in analysis.json alongside the holding-specific corrected state.

## What to do next

1. Audit the independent MULTI comparison at identical atom/atmosphere/rate/EOS/
   background settings. A carefully labeled fixed-population replay can isolate
   formal-transfer response from SE/population differences; b coefficients from
   different LTE normalizations must not be silently treated as identical physics.
2. Add qualified weak blends/background-opacity controls and an independently
   constrained local LSF/velocity check. Keep known solar parameters fixed; do not
   fit away a weak NLTE response by increasing temperature or abundance.
3. Validate the structure-preserving adapter on more held-out native points and
   refine with new physical calculations where its error is too large. The existing
  36 newly solved1D nodes are a small experimental spectrum/population grid,
   not a freshly solved atmosphere grid or a full3D result.
4. Select a stellar benchmark with a Gaia identity, independent stellar parameters
   **and oxygen abundance**, matched treatment/uncertainty and actual777nm spectra.
   Gaia identity or Fe/H alone is insufficient. Prefer qualified existing holdings;
   freeze inputs before comparing so parameter degeneracy cannot hide model error.
5. Follow the [adaptive/full3D/all-band roadmap](../../decisions/bride_adaptive_direction_rya1281.md).
   Source adapters, error-controlled refinement and backend capability contracts
   let us build beyond existing components. True horizontal3D transport, complete
   band opacity/chemistry and new atmosphere-generation physics remain distinct
   validation/compute tasks, not claims made by this1D experiment.

## Reproduction, artifacts and checks

- [All diagnostics and policy receipts](../../../data/results/rya1281/mac_parameter_checks_v1/analysis.json)
- [Experimental36-node grid index](../../../data/results/rya1281/mac_parameter_checks_v1/grid_index.json)
- [Compact native/broadened profile grid](../../../data/results/rya1281/mac_parameter_checks_v1/profile_grid.npz)
- [Fine-temperature responses](../../../data/output/rya1281/mac_parameter_checks_v1/fine_temperature_response.png)
- [Near-solar profiles](../../../data/output/rya1281/mac_parameter_checks_v1/fine_solar_profiles.png)
- [Parameter-drift profiles](../../../data/output/rya1281/mac_parameter_checks_v1/parameter_responses.png)
- [Broadening/alignment](../../../data/output/rya1281/mac_parameter_checks_v1/broadening_alignment.png)
- [Interpolation holdouts](../../../data/output/rya1281/mac_parameter_checks_v1/interpolation_holdouts.png)

Five source summaries and36 unchanged manifests are archived under the result
directory. Full original histories, absolute flux and populations remain in the
registered experiments.bride_mac_grid store; the index records their relative
paths and hashes.46 result artifacts have actual producers/invocations registered
in GENERATORS.yaml. Use fresh destinations; prior outputs are never overwritten.
The index intentionally retains the poor interpolation candidates as negative
controls; its explicit HOLD capability contract means it is not an admitted
inference/training grid. The interpolation API now requires an explicit coordinate
method rather than silently defaulting to the poor linear-height candidate.

RT/SE Python3.12.13: LW0.16.1, numpy2.5.2, scipy1.18.1, numba0.67.0. Project-side
policy/iSpec analysis records its separate environment and installed effects-source
hash. No new dependency was installed. Each phase preserves its generator commit,
plan fingerprint, source atmosphere/atom hashes, domain checks and convergence.
New geometry/pressure controls were chosen using native holdouts, not observed fit.

Local verification:163 LW/regression tests pass with666 known deprecation warnings;
14 project policy/effects/input checks pass with4 installed-iSpec warnings. Existing
holding-route and global-generator debt are not waived by those focused passes.
The Mac native baseline differs from the prior Sirius normalized flux by at most
4.21e-5 with different angular/frequency settings; this is a consistency check,
not an identical-settings backend-accuracy claim.

No production model bug is assigned from experimental interpolation failures;
the failed approximations are preserved as negative controls. Atom/reference
lineage, complete opacity and full uncertainty admission remain open. No raw/FITS,
frozen gold, published/adopted values, treatment token, merge or publication changed.
Refinement debt remains informational and unresolved. This push uses a scoped
[skip ci] marker to keep work off the Sirius self-hosted runner at Ryan's request;
it does not claim a green CI run or change workflow policy.

Attribution: Kitt Peak FTS data were produced by NSF/NOAO; the selected residual
irradiance atlas is credited to Kurucz2005. Lightweaver and the installed iSpec
effects retain their authorship/license provenance; no external solver source
or raw atom/grid file was vendored into this change.
