# RYA-1281: numerical refinement at thermal-axis endpoints

Predeclared 2026-11-04, before new RT. Scientific HOLD. Continue the accepted
six-case thermal experiment without changing its source plan or physical inputs.
Run thermal_T5572, thermal_T5822, thermal_g4338 and thermal_g4538, in that order,
serially on the Mac using the already committed run_bride_thermal_sweep.py
--refine-numerics. Reuse each exact producer certificate; original results stay
immutable. Rays3->5, qscale1->2, flux rays5->7 only; retain original1500iteration/
600second limits and all residual/conservation tolerances. No retries with relaxed
budgets or favorable case selection. Retain any refusal/failure.

Compare each refined run with its exact unrefined physical counterpart, including
complete source/binding/controls checks, saved independent R+C residual and
conservation, raw spectra, three O-aperture areas and maximum profile differences.
The solar nominal refinement is an existing control, not a new solve. This checks
combined sampling sensitivity at these endpoints, not individual error sources,
depth convergence, atmosphere interpolation or physical correctness. Numerical
differences will be reported directly; no new accuracy threshold invented after
seeing results. Continuum/observation choices and adopted solar parameters stay
unchanged. Full EOS, independent O reference/atom and continuum gates stay open.

New immutable Mac store root: rya1281-thermal-extremes-refined-v1. Full R/C/J
archives remain there if over the repository size ceiling, with verified receipts.
Base main f8c723a85f905425b1f7ae189f15a8763998ffc0; stacked base
4635c7421a7dee8671ffefe05d771c68da1cd24e; initial HEAD
d61e4ab9b8034d337b2bf2564688a4da4e46094e. No merge/production/publication.

## Conditional initialization improvement, declared before the first outcome

The legacy runner starts every thermal refinement from the old solar control,
including off-solar cases. If this bounded replay exhausts its budget, retain
that complete failed attempt and stop the legacy serial batch. The next bounded
attempt may instead initialize each refinement from its own already converged,
exact physical baseline. Validate case/source/certificate/local fields and
original controls, saved population/rate/history hashes, positive nuclei totals
and exact LTE state before using that guess. Re-solve with all original joint,
fresh-rate, conservation, iteration and time limits unchanged. This improves an
initial guess, not the atmosphere, atom or fitted parameters. Record the alternate
initialization and code lineage explicitly; do not call the failed replay passed.
Use a separate immutable root rya1281-thermal-extremes-matched-v1 and keep the
legacy attempt root intact. No increase of the600-second per-attempt budget.
The manifest's legacy initial_populations/report source_populations receipts
retain the original solar lineage. The new initialization receipt identifies
the actual matched starting populations. The solver also saves the actual
initial array in populations.npz; analysis requires exact array identity with
that case's archived converged baseline, not just a matching label.

## Independent light analysis while the serial RT batch runs

Use only the already archived thermal/xi profile products to compare finite
temperature5672->5822K, gravity4.338->4.538 and xi0.9->1.1km/s response directions.
Check exact shared nominal arrays, observation, broadening and source receipts.
Compare unit Euclidean directions over the union of the same three AIR O777
apertures, retaining both broadening/RV variants. Report signed cosines and
singular values of the unit-column profile matrix, alongside existing area
directions. This is unweighted geometry on a fixed grid: no noise model,
covariance, inversion, uncertainty estimate or observational identifiability
claim. No chosen favorable variant, new RT, continuum correction or fit. It
asks whether line shape carries distinctions lost by integrating three areas.

## Execution contract and replay

Mac science interpreter: /Users/ryanschmitt/bride_rya1118/.venv/bin/python,
Lightweaver0.16.1/NumPy2.5.2, one context thread. No installs, Sirius science or
CI load. Source host profile config/hosts/mac.env. The same-case runner commit
28aa7c30 precedes all four restarted solves. Coupled-control, context builder,
opacity/pressure/Stark kernels, certificates, atom and collision sources remain
unchanged. The initialization loader is separately hashed. No automatic fallback
on a rejected initial state and no relaxed convergence/iteration/time budget.

Actual invocation below; existing output paths refuse reuse. A reproduction must
choose a new immutable output directory, preserving these originals.

```sh
source config/hosts/mac.env
BRIDE_RT_PY=/Users/ryanschmitt/bride_rya1118/.venv/bin/python
BRIDE_CHECK_STORE=/Users/ryanschmitt/codex/bride_mac_checks
for thermal_endpoint in thermal_T5572 thermal_T5822 thermal_g4338 thermal_g4538; do
  "$BRIDE_RT_PY" scripts/run_bride_thermal_sweep.py \
    --plan config/experiments/rya1281_mac_thermal_slice.json \
    --case "$thermal_endpoint" --refine-numerics \
    --warm-start-run "$BRIDE_CHECK_STORE/rya1281-thermal-coupled-v1/$thermal_endpoint" \
    --coefficients "$BRIDE_CHECK_STORE/rya1281-thermal-coefficients-v1" \
    --case-directory "$BRIDE_CHECK_STORE/rya1281-solar-bound-v1/solar5772_bound_Fe_sampled" \
    --lines "$BRIDE_CHECK_STORE/rya1281-trace-fe-lines-v2/lines.json" \
    --pressure-source "$BRIDE_CHECK_STORE/rya1281-fe-pressure-source-v2/report.json" \
    --stark-source "$BRIDE_CHECK_STORE/rya1281-fe-stark-source-v2/report.json" \
    --output "$BRIDE_CHECK_STORE/rya1281-thermal-extremes-matched-v1/$thermal_endpoint" \
    > "/private/tmp/rya1281-matched-$thermal_endpoint.log" 2>&1 || exit $?
done
```

The first cold-start command used the same source arguments without
--warm-start-run and output rya1281-thermal-extremes-refined-v1/thermal_T5572;
exit1, fixed-population cap600s budget refusal after139 completed outer rows.
There is no accepted flux/report for that attempt. Its original manifest,
history and exception log are retained, not replaced by a successful label.
The remaining original-solar-start cases were not attempted; all four restarted
cases belong to the separately declared matched-initialization batch.

Analysis commands (use fresh output names for replay):

```sh
"$BRIDE_RT_PY" scripts/analyze_bride_thermal_refinement.py \
  --plan config/experiments/rya1281_mac_thermal_slice.json \
  --failed-attempt-log /private/tmp/rya1281-extreme-thermal_T5572.log \
  --failed-refinements thermal_g4338 --failure-log-dir /private/tmp \
  --output "$BRIDE_CHECK_STORE/rya1281-thermal-extremes-review-v1"
"$BRIDE_RT_PY" scripts/analyze_bride_profile_geometry.py \
  --output "$BRIDE_CHECK_STORE/rya1281-profile-geometry-v2"
```

The archived log can replace the temporary log in the first command. Full R/C/J
arrays resolve through experiments.bride_thermal_extreme_refinements; coefficients
through experiments.bride_thermal_coefficients. Original timeout root resolves
through experiments.bride_thermal_extreme_legacy. All three paths verified on Mac.
Profile geometryv1 remains intact; v2 corrects plot autoscaling and wavelength
tick formatting only, with byte-identical direction CSVs.

## Endpoint results and retained evidence

All four declared matched endpoints attempted. Three converge; one exhausts
the original budget. No failed point is dropped from the report or replaced
with a favorable spectrum. Existing nominal refinement is reused, not rerun.

| Case | Status | Iterations | Seconds | Independent fresh rate | 7771/7774/7775 area delta mA |
|---|---|---:|---:|---:|---|
| T5572 | converged |98|390.311|4.091227e-7|+0.0048073/+0.0045403/+0.0033531|
| T5822 | converged |178|459.824|3.435123e-7|-0.0025101/+0.0007394/+0.0024458|
| g4.338 | budget exhausted, NOT ACCEPTED |223|600s limit|not accepted|unqualified|
| g4.538 | converged |152|481.165|4.876575e-7|-0.0003267/+0.0017252/+0.0024403|
| nominal, earlier control |converged, reused|153|484.133, earlier|3.572884e-7|-0.0003975/+0.0019795/+0.0028804|

Three completed new replays total1331.300s (22.19min); two separate budget
refusals retained (legacy cold and matched g4.338). Completed-new process peak
RSS max1366327296bytes, including imports/setup; one context thread, host loads
in each report. Not full-grid/3D throughput. Upper-gravity standalone command
exit0; matched batch exit1 at the lower-gravity case as documented above.

Eight independent saved-R+C reconstructions pass for the four completed pairs
(three new refinements plus reused nominal, each with its old baseline).
New refined conservation errors<=1.111e-15. Exact saved initial-population arrays
match their own physical baselines. Original joint1e-6, conservation1e-10,
inner1e-8,1500iteration/600s controls are unchanged. Maximum area difference over
completed pairs0.0048073mA; maximum normalized profile difference9.51114e-5.
Sampling effects in these completed cases are much smaller than the several-mA
observational/reference discrepancies. Lower-gravity accuracy is unbounded by
this failed refinement, and individual quadrature errors, depth/cap/source,
interpolation and physical accuracy remain unqualified.

40 derived result artifacts plus2figures are byte-checked copies under
data/results/rya1281/thermal_extremes_v1 and data/output/rya1281/thermal_extremes_v1.
Actual generators/wrappers are registered in GENERATORS.yaml. Largest committed
copy1307211bytes; no size-hook bypass or renamed binary archive. Full new R/C/J
arrays remain in the registered matched root, with these archived hash receipts:

| Relative external file | Bytes | SHA256 |
|---|---:|---|
|thermal_T5572/rates.npz|6399558|c3dc2e25742e64da4b022e6931d0534ced4201895f9dcff43adda7a506d7e7d9|
|thermal_T5822/rates.npz|6382114|46b9071692962608f08b63ed43ba9fdd257db7e81d147a9d164af3e4f072d300|
|thermal_g4538/rates.npz|6383161|b2f9a88fd48bf71b4c9c89e90a99c97040fa0e13fec8f4672f3259b52cd4965a|

[Completed-pair numerical plot](../../../data/output/rya1281/thermal_extremes_v1/analysis/endpoint_refinement.png),
[finite profile geometry](../../../data/output/rya1281/thermal_extremes_v1/profile_geometry/profile_geometry.png),
[complete statuses/rates/receipts](../../../data/results/rya1281/thermal_extremes_v1/analysis/report.json).
Both plots visually inspected; no changed observation or chosen favorable case.

## Finite profile geometry results

Existing three-area directions have signed cosines T/g=-0.999386,
T/xi=+0.997084 and g/xi=-0.993800. In the full sampled O-aperture profiles after
the same canonical broadening, those become -0.988107/+0.586113/-0.478944.
The gravitational-shift variant gives -0.988115/+0.586631/-0.479516.
Shared nominal profiles and observations match exactly between source archives.
All600 selected pixels use the same0.005A AIR grid; no interpolation or new RT.

Smallest/largest singular-value ratios for unit response columns are0.000236758
for areas and0.0398744/0.0398403 for the two profile variants. This is geometric
conditioning after discarding physical amplitudes/parameter units, without
noise or nuisance columns, not a covariance or an observational rank decision.
Microturbulence has a different core/wing pattern; temperature/gravity remain
similar and opposed. Macroturbulence, rotation/LSF, abundance, continuum and
other missing physics may confound these distinctions. The source intervals
are finite150K/0.2dex/0.2km/s, not local derivatives or parameter uncertainties.
No input is chosen to fit Kitt or to minimize a residual.

Kitt remains AIR for comparison, internal RT vacuum; original observations,
shipped baseline and earlier ORCH synthesis/IAG continuum legs are preserved.
This phase compares intrinsic numerical replays and existing conditioned
profiles. It does not newly refit continuum/tellurics on the four refined
spectra, choose an anchor/reference or apply a second telluric correction.
Independent reference/continuum/LSF and full EOS/atom/opacity gates remain open.

## Software validation and remaining scope

50 targeted tests pass, exit0,6.82s on the held Mac interpreter:

```sh
"$BRIDE_RT_PY" -m pytest -q \
  tests/test_bride_profile_geometry_rya1281.py \
  tests/test_bride_thermal_initial_rya1281.py \
  tests/test_bride_thermal_refinement_rya1281.py \
  tests/test_bride_thermal_sweep_rya1281.py \
  tests/test_bride_xi_curve_rya1281.py tests/test_bride_xi_sweep_rya1281.py \
  tests/test_bride_solar_control_rya1281.py tests/test_bride_fe_stark_rya1281.py \
  tests/test_bride_fe_pressure_rya1281.py tests/test_bride_trace_fe_lines_rya1281.py \
  tests/test_bride_top_cap_rya1281.py tests/test_bride_observables_rya1281.py \
  tests/test_bride_reference_rya1279.py tests/test_rya1230_local_continuum.py
```

Controls reject wrong case/source/LTE state/nuclei/relaxed budgets, altered rate
operators and populations, disconnected zero-rate operators, false timeout
acceptance, and degenerate/nonfinite geometry inputs. New unit fixtures use
archived baseline/disposable data; no dependency on the transient timeout path.
Actual saved-initial-array equality and independent R+C/target-nuclei checks
are separate run validations, not extra stellar tests counted into the50.

Library audit exit1:193documents/127accounted/66unindexed/12declared off-Mac.
Earlier session had191/64; no new document is used as scientific evidence here.
Existing RYA-1277 indexing debt remains. Refinement debt report exit0:
10not-fired/5TBD/2discharged, unchanged; no element debt or product gate cleared.
Global path inventory has46ok/15broken entries, including other services/old
holdings. Required experiment resources resolve and exist under the Mac profile.
No resource fallback, new install/download, source-spectrum/grid/atom/gold write,
adopted abundance/product change, production, merge or publication.

Bug-ledger disposition: none discovered. The retained time-budget refusal is an
expected enforced limit; an incomplete attempt is never accepted as convergence.
The matched initializer is a bounded numerical starting-state option, not a
physical repair or a guarantee of a unique fixed point.

The matched serial batch exited1 at thermal_g4338:600s cap-radiation deadline,
223 completed outer rows, final joint maximum1.326199e-6 against1e-6 required.
Its small in-loop SE roundoff is not a fresh post-RT acceptance check. No accepted
refined spectrum/report exists for this case. The last thermal_g4538 command was
then launched separately with exactly the same declared arguments/limits; the
failed lower-gravity case was not retried with a larger budget or looser tolerance.
The history records zero invalid-update rejections in223rows. Its recurring
jumps are consistent with the declared Ng period, but this is an inference:
positive_stat_equil discards the returned IterationUpdate.ngAccelerated flag
and records only invalid/negative-update rejection. No paired acceleration
intervention was run here. Next numerical work should retain the actual Ng
flags and compare predeclared acceleration controls with the same physical
operator/tolerances, before attributing a cause or changing the production solver.

Next: qualify source-backed macroturbulence/rotation/FTS-LSF and continuum
nuisance responses before treating profile shapes as measured constraints;
resolve independent matched O reference/atom/EOS/opacity discrepancies. Broader
xi numerical extremes, composition and declared interactions remain owed.
Additional O diagnostics and C/N come before Fe abundance and a reference star.
Current static1D O777 VIS checks do not establish full3D, molecules or all bands.

Final QA: required resources and MIRROR targets, register structure/version pointer and whitespace pass. Result-generator guard retains exactly2436 pre-existing diagnostics, zero added/removed after intended artifacts are staged. Its initial pre-staging run correctly reported40 untracked archive entries; staging the intended40files resolves those transient notices without suppressing a rule. Scientific self-review checks case/source/frame identity, unchanged operator/controls, positive target-nuclei conservation, independent rate reconstruction, immutable archives, explicit failures and non-statistical interpretation. No material wrong-data or provenance ambiguity found; numerical/physics HOLDs remain explicit.
