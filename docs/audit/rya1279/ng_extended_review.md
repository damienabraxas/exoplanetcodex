# RYA-1281: predeclared longer enabled-Ng computation

Declared before RT,2026-11-04 local. One Mac thermal_g4338 replay with Ng enabled,
same verified starting population file/certificate/refined sampling/physical
operator, and native flags recorded. Increase only compute allowance600->1200s,
based on the measured two600s refusals. Preserve1500iterations, joint1e-6,
conservation1e-10, inner1e-8 and original source Ng parameters. No extra retry if
this bounded attempt fails; retain the failure. More time is not assumed to
guarantee convergence or physical accuracy. No fitted input or acceptance limit.

The original two600s attempts remain failed. Require an exact221-row shared
history with the earlier enabled control and real native call accounting. Only
after full convergence: independently reconstruct unscaled saved R+C, check O
nuclei/positive populations and actual initial-array identity, then compare
intrinsic O areas/profiles with the same unrefined physical baseline. No adopted
abundance, new continuum/telluric fit, production, merge or full3D/all-band claim.

Base HEAD4747ef72a5f4b2c6e522430743e057e03280f090; main refreshed unchanged
f8c723a8; stacked base4635c742, draftPR598. Same worktree/branch, no other edits.
Resource experiments.bride_ng_extended, new immutable rya1281-ng-extended-v1;
no source/old-result write. Original solver, recorder, context builder, collision,
pressure and opacity kernels remain unchanged.

## Implementation and reproduction

Declaration/runner commit48017bf2 precedes RT. Explicit --extended-budget is
restricted to the declared refined, matched-start, recorded-Ng case. The manifest
and report fingerprint the budget declaration and adapter. Unknown declaration
keys, wrong case/mode, noninteger or non-increased cost, undeclared extensions,
and altered residual/iteration controls refuse. No default allowance changes.

59 targeted tests pass in6.83s (exit0), including prior57 and2cost-only guards.
The native recorder's direct-backend arrays/history/retry tests remain included.
Actual RT source uses held Mac Python/LW0.16.1/NumPy2.5.2, one context thread;
no installs/downloads, Sirius science/CI or source-data changes. Existing
library193/127/66unindexed/12off-Mac and refinement10not-fired/5TBD debt remain.

Actual invocation, after source config/hosts/mac.env; these original directories
are immutable, so reproduction must choose a fresh registered output resource:

```sh
BRIDE_RT_PY=/Users/ryanschmitt/bride_rya1118/.venv/bin/python
BRIDE_CHECK_STORE=/Users/ryanschmitt/codex/bride_mac_checks
"$BRIDE_RT_PY" scripts/run_bride_thermal_sweep.py \
  --plan config/experiments/rya1281_mac_thermal_slice.json \
  --case thermal_g4338 --refine-numerics --ng-control record --extended-budget \
  --warm-start-run "$BRIDE_CHECK_STORE/rya1281-thermal-coupled-v1/thermal_g4338" \
  --coefficients "$BRIDE_CHECK_STORE/rya1281-thermal-coefficients-v1" \
  --case-directory "$BRIDE_CHECK_STORE/rya1281-solar-bound-v1/solar5772_bound_Fe_sampled" \
  --lines "$BRIDE_CHECK_STORE/rya1281-trace-fe-lines-v2/lines.json" \
  --pressure-source "$BRIDE_CHECK_STORE/rya1281-fe-pressure-source-v2/report.json" \
  --stark-source "$BRIDE_CHECK_STORE/rya1281-fe-stark-source-v2/report.json" \
  --output "$BRIDE_CHECK_STORE/rya1281-ng-extended-v1" \
  > /private/tmp/rya1281-ng-extended-run.log 2>&1
```

Analysis (no new RT; fresh output name required):

```sh
source config/hosts/mac.env
"$BRIDE_RT_PY" scripts/analyze_bride_ng_extended.py \
  --attempt-log /private/tmp/rya1281-ng-extended-run.log \
  --output "$BRIDE_CHECK_STORE/rya1281-ng-extended-review-v1"
```

The archived log may replace the temporary log. Analysis records exact history
and native-call prefix identities, full case/constructed-column/certificate and
original input provenance, and independently rebuilds the saved unscaled R+C
operator and target oxygen nuclei before numerical-spectrum metrics. A failed
attempt remains unaccepted and never receives fabricated flux/rates/area results.
Internal RT wavelengths remain vacuum; numerical spectral plots use AIR, with
original Kitt observation/shipped continuum and earlier ORCH controls unchanged.

Before inspecting any updated continuum result, retain both existing synthesis
and IAG reference legs on the newly converged spectrum, with the same comparator,
anchors, canonical post-effects and both RV variants. This is an analysis-only
stability check on the refined numerical spectrum, not a new continuum choice,
abundance/velocity/LSF fit or telluric correction. Preserve physical anchor and
independent-reference qualification gaps regardless of which residual is smaller.

## Results and interpretation

The single declared replay converges, exit0,233outer iterations,45actual native
Ng calls,630.712s measured setup/solve/output elapsed, peak process RSS1209335808
bytes. Recorder overhead0.04392s. One context thread, Mac host loads retained;
not full-grid or3D throughput. Entire221-row prior600s history and221native
calls match exactly. Twelve additional outer steps reach the original joint
criterion; neither the old timeout nor any accuracy criterion was relabelled.

Independent reconstruction of saved unscaled R+C gives residual3.5507904e-7,
oxygen-nuclei conservation error6.66134e-16; finite positive populations and
exact saved-initial-array identity with the physical baseline pass. Original
1500iterations, joint1e-6/conservation1e-10/inner1e-8 remain unchanged. This
closes this specific lower-gravity denser-sampling replay, not physical accuracy.

| AIR O aperture | Refined intrinsic area mA | Change from same unrefined physical case mA | Max normalized profile change |
|---|---:|---:|---:|
|7771.944|74.57846214|-0.00024548|7.55424e-5|
|7774.166|64.69511456|+0.00245097|3.15008e-5|
|7775.388|50.93607124|+0.00351838|2.24650e-5|

All four declared thermal/gravity endpoints now have a completed combined
angular/wavelength/flux-ray refinement at the original accuracy criteria, with
this lower-gravity point requiring the separately declared longer cost cap.
Existing nominal control retained. Maximum completed-pair area change across
these points remains0.0048073mA (cold endpoint), profile change9.51114e-5
(hot endpoint). No independent bound for individual numerical contributors,
all parameter points/xi extremes, depth/cap/source/interpolation or physics.

## Continuum stability, no favorable selection

Both existing continuum-reference analyses complete (exit0), retaining shipped
baseline/p95 envelope and both fixed RV variants. Reference code, comparison,
observations, nominal configuration and pre-clipping anchor geometry are identical
to their old unrefined-case controls; IAG source receipt also matches exactly.

| Reference leg | Guided observed7771/7774/7775 areas mA | Change from old same-case diagnostics mA |
|---|---|---|
|synthesis|74.57982009/65.89174905/46.82403801|+0.00000195/+0.00001226/-0.00000151|
|IAG|81.73581499/70.81716469/57.70338185|0/0/0|

Reference spread remains about7.156/4.925/10.879mA, far larger than numerical
refinement effects. This is diagnostic sensitivity, not an adopted error budget
or a choice of the apparently closer reference. First two synthesis and first
IAG pre-clipping selections remain entirely blue-side; recorded geometry is not
the final robust-fit anchor set. Complete coverage/numerical output does not
qualify continuum physics or full-opacity/reference parity. No renormalized
baseline replacement, second telluric correction, abundance/velocity/LSF fit.

Actual continuum commands, separate new output directories per reference:

```sh
for continuum_reference in synthesis iag; do
  /Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_kitt_residuals.py \
    --comparison data/results/rya1281/solar_bound_v1/kitt_peak/report.json \
    --wide-spectrum "$BRIDE_CHECK_STORE/rya1281-ng-extended-v1" \
    --continuum-reference "$continuum_reference" \
    --output "$BRIDE_CHECK_STORE/rya1281-ng-extended-continuum-v1/$continuum_reference"
done
```

15result artifacts/4visually inspected figures are byte-checked copies under
data/results/rya1281/ng_extended_v1 and data/output/rya1281/ng_extended_v1.
Actual writer/helper invocations registered in GENERATORS.yaml. Full6380614-byte
derived R/C/J archive remains immutable in experiments.bride_ng_extended/rates.npz,
SHA256 cc1796473ac55589186e093d7469c7824145daca8384dd869f0cb2709a8f67b9;
analysis/rates_external.json receipts it. Largest committed copy400167bytes;
no size-check bypass or renamed binary. Original600s refusals remain intact.

[Convergence history](../../../data/output/rya1281/ng_extended_v1/extended_history.png),
[intrinsic numerical spectrum](../../../data/output/rya1281/ng_extended_v1/extended_spectrum.png),
[synthesis continuum](../../../data/output/rya1281/ng_extended_v1/continuum/synthesis/continuum_controls.png),
[IAG continuum](../../../data/output/rya1281/ng_extended_v1/continuum/iag/continuum_controls.png).

Bug-ledger disposition: none discovered in this slice. Earlier setup/packaging
fixes and the baseline-reproduced historical ledger failure remain documented
in ng_control_review.md. No source/canonical parameter/adopted abundance/gold
or product changed. Default600s policy and numerical equations remain unchanged;
this opt-in resource allowance is not a universal performance policy.

Next: continuum-anchor/reference and source-qualified opacity/FTS broadening
controls, plus the independent matched O atom/EOS/reference gate. Numerical
sampling at these endpoints cannot explain the several-mA reference spread.
Broader composition/interactions, source-domain gaps, additional O diagnostics,
C/N before Fe and a reference star, full3D/all-band goals remain open. HOLD;
no adopted abundance, production, merge or publication. All owned jobs finished.

Final QA: MIRROR targets, registerv210 structure/version pointer, whitespace and intended artifact tracking pass. Generator guard retains exactly2436pre-existing diagnostics, zero added/removed. Scientific self-review checks unchanged physical kernels/source bindings, cost-only control difference, exact replay/native flags, saved starting arrays, independently reconstructed physical rates/target nuclei and explicit reference/anchor limitations. Existing gold/tracker/library debts are retained; no new bug or silent fallback discovered.
