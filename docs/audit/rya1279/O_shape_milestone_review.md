# RYA-1281 O777 diagnostic milestone: complete, model qualification HOLD

The declared O777 observational-shape milestone is COMPLETE. Removing sensitivity
to a common multiplicative level brings Kitt Peak and IAG into close relative
agreement, but the current O model remains shallower/weaker than both. That
conditional discrepancy survives both shoulder schemes, both RV controls, both
H normalizations and the seven declared broadening controls. It is not resolved
by the previously measured H or absolute-continuum effects.

This is a diagnostic milestone, not a solar abundance result or model admission.
The current native column is5750K/logg4.5, whereas canonical solar parameters
remain5772K/logg4.438/xi1.0. Missing opacity, atom/source, collision/cap/EOS and
independent reference gates also remain. The observations do not select a unique
physical cause or justify fitting abundance/width/RV to erase the discrepancy.

## Milestone contract and completed evidence

Plan `config/experiments/rya1281_O_shape_milestone.json` was committed before
running. Acceptance requires every declared case, numerical/source controls,
quantitative shape/sensitivity assessment and one consolidated scientific
handoff. `assess_bride_O_shape_milestone.py` independently verifies the complete
Cartesian identities and absence of duplicates in every output table.

|Evidence|Completed coverage|
|---|---|
|Current source spectra|Two coupled H cases plus four frozen-population controls, all source/report/flux receipts checked|
|Atlas views|Shipped Kitt, native IAG, actual pinned-ORCH Gaussian IAG; native/smoothed IAG are correlated|
|RV conventions|Zero and+0.633km/s; no fit or choice of a preferred shift|
|Triplet|AIR7771.944/7774.166/7775.388; no other bands or elements claimed|
|Shoulders|Fixed symmetric0.35–0.55A and0.65–0.85A offsets|
|Core/feature|Core±0.1A and finite feature±0.5A; not isolated-line EWs|
|Broadening|Nominal plus six one-at-a-time controls; no optimization or joint sweep|
|Sampling|Native-grid metrics; shared0.005/0.0025A grids for unweighted residual/sampling checks|
|Completeness|540metrics,3024comparisons,1008sampling rows,216interval records,12assessment contexts|
|Numerical identities|12exact nominal broadened-profile replays;1080multiplicative-scale canaries pass|

The shape statistic divides each profile by the mean of its two fixed shoulder
means and records core contrast, relative core area and finite relative feature
area. It is bookkeeping insensitive to a COMMON multiplicative scale. It is not
a fitted continuum, does not cancel wavelength-dependent normalization/blending,
and does not identify an isolated O transition or a formation depth.
Existing `feature_metrics` arithmetic is reused unchanged.

Observer centres follow each predeclared RV convention using the relativistic
Doppler factor. Observations themselves are not shifted/recentred. Both zero-RV
and gravity-shift windows are retained, as in the previous FeII profile controls.
Source AIR ownership remains unchanged; no vacuum/air confusion is introduced.

## Fixed broadening controls

Nominal:vmac3.8km/s,vsini1.8km/s,R500000,limb darkening0.6. Six one-at-a-time
controls use vmac3.3/4.3,vsini1.3/2.3,R400000/600000. These are fixed diagnostic
values, not calibrated error bars, best-fit values or a complete interaction
sweep. Every source population/RT spectrum is frozen. The exact same iSpec
effects source matches the prior receipt; all12nominal outputs replay the
archived processed profiles exactly. Edge padding is retained outside requested
diagnostic support and never used as a continuum/shape pixel.

Actual pinned de96ae51 continuum source is byte-checked, but no continuum fit
is performed in this phase. Its IAG helper reproduces the archived reference
exactly before the native/Gaussian views are split. Both holdings pass policy
guards; all native IAG flags in each±1.1Afeature input are exactly0. No second
telluric correction, new atlas smoothing choice, LSF calibration or deconvolution.

## Both atlases retain the same stronger relative feature

Plus-gravity, inner shoulders, coupled model-H benchmark:

|AIR centre(A)|Model core contrast|Kitt contrast|Native IAG contrast|Model relative feature area(mA)|Kitt area(mA)|Native IAG area(mA)|
|---|---:|---:|---:|---:|---:|---:|
|7771.944|0.229079|0.242669|0.243319|65.29767|69.77610|69.69174|
|7774.166|0.201215|0.215139|0.215462|56.38540|61.45247|61.36635|
|7775.388|0.161740|0.172605|0.173072|44.27730|47.53331|47.43643|

Kitt versus native-IAG relative feature areas differ only0.08436/0.08612/0.09688mA
here, while the nominal model deficits are4.47843/5.06707/3.25602mA against
Kitt and4.39406/4.98095/3.15913mA against native IAG. Gaussian-IAG results
remain close:69.68515/61.36346/47.43778mA. Repeated views are not independent
measurements or an empirical noise model.

Model-minus-Kitt contrast differences are-0.013591/-0.013924/-0.010865;
native-IAG-minus-Kitt differences are+0.000650/+0.000323/+0.000467.
The model therefore retains a core-strength/shape discrepancy after common
level sensitivity has been removed. The plot shows this directly; it is not
only a quoted continuum-dependent EW residual.

The H-normalization contrast responses are+0.000667/+0.000622/+0.000541,
small relative to the remaining model/atlas contrast differences. Frozen
population controls and all original source variants remain in the full tables.
No stronger H case is selected as physically preferred.

## What the full declared controls do and do not close

All12RV/shoulder/line assessment contexts retain observed contrasts AND relative
areas above the interval spanned by both coupled H models and all seven declared
broadening variants. None of the216individual source/atlas diagnostic intervals
contains its observed contrast or area. This is a finite-test result, not a bound
over all possible atmospheric/broadening/line-formation physics.

For plus-gravity/inner shoulders, even the most permissive coupled-model interval
falls below the least-strong of the three observed views by contrast
0.007724/0.008564/0.006436 and relative area4.14445/4.76442/2.98650mA.
Across BOTH RVs and BOTH shoulder schemes, minimum remaining gaps are
0.00284620contrast and1.60933mA. No variant is picked as an optimum; the full
interval is descriptive and not a probability/confidence interval.

Outer shoulders change relative areas and expose absorption absent from sparse
models. Maximum observed shoulder fractional range across the declared cases
is0.0104078. All24fixed-band candidate inventories are retained; e.g. plus-gravity
inner blue/red candidate counts are6/3,3/5,3/2across the triplet. Positions are
not blend identifications or predicted contributions. Nothing is dropped or
moved to improve agreement, and alternative shoulders are not used to select an
answer. Other O-triplet wings can contribute to these bands as well.

## Numerical and source checks

Multiplying whole diagnostic profiles by0.92and1.08changes contrast by at most
4.44089e-16and relative feature area by3.90799e-13mA, within prospective
512eps screens. This proves the specific multiplicative invariance, not general
continuum accuracy. Nominal broadening replay is exact in all12source/RV cases.

Native-grid metrics are primary. Resampling to0.005/0.0025Achanges model contrast
by at most5.38339e-5and model relative area0.00046023mA. Observed equivalents
are5.73399e-5and0.00241043mA. This measures interpolation bookkeeping only;
finer grids add no resolution or independent data, and do not tighten RT/SE or
physical uncertainty. Coupled source populations retain the previous1e-6joint controls
and independently audited rates/conservation; frozen-population controls remain
transfer-only and are not admitted SE solutions. No new RT occurred here.

All original source/observation/reader/model/plan hashes are protected. Main
remains7f023fd7, inspected PR615Ca notes; no merge/rebase. Library audit remains
212documents/127accounted/85unindexed/12off-Mac,exit1; no new paper used.
Actual atlas/source-lineage caveats remain as recorded in
[the same-pixel audit](H_reference_values_review.md).

## Reproduction and artifact package

Load `config/hosts/mac.env`, use project analysis virtualenv, OMP/OPENBLAS1:

```
python scripts/analyze_bride_O_shape_milestone.py --output <fresh-shape-run>
python scripts/assess_bride_O_shape_milestone.py --source <shape-run> --output <fresh-assessment>
```

Both exit0. Analysis whole-process elapsed4.20094s, one thread, excluding test
and handoff work. Ten original-name outputs byte-match external originals:
8run outputs and2assessment files, maximum1029903bytes. Actual writers and
invocations are registered. No raw atlas or binary source copy, split/rename,
5MiBguard bypass, Sirius, CI, installation, download or benchmark.

Tests:15project tests pass0.32s; the4FeII metric/Stark tests pass in the held-LW
environment1.06s. An initial mixed-environment invocation had6pass/1failure
because its Stark test imports LW, absent from the project environment. Routing
that test file to its existing LW environment resolves the dependency without
install/source/tolerance changes. No source/numerical defect was identified;
the analysis and all declared controls pass. Warnings from existing LW default
configuration/deprecated library imports remain visible, not suppressed.

## Milestone disposition and next physical experiment

|Gate|Current disposition|
|---|---|
|This declared observational diagnostic milestone|COMPLETE: all identities/cases and conditional assessment present|
|Frame/source/numerical interoperability|Checks pass; independent instrument/physical accuracy not inferred|
|Absolute continuum|Source-value/geometry dependence characterized, no adopted placement|
|Current column versus solar shapes|Does not match under declared controls; no solar-model/abundance admission|
|Atmosphere match|OPEN:5750/4.5native benchmark is not canonical5772/4.438solar column|
|O atom/collision/cap/EOS/independent MULTI|OPEN; no claim of unique physical cause|
|Real Fe/CN/provider completeness|OPEN; primary/EOS/native-state gates remain mandatory|
|C/N then Fe/reference star/adaptive3D/all bands|OPEN; provenance work may proceed, not abundance admission|

The next meaningful physical milestone is a MATCHED solar-column intervention,
not another continuum fitted to reduce these residuals. Predeclare the4point
matrix5750/4.5,5772/4.5,5750/4.438and5772/4.438using the SAME current O
atom/collision/background/xi and independent convergence/rate controls.
This separates temperature and gravity structure responses before attributing a
solar discrepancy to abundance or a new physical component.

Use the already audited nearest native cell5750–6000K/g4–4.5and explicit
log-positive/layer-thickness interpolation contract with all56depths; retain
source/model identities and interpolation caveats, no nearest/extrapolation
fallback or claimed independent interpolation accuracy. Keep model-H/ne owners
explicit through the existing LW constants/H-mass boundary. Never transplant
the exported5750native-network H/CN packet onto an interpolated solar column.
First replay a native5750baseline on the SAME field/rounding ownership as the
new columns, documenting any difference from the old native-export bindings.

Column mass was not validated by that interpolator. Any new cap policy must be
explicitly predeclared (e.g. retain the existing geometric cap for the matched
intervention), not obtained from an invented interpolated column mass. Use the
same saved fractional initialization for all fresh cases, original600second
per-case budgets, independent fresh C+R and conservation, fresh continua and
this complete shape package. No source/default/product/production/publication/
merge, physical-model or adopted abundance/continuum admission follows here.
RYA1281and parentRYA1278remain scientifically open.

Final MIRROR/version-pointer/sequence/staged-whitespace checks pass. Generator
invocation audit retains exactly2436historical diagnostics, zero added/removed.
Previously verified global ledger RYA587blank science_impact and tracker/gold
debts remain; no global-green claim or new implementation defect. Postcommit
stacked-base freshness and ordinary push hooks are required; no bypass.
