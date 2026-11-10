# RYA-1281 common-seed paired O NLTE restart

Scientific HOLD. This continues the [preserved first attempt](H_nlte_pair_review.md),
including its native-H timeout; no previous attempt or compute is erased.

## Predeclared intervention

The only numerical intervention is a closer common starting population field.
Both fresh cases use the completed v1 bare-background model-mass-H fractions,
renormalized separately to their respective O totals. Neither imports J/rates,
and native-H does not start from the newly completed v2 model-H state. The seed
is an initialization aid, not an independent reference or admission of one H/EOS.
Plan `config/experiments/rya1281_H_nlte_pair_restart.json` was committed before
execution. The runner gained an explicit `--plan` argument; its default retains
the original plan. Numerical helpers and thresholds were not changed.

The same 56-depth MARCS5750/logg4.5/xi1 column, T/Pe/ne/geometry/composition,
O8.70/41-level atom, source collision hypotheses, bare native backgrounds,
empty molecule table and experimental finite cap apply. No added Fe/CN opacity.
Controls remain5rays/q2/6796SE points, joint1e-6, nested radiation1e-8,
1500outer iterations and600seconds per case; flux7rays/4001AIR samples.
The diagnostics do not update chemistry, charge conservation or the atmosphere.
AIR-to-vacuum conversion occurs once; there is no observed-spectrum fitting,
broadening/RV, continuum adjustment or telluric correction in this experiment.

## Reproduction and evidence

Load `config/hosts/mac.env`; use the held Bride virtualenv, with
OMP_NUM_THREADS=1 and OPENBLAS_NUM_THREADS=1, executing cases serially:

```
python scripts/run_bride_H_nlte_pair.py --plan config/experiments/rya1281_H_nlte_pair_restart.json --case model_mass_H --output <fresh-pair>/model_mass_H
python scripts/run_bride_H_nlte_pair.py --plan config/experiments/rya1281_H_nlte_pair_restart.json --case native_network_H --output <fresh-pair>/native_network_H
python scripts/analyze_bride_H_nlte_pair.py --source <fresh-pair> --output <fresh-analysis>
```

The Context dry run passes. The model-H repeat passes3outer iterations in
17.223587040789425seconds, peak683458560bytes, one thread. Independent fresh
rate residual3.967907458149324e-7; conservation8.881784197001252e-16.
Its final continuum/flux complete. This shorter runtime measures initialization
cost, not improved physical accuracy or a preferred H budget.

The same five focused test files pass12tests in2.50s with the Mac profile.
Fresh reference-library audit:212documents/127accounted/85unindexed/12off-Mac,
exit1; no newly added source is used. Main remains7f023fd7, the previously
inspected PR615Ca notes; no merge/rebase. Serial Mac execution only; no Sirius/CI.

The first attempt used589.97seconds for model-H plus the native-H600second
timeout. Those costs remain part of this experiment's history, in addition to
the fresh restart times. No increased per-case budget or relaxed tolerance.

## Completed pair and intrinsic response

Native-H passes48outer iterations in92.27868958329782seconds, peak685801472bytes,
one thread. Independent fresh rate residual3.1753466052034624e-7;
conservation6.661338147750939e-16. Final population/J/cap changes are
6.536723900302631e-7/8.649511383891181e-7/5.3012182620716e-7.
Both fresh independent audits pass the original1e-6screen. The sum of fresh
run times is109.5023seconds; earlier first-attempt costs above remain included
in the history. Both600second limits are unchanged.

The committed analyzer exits0, rechecks source/manifest/hash/common-state identity,
starting fractions (max relative difference1.33227e-15), and independently rebuilds
per-class-summed-C plus transition-R net/gross residuals from saved arrays before
measuring profiles. The same fixed T/Pe/ne/geometry differs only in total-H
normalization, whose native/model ratio spans0.9919675919913356–1.0190607555015219.

| Fixed AIR centre (A), ±0.5A | Model-H area (mA) | Native-H area (mA) | Native minus model (mA) |
|---|---:|---:|---:|
|7771.944|68.77889666|68.97322329|0.19432663|
|7774.166|59.23292849|59.41000585|0.17707735|
|7775.388|46.28324511|46.43251495|0.14926984|

These are finite intrinsic residual-flux areas relative to each case's freshly
computed true continuum, not validated isolated EWs or inferred abundances.
Native-H slightly deepens the three cores: maximum absolute normalized-flux
difference0.0008426290723997365(0.0843percentage points). Maximum relative
continuum difference0.0023478579486202467(0.2348percent). The NLTE populations
respond: across all41levels/56depths the largest relative departure-coefficient
change is0.007881707180964925(0.7882percent), largest absolute departure change
0.10662785405934727, and largest absolute normalized-population-fraction change
2.299372549969484e-8. A tiny absolute fraction change can coexist with a larger
relative response in sparsely populated excited levels. No line contribution
function or atmosphere-average interpretation is assigned to these maxima.

Both plots were inspected: overlaid intrinsic profiles are close, the difference
panel localizes negative core changes, and the all-level departure map is strongest
in upper layers/excited states. All56depths are retained; no clipping of a difficult
level or plotting of only triplet levels. This is VIS triplet output, not an
all-band spectrum or a solar-data comparison. The full source atom's radiative
grid participates in SE, which does not establish all-band product validation.

For numerical context, direct comparison of the unchanged model-H v1/v2 outputs
gives max normalized-flux drift1.5886292281663827e-10, true-continuum relative
drift1.3988810110276972e-14, and respective window-area drifts
-2.9079683105948106e-8/-2.1822287032957632e-8/-1.5514270423899745e-8mA.
These were calculated with `window_integral` on the two original4001point flux
files. This one continuation is a useful restart-stability diagnostic, not an
independent accuracy bound, mesh convergence test or proof of unique convergence.
Both paired solutions still inherit joint1e-6accuracy and source-model uncertainty.

Twenty-two original-name artifacts are byte-verified against external originals:
two dry-run, seven per completed case and six analysis outputs. Maximum file size
is4239918bytes, below5MiB (no split/rename/guard bypass). The v1 timeout remains intact.

## What this establishes and next controlled task

The earlier LTE-fraction invariance does not extend to NLTE departures. Changing
H affects opacity and collision/SE/radiation together; this pair measures their
combined fixed-column response. It does not identify how much comes from the
continuum versus collision channels, establish a preferred H budget or reconcile
the native EOS. The source-defined normalization mismatch is measurable but this
signal alone does not establish an explanation for prior observed-fit discrepancies.

Next, predeclare a matched fixed-population transfer control: map the SAME
converged O fractions to each named H total, recompute each background/continuum/
radiation under unchanged geometry/cap, and hold those fractions fixed. Compare
that conditional transfer response with this fully coupled pair; any remainder
is an intervention-dependent SE feedback diagnostic, not a unique additive
physical decomposition. Require fresh fixed-population radiation convergence,
positive/finite flux, population identity and source guards. Do not demand SE
balance from deliberately frozen populations or label them converged SE.
Primary/EOS provider qualification precedes real CN insertion; no synthetic
native-purpose refusal is bypassed. Then return the admitted diagnostic spectra
to the pinned current-ORCH AIR/continuum/telluric comparison, without tuning the
observations or changing default products.

Primary atom lineage, independent MULTI departure discrepancy, collision/cap
hypotheses, molecular provenance/native provider, coherent EOS/chemistry,
C/N/Fe/reference star and adaptive3D/all-band qualifications remain open.
O8.70 is an unchanged input, not a derived solar abundance. No model/default/
product/production, merge or publication admission.

MIRROR/version-pointer/sequence/staged-whitespace checks pass. Generator invocation
audit retains exactly2436historical diagnostics, zero added/removed; this is not
a global-green claim. Existing RYA587blank science_impact and tracker/gold debts
remain as previously recorded. Postcommit stacked-base freshness and ordinary
push hooks are required; no bypass. No new numerical/source defect was identified.
