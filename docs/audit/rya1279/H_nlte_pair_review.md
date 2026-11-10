# RYA-1281 paired hydrogen-normalization O NLTE diagnostic

Scientific HOLD. This is a fixed-column sensitivity experiment, not an abundance
measurement or admission of either hydrogen/EOS convention.

## Inputs and controls

The two explicit H budgets from the preceding paired background/LTE experiment
are `model_mass_H` and `native_network_H`. Their ratio spans
0.9919675919913356–1.0190607555015219 across 56 depths. Both retain the native
MARCS 5750 K/logg4.5/xi1 column, T/Pe/ne/geometry/composition, O8.70, the same
41-level/256-line/22-continuum source atom and experimental geometric top cap.
This is not a 5772 K solar model or a self-consistent chemistry/EOS update.

Both start from the same saved canonical-control level fractions, rescaled to
their respective O totals. The prior Fe-background solution supplies initialization
only; its spectrum is not a matched baseline. No radiation/rates or sequential
warm start is imported. Bare native FastBackground and existing passive atoms
are used, with an explicit empty molecular population table and no added Fe/CN
line provider. Hydrogen/proton background continua remain active.

The predeclared plan and runner retain the previous numerical controls exactly:
5 angular rays, doubled line quadratures, 6796 SE wavelength rows, joint population/
rate/J/cap tolerance1e-6, nested radiation tolerance1e-8, max1500 outer iterations,
600 seconds per case, Ng with invalid-population rejection, and seven emergent
flux rays. Only the diagnostic-window endpoints augment the full O SE grid.
Emergent flux is sampled afterward at 4001 AIR points,7763–7783A/0.005A.
Canonical AIR-to-vacuum conversion occurs once before LW; no observation,
continuum fit, broadening/RV adjustment or telluric operation is involved.

Source TAUTSP/source-direction collision and author MN LTE-pair hypotheses are
unchanged and receipted. Their physical qualifications remain open. Installed
LW0.16.1 binary/source lineage limits from the preceding audit still apply.

## Reproduction

Load `config/hosts/mac.env`; use the held Bride virtualenv with
OMP_NUM_THREADS=1 and OPENBLAS_NUM_THREADS=1. Run serially:

```
python scripts/run_bride_H_nlte_pair.py --case model_mass_H --output <fresh-model-case>
python scripts/run_bride_H_nlte_pair.py --case native_network_H --output <fresh-native-case>
```

Both dry runs passed input, state, source and Context guards before real execution.
The comparison script requires two completed manifests and independently rebuilds
physical rate-flow checks before reporting any population/flux differences:

```
python scripts/analyze_bride_H_nlte_pair.py --source <pair-directory> --output <fresh-analysis>
```

It predeclares fixed ±0.5A triplet windows and all41levels/56depths. Window areas
are intrinsic residual-flux diagnostics, not validated isolated-line EWs or
abundance errors. Both cases inherit joint1e-6 numerical accuracy.

## Completed model-mass-H case

Passed in589.9710480831563seconds,282outer iterations, peak626524160bytes,
one thread. Fresh independently rebuilt class-C plus transition-R residual is
3.977245928878939e-7; O conservation error8.881784197001252e-16.
Final population/J/cap changes are7.99464292004337e-7,
8.000425464230864e-7 and8.000451693988869e-7. A fresh standalone true continuum
and4001point emergent spectrum completed within the same600second budget.
Nonmonotonic Ng iterations were retained; no tolerance or source adjustment.

Seven original-name artifacts are byte-identical to the external run directory,
including rates.npz4239894bytes (below5MiB); no file split, rename or guard bypass.

## Native-network-H timeout: no admitted pair

The serial native-network-H attempt exited1 at the original600second deadline
inside fixed-population cap radiation, with
`TimeoutError: Fixed-population cap radiation budget exceeded`.
The last complete record is iteration274: population change4.918023753026354e-4,
rate residual2.4399927082155957e-4, J change1.6082772845897125e-6,
cap change1.6082819133503577e-6 and conservation1.106004177131581e-12.
These are incomplete-iteration-history diagnostics, not a fresh final audit.
No final populations/rates/continuum/flux/report were emitted. Input and full
history are preserved unchanged; stderr remains at
`/private/tmp/rya1281-H-nlte-native-network.log`. The model case's stderr is
`/private/tmp/rya1281-H-nlte-model-mass.log`.

The paired analysis was deliberately not executed: its two-completed-manifest
precondition is unsatisfied. There is no valid paired departure, continuum or
triplet-area difference from this attempt. Neither runtime nor failure establishes
a preferred physical H normalization. This is a characterized budget/convergence
limit, not evidence of a physical defect or an excuse to relax tolerances.
Total attempted solver allocation was two600second cases; successful model run
589.97seconds plus the native timeout, excluding dry runs/tests/analysis setup.

Next bounded task: predeclare a new shared initialization from the completed
bare-background model-mass-H populations, renormalize those same fractions for
BOTH fresh cases, and retain all original controls/600second per-case limits.
Start neither case from the other case's newly converged state. Import no J/rates;
independently rebuild class-C/transition-R and conservation again. Preserve this
attempt and account for both attempts' compute. Only after both pass may the
already committed comparison script run. This closer initialization tests numerical
cost; it does not validate a new EOS or establish unique physical convergence.

Thirteen original-name result artifacts are archived: seven completed-model,
two incomplete-native, and four dry-run files. All copies were byte-verified.
The failed native run contains no admitted spectrum to compare against Kitt Peak.

## Verification and limits

The five focused test files for H initialization, rate audit, finite cap, native
spline and charge contract pass12tests in2.60s with the Mac resource profile.
An initial invocation without the profile had7pass/5explicit resource skips;
the full profiled run resolves those skips. Numerical kernels remain unchanged.

Fresh library audit:211documents,127accounted,84unindexed,12off-Mac (exit1).
Eight documents were added since the previous203-document audit. No new paper
was used here. The first lookup used an incorrect library path and failed;
the audit above uses the registered `refs.library` holding.

Main remains7f023fd7(PR615Ca notes), as already inspected; no merge/rebase.
Independent MULTI departure discrepancy,41-level atom lineage versus published
50/122-level references, primary molecular provenance, native provider parity,
chemistry/EOS, C/N/Fe, reference-star and adaptive3D/all-band work remain open.
No preferred normalization, changed model/default/product, measured solar
abundance, production, merge or publication admission follows from this test.

MIRROR, version-pointer, sequence and staged-whitespace checks pass. Generator
invocation audit retains exactly2436historical diagnostics, zero added/removed,
after staging the13new artifacts; this is not a globally green audit claim.
Existing global bug-ledger RYA587blank science_impact and tracker/gold annotation
debts remain as previously recorded. Ordinary push hooks and postcommit freshness
checks are required; no bypass. No numerical/source defect was established by
this budget-limited attempt, so no speculative bug-ledger correction is made.
