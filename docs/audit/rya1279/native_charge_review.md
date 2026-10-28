# Native MULTI interpolation and reverse charge-transfer investigation

Continuation of RYA-1281 from ae5e65c3. Worktree codex/rya1279, existing stacked
draft PR #598. Refreshed main is now 5d517cbc (Na PR #595 merged); the experimental
branch is not merged/rebased into that new scientific state. All new calculations
are bounded one-thread Mac diagnostics. Scientific admission remains **HOLD**.

## Main result

An explicit reverse charge-transfer hypothesis substantially improves the
independent reference comparison in the declared photospheric diagnostic subset,
without changing abundance or source coefficients. It also restores full-atom
collision-only LTE closure. This is strong evidence about a missing reverse path
under the stated LTE-reservoir assumptions, **not** proof of the original grid
producer's implementation and not full atmosphere/population parity.

At logtau500=-1.2621:

| Control | Triplet lower b | Triplet upper b | O II ground b |
|---|---:|---:|---:|
| Independent grid reference | 1.03018 | 0.84559 | 0.87156 |
| Native interpolation, RH directions | 1.00884 | 0.99283 | 2.93e-8 |
| Native interpolation, source directions | 0.96136 | 0.77540 | 10.1705 |
| Above with LTE-restoring reverse CH0 | 1.02693 | 0.84906 | 0.94050 |
| Reverse CH0 plus LTE H background lines | 1.07133 | 0.88919 | 0.99304 |

For the reverse-CH0 case, upper triplet maximum absolute b difference over
-3<=logtau500<=0 is **0.00539**, versus about 0.18 in the original RH-convention
comparison. This subset was declared in earlier work; it is not a calculated
contribution function and must not conceal the outer-layer disagreement.

## Correcting the source search and building a direct numerical oracle

The previous search was incomplete: it missed the lowercase gencol/rcoll/gencl
routines in the already-held mul23_opacu.f. The inspected 1999 GENCOL copy was
therefore an unnecessarily old comparison. The pinned held source is
[pymulti 4d6aa3c, mul23_opacu.f](https://github.com/MingjieJian/pymulti/blob/4d6aa3c39d634b38da772e469a5692435eba1354/pymulti/files/multi/mul23/source/mul23_opacu.f).
Its 2007 RCOLL/GENCL path interpolates log CE/CI/CH0 tables with TAUTSP/PPVALU,
and CH in ordinary coefficient space. It retains downward CH and input-directed
one-way CH0. The spline flexibility parameter is the source's 5.5.

scripts/build_bride_multi_spline_oracle.py verifies the pinned source bytes
against Git, extracts only TAUTSP, PPVALU and INTERV, expands the source PREC
double-precision mechanism and builds a small local library using already-held
GNU Fortran 16.1.0. No full MULTI installation, dependency installation or source
modification. The generated source, compiler command/version, routine hashes,
original source hash and library hash remain in an immutable build receipt.
The library is resolved through experiments.bride_multi_spline_oracle, with an
explicit environment override for new immutable rebuilds.

pipeline/bride_native_spline.py calls the actual compiled TAUTSP/PPVALU routines,
not a substitute SciPy spline. Tests cover affine laws, endpoint behavior,
log-linear exponential laws, units/directions, thermal balance and rejection
outside the supplied temperature domain. Native INTERV has SAVE state, so calls
are serialized. The control caches coefficients only for the verified unchanged
atmospheric temperature array. CE/CI/CH/CH0 values remain source values; native
negative-CH spline clipping is explicitly counted and was **zero** in these runs.

This is a native interpolation oracle with inspected rate formulas in the adapter,
not a full native MULTI transfer/EOS calculation. The exact modified executable
and input files that generated the May-2021 departure grid are still unidentified.

## Separating interpolation, direction and reversibility

Four predeclared full-atom NLTE solves retained native MARCS 5750/4.5, xi=1 km/s,
A(O)=8.70, every held level/line/continuum, fixed ne/nH, boundaries and numerical
gates. All four converged, with 197.37 s summed case wall time:

| Case | Iterations | Seconds | Intrinsic triplet window areas, mA |
|---|---:|---:|---|
| Native spline, RH directions | 433 | 29.75 | 61.31457 / 53.08074 / 41.84159 |
| Native spline, source directions | 661 | 43.92 | 71.87854 / 62.31535 / 49.28368 |
| Native spline, source directions + reverse CH0 | 437 | 32.47 | 68.77840 / 59.22948 / 46.27885 |
| Reverse CH0 + LTE/CRD H background lines | 1131 | 91.23 | 69.13147 / 59.53236 / 46.51256 |

The original RH/WENO control was 61.31473/53.08087/41.84167 mA. Thus native
interpolation alone is negligible at this node/tolerance; it does not explain
the weak triplet NLTE response. Earlier linear-log CH0/source-direction control
gave 71.82768/62.26735/49.24060, also close to its native-spline counterpart.
The large changes are associated with rate directions and the reverse-path
hypothesis, not with selecting a more convenient interpolation scheme.

All source temperatures lie inside their own coefficient grids. No extrapolation
or rate/abundance scaling. Joint population/J/physical-rate tolerance is unchanged
at 1e-6, particle conservation at 1e-10. Largest terminal physical-rate residual
is 1.80e-14. The source-direction run rejected one invalid Ng proposal via the
already-tested safeguard; the other three rejected none.

These are the original fixed +/-0.5 Angstrom apertures, not isolated full oxygen
EWs or an abundance fit. The independent-b replay was 69.12376/59.50498/46.46724
mA. Corrected Kitt Peak observed blended windows remain 77.55651/67.20039/53.79430.
No new abundance is inferred from any of these differences.

## What the reverse-process hypothesis means

For stored forward ion-pair rate q, the control retains forward q*n(H I) and
derives its reverse from the local LTE O population ratio. Under an LTE collider
reservoir this makes C_forward*n_lower_star = C_reverse*n_upper_star. No free
factor is chosen. This represents a reversible charge-transfer pair under those
assumptions; it is not a newly sourced mutual-neutralization table and does not
prove the grid producer used this prescription. The adapter explicitly rejects
this hypothesis when hydrogen populations differ from LTE or the source process
does not connect adjacent ion stages.

Full 41-level collision-only LTE diagnostics use the actual backend C operator:

| Operator | Maximum absolute C*nstar / gross flow |
|---|---:|
| Original RH interpretation | 0.99999999982 |
| Native spline, RH directions | 0.99999999982 |
| Native spline, source one-way directions | 0.96403 |
| Native spline plus LTE-restoring reverse CH0 | 8.00e-16 |

Closure is expected from the reversible construction; the full-atom test checks
implementation/units and the presence of all processes. It does not independently
certify the absolute rate coefficients. The improved reference departure pattern
is separate evidence. A declared non-LTE external collider reservoir would require
its own consistent forward/reverse populations; this control cannot be reused for
that case by assumption.

## Remaining outer-layer and radiation discrepancy

The photospheric improvement does not extend over the whole atmosphere. At the
top native point (logtau500=-4.9177), the reverse-CH0 model has lower/upper b of
0.02512/0.03186 versus reference 0.98848/1.28139. Its maximum upper difference
over all depths is 1.24953, and all-level/depth maximum difference is 7.86243
(O II ground at logtau500=-4.7205). These are material unresolved differences.

A verified fixed-population replay under the paired operator extends process
budgets to logtau500=-4.9 and -4.0. New case receipts verify the native control's
manifest, populations and parent summary before replay. Local paired populations
retain maximum fresh-radiation physical imbalance 4.67e-7; LTE and independent
external populations give 0.89269 and 0.38729 under this operator, respectively.
Foreign populations are not expected to be equilibria of unproven local physics.

At the top point for triplet lower level 5, the external-population replay's net
imbalance is -0.00708 of gross flow, with its largest negative channel the
277.248 nm bound-free transition (-0.00416), followed by triplet transitions.
The local paired solution instead balances CH, CE, 135.56/135.85 nm radiative
coupling, triplet loss and that photoionization channel. This identifies useful
next diagnostics; it does not uniquely assign a wrong cross-section or background.

The H-background test adds the native H_6_atom as fixed LTE populations with CRD
lines. Oxygen nstar remains unchanged, and tests confirm hydrogen stays LTE.
It is **not** a coupled H NLTE/PRD calculation or a full line-blanketing model.
It changes populations and spectra, but does not fix the top layers. Its maximum
upper difference is 0.05529 over the declared photospheric subset and 1.21552 over
all depths, worse than the no-H-line paired case in the photospheric subset.
Yet its triplet areas nearly match the independent-b replay. This demonstrates
why matching these spectra alone cannot validate the depth-dependent solution.

## Native reference readiness and provenance limits

The held source checkout and wrapper do not include the required Uppsala ABSDAT
opacity input. It was also not found in the Mac's checked engine/default-pymulti
locations or a bounded read-only lookup under Sirius engines and /srv/codex.
The holdings registry contains no corresponding input receipt. No substitute
ABSDAT/EOS/opacity data were fabricated, and no full native transfer comparison
is claimed. The wrapper also expects atom-specific run/restart inputs, which
would need a verified setup for this exact atmosphere, composition and grid.

The distribution provenance pins the current O41 file and May-2021 grid bytes;
it does not contain the producer executable hash, complete atomic-data snapshot,
background-opacity configuration or a native rate/J dump. Those remain the
critical independent-reference requirements. A source checkout with the same code
family is useful evidence, but not a substitute for that provenance.

## Reproduction, artifacts and checks

```sh
source config/hosts/mac.env
"$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/build_bride_multi_spline_oracle.py \
  --source-root /private/tmp/rya1279-pymulti-source \
  --output "$CODEX_WORK_ROOT/bride_inputs/rya1281/native-spline-v2"
nice -n 15 "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/run_bride_mac_checks.py \
  --plan config/experiments/rya1281_mac_native_spline_controls.json \
  --output "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-native-spline-v1"
nice -n 15 "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/run_bride_mac_checks.py \
  --plan config/experiments/rya1281_mac_reverse_charge_control.json \
  --output "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-native-reverse-v1"
nice -n 15 "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/run_bride_mac_checks.py \
  --plan config/experiments/rya1281_mac_h_background_control.json \
  --output "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-native-h-lines-v1"
```

Use new directory suffixes on reproduction; outputs cannot be overwritten. For a
new oracle build, set CODEX_BRIDE_MULTI_SPLINE_ORACLE to its new directory. Final
analysis is native-charge-analysis-v3; verified rate replay is
native-paired-rates-v2 with explicit --case-directory and --case-receipt (exact
invocations in GENERATORS.yaml). Earlier analyses/replays/build remain preserved.
The native library and complete rate/J arrays remain external with hashes.

All four new SE jobs and all replays finished. No Sirius science job, dependency
installation, raw/FITS/gold/adopted-value change, production emission, merge or
publication. New emergent spectrum validation remains the 777 nm window; internal
UV radiation samples and H lines are not all-band observed-spectrum validation.
Corrected Kitt Peak/provider continuum/telluric/standard-air handling is unchanged.
Earlier O6300/O6158/Ni outputs retain their explicitly different RH-convention scope.

175 focused LW/parser/numerical/reference/observable/control tests and 3 installed
iSpec policy/effects checks pass, with existing 666/4 dependency warnings. The
oracle tests exercise actual compiled routines. The fixed-physics grid assembler
also rejects diagnostic background/rate variants. Prior bibliography audit remains
red (186 documents, 125 accounted for), and prior stale-six-band/generated-surface
aging issues remain documented. No green CI claim; final push skips Sirius CI.

Next: obtain the exact producer implementation and native opacity/atmosphere/rate
receipts, then isolate the remaining 135/277 nm radiative/background differences
under matched physics. Retain reverse-charge closure as a clearly labelled
hypothesis until that evidence arrives. More solar elements still precede another
star; full3D, all bands and grid admission remain separate gates.
