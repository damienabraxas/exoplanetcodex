# RYA-1354: NSIC0 outcome and extension from core-relaxed n3

Start004e3c20f08e9e0094d2e4f9f362557a4ed0de67, original
baseadefcaf073ff70035954cf45a9a2283774ac7d06. Refreshed main
45901e7678ab6c8624fe659e552fba3f2681457d; independent work preserved.
Existing task branch/worktree/draft PR616. Exposed/nonblind, scientific HOLD.
No production gf/grade/abundance/canonical/reference/stellar changes.

## How the NSIC0 run ended

`staged_nsic0_failed_v11` preserves38native paths/one565338-bytepart.
Receipt2d07e4062e760808050424b2d7ab103aa5592788679fed5e7756771d0279ee13;
manifest560b28e4be9879e18ea4d31fcc56707f7d21c5d6c80a19c361b9f141c7740067.
All hashes/native inputs/baseline sources replay locally.

The fixed_2s_ground preconditioner fails at352outercycles, native exit1. Final
2p method2 reaches24/24inneriterations with1node instead of0; bracket width
13.6428101294 remains above0.125. This is a cap/node refusal, unlike earlier
final narrow-bracket refusals. Last complete351: maxresidual32.7, weighted
energy-74.900863347Hartree. All351complete cycles contain14regular orbital
updates and no extras. Removing extra updates did not resolve the instability.
No physical state/fresh restart/released phase/lifetime qualified.

The native loading baseline remains exact for1s, but2s differs by2.22e-16 inP.
Upstream IMPROV recovery calls global ORTHSC, which renormalizes same-kappa
functions without LFIX filtering. ORTHSC SHA256:
8d528d69fa31fa246ad4d77ab55f83a6d7dbed0ee396e7d7883529439fb4d1dd.
This explains why a loading-only control does not guarantee fixed-record identity
after later recovery; the exact bit-change timing is not established. The audit
records the mismatch and retains the exact guard. No epsilon/tolerance waiver
or admission of the failed state. The initial exact2s audit assumption was
corrected against actual records, rather than silently normalizing them away.

`rya1354_nsic0_failure.py` generates
`source_curation_v2/nsic0_failure_v11.json` with failure/cap/record/update evidence.

## Different initialization and active-layer path

The already completed `runs/lifetime_relax_n3_v3` supplies core-adapted initial
radial functions. Source receipt SHA256:
27e8e9595a40840832cc6c050782ab4d3ce58565a61779df0054d5605bce5179.
The historical receipt lacks an explicitn field; actual artifacts determine
the layer. Validation checks every selected file's bytes/hash, native strict
criterion/439grid/H0.03125, physical state labels, fresh restart input identity,
eight-electron/frozen1s/full2s-hole CSFs and exact nine-orbital n<=3inventory.
All3states and their verification stages pass. The historical2p- spectroscopic
selection limitation is retained; this source is not calibrated truth or a
claim of new component qualification. Final new states must requalify explicitly.

`lifetime_plan.v11.json` and Linear declare the following before launch:

1. Keep original n4full2s-SD/DC target CSFs/nucleus. Native rwfnestimate loads
   each qualified n3 wavefunction and initializes only missingn4 orbitals using
   the previously tested hydrogenic initializer on the same439grid.
2. All3states optimize only4*, freezing allnine n<=3functions. Require exact
   native-loaded frozen records, original strict residual/identity/order and
   freshNSIC0 verification. No lifetime from this restricted intermediate stage.
3. Only after all3outer-layer states qualify, vary2*/3*/4* with1s fixed and
   requalify all3states/fresh restarts before13channel final evidence.

NSIC0/NSOLV24/H0.03125/originalACCY5.9604644775390625e-8/damping0.5/1500cycles
unchanged. No failed-wave/observed input, automatic fallback or increased cap.
This tests a core-adapted initialization/active-layer path, not experimental fitting.
Native source unmodified; exact dependencies copied into run/drivers. Any source,
fixed-record, numerical or physical gate failure stops the whole campaign.

Command from resolved Sirius task root,8.8GiB available/127GiB free before launch:

```sh
PYTHONPATH=. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nohup nice -n 19 python3 evidence/staged_driver_v11/rya1354_lifetime_staged.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --core runs/lifetime_core_n4_v1 --n3-core-relaxed runs/lifetime_relax_n3_v3 --optimization-nsic 0 --out runs/lifetime_relax_n4_from_relaxed_n3_v7 > runs_lifetime_relax_n4_from_relaxed_n3_v7.stdout 2> runs_lifetime_relax_n4_from_relaxed_n3_v7.stderr < /dev/null &
```

Initial driver3030614/rmcdhf3030719, single-thread/nice19. Discover current
arguments before trusting PIDs. Watcher3030660,30second heartbeat:

```sh
PYTHONPATH=. nohup nice -n 19 python3 evidence/collector_driver_method3_v6/rya1354_collect_pending.py --root . --out evidence/from_n3_collection_v11 --campaign lifetime_relax_n4_from_relaxed_n3_v7 --watch > evidence/from_n3_collection_v11.stdout 2> evidence/from_n3_collection_v11.stderr < /dev/null &
```

Live status `evidence/from_n3_collection_v11/status.json`; run progress.json
after qualified stages only. Neither grants scientific admission.

## Preserved initialization and early native scope

`from_n3_checkpoint_v11`:36paths/one431528-bytepart. The posthoc checkpoint
captures closed initializer directories for all3states, plus an asynchronous
outer_n4_ground log prefix. It is NOT a terminal campaign receipt or accepted
state. Source/target receipts and every captured initializer file/hash replay.
Each previous.w exactly matches the qualified historical n3_verify wavefunction;
each target CSF exactly matches original n4core input. Native initializer stdin
matches the declared hydrogenic/fine-grid control. All16target orbitals present.

At2026-11-11T04:08:44UTC, snapshot printed13cycles/12complete; last complete
maxresidual9.51e-5/Eweighted-75.001292555Hartree, still unqualified. Native
rows confirm only seven4* functions vary, with seven regular updates percycle.
No n<=3function is selected for ordinary optimization. Exact fixed-record
qualification remains a later explicit gate; no accuracy improvement inferred.

### Closed restricted ground stages now independently replayed

By04:13:33UTC, outer_n4_ground optimization and its fresh verification had
completed. `from_n3_ground_partial_v11` preserves71native paths in one610872-bytepart
with an explicit posthoc closed-stage checkpoint, NOT a terminal campaign receipt.
`rya1354_ground_preconditioner.py` generates `ground_preconditioner_v11.json`.

| Restricted ground stage | Iterations | Largest native residual |
|---|---:|---:|
| Optimize4* with n<=3 fixed | 32 | 5.11e-8 |
| FreshNSIC0 verification | 1 | 1.77e-8 |

Both pass originalACCY5.96046448e-8 and native physical ground labels/grid.
Offline replay verifies actual controls/initial source/fresh wave, unchanged
target CSFs/nucleus, promoted wave/mixing bytes versus native outputs and exact
equality of allnine frozen n<=3records to each same-input loading baseline.
Only documented mcp.* angular caches remain outside the bundle, retained onSirius;
checkpoint hash/membership checks distinguish these from required artifacts.
This qualifies the **restricted numerical ground preconditioner**, not a final
fully relaxed state, completed campaign, lifetime or accuracy grade. Quintet,
triplet and all full-release stages remain pending at this checkpoint.

## Verification and resume

Tests cover full NSIC0 failure/cap replay and record mismatch, historical n3
source qualification/tamper refusal, unchanged default/fresh controls, outer-only
scope and frozen prior-layer records, actual three-state source/target/native
initializers and seven-update native snapshot. Initial source validation was
corrected for native CSF whitespace; no physics launched from a failed source.

```sh
python scripts/rya1354_nsic0_failure.py --out data/audit/rya1354/source_curation_v2/nsic0_failure_v11.json
python -m pytest tests/test_rya1354_native_parts.py tests/test_rya1354_lifetime_path.py tests/test_rya1354_strict_scf.py tests/test_rya1354_decay_inventory.py tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

Existing numerical-failure ledger row amended with NSIC0 outcome and recovery
semantics; no claim that instability or orthogonalization is fixed. Inspect native
terminal receipt, watcher and progress before proceeding. On termination collect
stable full outputs into a NEW local directory and audit every initializer,
baseline/frozen record/identity/fresh restart/final channel or preserve failure
before another declared control. No intermediate lifetime, calibrated grade,
adoption, merge, publication or Done. Larger-space/radial/orbital and independent
relativistic convergence, experimental curation/calibration/covariance and human
review remain owed. Merged1378 research queue remains as reconciled in v9,
not certified experimental or blind calibration truth.

Final focused regressions: **105passed in41.41s**, including closed restricted
ground/fresh/frozen/native-output replay; existing MacPython3.12.13
metadata/native replay only; native SiriusPython3.14.4/pinned toolchains.
Deterministic failure replay and native source/initializer checks pass. Scientific
self-review checked failure versus qualification, conservative fixed-record
guard, source model/grid/CSF/wave identity, no intermediate admission and phase
ordering. Original solver source clean; protected production paths unchanged
against task base. Whitespace/normal freshness pass with existing informational
stale-tracker notice; no state surfaces change. Library audit retains existing
civis2013/mallinson2024 vocabulary failures before traversal. Initial test command
had a nonexistent filename and ran no tests; the corrected full command above
passed. Final SHA, push/post-commit guard, actual check status and live resume
point are recorded in Linear. No hosted-CI/full-suite/whole-library pass inferred;
documented draft pre-push exception remains scoped to unrelated merged1352/1353
SEQUENCE debt, with file-size checks enforced.
