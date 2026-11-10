# RYA-1354: adaptive extra-update and native loading controls

Start32acb2ea473d239323b75890a309e0f7b0a48455; original task
baseadefcaf073ff70035954cf45a9a2283774ac7d06. Refreshed main
45901e7678ab6c8624fe659e552fba3f2681457d. Existing task branch/worktree
and draft PR616 retained. Exposed/nonblind, scientific HOLD. No production
atomic values, grades, abundance, canonical/reference or stellar changes.

## Staged preconditioning failed before any qualified state

`staged_failed_v10` preserves24native paths/one419784-bytepart. Terminal receipt
015df0fa8505eb0c0a6949ea6d8fbd8b7623def800e188a5543dada016141e27;
manifest6f864f3dd0d93136b83eaa407d507e0656721e63b423c95a1d679ba9e3fa310c.
Native fixed_2s_ground fails at319outercycles. Method1 counts5nodes instead of0
for2p-; fallback method2 counts2instead of0, inner9/24, bracket width0.107468043558
below native0.125. Last complete cycle318 maxresidual2.61, weighted energy
-403.93909669Hartree. No physical preconditioner, fresh restart, released phase
or lifetime qualified. Fixing2s shifted the failure to2p- without curing it.

`rya1354_staged_failure.py` deterministically replays the native receipt/input,
failed attempts and last complete cycle. It also joins the loading control below.
No divergent output becomes an accepted prediction or uncertainty estimate.

## NSIC runs inside every cycle and can adapt

Pinned SCF.f90:154-187 places the extra MAXARR/IMPROV loop inside every outer
cycle, despite GETSCD wording about the first cycle. IMPROV:148 can decrement
NSIC when DEL1<0.0001 andNSIC>1. Thus20 is an initial allowance, not exactly
20forever. Native failed-run cycles1/2/3 contain34/23/23 update rows:14regular
plus20/9/9extra. This corrects the earlier shorthand of20extra everycycle.
It establishes implementation semantics, not that extra updates caused collapse.

| Inspected source | SHA256 |
|---|---|
| scf.f90 | 1873d0794f3a90614f2dd51a96960bd1b9d526ddc69cca5272f001a67d870ad3 |
| improv.f90 | 7695750e89a4572bfb639a2348d1a30e40f92ac127af9f9dbdcd06f011494679 |
| lodrwf.f90 | 413f6fc327a4311ec8ca16e504569e63d0a675fa8ae65c4c268a9c4805703069 |
| intrpq.f90 | 221345ba5ef75bc5bb38f1aac674f5f776f60443865c8c5a61d31a5622bc9594 |

## Native loading baseline fixes a potential false rejection

Raw pre-load fixed1s/2s records differ from saved output by normalization
roundoff: maximumP changes2.22e-16/1.11e-16, PZ change about1.78e-15.
LODRWF calls INTRPQ, which renormalizes via1/sqrt(RINT). The old raw-equality
guard would reject legitimately fixed records even after successful convergence.
It was never reached in the failed campaign, which stopped inside the solver.

A separate all-fixed native same-input/grid control reproduces those saved
fixed1s/2s records EXACTLY. `fixed_load_control_v10` preserves22paths in
one314869-bytepart, including original inputs/rejected output, exact source
drivers, native commands/executable hashes and control receipt. This proves
loading roundoff separately from orbital optimization; no state or lifetime is
admitted by the control.

The runtime now builds that same native loading baseline before **each**
optimization and verification, then compares frozen records exactly to it.
No epsilon, orbital criterion change or tolerance widening. Actual changed2s
records still fail. Baselines use the current incoming wavefunction, so fresh
restarts also account for their own native loading normalization. They run in
separate directories and do not replace the optimizer input.

Native control from task root (completed output must not be overwritten):

```sh
PYTHONPATH=. nice -n 19 python3 evidence/staged_driver_v10/rya1354_fixed_load_control.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --failed-run runs/lifetime_relax_n4_staged_v5 --out validation/fixed_load_control_v10
```

## Declared NSIC0 matched staged run

`lifetime_plan.v10.json` and Linear declare this before launch. Same v9 staged
path, original qualified seeds/full2s-SD/DC n4 CSFs, explicit both2p components,
NSOLV24/H0.03125/N439/originalACCY5.9604644775390625e-8/damping0.5/1500cycles.
Only the initial optimization extra-update allowance changes20to0. Each cycle
has a regular pass without additional least-self-consistent updates. Fresh
verification remains0. Other native method/node/bracket behavior is unchanged.

All3 fixed2s preconditioners must pass strict residual/identity/order/frozen-record
and fresh verification before2s release; all3released states must pass again
before13channel final evidence. No failed-wave/observed inputs, fallback or
preconditioner lifetime. Any failure stops the whole campaign. Original solver
source unmodified; exact driver dependencies copied into run/drivers.

```sh
PYTHONPATH=. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nohup nice -n 19 python3 evidence/staged_driver_v10/rya1354_lifetime_staged.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --core runs/lifetime_core_n4_v1 --optimization-nsic 0 --out runs/lifetime_relax_n4_staged_nsic0_v6 > runs_lifetime_relax_n4_staged_nsic0_v6.stdout 2> runs_lifetime_relax_n4_staged_nsic0_v6.stderr < /dev/null &
```

Initial driver2956988/rmcdhf2957109, one thread/nice19; about8.8GiB available
before launch. Discover live arguments before trusting PIDs. Watcher2957035:

```sh
PYTHONPATH=. nohup nice -n 19 python3 evidence/collector_driver_method3_v6/rya1354_collect_pending.py --root . --out evidence/staged_nsic0_collection_v10 --campaign lifetime_relax_n4_staged_nsic0_v6 --watch > evidence/staged_nsic0_collection_v10.stdout 2> evidence/staged_nsic0_collection_v10.stderr < /dev/null &
```

Live heartbeat `evidence/staged_nsic0_collection_v10/status.json`; per-qualified
stage progress in the run's progress.json. Neither grants admission. The saved
`staged_nsic0_logs_v10` is a timestamped asynchronous log prefix, not a terminal
receipt or qualified state. Actual native inputs and completed-cycle rows verify
NSIC0 and14regular updates without extras. All part/member/history hashes replay.

## Validation and next step

Existing MacPython3.12.13 runs metadata/native replay; native controls and solver
use SiriusPython3.14.4/pinned toolchains. Replay commands:

```sh
python scripts/rya1354_staged_failure.py --out data/audit/rya1354/source_curation_v2/staged_failure_v10.json
python -m pytest tests/test_rya1354_native_parts.py tests/test_rya1354_lifetime_path.py tests/test_rya1354_strict_scf.py tests/test_rya1354_decay_inventory.py tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

The first update-count assertion assumed a constant20 and failed; native
IMPROV/source and actual34/23/23 rows resolved the incorrect test expectation.
Regression coverage now preserves adaptive semantics, native loaded baseline
versus raw false rejection, unchanged default/fresh controls and actual NSIC0
row counts. Historical summaries remain unchanged. Existing numerical-failure
ledger row amended; one new fixed-record guard defect fixed with exact native
baseline, not a solver/source or physical-accuracy claim.

Inspect terminal receipt and live watcher/progress first on resume. On termination,
copy stable full bundle into a NEW local directory and audit every baseline,
phase/identity/frozen record/fresh restart/final channel, or retain explicit
failure before another declared control. Larger-space/radial/orbital and independent
relativistic convergence, experimental curation/calibration/covariance and human
review remain owed. Merged1378 queue remains as reconciled in [v9](staged_relaxation_v9.md),
not certified calibration truth. No grade, adoption, merge, publication or Done.

Final focused regressions: **100passed in34.64s**. Deterministic failure/native
load-control replay passed. Scientific self-review checked native adaptive-loop
semantics, normalized versus raw record identity, unchanged criterion/defaults,
per-stage baseline isolation, phase/physical admission and no intermediate
lifetime. Original solver source clean and protected production paths unchanged
against task base. Whitespace and normal freshness pass with existing informational
stale-tracker notice; no state surfaces change. Library audit retains existing
civis2013/mallinson2024 vocabulary failures before traversal. Final SHA, push,
post-commit guard, actual check status and owned-job resume point are in Linear;
no hosted-CI/full-suite/whole-library pass inferred. Draft pre-push exception
remains scoped to unrelated merged1352/1353 SEQUENCE debt; size checks enforced.
