# RYA-1354: terminal node failures and staged relaxation

Continuation starts at42227155fc8bd0f0aade722bd3b9cd199e2b6c48; original
baseadefcaf073ff70035954cf45a9a2283774ac7d06. Refreshed main
45901e7678ab6c8624fe659e552fba3f2681457d contains the independent RYA-1378
landing. Existing task branch/worktree/draft PR616 retained. Scientific HOLD,
exposed/nonblind. No production gf, grade, abundance, canonical/reference or
stellar changes; source/other-worker outputs preserved.

## Both previous controls failed; larger budgets were insufficient

Both collectors completed; full terminal outputs now replay locally, separately
from the historical running snapshots. Each contains18native paths.

| Control | Local bundle | Compressed bytes | Final outer cycle |
|---|---|---:|---:|
| NSOLV24, previous component selection | nodebudget_failed_v9 | 340907 | 194 |
| NSOLV24, explicit both2p | pcomponents_failed_v9 | 355289 | 206 |

Both remain FAILED_RETAINED with native exit1, no optimized ground label,
fresh verification or relaxed lifetime. The first193complete cycles of the
nodebudget run reproduce the original default-budget trajectory exactly.
Last complete maximum residuals62.4 and619 are far above originalACCY.
Unchanged qualified initializer/CSF/nucleus bytes and actual native stdin verified.

| Control / attempt | Inner cycles/limit | Bracket width | Width floor | Nodes/target |
|---|---:|---:|---:|---:|
| nodebudget method1 | 3/24 | 0.095302707465 | 0.125 | 3/1 |
| nodebudget method2 | 9/24 | 0.088909310206 | 0.125 | 3/1 |
| both2p method1 | 24/24 | 0.924328225400 | 0.125 | 3/1 |
| both2p method2 | 13/24 | 0.069725811336 | 0.125 | 3/1 |

Pinned SOLVE rejects when inner iteration allowance is exhausted OR the native
energy-parameter bracket is narrower than1/n^3. Both final method2 refusals
match the bracket guard before24iterations, still counting3nodes instead of1
for2s. The both2p method1 attempt did exhaust its allowance; the distinction is
retained. Increasing the budget again is not assumed to address the instability.
This identifies the numerical refusal, not the underlying physical cause.

`source_curation_v2/node_failures_v9.json` is deterministic output of
`rya1354_node_failure_summary.py`, including original receipts/part hashes,
failed attempts and separate last-complete versus unfinished cycles. No failed
intermediate value is adopted or treated as an uncertainty estimate.

## Declared staged optimization path

`lifetime_plan.v9.json` and Linear declare this before launch. The final physical
model and full2s-SD/DC CSF space remain unchanged. Freeze the **2s radial function**
temporarily, not its configurations or CI coefficients:

1. All3states: vary2p*/3*/4*, fix1s/2s; explicit2p-/2p (+3s odd) spectroscopic
   selection. Require strict original orbital criterion, exact frozen native
   radial records, physical identities and freshNSIC0 verification.
2. Only after all3preconditioners pass, release2s: vary2*/3*/4*, fix1s; optimize
   and fresh-verify all3states again under original criteria.
3. Only final released states may generate13channel transition evidence for
   offline review. The restricted intermediate model supplies no lifetime.

Same qualified original n4 core initializers, NSOLV24, H0.03125/N439,
ACCY5.9604644775390625e-8, damping0.5,1500cycles. No failed-wave, observed-energy,
strength or lifetime inputs. No automatic parameter fallback or criterion waiver.
Any numerical refusal, fixed-record/CSF/nucleus mismatch, incorrect state/order
or failed verification stops the whole campaign and preserves outputs.
Native source is unmodified; code dependencies are copied into the run's drivers
directory. Progress records are not terminal receipts or scientific acceptance.

Sirius command from task root,9.1GiB available/127GiB disk free before launch:

```sh
PYTHONPATH=. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nohup nice -n 19 python3 evidence/staged_driver_v9/rya1354_lifetime_staged.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --core runs/lifetime_core_n4_v1 --out runs/lifetime_relax_n4_staged_v5 > runs_lifetime_relax_n4_staged_v5.stdout 2> runs_lifetime_relax_n4_staged_v5.stderr < /dev/null &
```

Initial driver2818729/rmcdhf2818804, single-thread/nice19. Watcher2818782:

```sh
PYTHONPATH=. nohup nice -n 19 python3 evidence/collector_driver_method3_v6/rya1354_collect_pending.py --root . --out evidence/staged_collection_v9 --campaign lifetime_relax_n4_staged_v5 --watch > evidence/staged_collection_v9.stdout 2> evidence/staged_collection_v9.stderr < /dev/null &
```

Heartbeat `evidence/staged_collection_v9/status.json`, native progress
`runs/lifetime_relax_n4_staged_v5/progress.json` after qualified stages.
Watcher archives terminal success/failure without launching physics or granting
admission. Discover live arguments before trusting old PIDs.

`staged_logs_v9` preserves5files/one12795-bytepart: asynchronous log-prefix
snapshot at15:02:02UTC, fixed_2s_ground cycle11, last maximum residual0.000284.
Actual native stdin and14varied-orbital rows confirm1s/2s absent and both2p
METHOD1. No stage is qualified. Snapshot files/hashes/history replay; this is
not a terminal receipt. Capture command (existing output must not be overwritten):

```sh
python3 rya1354_scf_history_v9.py --run runs/lifetime_relax_n4_staged_v5 --stage fixed_2s_ground --out evidence/staged_log_snapshot_v9
```

## Merged RYA-1378 queue reconciled

PR621 is merged at45901e76. Directly inspected merged generator/README and
all27CSV/JSON pairs; row/held-reason/flag/priority/fit-evidence counts agree.
`source_curation_v2/gap_register_intake_v9.json` records per-file/source/canonical
hashes and pointers, not another parameter authority or copied gf values.

| Inspected population | Candidate rows | Priority1 |
|---|---:|---:|
| All27elements | 126399 | 20117 |
| Cr | 17817 | 2939 |
| O | 327 | 2 |

This supersedes the preliminary untracked17772-row Cr receipt as the current
queue; old receipts remain historical. Population/rules differ, not an inferred
data defect. UNRATED_GF means no primary lab attachment; held_reason distinguishes
actually UNGRADED_GF from GRADED_NOT_LAB. Priority1 means no graded species/band
route, not poor atomic data. Fit evidence requires post-01f51846 continuum-gated
runs. Grade reliability, empirical line quality and research priority stay distinct.

Opriority1 rows: neutral gf_165007/pk_1cdd19e1b14b at3692.395Aair, OII
gf_109001/pk_7d85b7952f40 at9199.816Aair. Canonical records lack full upper/lower/J
mapping and attached primary lab uncertainty. They are identity/source leads,
not certified calibration truth. No OII expansion of the neutral pilot is implied.
Observability, primary dependence/uncertainty, model feasibility and exposure
review remain owed before later Li/P/S/Cr or oxygen-line tests. No sibling work
or merged canonical parameters changed.

## Verification and resume

Terminal native-failure replay checks actual inputs/initializers and proves the
pre-cap bracket refusals. Tests preserve historical controls, distinguish cap
versus bracket failures, enforce fixed radial records and prove a failed
preconditioner cannot release2s or complete a phase. Actual staged snapshot
tests verify the requested native scope and unaccepted-state metadata.

```sh
python scripts/rya1354_node_failure_summary.py --out data/audit/rya1354/source_curation_v2/node_failures_v9.json
python -m pytest tests/test_rya1354_native_parts.py tests/test_rya1354_lifetime_path.py tests/test_rya1354_strict_scf.py tests/test_rya1354_decay_inventory.py tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

Existing numerical-failure bug-ledger row amended with both terminal outcomes;
underlying instability remains characterized, not fixed. Read native receipt
and live watcher/progress before continuing. On termination, collect full stable
evidence into a NEW local directory and audit each phase/identity/frozen record,
fresh restart and final channels, or retain explicit failure before another
declared control. No intermediate lifetime, calibrated grade, merge, publication,
adoption or Done. Larger-space/radial/orbital and independent-relativistic
convergence plus curated experimental calibration/covariance and human review
remain owed.

Final focused regressions: **96passed in36.09s**, existing MacPython3.12.13
metadata/native replay only; SiriusPython3.14.4/pinned toolchains for native
calculation. Terminal audit exit0 and deterministic replay; no full-suite/hosted
CI pass inferred. Scientific self-review checked actual native refusal guards,
unchanged initializer/CSFs, staged source/frozen-record/phase admission, no
intermediate lifetime and merged research-queue/calibration distinction. Original
solver source clean; protected production paths unchanged against original base.
Whitespace/normal freshness pass with existing informational stale-tracker notice;
no state surfaces change. Bibliography audit retains existing civis2013/mallinson2024
vocabulary failures before traversal. Final SHA/push/actual check status and live
owned-stage observation are in Linear. Draft push retains the documented exception
for unrelated merged1352/1353 SEQUENCE debt; file-size checks remain enforced.
