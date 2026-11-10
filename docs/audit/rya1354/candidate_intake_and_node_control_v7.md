# RYA-1354: Cr candidate intake and oxygen node-checked continuation

Startab666ac3b2b7e0641014e9c219ace8cda9d38bc1, original task
baseadefcaf073ff70035954cf45a9a2283774ac7d06, refreshed main
10f0dc955803253570a944ca37012654e80ead9f. Existing RYA-1354 worktree,
branch and draft PR616. Scientific HOLD; exposed/nonblind. Production atomic
data, grades, reference lines, abundances and other workers remain unchanged.

## RYA-1378 intake without replacing its register

Ryan requested checking RYA-1378 while continuing oxygen. Read that issue and
its comment, then inspected actual RYA-728 CSV/generator read-only. The sibling
worktree has active uncommitted changes; CSV/generator are still untracked and
are not committed by sibling HEAD31a2cee42a074e96041735e5dfe64ad63b4a85bf.
No sibling files modified or duplicated into a competing register.

`source_curation_v2/cr_candidate_intake_v7.json` records source/generator hashes,
available columns, measured counts, source scope and intake requirements. CSV
SHA256b5f3de7b9be911ca95e8afd086312662aa407ad0f76d338aa7de07be05cabd4a.

| Measured item | Rows |
|---|---:|
| Total candidate rows | 17772 |
| UNRATED_GF | 17770 |
| MISSING_GF | 2 |
| PROVENANCE_GAP | 3909 |
| GF_DISCREPANCY | 7 |
| SOLAR_DISCREPANCY | 10 |
| Priority1 CrI H | 65 |
| Priority1 CrI K | 67 |
| Priority1 CrII red-optical | 2809 |

Flags coexist. The source comment says UNRATED on all; actual two missing
canonical rows have MISSING_GF and no physical_id instead. No duplicate
nonempty physical_id found. This receipt preserves the discrepancy; it does
not rewrite the producer or its definitions.

These are proposed research targets, not known-grade experimental calibration
controls. Missing grade does not mean poor atomic data. Full lower/upper/J
identity and source uncertainty/dependence are not supplied by this preliminary
table. Observability, blends/tellurics and theoretical feasibility remain owed.
Reconcile the reviewed upstream register before selecting specific calculations;
preserve CrII diagnostic-only restrictions and distinct atomic reliability,
empirical line quality and research-priority axes. Existing-source processing
does not attest blindness. Oxygen qualification continues independently.

## Rejected method3 control

`lifetime_method3_failed_v7` preserves18native paths in one404421-bytepart.
Manifest SHA256af7bfbe523a87c62b01d2926ed360736cb7b33e2c3def6be884951c71e69d4c0;
terminal native receipt238e16bfbd8c31c32524179744a1e0b3462b76d330beed9cd0f1fd1375dea5d1.
All part/member/original-receipt hashes replay locally. Native ground iteration307
exits1: IN refuses radial extent beyond439, warning of possible continuum state.
Final printed residuals2s8.23/2p-0.57/2p47.8 exceed originalACCY5.96046448e-8;
energies also diverged. No ground label, fresh restart or lifetime qualified.
This is not repaired by relabeling the failed output or enlarging the grid silently.

`rya1354_method3_failure.py` replays the failure and refuses a forged completed
status. `source_curation_v2/method3_failure_v7.json` retains diagnostics and null
accepted lifetime/uncertainty. Previous completed MR results remain as audited
in [v6](terminal_continuation_v6.md); that agreement/model sensitivity is not
calibrated accuracy.

## Declared next diagnostic: native inner-iteration floor

Before launch, `lifetime_plan.v7.json` and Linear declare original node-checked
methods, NOmethod3, with NSOLV24. Pinned getscd.f90 calls NSOLV a factor;
actual solve.f90:93 computes `NLOOPS=MAX(NSOLV,3*NP(J))`. This is an iteration
floor, not24*n. Ground2s limit6 becomes24, while native node targets and
bracket-width refusal remain unchanged. Original method1 failed on narrow bounds
before its cap, so a larger allowance is not assumed to cure it. No model tuning
to literature or automatic fallback.

Same qualified full2s-SD/DC n4 initializers; all3states vary2/3/4, fix1s;
H0.03125/N439/originalACCY/damping0.5/1500outercycles unchanged. Same strict
orbital and physical-state checks; fresh NSIC0 restart with node-checked methods
and declared allowance required before any lifetime. No failed-wave initializer.
Native solver source remains unmodified. Launch resources9.3GiB available,
127GiB disk free; no other active atomic solver observed. Single-thread/nice19.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nohup nice -n 19 python3 rya1354_lifetime_relax_v7.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --core runs/lifetime_core_n4_v1 --n 4 --damping 0.5 --cycles 1500 --radial-iteration-floor 24 --out runs/lifetime_relax_n4_nodebudget_v3 > runs_lifetime_relax_n4_nodebudget_v3.stdout 2> runs_lifetime_relax_n4_nodebudget_v3.stderr < /dev/null &
```

Initial driver2667076/rmcdhf2667104. Discover live arguments before trusting PIDs.
Tested watcher, reusing preserved collector code with explicitly named campaign:

```sh
PYTHONPATH=. nohup nice -n 19 python3 evidence/collector_driver_method3_v6/rya1354_collect_pending.py --root . --out evidence/nodebudget_collection_v7 --campaign lifetime_relax_n4_nodebudget_v3 --watch > evidence/nodebudget_collection_v7.stdout 2> evidence/nodebudget_collection_v7.stderr < /dev/null &
```

Heartbeat `evidence/nodebudget_collection_v7/status.json`; it archives terminal
success/failure, launches no physics and grants no scientific admission.
Inspect native terminal receipt first; copy stable evidence to a new local
directory. Keep failures, require all3states/fresh restarts before13channel
lifetime. Larger MR/radial/orbital convergence, independent relativistic checks,
curated experimental calibration/covariance and human scientific review remain
owed. No merge, publication, grade, adoption or Done.

## Verification and handoff

At05:52UTC the new ground optimization reached printediteration8; no terminal
receipt or accepted state. Native stdin explicitly sets NSOLV24. The saved
`source_curation_v2/nodebudget_pending_snapshot_v7.json` is a timestamped
observation, not proof of current progress/completion. Watcher status above
remains the live collection observation; native receipt governs its terminal state.

Mac existing Python3.12.13 metadata/native replay only; Sirius Python3.14.4
and pinned native toolchains. Focused regression command:

```sh
python -m pytest tests/test_rya1354_native_parts.py tests/test_rya1354_lifetime_path.py tests/test_rya1354_strict_scf.py tests/test_rya1354_decay_inventory.py tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
python scripts/rya1354_method3_failure.py --bundle data/audit/rya1354/lifetime_method3_failed_v7 --out data/audit/rya1354/source_curation_v2/method3_failure_v7.json
```

**88passed in36.18s**; failure summary replay exit0. Tests discriminate forged
completed status and prove the sole changed native control is the inner-limit
prompt for all3states and fresh restarts. Historical v4/v5/v6 scientific outputs
replay unchanged. Source part/member/receipt identity and failure non-promotion
verified. No further scientific acceptance inferred from passing tests.

Scientific self-review checked preserved model/wave/nucleus identities, native
node/bracket/ACCY semantics, explicit source-code versus prompt discrepancy,
separate failed states and research-candidate/calibration-truth distinctions.
No new wrong-data/provenance defect found. Existing numerical-failure bug-ledger
row amended with method3 divergence; no duplicate row or claim of a fixed solver.
Protected production-path diff and original GRASP-source cleanliness pass.
Whitespace/normal register freshness pass with existing informational stale
tracker notice; no state surfaces changed. Bibliography audit still exits1
before traversal on pre-existing civis2013/mallinson2024 vocabulary errors.
Final SHA, post-commit guard, draft push/check status and live resume point are
recorded in Linear. Keep previous documented draft pre-push exception scoped to
unrelated merged1352/1353 SEQUENCE debt; no file-size-hook bypass.
