# RYA-1354 terminal MR audit and numerical continuation

Continuation starts at143b6d63f0ddf1a37512e2199cfe91ba63241759, original
task base adefcaf073ff70035954cf45a9a2283774ac7d06. Refreshed main
10f0dc955803253570a944ca37012654e80ead9f is independent ongoing work.
Existing task branch/worktree and draft PR616 retained. NONBLIND; scientific
HOLD. No production values, grades, products or canonical/reference lines change.

## Evidence collection failure and repair

Both original campaigns ended. MR completed all four cases. Separate n4
relaxation failed ground2s: methods1/2 counted3nodes versus1required, then
IMPROV/Convergence not obtained, native exit1. No ground or dependent lifetime
qualified. The v5 collector archived that failure but crashed on one MR native
`lower.upper.+2T` file larger than the5MiB compressed ceiling. Its waiting
heartbeat became stale; that heartbeat never supersedes native receipts.

The v2 bundle format fragments an oversized individual file into contiguous raw
chunks, each with offset/bytes/SHA256. Archive parts and stored members are
hashed, then original paths are reconstructed and checked against original
hashes and receipts. No filesystem extraction. Old v1 readers/generators and
historical v4/v5 summaries remain unchanged. `rya1354_native_parts.py` reads
both formats; `rya1354_terminal_summary.py` audits the new terminal evidence.
The collector writes an explicit failure heartbeat and re-raises errors.
All original raw outputs and partial v5 collection remain on Sirius. Only
documented mcp.* angular caches and explicitly hashed regenerable rci.res
matrix caches stay outside the bundle, retained on Sirius.

New driver copies: `evidence/collector_driver_v6`; output:
`evidence/pending_collection_v6`. Command from the resolved Sirius task root:

```sh
PYTHONPATH=. python3 evidence/collector_driver_v6/rya1354_collect_pending.py --root . --out evidence/pending_collection_v6
```

Full MR audit requires unchanged qualified n4 wave/nucleus inputs, exact prior
SR/DC reproduction, actual CSFs equal to the independently generated MR
inventory, identical core-MR DC/Breit CSFs, physical states/downward order and
all13decay channels. A failed relaxation is retained separately with a null
accepted lifetime. Gauge agreement/model shifts never define sigma or grades.

Collection completed at00:56:45UTC on2026-11-10. Local reconstruction,
original receipt hashes and complete model/identity/channel replay passed.
`lifetime_mr_complete_v6` preserves490paths/344931569original bytes in35parts
totaling99400208compressed bytes; maximum part5205280bytes, below5MiB.
Six individual transition-output files require fragmentation. The original
failed relaxation is18paths/one340369-bytepart in `lifetime_relax_failed_v6`.
Part/member/fragment/original SHA256 and explicit cache receipts all replay.
Native RCI stdin matches each declared Hamiltonian and requested roots.
`source_curation_v2/terminal_collection_v6.json` preserves terminal collection
status; `terminal_summary_v6.json` is the deterministic scientific audit.

| Fixed-orbital n4 model | B lifetime(ns) | C lifetime(ns) |
|---|---:|---:|
| SR/DC | 2.500065889 | 2.180712938 |
| VV-MR/DC | 2.459565977 | 2.160236577 |
| full2s-MR/DC | 2.295118548 | 1.733118149 |
| full2s-MR/staticBreit | 2.295984146 | 1.732686510 |

The same-space SR control reproduces the previous calculation exactly. Core
configuration expansion materially changes the lifetime; substantial B/C
separation remains. Neither a closer individual gauge to the exposed Smith1971
value nor the small static-Breit shift calibrates uncertainty or accuracy.

## Explicit optimization control and mandatory verification

Declared before launch in `config/atomic/rya1354/lifetime_plan.v6.json` and Linear.
Native source inspection: getscd.f90 selects method1..4 with four subshell
lists. solve.f90 METHOD>2 skips methods1/2 node-based energy correction;
improv.f90 retains upstream method2 recovery on integration failure. Method3
alone is therefore not node qualification. Use method3 for2s optimization
only, then restore original default node-checked methods for a fresh NSIC0
restart. Original strict orbital criterion and physical identities must pass
for every state before any lifetime. No automatic wrapper fallback.

Same qualified `runs/lifetime_core_n4_v1` initializers, full2s-SD/DC CSFs,
allthree states; vary2/3/4 and fix1s. Same gridH0.03125/N439 and original
ACCY5.9604644775390625e-8, damping0.5,1500cycles, optimizationNSIC20.
No failed-state initializer, observed energies/strengths, solver-source patch
or fitting to experiment. New output directory; preserve every failed stage.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nohup nice -n 19 python3 rya1354_lifetime_relax_v6.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --core runs/lifetime_core_n4_v1 --n 4 --damping 0.5 --cycles 1500 --two-s-method3 --out runs/lifetime_relax_n4_method3_v2 > runs_lifetime_relax_n4_method3_v2.stdout 2> runs_lifetime_relax_n4_method3_v2.stderr < /dev/null &
```

Initial driver2604334/rmcdhf2604387; discover current arguments/status rather
than trusting stale PIDs. About9.5GiB memory available/127GiB disk free before
launch, single-thread/nice19. Other workers untouched. Read terminal receipt
and logs; process execution/SCF stopping alone is not physical validation.

At00:57UTC ground optimization reached printediteration26, unqualified.
Detached collector2606014 watches only this explicitly named campaign:

```sh
PYTHONPATH=. nohup nice -n 19 python3 evidence/collector_driver_method3_v6/rya1354_collect_pending.py --root . --out evidence/method3_collection_v6 --campaign lifetime_relax_n4_method3_v2 --watch > evidence/method3_collection_v6.stdout 2> evidence/method3_collection_v6.stderr < /dev/null &
```

Its tested30second heartbeat is `evidence/method3_collection_v6/status.json`.
The saved `source_curation_v2/method3_pending_snapshot_v6.json` is explicitly
an observation at its timestamp, not a current-state or completion claim.
The watcher archives terminal success/failure; it launches no science and
makes no acceptance decision. Named campaigns cannot escape task run paths.

Inspected pinned source SHA256:

| File | SHA256 |
|---|---|
| getscd.f90 | 67ff53450acf233f21bcb9191eb350ec559259967585136044f291210e96e96b |
| solve.f90 | 2d0ea3e544f0ba914701da089b5d043bd674bfc274bc6f3afb64e0ea6ce2c85f |
| improv.f90 | 7695750e89a4572bfb639a2348d1a30e40f92ac127af9f9dbdcd06f011494679 |

## Verification and remaining work

Mac existing Python3.12.13 only replays evidence/metadata; native computations
remain Sirius Python3.14.4/pinned toolchains. Fragment tests discriminate
missing, corrupt, duplicate, unclaimed and reordered chunks; old v1 evidence
replays. Forced collection error replaces stale heartbeat while propagating
failure. Native-method input test proves only optimization method selection
changes; fresh verification stdin is identical to defaults.

Focused regression command: **86passed in35.01s**, including full terminal
bundle replay, exact original native stdin preservation, misdeclared native
Hamiltonian refusal and explicitly scoped watcher tests:

```sh
python -m pytest tests/test_rya1354_native_parts.py tests/test_rya1354_lifetime_path.py tests/test_rya1354_strict_scf.py tests/test_rya1354_decay_inventory.py tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

Two bug-ledger findings: oversized-evidence/stale-status defect fixed; n4 ground
node failure characterized, physical root cause unresolved. Refinement-debt
report21rows/10not-fired/5no-ticket, none cleared. Library audit exits1 before
traversal on existing civis2013/mallinson2024 vocabulary errors.

Self-review checked native receipt/fragment/part identity, original fixed-wave
and nucleus bytes, actual MR CSFs, matched Hamiltonian controls, complete
operator/channel inventory, failed-state separation and mandatory default-method
verification. No unresolved wrong-data/provenance defect found in admitted
exploratory MR results. No physical accuracy/convergence sign-off inferred.
Original GRASP source clean; protected production paths unchanged against task
base. Whitespace and normal register freshness pass (existing informational
stale-tracker notice retained); no state surface changes. Previous hosted check
on143b6d63 was CANCELLED. Final post-commit PR-mode guard and push are recorded
in the Linear handoff. Draft push retains the documented pre-push exception for
unrelated merged1352/1353 SEQUENCE debt; file-size ceiling is enforced.

Next: monitor the separately declared optimization control through mandatory
fresh restarts. Inspect its terminal receipt and collector status, copy stable
evidence into a new directory and audit successful states or retain failures
before considering further numerical or physical controls.
Larger MR/active-space and radial/orbital convergence, independently implemented
relativistic checks, primary experimental independence/uncertainty/covariance
curation and human scientific review remain owed. Li/P/S are preparatory leads,
not scored or certified blind benchmarks. Keep PR draft and ticket In Progress;
no merge, adoption, calibrated theoretical grade or Done from this checkpoint.
