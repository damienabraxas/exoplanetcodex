# Layerwise n5 continuation and native rate decomposition

The declared n5 frozen-prior-layer sensitivity is running on Sirius in
`runs/lifetime_frozen_prior_n5_v1`. No n5 state or lifetime is qualified at this
checkpoint. All prior full-release failures remain rejected. Linear declaration
`9218e4a5-e9a1-4ac1-8201-93d434ef1f60` preceded execution.

`lifetime_plan.v13.json` defines the test: use the qualified frozen-layer n4
states, revalidate their original six source stages and all source receipt
hashes, retain the existing full2s-SD/DC n5 configuration spaces, initialize only
the nine missing n5 functions natively, and optimize 5* only. All sixteen n<=4
records must remain exactly equal to a same-input native all-fixed loading
baseline. Every state needs original H=.03125/N439/ACCY, physical labels,
target CSF/nucleus identity and strict residual/fresh-restart qualification.
NSIC0, damping .5, NSOLV24 and 1500-cycle ceiling remain explicit. Historical n3
component-constraint limitations remain; this is not full core relaxation.

Only after all six stages and physical odd-state ordering pass may the driver
calculate the thirteen channels with both gauges. Any failure stops the campaign
and preserves its outputs. No failed-wave reuse or automatic parameter fallback.
The active calculation runs single-thread, nice19, without changing native source.

## Additional finding

`scripts/rya1354_layer_change.py` independently replays native n3 and frozen-n4
rate tables and checks A proportional to energy cubed times line strength for
each dominant ground-triplet E1 component and both gauges. The native
energy-cubed ratios are .95638-.95669; B strengths increase 3.58-3.62%, while C
strengths decrease 5.20-5.27%. Corresponding rate ratios are .99058-.99129 and
.90594-.90691. Thus energy changes alone do not explain the differing gauge
response. This is a theory-to-theory decomposition, not a causal separation of
configuration-space and optimized-orbital effects or an uncertainty estimate.
Machine record: `source_curation_v2/layer_change_v13.json`.

## Validation and exact resume

The focused suite passed 113 tests in 54.53s, including initial source tests.
The final four-test extension suite passed in 8.51s, adding decomposition
corruption rejection and terminal frozen-status collection without scientific
admission. Production/reference/catalog comparison against original task base
and whitespace checks pass. The initial lossless archive replays 65 paths in a
321216-byte part, including complete static ground initialization, executed
drivers, plan and source/target receipt hashes. It contains no SCF qualification.

Solver driver:
`evidence/frozen_driver_v13/rya1354_frozen_extend.py`.
Launch logs: `evidence/frozen_driver_v13/launch.{stdout,stderr}`.
Observe `runs/lifetime_frozen_prior_n5_v1/progress.json` and the eventual
`receipt.json`; absence of progress is not a qualified stage. During handoff the
native ground all-fixed loading control was running, before outer optimization.
The owned driver PID at launch was 3371026, collector 3371111; verify live
arguments before using PIDs. Do not relaunch or overwrite directories.

Independent watcher:

```sh
python3 evidence/frozen_driver_v13/rya1354_collect_pending.py --root /mnt/codex-data/codex/rya1354 --out /mnt/codex-data/codex/rya1354/evidence/frozen_collection_v13 --campaign lifetime_frozen_prior_n5_v1 --watch
```

This command is already running; read its heartbeat
`evidence/frozen_collection_v13/status.json`. Terminal outputs will be archived
under `evidence/frozen_collection_v13/lifetime_frozen_prior_n5_v1`. The collector
now recognizes the distinct frozen-model terminal status while preserving HOLD.
Retrieve all archive parts only after completion, replay hashes and independently
audit all six native stages and thirteen channels before reporting any n5 rates.
Compare against the declared n4 model without selecting gauges or observed
energies. EXPOSED_NONBLIND; adoption false; accepted lifetime, uncertainty,
covariance and calibrated sigma null. No grade, canonical product/reference/star,
publication, merge or Done transition.
