# Joint-outer timeout and matched qualified-orbital CI continuation

The joint4*/5* n5 attempt is **FAILED_RETAINED** after its43200-second ground
limit. All459 complete cycles exceed original ACCY5.960464477539063e-8;
the best maximum residual is5.95e-5 at49 and the last complete maximum is1.87e-4
at459. Cycle460 is unfinished. The trajectory oscillates with excursions;
longer wall time alone is not supported as a convergence remedy. No state,
fresh restart, odd-state calculation or lifetime qualified. The physical cause
remains unresolved; this timeout is distinct from earlier radial-node failures.

Native receipt SHA256:
`5c78b4cda419125225ad01246a43c8552d3f82fa0444c8e7427f38223116c456`.
Manifest SHA256:
`ed2f8bab59280a27c31f401464cbacbdcb67fbbf5d64ac5f9e55bfa03e380ce6`.
`joint_failed_v16` preserves89 paths in one1018658-byte lossless part.
Independent `rya1354_joint_failure.py` verifies qualified n5 source, raw controls,
loading input identity and complete/unfinished cycle separation. All9 fixed n<=3
records in retained ORBOUT match the native loading baseline; exact file-to-cycle
timing is not asserted. This does not admit the failed wave. Machine report:
`source_curation_v2/joint_failure_v16.json`. Existing numerical bug-ledger row is
amended with the timeout and remains characterized, not resolved.

## Next physical comparison, declared before execution

Linear comment`4fdf2457-f08e-455e-8935-a3acdf5db888` and plan
`lifetime_plan.v16.json` declare matched native RCI on the **verified sequential
n5** states. It uses no failed joint/core-release wave and performs no orbital
optimization. Both DC and static-limit Breit use identical full2s-SD n5 CSFs,
orbital files and nucleus, all3 states and13 channels per case.

DC first measures representation/solver consistency against prior SCF rates;
native internal RCI loading/interpolation is retained and exact equality is not
assumed. Static Breit then tests that Hamiltonian contribution at matched inputs.
Both gauges remain, with no observed-energy substitution or experimental tuning.
This is the same implementation, not an independent relativistic method or an
accuracy calibration. Native stationary-nucleus warnings remain part of evidence.
Historical n3 spectroscopic-component limitation and rejected full release remain.

The source helper replays complete n4/n5 provenance and stage/channel checks
before creating outputs. Any semantic refusal, changed source C/W/nucleus,
missing mixing, wrong physical state/order or missing transition stops with all
outputs retained. Source C/W/nucleus equality is checked again after native CI.

## Validation and live handoff

Focused suite126passed85.06s. An additional targeted startup-refusal test passed
after the terminal CI auditor was prepared. Tests replay full joint timeout,
reject relabelling it successful, and match DC/Breit/root controls byte-for-byte
against retained native n5 inputs. `rya1354_qualified_ci_audit.py` prepares full
terminal replay: source provenance, matched state inputs, physical labels,
mixing, exact Hamiltonian/root/operator controls, all13 channels each, and native
SCF/DC/Breit E1 decomposition. Its successful new CI path awaits terminal output;
no completed CI result is claimed now. Startup archives cannot be admitted.

`qualified_ci_initial_v16` replays67 paths in251162bytes. Static ground
C/W/nucleus/labels exactly match the qualified n5 source; executed drivers and
plan are preserved. No partial native rate is presented as a result. All new
archive files stay below5MiB. Production/reference/catalog/star comparison
against task base and whitespace/register checks pass. Standing bibliography
audit still exits1 on pre-existing civis2013/mallinson2024 vocabulary errors before
library traversal; no library-success or remote-CI-pass claim.

Sirius root:`/mnt/codex-data/codex/rya1354`. Run:
`runs/lifetime_qualified_ci_n5_v1`; immutable drivers:
`evidence/qualified_ci_driver_v16`; source:`rya1354_qualified_ci.py`.
Driver PID3707710 and collector3707794 were verified at launch; check actual live
arguments before acting. Native DC ground mixing exists and the remaining states
are running; no full DC/Breit case completion is claimed at this checkpoint.

Already running, single-thread/nice19:

```sh
python3 evidence/qualified_ci_driver_v16/rya1354_qualified_ci.py --source runs/lifetime_frozen_prior_n5_v1 --n4-source runs/lifetime_frozen_prior_n4_v1 --parent runs/lifetime_relax_n4_from_relaxed_n3_v7 --n4-core runs/lifetime_core_n4_v1 --core runs/lifetime_core_n5_v1 --grasp-bin toolchains/grasp/build/bin --out runs/lifetime_qualified_ci_n5_v1
```

Read launch.stdout/stderr in that driver directory, run progress.json/eventual
receipt.json and `evidence/qualified_ci_collection_v16/status.json`. The existing
collector watches lifetime_qualified_ci_n5_v1 and archives terminal success or
failure to qualified_ci_collection_v16/lifetime_qualified_ci_n5_v1. It retains
hashed rci.res caches on Sirius and preserves all other non-angular artifacts.
Do not relaunch/overwrite. Retrieve every part after terminal completion, then:

```sh
/Users/ryanschmitt/codex/venv/bin/python scripts/rya1354_qualified_ci_audit.py --bundle <new-local-terminal-bundle> --out <new-audit-json>
```

If failed, retain and audit failure instead of calling the successful result
path. Reviewed CI changes are model sensitivities, not calibrated uncertainty.
EXPOSED_NONBLIND; scientific HOLD; adoption false; accepted lifetime, uncertainty,
covariance and calibrated sigma null. No canonical line/gf/grade/product/reference/
star, source spectrum, frozen gold, publication, merge or Done change.
