# Verified n5 completion and joint outer-layer continuation

The sequential-layer n5 model completed and independently replayed all six
strict/fresh states and all thirteen decay channels. Native terminal receipt
SHA256 is `212310c7215b9948bb557f50dc3f380727c7416e66bd6bfc6085105c9bd4b07a`;
manifest SHA256 is
`c58f005c735bb967b71e3ca01ad237589dd829478ed6c1714e9b7ff75aab62eb`.
The terminal receipt supersedes its earlier progress file's RUNNING field.

| State | Optimization cycles / maximum residual | Fresh cycles / maximum residual |
|---|---|---|
| Ground | 44 / 5.53e-8 | 1 / 2.98e-8 |
| Quintet | 40 / 5.79e-8 | 1 / 2.96e-8 |
| Triplet | 48 / 5.28e-8 | 1 / 3.72e-8 |

Original ACCY is5.960464477539063e-8 on H=.03125,N439. Independent replay
verifies physical labels, controls, exact target CSF/nucleus, source/fresh inputs,
25-orbital inventory, all sixteen fixed n<=4 records relative to native loaded
baselines, odd-state ordering, transition inputs/operators and complete decay.

| Model | Gauge B lifetime (ns) | Gauge C lifetime (ns) |
|---|---|---|
| Frozen-prior-layer n4 | 2.402393349410827 | 1.6803476177190935 |
| Sequential-layer n5 | 2.2844837477524753 | 1.7132068154876408 |

The n5 rates are437735659.5264999 and583700689.8173959 s^-1. Lifetimes change
about-4.9% and+2.0%; substantial gauge dependence persists. Native E1 decomposition
for each dominant ground component verifies energy-cubed factors1.01535-1.01538,
B strength factors1.03550-1.03581 and C factors.96553-.96621. Energy changes alone
do not explain the different gauge response. These are exploratory sensitivities,
not a single-cause decomposition of configuration/orbital effects or error bars.
No gauge is selected against experiment; no observed energies enter the model.

`scripts/rya1354_n5_terminal.py` reproduces
`source_curation_v2/n5_terminal_v15.json`. Full source evidence is
`n5_complete_v15`:347 native paths, seven lossless parts totaling19704749bytes,
each below5MiB. The independent auditor shares the same byte-based checks with
native initializer validation; archive hash replay precedes semantic checks.
New tests reproduce the full native model and reject substituted loading and
transition waves. The native source-initializer test independently replays all
provenance and refuses changed mixing coefficients.

The focused suite passed123tests in75.65s, including successful native terminal
replay, source initialization and corruption refusals. Startup archive replay
also verifies exact qualified n5 ground C/W/nucleus and input aliases. Whitespace,
production/reference/catalog/star comparison against original task base, manual
archive-size checks and register freshness pass. Standing bibliography audit still
fails on pre-existing civis2013/mallinson2024 verification vocabulary before
library traversal; no library-success or remote-CI claim is made. Native source
remains clean; prior owned campaigns finished before this new run.

## Next declared test, now running

Before execution, Linear comment`dd38f202-e6cf-495f-b059-e1f89b8a28e3` declared
the joint outer-layer sensitivity. It starts only from verified n5 *_verify
states, retains identical n5 full2s-SD/DC CSFs/nucleus, and jointly varies4*,5*
with all nine n<=3 records fixed against same-input loading baselines. This tests
dependence on freezing n4 during the sequential n5 step. It does not release the
unstable inner core or use failed waves. Plan:`lifetime_plan.v15.json`.

Original grid/ACCY, NSIC0, damping.5, NSOLV24, explicit physical2p pair,
1500-cycle ceiling and fresh restart remain. All six states and odd ordering
must pass before thirteen channels. Any refusal or changed fixed record stops
the campaign and retains all outputs, without automatic fallback. Historical n3
component-constraint limitation and failed full core release remain explicit.

Sirius run:`runs/lifetime_joint45_n5_v1`; immutable drivers:
`evidence/joint_driver_v15`; source driver:`rya1354_joint_outer.py`.
Native driver PID3543975 and collector3544054 were verified at launch; check
actual live arguments before acting. At handoff ground angular preparation is
active; no joint model state or lifetime is qualified. Startup archive
`joint_initial_v15` replays65 paths in246791bytes: source inputs, executed
drivers, plan and source-receipt identity, with no SCF qualification.

Existing native command, already running:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nohup nice -n 19 python3 evidence/joint_driver_v15/rya1354_joint_outer.py --source runs/lifetime_frozen_prior_n5_v1 --n4-source runs/lifetime_frozen_prior_n4_v1 --parent runs/lifetime_relax_n4_from_relaxed_n3_v7 --n4-core runs/lifetime_core_n4_v1 --core runs/lifetime_core_n5_v1 --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --out runs/lifetime_joint45_n5_v1
```

Read launch.stdout/stderr in the immutable driver folder, the run's progress.json
and eventual receipt.json, and `evidence/joint_collection_v15/status.json`.
The independent collector already watches lifetime_joint45_n5_v1 and will archive
either terminal result under joint_collection_v15/lifetime_joint45_n5_v1. Do not
relaunch or overwrite. Retrieve every terminal part and independently replay joint
controls and nine fixed records before reporting any joint-model rates; the
sequential-n5 auditor's sixteen-fixed-record scope does not apply to this new model.

EXPOSED_NONBLIND; scientific HOLD; adoption false; accepted lifetime, uncertainty,
covariance and calibrated sigma null. Full relaxation, independently implemented
relativistic validation, calibrated accuracy and known-grade/blind tests remain
owed. No canonical line/gf/product/reference/star, source spectrum, frozen gold,
grade, publication, merge or Done transition.
