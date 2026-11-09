# RYA-1354: n4 multireference and orbital-relaxation continuation

Start: b04ed3d71fb8220d9858cea008adefd7332136cd. Original task base:
adefcaf073ff70035954cf45a9a2283774ac7d06. Refreshed main:
078cf7dedfa961452c8bd9af43c4bf83c5463fdf (independent Al/Si landing).
Worktree/branch and draft PR616 remain the existing RYA-1354 research branch.
No production products, adopted atomic values, grades or stellar results change.
All comparisons remain exposed/nonblind. Adoption and calibrated uncertainty HOLD.

Ryan authorized continuing toward a validated theoretical lab and later tests on
other elements, and requested ongoing ticket updates. Linear records the plan,
launch, native inventory, intermediate findings and failures chronologically.
RYA-1355 retains its oxygen-pilot/grade-review dependency; this continuation does
not declare that gate passed. The prior n5/core results remain in
[lifetime_continuation_v4.md](lifetime_continuation_v4.md).

## Experiments declared before their results

`config/atomic/rya1354/lifetime_plan.v5.json` owns these exploratory choices.

1. Extend the n3 state-specific core-orbital relaxation to n4. Initialize from the
   completed n4 full2s-SD/DC states; reoptimize n2/n3/n4, fix1s. All three states
   receive equal radial scope. Same439grid/H0.03125/originalACCY5.96046448e-8,
   damping0.5,1500cycles, native strict orbital stopping and freshNSIC0 restart
   before any lifetime. Every failed stage is retained. Same-CSF continuation is
   explicit, with byte-identical CSF/nucleus required; no automatic fallback.
2. Generate a bounded n4 multireference inventory before diagonalization. The
   reference families are an explicit n<=4 subset of Li2023 Table1 MCDHF families:
   even2s2(2p4,2p3 3p,2p3 4p,2p3 4f); odd2s2 2p3(3s,4s,3d,4d), all1s2.
   Compare SR/MR and VV/full2s-SD. These are reduced models, not reproduction of
   the published much larger calculation. No observed energy/strength is input.
3. Using unchanged qualified n4 VV orbitals, compute identical SR/DC control,
   VV-MR/DC, full2s-MR/DC and identical core-MR CSFs plus static-limit Breit.
   Generate each state's proper J blocks, check physical labels/downward ordering,
   retain all13 decay channels. This isolates configuration/Hamiltonian effects;
   MR orbitals are not reoptimized in this experiment. Relaxation remains separate.

Sirius is the native runner, single-thread/nice19. Initial resources8.0GiB available,
133GiB free data disk; before MR launch9.5GiB available. Existing workers are untouched.
No GitHub-CI-driven scientific space/resource cap is introduced.

## Completed native inventory

| Parity/scope | Single reference CSFs | Multireference CSFs | Added CSFs |
|---|---:|---:|---:|
| even VV | 801 | 4331 | 3530 |
| even full2s | 3459 | 29185 | 25726 |
| odd J1,2 VV | 1159 | 3722 | 2563 |
| odd J1,2 full2s | 6337 | 25382 | 19045 |

Every MR set is a strict superset of its matching SR set, checked from actual
native three-line CSF records. Each record has8electrons, closed1s, expected2s
constraints and parity; no duplicate CSFs. Largest MR core blocks are14432even
and14490odd. This inventory alone produces no atomic prediction or grade.

`data/audit/rya1354/mr_inventory_v5` contains59native paths in one360355-byte
archive part. The manifest checks part/member hashes and original receipts.
`source_curation_v2/mr_summary_v5.json` replays it deterministically. Its optional full-campaign
lifetime fields remain null until terminal native campaigns are separately
archived and verified. A separate completed-subcase field reports only SR/DC and
VV-MR/DC, from237 stable native paths in three6442873-byte-total parts in
`mr_partial_v5`. This archive explicitly records a posthoc completed-subcase
checkpoint at18:25:43UTC, not an original final campaign receipt. All checkpoint
file hashes/bytes and explicit omitted cache declarations replay successfully. `rya1354_mr_summary.py` refuses incomplete campaigns,
requires fixed-wave/nucleus identity, identical core-MR DC/Breit CSFs and exact
prior SR/DC reproduction before admitting an exploratory MR lifetime.

## Native commands and monitoring

Run from the existing task-owned Sirius workspace resolved from the work root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_lifetime_relax_v5.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --core runs/lifetime_core_n4_v1 --n 4 --damping 0.5 --cycles 1500 --out runs/lifetime_relax_n4_v1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_mr_inventory.py --grasp-bin toolchains/grasp/build/bin --out runs/mr_inventory_n4_v1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_lifetime_core_v5.py --grasp-bin toolchains/grasp/build/bin --ladder runs/lifetime_ladder_v3 --n 4 --multireference --out runs/lifetime_mr_n4_v1
python3 rya1354_pack_native.py --root . --out evidence/mr_inventory_v5 runs/mr_inventory_n4_v1
```

Version-suffixed remote drivers are byte-identical copies of the versioned local
relax/core drivers, preserving previous remote drivers and all output directories.
Each native stage retains its exact stdin/stdout/stderr and executable hash.
The initial relaxation driver PID is2493107; initial ground rmcdhf PID2493136.
PIDs identify this launch only; check live process arguments, never kill by a stale
PID. Progress is in each run's native logs and the MR `progress.json`; final
`receipt.json` must be inspected before accepting a completed campaign.

First completed MR subcases: SR/DC B/C2.500065889/2.180712938ns exactly reproduces
prior n4. VV-MR/DC B/C2.459565977/2.160236577ns. These are complete13-channel
subcases within an ongoing larger campaign, not its final core-MR result.
The modest shift does not establish convergence or calibrated accuracy.

## Cross-element source preparation and source-access failures

`source_curation_v2/cross_element_sources_v5.json` records the Li lifetime source
lead and its explicit scope/uncertainty/access limits. Publisher:
[McAlexander et al.1996](https://doi.org/10.1103/PhysRevA.54.R5).
The existing Li/P/S shortlist remains a curation shortlist, not a blind experiment
or a certified existing grade. No numerical benchmark is scored here.

Institutional full-text extraction could be inspected, but raw PDF download
returned HTTP202/zero bytes with a WAF challenge; screenshot retrieval also failed.
No preserved-PDF hash or distinct-curator receipt is invented. A guessed author
page path was inaccessible, and the author host did not resolve. Oxygen Golly2003
and Baclawski2008 remain relative-strength/source leads without inspected original
numeric tables. No relative multiplet strength becomes an independent absolute gf.
Mac `pdftotext` was absent; existing Python/PyMuPDF inspected the already-held
Li2023 methods/Table1 without installing dependencies. No source spectrum or
frozen gold was touched. Source-page audit still fails before traversal on existing
civis2013/mallinson2024 vocabulary; no successful regeneration is claimed.

## Verification and continuation

New tests reject wrong-layer or tampered native core initializers and prove n-aware
radial scope, actual-native CSF parity/electron/duplicate checks, strict MR set
inclusion, receipt-count tampering and inventory-only null predictions. Existing
lifetime/regression tests continue to cover original n3 behavior and old summaries.
The informational refinement-debt report remains21rows/10not-fired/5no-ticket;
none discharged. No new numerical defect was discovered by the completed inventory;
new solver outcomes remain pending until inspected. Final focused regressions passed **77tests in32.52s**, including completed-subcase
archive replay, wrong-wave/wrong-MR-space refusal and failed-campaign collection.
Mac Python3.12.13 runs metadata/arithmetic/native-output replay only; native
solvers remain on Sirius Python3.14.4 with the pinned toolchains from v4.

Local replay:

```sh
python scripts/rya1354_mr_summary.py
python -m pytest tests/test_rya1354_lifetime_path.py tests/test_rya1354_strict_scf.py tests/test_rya1354_decay_inventory.py tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

Exact next action: inspect terminal receipts for `lifetime_mr_n4_v1` and
`lifetime_relax_n4_v1`, preserving failures. Archive only stable completed outputs
with explicit `--omit-rci-restart-cache` policy where necessary; preserve cache
hashes/bytes on Sirius. Run optional `rya1354_mr_summary.py --mr-bundle ...
--relax-bundle ...` only for verified complete campaigns. Failed relaxation cannot
supply an accepted lifetime; declare any continuation and retain the same criterion.
Continue larger MR/AS and targeted orbital strategies, independently implemented
relativistic checks, and primary experimental dependence/covariance review. No
choice of gauge/model based on matching experiment; no production adoption,
publication, merge, calibrated theoretical grade or Done status from this checkpoint.


### Active-work checkpoint

At18:31UTC, n4 ground relaxation was at printediteration64, still unqualified;
core-MR/DC ground finished and the quintet stage existed. The core-MR campaign
and separate relaxation were both still running. Their pending status is explicitly
preserved in `source_curation_v2/pending_collection_snapshot_v5.json`.

Detached collector PID2500224 runs the tested `rya1354_collect_pending.py --root .
--out evidence/pending_collection_v5 --watch` on Sirius. It checks every30seconds,
retains success/FAILED_RETAINED distinctly, refuses changed terminal receipts,
and losslessly packages terminal task outputs into size-bounded parts. Heartbeat:
`evidence/pending_collection_v5/status.json`; its stdout/stderr are the sibling
`.stdout`/`.stderr` files. It launches no physics, changes no thresholds, makes no
scientific acceptance decision and sends no messages. A missing/incomplete receipt
remains pending. If a solver is externally killed before its finally receipt,
inspect the retained outputs rather than manufacturing a terminal result.

The collector's generated bundles still need physical/model review and local replay.
Read its latest heartbeat and both native receipts first on resume; copy stable
bundles, run the MR summary optional arguments only for qualified complete campaigns,
and post all findings/failures to RYA-1354. The active calculations have not been
stopped or mislabeled complete at this checkpoint.


Final self-review: native model identity, fixed-wave/nucleus preservation,
MR-CSF equality to the separately generated inventory, full13-channel accounting,
strict stopping/fresh-restart admission and null uncertainty/grade fields checked.
No unresolved wrong-data/provenance finding in the completed subcases; no new
run-bug-ledger row warranted. Original GRASP source worktree remains clean.
Production-path diff against original task base exits0. Whitespace check passes.
Register freshness passes with the existing informational stale-tracker notice;
post-commit PR-mode check remains required. Push uses the previously documented
draft --no-verify pre-push exception for unrelated merged RYA-1352/1353 sequence
entries; no file-size-hook bypass. Each new archive part is below5MiB.
