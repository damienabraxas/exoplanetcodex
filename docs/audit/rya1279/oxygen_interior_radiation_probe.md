# Bounded actual interior-radiation probe

RYA-1371, basecc98a92d, worktreecodex/rya1371, branch
ryandamienschmitt/rya-1371-oxygen-reference-provider; main4cdc449a.
Latest Linear handoff checked; scientific HOLD, draftPR620. Repository intake,
briefing, scientific-review and state-register procedures continue this session.

## Predeclared selection and run contract

Use the existing matched joint coarse17263/fine19963actual radiation fields.
For EACH22BF/56depth/upward-downward term, compare actual fine J against the
coarse field's linear interpolation at fine nodes. Compute positive absolute
field-difference integrands using actual qualified source Hermite alpha and
held-LW thermal/LTE factors. Spontaneous terms cancel in this field difference.
Integrate each fine panel by trapezoid, normalize each direction/depth by its
full fine compiled BF rate, and maximize over directions/depths/continua.
This is a source/interpolant indicator, NOT native-residual or abundance fitting.
Stable score ties follow wavelength order. Select8highest positive fine panels;
add4interior Gaussian nodes each,32total. Preserve all19963ancestors/endpoints,
all41levels/256BB/22BF and original20000point cap:19995nodes,5slots unused.
No larger packet, point removal, wavelength tiling or tolerance/resource change.

First test selection, ties, all-node retention, positive quadrature weights,
invalid data/budgets and cap refusal. Dry producer binds actual parent plans,
radiation/rate files, native-state receipts, complete physics controls and source
oracle. Loader revalidates those identities, actual native/declared LTE states
and exact selector/grid/weight reconstruction before context allocation. Fresh
baseline must reproduce all8previous fine joint arrays exactly. Then run the
19995node diagnostic under unchanged one-thread600s controls; actual n/nStar/C,
background emission, source coefficients, alpha and boundary conventions remain.

Compare each new actual J against original fine linear interpolation at all56
depths, report all32samples and inherited-node differences; no selected depth or
weak level removal. Independently reconstruct total chi/thermal eta and all
radiative coefficients from actual saved J under prior arithmetic screens.
Any Gaussian weights used later for a sampled-panel field-response estimate
are a shadow diagnostic only; compiled quadrature and coupled runner unchanged.
No numerical convergence, between-node error bound, absolute EOS/partition,
observed solar or model admission follows from32new samples.

## Executed milestone (v298, 2026-11-16)

Dry selection completed0.669645s. Selected parent intervals2792,2793,2737,2736,
2787,2786,1932,1933 in descending source-indicator order; exact ties retain
wavelength order. All32new points lie105.6409007--130.4697080nm VACUUM. The
selection was fixed before the new RT calculation. No native residual enters
selection. Actual grid19995 retains every19963ancestor and all endpoints;
5slots remain under original20000cap, no tiling or point removal.

Fresh baseline reproduces all8preceding joint rate/state/J arrays exactly;
its total chi/thermal-eta/background arrays also reproduce the saved parent
exactly. Baseline6.59140s and probe6.64177s, each9radiation iterations, excluding
preflight. Original one-thread600s controls, actual n/nStar/C, native reference
R/C, alpha, widths, emission/ghost and source forward rates unchanged. No SE
update, backend build, new memory measurement or production/default change.

### Direct radiation measurement

Actual new J differs from parent linear interpolation by up to
**0.0254532324906587 relative (2.545323percent)** at130.231096151674nm VACUUM,
top-depth0, positive actual/interpolated difference. Selected-sample RMS relative
0.00159192366 is a descriptive statistic over32x56cells, not an error bar or
global sampling uncertainty. Per-selected-panel max differences are
6.51705e-6,0.0150881,0.0248653,2.38031e-6,1.08830e-6,0.0254532,
1.57886e-4,1.57005e-4. All32samples/all56depths are retained in CSV.

Every inherited J value is bit-identical. Inherited total chi/thermal eta and
background chi/eta/scattering are also exact. This isolates wavelength sampling
under the same prescribed physical state. It directly demonstrates that the
small Gaussian-order differences of the *assumed* linear-J interpolant did not
establish actual interior-field accuracy. These are selected source-indicator
panels, not a global UV validation or comparison with observed UV spectra.

Full compiled rate-grid refinement max1.99807917e-4. Overall physical balance
residual only changes2.72286197895160e-4 ->2.72292163791942e-4, both level41/top0;
a nearly unchanged global maximum does not certify local J convergence.
Independent reconstruction of all256BB/22BF/31136positive directional-depth
rates matches3.07532e-14 relative; total chi/thermal eta matches4.44089e-16
across every wavelength/depth cell. Existing arithmetic screens retained.

A separate Gauss4-weighted(actualJ-linearJ) estimate on ONLY the8sampled parent
panels changes some BF coefficients by2.14850265e-4 of original rates; spontaneous
terms cancel, BB/C untouched. Whole-graph shadow residual2.66945671826650e-4 is
recorded without adoption or fit selection. It is not a complete integration,
global error bound, replacement solver operator or unique attribution of the
compiled-grid rate change, which also changes quadrature weights. Every measured
and shadow outcome remains separate.

## Failure, tests and provenance

The first DRY preparation failed BEFORE output creation: require_codex_path
incorrectly demanded that the new output resource already exist. Actual failed
driver bytes and an explicitly authored note are preserved. Corrected write-only
destination validation uses codex_path and refuses mismatched or existing targets;
all READ inputs still require existing registered paths. No dummy input/fallback
or source guard weakening. No native/RT call on failed preflight. Canonical
run-bug ledger records the fixed path-resolution defect; successful preparation
and new overwrite/location-refusal tests verify it.

143issue-focused tests pass2.99s (four new selector/semantic-loader/output tests).
Loader tests reject removed ancestors even when packet checksums are refreshed,
changed parent physics and changed source provenance. Initial focused8tests
pass0.14s, output-fix4tests pass0.13s. All actual RT/audit/comparison/archive calls
exit0. Existing LW default-configuration warning remains visible; no benchmark,
install, download, Sirius/hostedCI or observational/product publication.

Evidence `data/results/rya1371/interior_radiation_probe_v1`:526raw receipts /
48compact exact copies plus archive receipt, unchanged5MiBcap. Immutable dry
probe packet/plan and producing source snapshots retained; comparison component
snapshots explicitly after-run and checked against actual producing receipts.
Raw arrays/private alpha-oracle requests remain external receipt-only. Native
source input tables untouched. Failed preflight is retained separately from
successful plan and actual runs. Config probe plan byte-matches its generated
probe_plan.json; no manually invented input manifest.

### Reproduce

```sh
source config/hosts/mac.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LW="$CODEX_BRIDE_EVAL/.venv/bin/python"
RAW="$CODEX_WORK_ROOT/bride_mac_checks"
"$LW" scripts/prepare_oxygen_radiation_probe.py --fine-plan config/experiments/rya1371_hop_joint_lte.json --coarse-plan config/experiments/rya1371_hop_joint_lte_midpoint.json --fine "$RAW/rya1371-joint-lte-replay-v1" --coarse "$RAW/rya1371-joint-lte-midpoint-v1" --resource experiments.joint_radiation_probe_grid --output "$RAW/rya1371-radiation-probe-plan-v1"
"$LW" scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_joint_lte.json --output "$RAW/rya1371-radiation-probe-baseline-v1"
"$LW" scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_interior_probe.json --output "$RAW/rya1371-interior-radiation-probe-v1"
"$LW" scripts/audit_oxygen_integrated_thermal_rates.py --plan config/experiments/rya1371_hop_interior_probe.json --replay "$RAW/rya1371-interior-radiation-probe-v1" --output "$RAW/rya1371-interior-radiation-probe-audit-v1" --absolute-replay --profile-contract --coefficient-replay
"$LW" scripts/compare_oxygen_radiation_probe.py --plan config/experiments/rya1371_hop_interior_probe.json --grid "$RAW/rya1371-radiation-probe-plan-v1" --baseline "$RAW/rya1371-radiation-probe-baseline-v1" --probe "$RAW/rya1371-interior-radiation-probe-v1" --audit "$RAW/rya1371-interior-radiation-probe-audit-v1" --output "$RAW/rya1371-interior-radiation-comparison-v1"
"$LW" scripts/archive_oxygen_reference_checks.py --raw-root "$RAW" --output data/results/rya1371/interior_radiation_probe_v1 --radiation-probe-supplement
"$LW" -m pytest -q tests/test_*rya1371.py
```

Fresh output versions required. Re-register a new versioned write resource for
new preparation; existing successful paths are never overwritten. Preserve the
generated config/source-bound packet identities before launch. Predeclaration
97d4da27, output-fix/generated-plan4f35bbb1. Basecc98a92d; existing draftPR620.

Self-review: no source state/data substitution or residual-selected samples;
complete source graph, all ancestor nodes and original cap preserved. Semantic
loader and immutable-output refusals tested. Actual coefficient/rate reconstruction
and unchanged inherited arrays isolate this numerical probe; selected panels do
not establish global convergence. Absolute EOS/partition, finite-CRD/BB integration,
coupled action/conservation/convergence/restart and canonical-solar observed AIR/
intensity-flux/continuum/telluric/fullband/3D remain OPEN, Fe-next/CN/1374 unchanged.

Next bounded step may use one final4point nested probe under the5remaining slots,
selected from actual rate-normalized probe-field indicators, while retaining all
19995nodes. A new mesh/compute strategy needs an explicit scope and its own
normalization/coupling validation; do not remove ancestors, silently tile or
raise the packet cap to claim convergence. Absolute EOS/partition qualification can proceed independently of this sampling
limit.

Verification close: all25834raw receipts/1181compact copies acrossRYA1371
verify by hash/length/exact bytes. Staged generator error SET unchanged5120,
zero added/removed. New registered probe resource now exists:245paths present,
18pre-existing broken. Fresh library audit exits1 with215documents/132indexed/
83unindexed/12absent bibliography paths, unchanged. Registry pointer/params mirror,
stagedwhitespace/reconciler pass; tracker/gold notice retained. Bug-ledger schema
valid113rows; new resolver defect fixed at4f35bbb1 and regression-tested.
Actual CLI overwrite and conflicting archive flags refuse before output, with
retained packet/report hashes unchanged. No global green claim. Generated config
and original probe_plan.json match byte for byte; all other worker files retained.
