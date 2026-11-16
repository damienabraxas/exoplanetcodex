# Full oxygen LTE population and inverse-input milestone

RYA-1371 continuation from53055422. Explicitly read-only prerequisite for a
common active thermal convention; no coupled solve or physical admission.
Earlier source reconstruction was same-stage only. This audit covers all41
levels, both OI/OII stages and56depths, including the disconnected excited OII
level41. Source `atom_mod.lte_pop` consumes FP32 e_kT, lsaha and Debye lowering,
with FP64 consumed level energies/statistical weights and exponentiation.
`composition_mod.set_aux_data` owns lsaha and lowering. Source constants and
actual DEBYE_SHIELD/LEVEL_PARTF flags are pinned. No parameter is inferred by
fitting nStar. The1e-12 arithmetic screen is inherited from the prior same-stage
audit and declared before this execution; it is not a physical error budget.

Independently reconstruct all normalized source populations and ground ratios;
compare a separate held-LW log formula against actual held `lte_pops`. Retain
native FP32, native FP64, no-shield negative control and held-LW alternatives.
Keep the actual total oxygen density at every depth. This qualifies finite-atom
shape ONLY: absolute element density, external partition-function normalization
(LEVEL_PARTF=F), chemistry and EOS remain open. No substituted LTE population
is supplied to a solver.

Audit fixed collision reciprocity under each diagnostic population, and a
clearly algebraic reverse-ratio projection preserving source upward coefficients.
The latter tests ownership coupling, not independent microscopic inverse physics
and never becomes solver input. No missing level41 collisions are invented.
All prior files, model/production defaults, observations/gold, source tables,
resource constraints and scientific HOLD remain unchanged.

Before v2 execution: extend the same read-only audit to actual candidate BB
line-centre thermal identities under held-LW Einstein/Planck conventions, and
the inverse nStar factors for every256BB/22BF source link. Retain v1. These
canaries distinguish population ownership; they do not assert finite-profile
Kirchhoff equilibrium or qualify an opacity/rate replacement.

## Milestone achieved: complete finite-atom shape and inverse-input ownership

Actual source-shape reconstruction covers41levels/both OI-OII stages/56depths:
maximum population relative difference1.0769163338864e-14 and ground-ratio
difference1.0325074129014e-14. Independent held-LW log reconstruction agrees
with actual held `atomic_set.lte_pops` to3.5527136788005e-15. Both are below the
predeclared1e-12 arithmetic screen. There is no fitted Saha/ionization parameter.
Source Debye lowering spans9.254e-5--0.0229416eV; FP32 lsaha13.5795--24.1649.

| Read-only population convention | Max relative population change from native | Fixed-C reciprocity defect |
|---|---:|---:|
| Source native FP32 | 1.07692e-14 | 7.39760e-15 |
| Same constants/energies, FP64 auxiliaries | 2.79548e-6 | 1.21455e-6 |
| Source native, shielding disabled (negative control) | 0.0260153 | 0.0136240 |
| Held-LW including its Debye lowering | 0.000612918 | 0.000249136 |

All43680collision pair-depth edges remain included; level41still has no
collision edges. Algebraic inverse projection preserves every upward
coefficient and restores reciprocity below1.10e-16 for all candidates, as its
construction requires. This proves the ownership dependency; it is not
independent physics validation or permission to change source inverse rates.

Actual held-LW full-LTE populations satisfy all256BB line-centre thermal
canaries to5.49167e-15 under held-LW Einstein/Planck conventions, versus
2.25475e-4 using the source-native populations. Their inverse LTE factors differ
from imported native by up to4.50847e-4 for BB and5.45321e-4 for all22BF links.
This measures input consistency at line centre; it does not eliminate the
finite-CRD-profile thermodynamic residual, establish frequency-integrated
closure, or qualify a modified opacity/rate operator. Background/LTE/collisions
and active emission/inverse terms must be coordinated before another coupled
solver experiment. No candidate is chosen by a residual or abundance.

Crucial remaining limit: measured total oxygen per depth is retained. Actual
LEVEL_PARTF=F uses external partition functions; common normalization cancels
in the reconstructed ratios. Absolute element density, partition/EOS/chemistry
normalization and molecular/observational/solar/fullband validation remain OPEN.
This milestone closes full finite-level LTE *shape* reconstruction only.

## Reproduce and evidence

```sh
source config/hosts/mac.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LW="$CODEX_BRIDE_EVAL/.venv/bin/python"
RAW="$CODEX_WORK_ROOT/bride_mac_checks"
"$LW" scripts/audit_oxygen_full_lte_contract.py --plan config/experiments/rya1371_hop_grid_refined2_native_alpha.json --replay "$RAW/rya1371-native-alpha-joint-v1" --output "$RAW/rya1371-full-lte-contract-v2"
"$LW" scripts/archive_oxygen_reference_checks.py --raw-root "$RAW" --output data/results/rya1371/full_lte_v1 --full-lte-supplement
"$LW" -m pytest -q tests/test_*rya1371.py
```

Fresh output paths required. v1 source-shape audit is retained. First v2
extension failed before output creation: builtin abs(list) in radiative summary
raised TypeError. Failed driver bytes and an explicitly authored failure note
were preserved before fixing np.abs aggregation; successful v2 retains actual
producer source snapshots. No RT/SE, input write or scientific output from the
failure. Canonical run-bug ledger records this reporting defect as fixed.
Diagnostic NPZs remain external, receipt-only; compact CSV/report/source copies
are exact. `full_lte_v1` adds12raw receipts/10compact copies plus archive receipt,
under unchanged5MiBcap.131issue-focused tests pass2.57s (four new tests), including
independent two-stage analytic Saha/normalization, negative shielding/precision
controls, invalid-state refusals and actual-library parity. No new physics run.

Worktree codex/rya1371 / branch ryandamienschmitt/rya-1371-oxygen-reference-provider,
base53055422, predeclaration55f5df8e; draftPR620. Main refreshed to4cdc449a/PR625:
RYA728Cr theoretical-line-set/result additions and line-set/grade descriptors;
not imported or exercised by these pinned native-atom diagnostics. Existing
production products untouched. Earlier83unindexed-library/18broken-path/global
generator debt remains separately reported. No installs/benchmarks/Sirius/CI.

Self-review: ratio arithmetic and units tested; source auxiliaries calculated
from actual consumed T/ne and source constants/flags, never inverted from nStar;
all source levels/transitions retained; normalization limitation explicit;
no solver population/collision or default change. Fixed-C failure and algebraic
projection are separate diagnostic statements. Scientific HOLD remains.

Next exact milestone: qualify a joint fixed-population active-opacity/emissivity
and rate operator with the declared population/inverse convention; first retain
full reproduction and finite-profile CRD controls. Any absolute EOS/partition
normalization claim requires independent corresponding evidence. Then numerical
refinement, a qualified coupled solve/restart and observed canonical-solar
validation remain. Fe-next atom/resource/FeII and C/N/1374 gates unchanged.

Verification close: all17432raw receipts/1081compact copies across RYA1371
verify by length/hash/exact bytes. Generator diagnostic set unchanged5120errors,
zero added/removed. Registerv295 pointer, params mirror, staged whitespace,
reconciler and freshness pass; existing tracker/gold notice remains. Bug-ledger
CSV schema intact112rows, including the new fixed reporting defect. No other
worker files or source/production inputs changed. Existing draftPR620 receives
this checkpoint; no merge or publication.
