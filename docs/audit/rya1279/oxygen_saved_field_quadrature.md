# Joint oxygen saved-field quadrature budget

RYA-1371 continuation from00579232, branch
ryandamienschmitt/rya-1371-oxygen-reference-provider, worktreecodex/rya1371,
main4cdc449a, draftPR620. Scientific HOLD; latest Linear handoff unchanged.
Task-intake/implementation-brief/scientific-review/state-register procedures
continue; no new acquisition, model or observational production run.

## Predeclared bounded milestone

Audit BOTH existing joint17263/19963-node radiation fields. Preserve all41levels,
56depths,22BF/256BB, physical populations/LTE/C and source alpha/thermal factors.
Use physical angle-weighted J_nu linearly interpolated in VACUUM wavelength;
this is a declared saved-field representation, not new radiation evaluations.
Integrate its22BF rates at Gaussian orders4/8/16 on panels split at every saved
J node and raw source-alpha knot. Actual qualified compiled frequency-Hermite
FP32 alpha is queried, not replaced by a Python interpolation. Preserve actual
angular-weight sum in spontaneous terms. Native alpha queries remain bounded
by20000; stream256panels per request (maximum4096nodes). No larger HOP/RT/SE
packet, field/population update, backend change or guard bypass. Original
one-thread600s deadline retained; raw requests private/receipt-only.

Before Gaussian comparisons, source-alpha/explicit wavelength trapezoid must
reconstruct every saved compiled BF directional rate within the existing
16*maxnodes*5rays*2directions FP64 positive-summation screen. Every BF/depth
included. Independent tests verify positive Gaussian weights and exact known
polynomial moments, interpolation/depth axes, no extrapolation and size refusal.
Gaussian agreement is reported numerically, never promoted into physical error
or convergence policy. FP32 alpha rounding remains part of the queried function.

Separate the observed BF refinement into the exact telescoping identity:
Tfine-Tcoarse=(Gfine-Gcoarse)+[(Tfine-Gfine)-(Tcoarse-Gcoarse)]. G16 integrates
each saved-J representation. Report roundoff closure relative to summed positive
coefficients under128eps, plus changes separately; maxima are not additive.
Keep whole-graph shadow matrices with BB untouched for every order/grid. These
are NOT a modified physical RT/SE operator or a selected best-fitting result.
Unknown between-node radiation accuracy, finite-CRD effects, absolute EOS and
coupled/observational admission remain open regardless of Gaussian agreement.

## Executed milestone (v297, 2026-11-16)

Completed33.3056s,2339serial qualified-oracle calls, largest actual query10043
nodes (saved-node reconstruction; Gaussian batches at most4096), all below the
unchanged20000query cap/600s deadline. No new RT/SE/physics run or backend build.
Both41level/56depth/22BF graphs retained. Independent source-alpha trapezoids
reconstruct compiled BF coefficients to3.04201e-14 fine/3.41949e-14 coarse,
under the existing positive-summation arithmetic screens.

| Saved field | Nodes | Max G8/G4 relative difference | Max G16/G8 relative difference | Max G16/trapezoid relative change | Whole-graph G16 shadow residual |
|---|---:|---:|---:|---:|---:|
| Joint fine | 19963 | 3.98946e-9 | 2.14107e-9 | 1.19737574e-4 | 2.98009638342169e-4 (level41/top0) |
| Joint coarse | 17263 | 7.74200e-9 | 3.48760e-9 | 4.78831287e-4 | 3.82634102601162e-4 (level41/depth1) |

All orders4/8/16 retained. Small order differences concern this declared
linear-J/source-Hermite-alpha interpolant, including FP32 alpha quantization;
they do not bound the unknown true radiation between saved samples. Neither
Gaussian order nor physical convention is chosen by native agreement. Shadow
residuals WORSEN from compiled joint2.72286e-4 fine/2.81436e-4 coarse, and remain
inadmissible as a new operator or model result. BB rates and C untouched.

Largest fine/downward coefficient correction is4->41/top0,-1.19737574e-4;
upward max5->40/top0,-7.93769e-5. Coarse counterparts are4->41/top0,-4.78831e-4
and5->40/top0,-3.17449e-4. Excited OII4->41support has1855fine nodes versus1357
coarse nodes. At top0, its upward corrections are-6.82955e-5 fine/-2.73170e-4
coarse. The last depth/weak ion level is not dropped.

### Exact matched-field budget

Inherited coarse J nodes are bit-identical in the fine field; n/nStar/C and
native reference matrices are exact across both saved cases. Full BF coefficient
refinement max3.74341029e-4 (the earlier whole-R max3.75128e-4 includes BB).
G16 fine/coarse represented-field change max1.71579136e-4; change in trapezoid
integration defect normalized by coarse rates max3.59136715e-4. The exact
per-coefficient telescoping identity closes with zero measured normalized
roundoff in this implementation. These maxima occur at different cells and
are NOT added as error bars or unique physical attributions.

This separates a measured numerical integration response from a measured change
of the represented saved field. It does not certify either represented field,
complete radiative quadrature convergence, finite-CRD line integration, full
thermodynamic consistency or a self-consistent coupled solution.

## Reproduction, evidence and review

```sh
source config/hosts/mac.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LW="$CODEX_BRIDE_EVAL/.venv/bin/python"
RAW="$CODEX_WORK_ROOT/bride_mac_checks"
"$LW" scripts/audit_oxygen_saved_field_quadrature.py --plan config/experiments/rya1371_hop_joint_lte.json --coarse-plan config/experiments/rya1371_hop_joint_lte_midpoint.json --fine "$RAW/rya1371-joint-lte-replay-v1" --coarse "$RAW/rya1371-joint-lte-midpoint-v1" --output "$RAW/rya1371-saved-field-quadrature-v1"
"$LW" scripts/archive_oxygen_reference_checks.py --raw-root "$RAW" --output data/results/rya1371/saved_field_quadrature_v1 --saved-field-quadrature-supplement
"$LW" -m pytest -q tests/test_*rya1371.py
```

Fresh outputs required. Raw BF matrices/native oracle requests and outputs remain
external/receipt-only. `saved_field_quadrature_v1` adds7023raw receipts/5compact
exact copies plus archive receipt under unchanged5MiBcap. Actual driver/helper
source snapshots preserved before oracle execution; all inputs checked again
after integration.139issue-focused tests pass2.99s (six new Gaussian/interpolant
controls); focused initial8tests pass0.15s. All actual commands exit0; no failed
run, fallback or new run-bug-ledger defect. Existing held-LW configuration warning
remains visible; no install/benchmark, Sirius/hostedCI or production change.

Predeclaration6d3b04f4, base00579232; worktree/branch and main above. Registerv297
records the budget, not numerical/model admission. Existing83unindexed-library,
18broken-path and global generator/tracker/gold debt remain separately reported;
no acquisition or invented metadata. Self-review confirms all BF/depths, actual
angular normalization, source-alpha ownership, positive quadrature weights,
no extrapolation and bounded requests; both grids/orders retained without tuning.

Next exact scope: define and validate a bounded actual-radiation probe at new
interior wavelengths chosen from source/interpolant error indicators, retaining
all inherited nodes and the20000point cap (the19963grid leaves37slots). Compare
new sampled J with its declared interpolant before using a replacement rate
quadrature. No larger HOP packet, implicit wavelength tiling or coupled update
is licensed by this frozen-field audit. Absolute oxygen/EOS/partition,
coupled action/conservation/convergence/restart, observed canonical-solar AIR/
intensity-flux/continuum/telluric/fullband/3D and Fe-next/CN/1374 gates persist.

Verification close: all25308raw receipts/1133compact copies acrossRYA1371
verify by hash/length/exact bytes. Staged generator diagnostic SET unchanged5120,
zero added/removed; Mac paths244present/18existing broken. Fresh library audit
exits1 with215documents/132indexed/83unindexed and12absent bibliography paths,
unchanged. Held thermal constants exactly match BOTH saved producing input
receipts; conflicting archive supplements refuse before output. Registerv297
pointer/params mirror/stagedwhitespace/reconciler pass. Existing tracker/gold
regeneration notice and historical immutable native-log whitespace persist;
no global green claim. No new run-bug-ledger defect or other worker changes.
