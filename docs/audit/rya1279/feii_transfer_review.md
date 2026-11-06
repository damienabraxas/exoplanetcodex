# RYA-1281: fixed-O natural FeII transfer differential

Predeclared2026-11-05 before evaluating new transfer. Startd19af347, refreshed
main9a74d843 unchanged, stacked base4635c742/draftPR598. Mac only, one thread.
User authorizes continued bounded oxygen controls; no product/adoption/merge.

Manifest config/experiments/rya1281_feii_transfer.json selects gf_079010 after
source and LTE strength qualification. Reconstruct exact source row at runtime,
require both frozen strength screens to pass, compare current main target row,
refuse ambiguity/isotopes/source changes. Add only its natural-Voigt LTE FeII
bound-bound opacity/emission, same transfer with existing41FeI pressure/Cowley
and bound background unchanged. No pressure/Stark guessed for FeII; no full-opacity
claim. No source/width/gf/abundance fit, no O/charge/SE feedback.

Both canonical/lower-g columns use saved O populations, same physical inputs,
rays/qscale/flux quadrature and dense source grid. Empty wrapper must reproduce
existing background chi/eta/sca exactly. Each baseline transfer must reproduce
saved normalized flux within1e-6, the original outer criterion, before interpreting
an extra-line result. No relaxed thresholds or best-case selection. Paired fixed
radiation/cap convergence1e-8, max1000,180s per leg; all old failures retained.

Reuse the exact hash-verified saved standalone continuum Fnu for both legs; no
refit or re-normalization. Both continuum references and both fixed RV controls
retained in post-analysis. Historical3percent continuum source remains frozen;
main's15percent model-guided rule is not exercised or represented by these fits.
All profiles in AIR; transfer internally vacuum; canonical gravitational shift
only in observed comparison. No additional telluric correction.

Two initial real-backend tests pass: empty-wrapper identity, FeII-stage-specific
coefficients/thermal emission, unchanged scattering and outside-support arrays,
deepcopy, cached-field refusal and invalid/unsupported line data. Local-opacity
ratios from the previous screen do not predict flux changes in this test.

## Transfer result and direct/reference separation

Plan/provider/runner committed before evaluation at4970270a. Both paired cases
complete. Empty wrapper chi/eta/sca exactly equals inherited background. Both zero
normalized spectra have max absolute difference0 from saved frozen spectra, and
the entire zero flux CSV fingerprint also equals its frozen source. The originally
declared1e-6guard was not relaxed; observed agreement is exact.

Canonical zero/extra fixed-radiation solves each take65iterations, final J change
8.2973e-9/cap1.1184e-9; lower-g each67iterations, J8.6213e-9/cap1.5218e-9.
These are radiation/cap convergence checks with O populations deliberately fixed,
not fresh joint O statistical-equilibrium or abundance certificates. Saved O and
continuum arrays remain exactly unchanged. No new R+C equilibrium claimed.

The source line is gf_079010, rest AIR7777.099A (vacuum777.923920376nm),
excitation3.199eV/loggf-5.33/natural gamma2.95e8s^-1. Canonical/GES RU-tag gf
retained; physical gf accuracy is not newly qualified. Observed feature appears
near7777.12A after the existing canonical gravitational shift. It changes the
intrinsic wider-region model by at most0.01348195/0.01489724 normalized flux.
Those1.35-1.49%values belong to the added FeII neighborhood, not the O apertures.

Direct O-window model-area changes in mA, with canonical broadening plus fixed
gravitational shift; the no-shift broadened variant is separately retained:

| Frozen case |7771.944|7774.166|7775.388|
|---|---:|---:|---:|
|canonical|0.000004561|0.000014811|0.000048393|
|lower-g|0.000005036|0.000016345|0.000053393|

Across both broadened/RV variants the maximum is0.000054737mA; largest intrinsic
O-window flux change1.1734e-7. Added natural-only wings therefore supply negligible
DIRECT O-window absorption in this declared fixed-O control. This is not a bound
on pressure wings or other source physics omitted from it.

The much larger change is through the observed-continuum reference:

| Case |7775synthesis level before/after|Synthesis-guided observed area before/after(mA)|Delta(mA)|IAG-guided area(mA)|Remaining IAG-synthesis gap(mA)|
|---|---|---|---:|---:|---:|
|canonical|0.992636942/0.994747090|46.776057/48.797953|+2.021896|57.703382|8.905428|
|lower-g|0.992686907/0.994865972|46.824038/48.911594|+2.087556|57.703382|8.791788|

The IAG-reference fitted levels/supports and observed areas are exactly unchanged for every line, as expected with the
same observation/atlas. Synthesis-guided7771/7774changes are only~3-6e-6mA.
No preferred reference, continuum or abundance is selected from improved agreement.
The before/after reference gap is a diagnostic spread, not an uncertainty budget.
Do not infer this fixes the independent O atom/EOS/reference/3D discrepancies.

## Actual anchor replay

Separate trace replays all24native records and areas exactly across2cases x2legs
x2references x3lines, using the same source/hash-checked observation and arrays.
At7775, initial selected sets have74/78symmetric-difference pixels (37/39replaced
out of90); actual coefficient supports90->83/86, symmetric differences73/76.
Actual blue/red support changes31/59->29/54 and30/60->30/56. Both final fits
stop ALL_RETAINED on pass4, reported counts equal actual support. The reference
and selection/clipping both change; no selector-only causal fraction is inferred
without an additional fixed-support/crossed decomposition.

7771/7774synthesis masks do not change. Existing7774fifth-pass73coefficient/70
survivor discrepancy remains characterized; no shared-estimator fix or postclip
refit adopted. IAG selectors/supports unchanged. Historical3percent native source
is checked by hash; these fits do not exercise current main's15percent rule.
Both earlier no-gravity/plus-gravity residual controls and every numerical refusal
are retained. All measured continuum fits in this control pass their historical
numerical bound; that does not sign off physical continuum quality.

## Reproduction, compute and archived evidence

```sh
source config/hosts/mac.env
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/run_bride_feii_transfer.py --case canonical_refined --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-transfer-v1/canonical_refined
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/run_bride_feii_transfer.py --case lower_gravity_refined --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-transfer-v1/lower_gravity_refined
/Users/ryanschmitt/codex/venv/bin/python scripts/analyze_bride_feii_transfer.py --source /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-transfer-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-analysis-v1
/Users/ryanschmitt/codex/venv/bin/python scripts/trace_bride_feii_anchors.py --transfer /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-transfer-v1 --analysis /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-analysis-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-anchor-trace-v1
```

All exit0. Immutable experiments.bride_feii_transfer/analysis/anchors originals
remain; reruns require fresh output directories. Canonical zero/extra legs take
17.261/18.018s; lower-g17.550/18.419s, all below180s per leg. Total case including
setup37.558/37.972s; peak process RSS1.072/1.073GB. Serial one-thread Mac runs,
no nice-priority claim, no Sirius/installation/download/benchmark or CI work.
No original failed output, source/canonical parameter/gf/gold/product altered.

23result artifacts and2visually inspected response figures copied byte-for-byte
under data/results/rya1281/feii_transfer_v1 and data/output/rya1281/feii_transfer_v1.
Largest result507694bytes; no raw external spectral flux copied. The8individual
helper comparison PNGs/logs also remain in the registered external analysis
store; they are not duplicated in the archive. Actual root/leg/helper/trace writer
invocations registered in GENERATORS.yaml, no wrapper masquerading as the writer.

[Canonical response](../../../data/output/rya1281/feii_transfer_v1/canonical_refined_response.png),
[lower-g response](../../../data/output/rya1281/feii_transfer_v1/lower_gravity_refined_response.png).

77tests pass in7.99s,292existing backend printUpdate warnings. Exact expanded
command: previous75-test suite invocation in lte_opacity_screen_review.md
plus tests/test_bride_feii_control_rya1281.py; OMP_NUM_THREADS=1 and
OPENBLAS_NUM_THREADS=1. New real-backend tests prove empty identity, FeII stage1
rather than FeI stage0 coefficients, thermal emission, unchanged scattering and
outside-support arrays, deepcopy, cache refusal and invalid input guards. Runtime
analysis also verifies the O atom hash against the previous coefficient screen
before interpreting outputs. No new bug or silent fallback discovered; previous
characterized/fixed ledger entries retained, no empty placeholder row.

Self-review: source identity/qualifications and main target row rechecked; source
line added once only, native overlap refused, same transfer and inherited FeI
provider unchanged. Fixed O/continuum and exact baseline prevent source/iteration
changes being spent on the line intervention. Max wider-region change cannot be
mislabelled O absorption. Selection/reference effects are retained separately
from direct model-window areas; all references/RVs and historical rule identity
remain explicit. Fixed-radiation convergence is not NLTE physical admission.

Next: compare the new FeII neighborhood itself to both atlases without fitted
gf/RV/width, perform fixed-support/crossed reference accounting, and qualify
source FeII pressure/Stark conventions before physical profile/adoption claims.
Only then assess coupled-O feedback with independent fresh physical rates. Keep
source-qualified missing wings/molecules/continuum and independent O atom/EOS/
reference, wider xi/T/g/composition interactions and additional O diagnostics open.
Previous parameter sweeps retain their historical line subset and are not silently
promoted to this new opacity. C/N then Fe then a reference star, adaptive3D/all
bands remain owed. Scientific HOLD; RYA-1281 In Progress / RYA-1278 In Review.
No native/default/production, abundance/adoption, merge or publication change.

Final repository checks: MIRROR targets/v216 pointer/whitespace/main landing reconciliation pass. Generator diagnostics remain2436exactly, zero added/removed versus prior LTE screen. Fresh library194documents/127accounted/67unindexed/12off-Mac registered paths. Historical RYA-587 ledger/gold/tracker debts remain documented; no green CI claim. No new bug or silent fallback; no empty ledger placeholder. No owned compute jobs remain.
