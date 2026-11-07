# Coupled oxygen response to the FeII neighbor

2026-11-07, RYA1281/1278, scientific HOLD. Historical worktree
`/Users/ryanschmitt/codex/rya1279`, branch
`ryandamienschmitt/rya-1279-1280-bride-oxygen-gates`, stacked base4635c742.
Starting evidence d7b834bf; latest ORCH pulled separately at main de96ae51 and
all48native continuum comparisons already exact. No other worktree changes.

## Predeclared scope

`config/experiments/rya1281_feii_coupled.json`: canonical column only,
three modes natural/cowley/provided, identical saved O populations as initial
guess for every leg. Fresh O SE plus radiation, unchanged atmosphere/EOS/ne/H/He,
bound background and41FeI. No charge/temperature/hydrodynamic feedback. Original
600s/1500iterations/1e-6outer convergence; original5rays/qscale2/7flux rays and
1e-8inner radiation. No compute extension or tolerance relaxation. Serial Mac,
one thread; pre-run load2.52/2.63/2.44. No Sirius/CI or priority claim.

Use unchanged source-matched collision spline and author-MN receipts. Require
joint population/radiation/cap/physical-rate convergence, then evaluate individual
collision classes afresh and compare to native C; independently reconstruct
transition R plus class C net/gross flows and conservation without accelerated
Gamma. Residual must remain<1e-6, conservation<=1e-10. Fixed O transfer is not
accepted as a new SE result. Each new final O state gets its own standalone
continuum, rather than reusing an old continuum after changing populations.

Source/canonical/main target, native widths and pinned source receipts, atmosphere,
saved populations, inherited provider/cap/solver and binding identities checked.
Natural and both ionic Stark prescriptions retained without fit/preference.
Preserve old outputs. Populations/rates/J/history stay in immutable registered
Mac store; no binary rename or5MiBguard bypass. Compact reports/flux/history can
be copied with actual generators. Downstream continuum uses verified latest1315
implementation; both references/RVs retained when comparing spectra.

No new adopted abundance/opacity/default/product, production/merge/publication,
full3D or all-band claim. Broader source/model/EOS/reference and missing neighboring
opacity gates remain OPEN. Precommit plan/runner before observational evaluation.

Exact runner: held Bride interpreter after `source config/hosts/mac.env`, with
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`,
`python scripts/run_bride_feii_coupled.py --case canonical_refined --mode natural --output <fresh-directory>/natural`;
repeat serially for cowley and provided. Do not reuse final populations between modes.

## Accepted coupled runs

Predeclared code/plan8f7e84c2; all3commands exit0, no retry, budget extension or
tolerance change. Every run converges after3outer iterations. Histories retain
joint population/J/cap/rate/conservation checks; fresh radiation and independent
physical-flow checks pass. Native C is freshly evaluated; separately evaluated
CE/CI/CH/CH0 class sum matches it within256machine eps. Transition R plus class C
flows do not use accelerated Gamma. All56native depths and41O levels retained.

| Mode | Fresh independent R+C residual | Conservation error | Time(s) | Max normalized-flux change vs matched frozen |
|---|---:|---:|---:|---:|
| Natural |3.5646523443e-7|6.6613e-16|32.8862|1.0284643e-8|
| W+Cowley |3.5646516515e-7|6.6613e-16|32.5462|1.1244467e-8|
| W+supplied |3.5646516526e-7|6.6613e-16|32.8800|1.1244931e-8|

Identical initial arrays verified across modes. Max relative population movement
from original initial state is2.15105e-6 in each run; this includes the common
three further convergence steps. Width-versus-natural final population differences
are2.57015e-8Cowley/2.57740e-8supplied over all levels/depths. Those differences
are below the retained1e-6joint stopping tolerance; they do **not** establish a
physical no-feedback limit or tighter abundance accuracy. Final fresh independent
residual~3.565e-7 is deliberately reported alongside the tiny spectral changes.
Original opacity, source atom/collision hypotheses and EOS validity are not
validated by this restricted-model rate balance.

Each state gets a freshly calculated standalone continuum. Continuum changes
relative to the corresponding frozen continuum are<=2.5e-14, not forcibly zero.
O nStar and particle totals stay unchanged; the atmosphere/ne/H/He/Fe background
is fixed. This is radiation plus oxygen SE, **not** charge/EOS/temperature or
hydrodynamic feedback. No adjustment to abundances, widths, gf or source inputs.
Peak process memory1.358GB, one-thread jobs strictly serial; all jobs finished.

## Current-ORCH spectrum/continuum analysis

`scripts/analyze_bride_feii_coupled.py --source <registered-feii-coupled> --output
<fresh-registered-feii-coupled-analysis>` exits0 under one-thread project Mac
interpreter. Actual latest1315/native continuum is loaded from the pinned full
checkout de96ae51, source bytes checked against Git, model/envelope bounds checked
(.15/.03). Kitt product, actual IAG reader/path and source hashes checked. Same
canonical broadening, both fixed RV conventions(no-gravity/plus-gravity), both
references, all3O lines.36paired diagnostic rows;72current-native fits all applied.
Repeated reference/RV views are correlated, not36independent validations.

| Mode | Max absolute model O-window change(mA) | Max observed continuum-guided area change(mA) |
|---|---:|---:|
| Natural |2.3851492e-6|2.7826019e-8|
| W+Cowley |2.6505791e-6|3.0051751e-8|
| W+supplied |2.6508785e-6|3.0053158e-8|

Changes are **coupled minus matched frozen** for each FeII prescription, not the
whole added-line effect. IAG-guided observed areas remain exactly unchanged;
model residuals include the tiny flux changes. Synthesis-reference changes are
tiny too. The previous~8.8-9.0mA reference spread and partial-opacity/structured
shoulder/one-sided-anchor issues remain. No reference preferred and no fitted
continuum adopted. Window integrals are blended diagnostic areas, not isolated
O EWs. AIR comparison/internal vacuum and provider-corrected telluric handling
retained, no second correction or source conditioning.

## Evidence and validation

Registered immutable resources `experiments.bride_feii_coupled` and
`experiments.bride_feii_coupled_analysis`. Source runs/logs/populations/full R/J,
histories and independent rates remain there.20compact result copies at
`data/results/rya1281/feii_coupled_v1` byte-match originals; largest888297bytes.
The original~6.1MiB`rates.npz` files stay external with their original names and
hashes in each coupled report; no rename, exemption or5MiBguard bypass. Independent
rate arrays, populations and3-step histories are included in review copies.
Actual writer/invocation registered. No raw observational/source data copied or
prior result overwritten.

Preflight11tests pass3.87s. Final focused suite23passed/3existing dependency
warnings4.52s, exit0: FeII natural/width providers, physical-rate accounting,
charge contract/native spline, continuum/crossed controls and solar context.
Exact log `/private/tmp/rya1281-feii-coupled-final-tests.log`. Commands use held
Bride interpreter/one thread. Independent runtime physics/source checks are
additional evidence, not counted as extra unit tests. Full earlier87suite not
rerun or claimed anew; no changed shared solver kernel. No new bug or silent
fallback discovered; prior ledger dispositions retained.

Latest origin fetch unchanged at main de96ae51. Fresh library audit retains
195documents/127accounted/68unindexed/12off-Mac declarations, exit1. No new claim
depends on unidentified literature. Existing global-ledger RYA587blank science
impact/tracker/gold/refinement debts retained. Statev221/SEQUENCE record this
restricted coupled check, not scientific admission. No Sirius/CI/benchmark,
installation/download or nice-priority claim.

Self-review checked same initial arrays, source/native identity, active versus
fixed contexts, neutral-collider cache invariants, unchanged nStar/totals,
independent physical C+R rather than accelerated Gamma, fresh continuum,
two separate Stark modes, both references/RVs, currentORCH source and numerical
tolerance limitations. No new blocking implementation finding. Scientific HOLD.

Artifact guards: generator check exit1 with exactly2436existing diagnostic
strings, zero added/removed versus preceding phase; no new registration or
invocation failures. MIRROR targets, version pointer, main-landing reconciliation
and whitespace checks pass. Postcommit freshness and standard pre-push guards
run before handoff; no bypass. Historical global ledger1failed/1passed remains
documented from prior check, not rerun or represented as green here.

## Next bounded direction

Prioritize missing neighboring opacity and continuum information, plus other
oxygen indicators, rather than fitting a width/abundance to close the remaining
reference spread. A tighter paired convergence study is required before claiming
feedback accuracy below the current1e-6floor; do not rebrand these tiny differences
as resolved physics. Additional O diagnostics and qualified solar C/N precede Fe
abundance/reference star; broader xi/T/g/composition interactions, full source/EOS/
independent-reference closure and actual adaptive3D/all-band grids remain OPEN.
RYA1281 In Progress/RYA1278 In Review, draftPR598; no adoption/default/product/
production/merge/publication change.
