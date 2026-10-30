# Bride continuation: RYA-1279 / RYA-1280

## 2026-10-28 native LTE / opacity continuation

Actual pinned LTEPOP/GETWRD plus BLOCK DATA compiled on the Mac; normalized
fractions compared across all41 levels/56 depths at identical T/ne. Maximum
relative fraction difference0.061374%; top O II0.049668%. This excludes a large
matched-input formula discrepancy, not producer ne/EOS/normalization differences.
Final immutable audit bride_mac_checks/rya1281-native-lte-v2; v1 retained.

Exact independent MARCS opacity companion acquired through the official model
search/download form after holdings checks and bounded read-only Sirius metadata
queries. No new Sirius science. T/Pe/rho,92 abundances and xi match printed .mod
exactly; max logtau difference5.44513e-5 retained. Source/official reader held at
experiments.bride_marcs_opacity (marcs-opacity-v2). Wavelength medium UNVERIFIED;
no wavelength-aligned opacity comparison, coefficient replacement or SE/RT rerun.
Native MULTI ABSDAT and producer ne/nstar/J/rates still missing.44 focused tests
pass; initial analytic-literal test failure and mistyped test command are explained
in native_lte_opacity_review.md. Four review artifacts registered. Registerv186
and SEQUENCE preserve HOLD; no new verified run bug, existing collision/boundary
rows remain OPEN. Main dbda7faa/PR605 scoped by diff, not merged into this branch.

Standing bibliography audit:188documents/125accounted,63 unindexed; unchanged
pre-existing library debt, not a clean-library claim. Resume with the OPA wavelength
source contract, then independent continuum absorption/scattering comparison with
active O counted explicitly. Verify native producer inputs before attributing the
O II mismatch. Validate column interpolation before capped solar reruns; O then
C/N, Fe, another star. PR598 remains draft; final checkpoint requires[skip ci].

Base: 4635c7421a7dee8671ffefe05d771c68da1cd24e (unmerged parent PR #597).
Worktree: codex/rya1279. Branch: ryandamienschmitt/rya-1279-1280-bride-oxygen-gates.
Ryan authorized continuing through the physics gates toward oxygen results and
subsequently said to proceed without waiting for CI. That supersedes an additional
permission checkpoint for the bounded experiment; scientific validation still applies.
No merge/publication authorization was given.

## Verified work

- Held Gerber `atom.o41f` receipt and MD5 match. This candidate has 41 levels,
  256 lines, 22 continua and 1106 collision tables. All triplet components exist.
  No level/rate reduction was made. Full source remains on Sirius; local small-file
  inspection is evidence only, never a production compute-input fallback.
- Found and fixed a converter defect: only the first TEMP grid was retained;
  976 later rate sets were truncated/misassigned. Each collision now owns its active
  grid: CE 130x7 + 611x10, CH 302x10, CI 39x10, CH0 24x10. All source rates and
  their own grids compare exactly after the documented cm3-to-m3 conversion.
- Found and fixed discarded damping coefficients: triplet packed ABO 453.234 had
  become Unsold 1. Native conversion now keeps sigma=453, alpha=.234, radiative
  gamma=3.89e7. Supported subset is zero/standard Unsold/packed ABO; unvalidated
  negative encodings, nonzero Stark and unsupported profiles/rates fail closed.
- Parser rejects unknown/truncated/extra rate data, bad grids, missing END and
  fixed radiative rates rather than skipping them. No existing Fe result is rerun.
- Synthetic actual-backend checks: continuum slab/thin limit with explicit
  boundaries, 3 formal solvers, depth/ray cases; actual CE/CI+SE vs independently
  computed Saha-Boltzmann with matched constant vintage and Debye lowering;
  alternate initial populations and rate scales; coupled CRD line+SE thermal
  enclosure. Deliberate wrong directional rates converge/conserve but fail the
  independent population test. Physical R+C diagnostics are distinct from MALI
  Gamma, which must not be mistaken for the physical SE operator.

## Numerical experiments retained honestly

Initial roundoff-only slab criterion failed four cases (largest absolute error
1.17e-7 for cubic Bezier). Inspection identified deliberate polynomial branches,
not a source/grid-tuning opportunity. Subsequent acceptance bounds are derived
from alternating-series remainders and rational conditioning in inspected 0.16.1
`LwInternal.hpp`, `FormalScalar.cpp`, `Bezier.hpp`. The tests do not choose a
constant from the observed discrepancies. These are bounded homogeneous controls,
not a universal spectral/abundance error budget.

The first LTE comparison used modern CODATA and omitted default Debye lowering;
that was a mismatched test reference. Matching the declared physics/constants is
required for an independent comparison. The corrected analytic reference shares
only constants with LW, never its LTE/population implementation. A hot synthetic
ion limit has a conservation-augmented condition number about 1.39e8; componentwise
roundoff assertions are inappropriate. Normwise backward error and an explicit
64*eps*kappa forward bound are tested, while component residuals remain reported.
Equation row equilibration preserves the homogeneous system; it changes no rates.
The coupled radiative control needed 51–54 iterations, exceeding an initial 20
iteration diagnostic budget. The bounded cap is now 80; population and J changes
must both reach 1e-11, with independently tested 1e-9 thermal-population accuracy.
No tolerance or physical input was changed toward an abundance target.

Exact vacuum (zero chi, undefined eta/chi source) is not labelled a passed backend
case. Positive-opacity thin-limit cases reach 1e-10 optical depth. Full-3D coupling,
nonuniform formal-transfer accuracy and real-atom resolution convergence are still
separate work; these controls do not certify a full-3D solver or a solar abundance.

## Provenance and baseline still being resolved

The distribution's literature attribution points to Gerber2023 and Bergemann2021.
[Bergemann2021 §3.1/Table1](https://arxiv.org/html/2109.01143v1) describes 122 levels,
not the held 41, and documents ABO and two alternative gf scales. No equivalence
of the reduced artifact to that paper or to Amarsi2018 is asserted. The historical
Gerber gateway anchor is not an identical-input independent numerical benchmark.
The raw source is not redistributed; TSFitPy's GPL does not supply a data license.

[Nahar2022](https://doi.org/10.1093/mnrasl/slac015) disputes photoionization data shown
in Bergemann2021. This is a primary methodological warning, not proof that the held
41-level artifact contains those disputed data. It must be traced rather than
silently inheriting or rejecting a model by author label. No cross-section altered.

Next: extract the already-held independent MULTI departure reference at the exact
MARCS benchmark node; freeze its atmosphere/atom/abundance/opacity assumptions and
compare, then run the bounded STAGGER O experiment if applicable scientific gates
are satisfied. MARCS 5750/4.5 is a benchmark node, not the canonical Sun's stellar
parameters and not a substitute for the STAGGER cube.

## Tests / operational evidence

Held Bride interpreter, from checkout:

```
PYTHONPATH=. /Users/ryanschmitt/bride_rya1118/.venv/bin/python -m pytest -q -p no:cacheprovider tests/test_multi_atom_rya1279.py tests/test_bride_numerics_rya1280.py tests/test_multi_atom_rya1119.py tests/test_multi3d_read_rya1118.py tests/test_bride_foundation_rya1278.py --tb=short
```

138 passed in 2.94s, exit 0. LW emitted 228 known stat_equil deprecation warnings.
Project interpreter (without LW), new parser/numerical plus reader/foundation tests:
40 passed, 7 skipped, exit 0. No full-suite/CI claim. Whitespace check exit 0.

A dedicated Sirius environment is being installed under the ticket work root using
LW0.16.1/numpy2.5.2/scipy1.18.1/numba0.67.0. No pre-existing environment modified.
Host load was 10.51/11.82/10.09 with 21 users; no heavy science job started. Future
jobs must be bounded, low priority, single-threaded initially, and separately logged.
Original FITS, frozen gold, Fe/Al products and published values remain untouched.


Further checkpoint: the explicit source Q0=3 for the triplet is now preserved
instead of guessing qCore=qmax/4. Q0=0 retains a documented automatic numerical
prescription; native frequency-grid equivalence is not asserted. The O grid has
an abundance axis; reference intake requests the full atom's stored abundance
exactly and refuses nearest-node/nearest-abundance substitution. The new reader
checks binary record identity/dimensions, uses (level,depth) ordering and checks
MARCS row identities and SI conversions on independent non-square synthetic data.
Source values are never adjusted to improve a comparison.

Sirius checkpoint: isolated environment installed successfully; pinned LW0.16.1
synthetic tests at 4b3da23f passed 95/95 in 38.34s, with the same 228 deprecation
warnings. Exact MULTI reference intake at 00814f76 succeeded: model 5750/4.5,
A(O)=8.70, 56 depths, 41 levels; max logtau500 mismatch to printed MARCS structure
1.89e-8 (Rosseland differs by 0.226 dex, so those axes are not interchangeable).
No solver had run at this checkpoint. Existing STAGGER prov.json was located;
its acquisition identity is useful, but its raw-EOS/byte-count claims are older
than the corrected reader/header audit and cannot override that evidence.

The bounded MARCS model-under-test runner now prepares a manifest by default.
It validates all source hashes, uses all MARCS header abundances with the exact
reference trace O abundance, explicit source gravity/units/boundaries, native
passive atoms, fixed ne, no molecular opacity, one thread, physical R+C diagnostics,
a 900-second budget, population/J/rate convergence 1e-6 and particle conservation
1e-10. These are declared numerical experiment controls, not abundance-admission
thresholds. It saves independent-reference departures and absolute triplet
intensities in W m^-2 sr^-1 Hz^-1, plus an isolated LTE comparator. It must label
background/EOS/model differences honestly; a converged experiment is not a passed
paper replication or solar-abundance product. The user has authorized continuing
these validation experiments; no new merge/publish permission is inferred.

## First real-input attempt and resolved format difference

At 85820b0b the prepared exact-node manifest passed source-hash checks but the
first run failed during wavelength-grid construction, before RT/SE. Source
continuum 4→41 has an energy-derived edge 82.7947 nm while its explicit table
runs 94.3253–106.4179 nm. Standard LW ExplicitContinuum clipped all its points
and raised IndexError. The original failed outputs/log remain intact under
Sirius ticket root `benchmark-q1-rays3-v1`; no fabricated spectrum was emitted.

The original MULTI2.3 source was located in the public pymulti distribution at
commit `4d6aa3c39d634b38da772e469a5692435eba1354`,
`pymulti/files/multi/mul23/source/mul23_sub.f`. ATOM 540–572 and FREQC 2325–2376
explicitly consume the supplied absolute wavelength table for QMAX<0, without
filtering it to the energy-derived threshold. The radiative Milne factor at
1010–1019 uses the actual photon frequency and LTE population ratio. Thus a
faithful legacy-model conversion must preserve the explicit table separately
from energy normalization, not truncate/relabel/shift its source data.

`MultiExplicitContinuum` implements that named convention, retains the energy
lambda0 and exposes both energy/table edges for audit. All 22 source continua
and 3778 base frequency points now construct. This is NOT a physical validation
of the source atom or an explanation of the level/table discrepancy; that caveat
remains. The converter does not manufacture cross-sections or alter level/rate
values. A synthetic conflicting-edge fixture verifies retention and zero opacity
outside the source table. Dry-run now checks native wavelength construction and
cross-sections before PREPARED, and records solver_started separately from setup.

Original-source permalink:
https://github.com/MingjieJian/pymulti/blob/4d6aa3c39d634b38da772e469a5692435eba1354/pymulti/files/multi/mul23/source/mul23_sub.f

## Real-input numerical controls after continuum conversion

At 33d95777 Sirius ran all three new test modules: 105 passed, 228 upstream
deprecation warnings, 13.33 seconds (`numerical-validation-v3.log`). The real
`benchmark-q1-rays3-v2` experiment then reached 300 iterations in 88.75 seconds.
Population and J changes were ~1e-13, but physical R+C residual stayed at
0.00531601279942. It correctly FAILED; no accepted spectrum was emitted.
Increasing iteration count or loosening the threshold is not a resolution.

At d7ee36cf the independent 100-iteration repeat preserved failed-population,
collision, radiative-rate and scaled-MALI matrices in
`benchmark-q1-rays3-v3-diagnostic/unconverged_diagnostics.npz`. The largest
physical imbalance is zero-based level40/depth54; scaled-MALI residual is
7.99e-16 while the physical residual remains 0.005316. Physical matrix column
conservation is 5.20e-16. This isolates a difference between the accelerated
operator fixed point and physical balance for this real model, without yet
assigning its cause to library, atom, or converter. No upstream defect is claimed.

Commit 14b79b90 adds an explicit `--lambda-iterate` control using LW's native
unaccelerated iteration. The otherwise identical bounded experiment is running
under `benchmark-q1-rays3-v4-lambda`. At iteration22 its physical residual is
1.81e-14; population/J changes still require convergence. All source arrays and
the original failed artifacts remain unchanged. Exact command: same reference,
one thread/nice19/timeout1000s, --run --lambda-iterate, default 300 iterations,
900-second internal budget, joint numerical threshold1e-6.

## Converged experimental oxygen outputs

See [oxygen_result_review.md](oxygen_result_review.md) for the final reviewable
results and links. Native pure-lambda+Ng runs v6/v7/v8 converged (453/832/1001
iterations;136.19/639.46/276.47s). The different initial populations and tightened
grid/tolerance controls change emergent intensity by4.20e-7 and1.91e-5 relative
to v6. Independent MULTI departures still differ by~0.18 in a diagnostic window;
no solar abundance or physical admission is claimed. Source model input8.70 was
not fitted. The unextrapolated v5 control failed its3000-iteration cap at2.06e-6
population/J change after816.47s; all failed evidence remains intact.

Ng options are explicit: order3, period5, delay10, threshold0.05, lower threshold0.
The same options and unaccelerated/accelerated transfer controls recover an
independent synthetic thermal solution from three initializations. Initial
pure-lambda synthetic controls failed the old80-iteration cap while still
approaching the exact solution; cap200, unchanged1e-11 iteration and1e-9 population
accuracy thresholds, passes. No physical parameter was tuned. Current focused
suite111passed on Sirius(12.52s),144passed across all five local LW modules,
and38passed8skipped without LW.666 known upstream deprecation warnings.

Original numeric arrays/rate matrices remain on Sirius; review copies include
CSV profiles, complete histories, immutable manifests and generated comparison.
Every committed result is registered to its actual producer/invocation. All
experiments have finished; no owned O solver is left running. State registerv176
records these experimental outputs and the open scientific gates.

Final result-generator guard:2436 existing violations, the same against main's
manifest/tracked tree; zero added violations. This is not a global green claim.
All nine copied experiment files are hash-identical to the comparison receipt;
all ten new outputs are registered to their actual producers. Register freshness,
source-page reproducibility and whitespace checks pass.

## First Kitt Peak observational diagnostic (Ryan-directed)

At9ddae7bd the shared solve/replay constructor and lightweight canonical Kitt
reader support fixed-population disk-flux/continuum calculation.158focused local
tests pass;20targeted Sirius tests pass. Original source histories remain intact.
The run completed in43.15s under `kitt-peak-fixed-v1` with no SE/abundance fit.
Disk-centre replay differs3.33e-15 from the original output;5/10-ray normalized
flux differs8.27e-6. Model±0.5Å window areas61.31/53.08/41.84mÅ vs corrected
Kitt Peak2005 areas77.56/67.20/53.79mÅ. Interpretation and limitations are in
[kitt_peak_result_review.md](kitt_peak_result_review.md); no canonical O value moved.

Extended project route tests:13passed20skipped, one pre-existing failure of the
old six-band1984 coverage assumption. Executing the original4542033b loader under
the same fixture setup reproduces its5000Å return (thus the old test's failure).
Updated the now-brittle inline-conversion test to exercise shared-reader behavior;
the data/measurement route itself retains identical conversion and normalization.
No new physics bug was identified; no fake bug-ledger row added. Registerv177 and
SEQUENCE record the new observation diagnostic and open scientific controls.

## Ryan-directed Mac controls and10K resolution

Ryan moved science compute to the Mac to avoid Sirius contention, requested
temperature/xi/other parameter responses and10K sampling, and reiterated the
adaptive/full3D/all-band/own-grid objective. All36 new1D solves ran serially on
Darwin arm64/one thread/nice15, with no new Sirius science job or dependencies.
Native Teff/logg/mixture controls, xi pairs/wide probes, abundance probes and depth
controls are completed; canonical iSpec broadening and frame controls are separate.

See [mac_controls_review.md](mac_controls_review.md) for the full case/evidence
matrix. Native held-out atmosphere tests uncovered large error in linear pressure/
height interpolation. Density-retaining, positive-layer coordinates reduced tested
temperature/gravity area error to~0.2–0.4% without using Kitt Peak as the target.
Five fresh NLTE solves at5752/5762/5772/5782/5792K resolve the requested10K sequence
with that prototype; this is not universal10K physical accuracy. The experimental
5772/4.438/O8.70 profile remains18–19% below observed window absorption.

All comparisons use air Angstrom; delivered Kitt Peak2005 vacuum nm is converted
once, engine frequencies remain vacuum, and solar gravitational redshift is a
separate control. Canonical ORCH policy serves the corrected pre-normalized atlas
unchanged. No continuum/RV/abundance/width optimizer or second telluric correction.

46 registered artifacts archive five source summaries,36 manifests, compact full
profile vectors and an experimental grid index with population hashes/relative
paths at experiments.bride_mac_grid. Full raw experiment outputs remain on Mac.
Generator/plan/atom/atmosphere identities are preserved.163LW/regression tests and
14project-side policy/effects checks pass. Source environment and analysis iSpec
hash/versions are stamped. All owned science jobs have finished.

The [adaptive roadmap](../../decisions/bride_adaptive_direction_rya1281.md) keeps
typed grid/atom/backend capabilities, independent held-out error controls, own
NLTE/atmosphere/3D grid-generation levels and a later qualified Gaia-identified
stellar benchmark explicit. Current capability remains experimental1D O777,
not a solved universal atmosphere/full3D/all-band model. State registerv178.

## 2026-10-27 population-audit checkpoint

Continued at d82e30e1, refreshed main remains 0ade2170. Four serial one-thread
Mac collision controls converged under the unchanged joint gate; matched
LTE/NLTE/external-departure replay and O6300/O6158 diagnostics completed.
See population_audit_review.md for exact commands, measurements, preserved
development runs, scope, failures and next work. No source data, adopted
abundance, production gate, merge or publication changed. Full outputs remain
in registered bride_mac_grid; compact review copies and plots are committed.
Register v179 records the population finding and blend work still owed.

## 2026-10-27 rate-convention and mixed-opacity checkpoint

Continuation at e75e2b94, refreshed main0ade2170. Physical process budgets
exposed CH/CH0 source-convention ambiguity; exact native grid producer remains
unknown. New explicit adapters and safeguarded Ng controls preserve all four
failed attempts and two converged controls without changing default converter.
Trace LTE Ni is added as actual shared opacity with fixed oxygen populations;
pressure/isotope/other blends remain outside prototype. See
rate_conventions_review.md, rate_conventions_v1 artifacts and open
collision_convention bug-ledger row. Registerv180 records these findings.
Next priority is producer-build/rate-dump parity, then qualified blends and
additional solar elements before another star. No science admission or merge.

## 2026-10-28 native interpolation and reverse-charge checkpoint

Continued from ae5e65c3; refreshed main5d517cbc includes Na PR595 and is
not merged into this stacked branch. Found previously missed lowercase
RCOLL/GENCL in held source4d6aa3c; built a direct DP TAUTSP/PPVALU oracle
with existing Mac gfortran, no source/dependency changes. Four native/paired/
LTE-H controls converged; reverse-charge hypothesis improves photospheric
reference populations and full-atom LTE closure, but outer layers still differ.
See native_charge_review.md for exact controls, bounds, original outputs,
verified rate replay, unavailable native opacity inputs and next steps.
Registerv181 and existing open collision_convention bug record the finding.
No production/default converter/gold/adopted value/merge or publication changed.

## 2026-10-28 reverse-reservoir contract start

Both tickets have consolidated findings and the five-step next-work list. Ryan
sets O, then C/N, then Fe before another star. The first bounded Mac reservoir
probe is complete; see charge_reservoir_review.md and the committed CLI.
Original reports remain in bride_mac_checks/rya1281-charge-reservoir-v1/v2.
No RT/SE solve, source/rate change or admission. Next: coefficient provenance and
matched reverse-channel/statistical-weight prescription, then outer UV controls.
Linear get_agent_context remains advertised but returns Tool not found.

## 2026-10-28 sourced MN / UV checkpoint

New report: author_charge_uv_review.md. All240 raw CH0 entries follow the
published MN reverse/Boltzmann candidate within0.48%; source/producer identity
is still not equivalent. Two predeclared Mac O41 solves converge: author MN
568iterations/43.0445s and47-line native LTE UV background396/250.2854s.
The surface mismatch remains; the held MULTI TRANSP/TRANF code explicitly uses
finite top optical depth, unlike our first-point ZeroRadiation condition.
No finite-depth cap has yet been implemented. No adopted abundance/default
mapping, source atom or complete author matrix changed or redistributed.

Original source/code/run receipts and diagnostics remain immutable under
CODEX_WORK_ROOT/bride_inputs/rya1281/barklem-oxygen-cd2f54ed and
CODEX_WORK_ROOT/bride_mac_checks/{rya1281-author-charge-v1,
rya1281-uv-background-v1,rya1281-charge-contract-v3,
rya1281-author-analysis-v2,rya1281-photoionization-v2}. Earlier versions preserved.
18 registered review artifacts and3 inspected figures archived in author_charge_v1.
All owned science jobs finished; no Sirius science jobs. Missing source ABSDAT
and exact producer run files remain a parallel prerequisite, not substituted.

Validation:167 selected LW/parser/numerical/reference/observable/control tests
pass (666 upstream deprecation warnings);3 installed iSpec checks pass (4 warnings).
After the final held-resource fixture addition, the2 new tests were rerun and
pass. Full fixed-physics grid CLI refuses the experimental UV run as expected;
the first canary invocation lacked the Mac host profile and failed before that
check, then the correctly profiled canary reached the intended rejection.
Sources page212 rows matches the bibliography. Standing library audit remains
186documents/125accounted, with pre-existing unindexed papers. Generator guard
remains2436diagnostics, exact same set as the preceding branch checkpoint:
0added/0removed.18new artifact entries and their invocation flags are admitted
once staged; no global-green claim. State Register v184/freshness/whitespace pass.

Current main fetched at9d3ebbea (PRs599/600/601 now merged); SEQUENCE records those
landings from merged diffs and ticket reports. No rebase/merge of main performed;
branch is still the existing stack on parent4635c742/PR597. Do not transfer older
main failure counts to this new tree. No full-suite or CI-green claim; final
checkpoint uses[skip ci] to avoid Sirius runner work. The existing open
collision_convention bug row was amended, preserving unrelated mixed line endings.

Linear get_agent_context and record_agent_learning are advertised but both return
Tool not found; verified findings are retained in the repo and issue comments.

Exact resume: verify MARCS RHOX/top optical-depth units, retain those native fields
at the atmosphere boundary, and predeclare a finite-depth top-cap diagnostic
using actual top opacity/source values. Test zero-depth, thick and angular limits
before a matched O41 run. Do not tune thickness/temperature or inject arbitrary
UV radiation. Then compare full-depth populations and fixed spectra, not only
photospheric agreement. Solar CNO milestone precedes Fe and another star.

## 2026-10-28 finite-column cap checkpoint

See top_cap_review.md. Native RHOX is retained as kg/m2; native column/density
sets effective84362.76027m cap with actual opacity/emission/scattering. Defaults
remain ZeroRadiation; unsupported/interpolated columns cannot borrow a native
anchor. Atmosphere interpolation now withholds unvalidated column metadata.
No mesh/temperature/opacity multiplier is fitted.

Original single-refresh plan (code5ddb7571) has2failed1500-step attempts, preserved
in bride_mac_checks/rya1281-top-cap-v1. Separate predeclared nested J/cap solve
(code3aca3092) converges both controls in291/836outer steps,199.8076/487.5533s:
source top upper b1.27969 vs reference1.28139, all-depth upper errors0.02620/0.03675.
Fresh radiation/R+C replays give4.646e-7/2.008e-7 physical imbalance. O II still
fails reference population parity. Visible triplet window areas nearly unchanged.
No new solar abundance, all-band validation or other-star claim.

Original stores: bride_mac_checks/rya1281-top-cap-nested-v1,
rya1281-top-cap-analysis-v2, rya1281-top-cap-author-rates-v1,
rya1281-top-cap-raw-rates-v1, rya1281-top-cap-photoionization-v1. Earlier analysis
and failed cases retained.28registered review artifacts and3figures archived in
top_cap_v1; full rate arrays external. All owned science jobs finished; no new
Sirius science jobs or dependency installs. Source constructor/adapter receipts
preserve exact code identities. Existing collision_convention row remains open;
a distinct boundary_contract row records this newly controlled geometry issue.

Validation:171focused science tests pass (666 upstream deprecation warnings),
3installed iSpec checks pass (4 warnings). Actual fixed-physics grid rejection,
source-page consistency, generator diagnostic comparison and register freshness
are checked at handoff. Standing library audit is now188documents/125accounted;
new unindexed Beichman2025/Sanghi2025 papers are unrelated to this experiment.
No global library-green, current-main suite-green or CI-green claim.

Main fetched42987358: PR602 uncertainty follow-up, PR603 withdraws97IAG VIS
products, PR604 runbook lessons. SEQUENCE sources those merged changes and the
RYA-1310 close; none merged/rebased into this stacked Bride branch. Kitt Peak
source/policy and other frozen experiment inputs unchanged. State Registerv185
preserves scientific HOLD and exact source/producer/normalization gates.

Exact resume: establish native producer depth mode and actual LTE populations,
electron density and opacity/radiation references before interpreting the O II
b discrepancy as an absolute ion-fraction error. Held LTEPOP includes the usual
0.5 electron phase-space factor; no factor-of-two defect was established. Then
check boundary normalization/background sensitivity without fitting factors;
validate native-column interpolation before capped solar-parameter reruns.
Continue O then C/N (CNO milestone), then Fe, then another star.

Handoff checks completed: fixed-physics grid CLI rejects the cap run with the
intended diagnostic-physics message; source page212 rows matches bibliography;
28new registered artifacts and invocation flags validate once staged. Global
generator diagnostic set is exactly the preceding branch's2436 (0added/0removed).
Whitespace and registerv185 freshness pass; only one new boundary_contract ledger
row was appended, preserving historical mixed line endings. No owned scientific
process remains running. Final checkpoint retains[skip ci] for the Sirius runner.

Native-LTE checkpoint validation:44 focused tests pass. New registered receipts
pass generator/invocation accounting; global generator guard remains exactly2436
diagnostics (0added/0removed versus prior branch). Results-ledger guard retains
3documented-owed Ba/Sc failures,0undocumented. Profiled path inventory has35present
and15unrelated absent holdings; both new experiment paths resolve. The initial
unprofiled path inventory is not an admission check. Register mirrors, pointer,
freshness, SEQUENCE reconciler and whitespace checks pass. No full-suite or CI
green claim. Existing collision_convention/boundary_contract run-bug rows remain
open; this checkpoint adds no verified science defect or new bug row.

## 2026-10-29 conditional continuum checkpoint

Four final Mac LTE captures (one thread,1.40s combined) at822 exact source samples
x56depths; both OPA coordinate hypotheses and matched O8.66 versus benchmark8.70.
No SE/population update or solver-opacity substitution. Local/MARCS absorption
at2769.76A top0.267519; frame effect0.067649%, O effect0.0000931%. Si/Mg contribute
40.976%/20.108% of local absorption there. FUV scattering is H resonant Rayleigh
dominated; independent TS source moves H-wing opacity between absorption and
scattering around1400A, so no missing-opacity or producer-error claim follows.
OPA frame remains UNVERIFIED. Detailed evidence:marcs_continuum_review.md.

Final external audit rya1281-marcs-continuum-v4 and analysis-v2 are immutable;
prior versions retained. Four review artifacts plus inspected figure archived;
full arrays remain external with hashes. Backend binary/passive factories and
inspected backend/TS source receipts included.120 focused tests pass,666 known
printUpdate deprecations. Statev187/SEQUENCE updated; defaults/gold/adopted values
unchanged; no new verified run bug, collision/boundary rows remain open. Main
dbda7faa fetched unchanged. Exact resume:component-resolved independent continuum
on the Mac at matched T/Pe/composition, prioritizing Si/Mg/Fe/H- and hydrogen-wing
partition, before a source-supported UV-background control and fresh OII/rate
checks. Continue O then C/N, Fe, another star. Final HEAD must retain[skip ci].

Continuum checkpoint guards: new artifact/invocation entries pass; generator
diagnostic set exactly matches prior2436 (0added/0removed). Source page212rows
regenerates/checks. Library gained Lawler2013 TiI PDF; first-page title/authors/DOI
verified against existing lawler2013 row, local_file linked and metadata_verified
set. Initial free-text verified label was rejected by schema and corrected before
regeneration. Library now189documents/126accounted,63pre-existing unindexed remain.
Pre-existing results-ledger Ba/Sc documented-owed failures remain unchanged.
Held Mac babsma_lu is present and ARM64 under engines.ispec; per-component DETABS
printout exists in inspected source, but no new TS computation was started here.

## 2026-10-29 native Fe I budget checkpoint

Protected-source TS source21269948 baseline/reporting derivatives built on Mac.
Original source/DATA/binaries untouched. Final build-v4, run-v8, analysis-v3;
all prior attempts retained. Exact native REAL32+E14.7 T/Pe conversion checked;
max density change0.0666% retained, no Pe recovery. Four cases0.83s total, one
thread. Reporting preserves byte-identical serialized opacity at both OPA frame
hypotheses. FeI84.37x local at277nm top while neutral Fe/kg1.0384x; netopacity
perneutral81.25x includes level weighting/cross-sections, not a single-level
cross-section measurement. TS absorption matches MARCS within0.12% there.
Post-DETABS H/NH additions explicitly accounted. Full evidence in
ts_fe_continuum_review.md and5review artifacts under ts_budget_v1.128tests pass,
666known deprecations. Report-format-capacity bug characterized; distinct
background_contract row OPEN. No opacity replacement/O NLTE run/adoption yet.

SEQUENCE pre-push reconciliation initially caught newly merged RYA-1317/PR606
(f8c723a8). Read merged diff and ticket00c7a562; recorded its CNO per-group
Asplund/sigma-stamp scope only. No main merge or import of product state here.
Registerv188 reflects only this experimental Fe background diagnosis.
Generator guard remains2436diagnostics,0added/removed against prior branch;
library189documents/126accounted,63pre-existing unindexed. No global green
claim. Next: source-derived FeI packet over the O lower/upper triplet NUV
photoionization supports; validate subtraction/Kirchhoff and EOS normalization
controls, then bounded Mac O UV J/population/rate comparisons. No fitted opacity
multiplier. Source frame/producer/full3D/all-band gates remain open.

## Source-derived FeI NUV control setup (continuing)

Native TS NUV baseline/reporting packet run rya1281-ts-nuv-v1 completed at
1601samples,200.0648-360.1027nm (native100sample extension retained). FeI-only
packet staged at experiments.bride_fe_packet/fe-nuv-packet-v1; source native
detail/opacities/build are hash-bound. Opt-in provider replaces old FeI bf,
preserves scattering and LTE emissivity, refuses changed T/ne/nH/depth, and
exposes raw_specific versus local_neutral controls. Actual-backend replacement,
identity/no-double-count, emission, scattering, repeated recomputation, deepcopy
and10K-node-refusal tests pass. Runner/rate/photoionization receipts preserve
Fe packet identity; fixed-physics grid refuses this intervention. Predeclared
plan rya1281_mac_fe_background_controls.json keeps1500iterations/600s per case,
outer1e-6 and inner radiation1e-8, source-MN collisions/native cap unchanged.
No new O solve at this implementation checkpoint; scientific HOLD remains.

## Completed FeI source-packet oxygen controls

Native NUV packet1601samples covers200.0648-360.1027nm; only FeI bf replaced,
old FeI subtracted, scattering unchanged and LTE emission preserved. Three Mac
runs (neutral/raw/neutral+packet-grid) converge in358/296/362outer steps and
197.38/174.39/239.96s at unchanged1500step/600s,1e-6outer,1e-8inner limits.
Top OII b16.252baseline becomes10.0588/9.9145/9.9645 vs reference8.1681.
All-level max b discrepancy8.084 ->1.891/1.746/1.796; upper-triplet discrepancy
slightly worsens, retained honestly. Grid effect on OII0.94%, area<=0.00184mA.
Fresh own rate imbalances4.54e-7/3.56e-7/4.82e-7; external-b/local-LTE remains
non-equilibrium0.0931/0.0857/0.0883. Independent lower-triplet top photo rate
16249 ->8241/8066/8111s^-1; J/B64.62 ->about32. Source-supported physical effect,
not a fitted multiplier or complete native producer match.

Dedicated Kitt Peak audit preserves corrected/pre-normalized observation and
AIR comparison; model internal frequency remains vacuum. Fixed canonical
broadening and separate deterministic gravitational-shift variant, no optimizer.
Sampled model areas68.893/59.336/46.368mA versus77.557/67.200/53.794observed.
Fe changes only0.10-0.15mA. Remaining native5750/4.5-vs-Sun, blends/continuum,
velocity and approximate LSF limitations are not solved by the UV correction.
Initial plot autoscale included off-view convolution edges; finalv2 plot is
readable and every metric/profile hash is unchanged.

Daemon restart lost the PTY handle of the sampled run, but fresh history writes
proved it continued; no duplicate solve started. All final cases/replays finished.
Nice19 was denied and first solver continued ordinary priority; all compute
one-thread Mac, no Sirius science. Initial iSpec tests used wrong isolated venv
(no pandas);3checks pass in project venv without workaround.131science tests pass.
Implementation commits85346551/06d31eb2/609a7d27; final source runs neutral-v1,
raw-v1,sampled-v1; comparison-v3,three rate/photo-v1, Kitt-v2 immutable externally.
47registered artifacts and4inspected figures in fe_background_v1. Registerv189
records qualified findings; existing background_contract annotated but remains
OPEN alongside boundary/collision contracts. No new run bug in this source-packet
phase, no source atom/default/gold/adopted/product change. Exact next: refine
numerical/background/producer contracts, validate native-column interpolation
and regenerate source opacities before canonical-Sun/10K controls; resolve visible
blends/LSF/velocity under frozen observed continuum/tellurics. O then C/N,Fe,star.

Final Fe-control handoff checks:131science tests and3project-iSpec checks pass
on final code, with666/4known deprecations. Fixed-physics CLI actually rejects the
Fe experiment with its diagnostic-physics message. Source page212rows checks;
library189documents/126accounted,63pre-existing unindexed. Result-generator
diagnostic set matches prior2436 exactly (0added/0removed); new47artifact and
invocation entries pass. Register mirrors/pointer/freshness and SEQUENCE checks
pass. Global results-ledger Ba/Sc documented debt remains; no full-suite, library
or CI-green claim. Final checkpoint must include[skip ci] for Sirius runners.

## 2026-10-30 column/interpolation checkpoint

Pure Mac input audit, no RT/SE. Nine native5500/5750/6000,g4/4.5/5models;
fouroutercorners and five held-out edge/central targets.45candidates retain
all existing thermo modes and linear/log/integrated-log-top column options.
Native rho*dz closure0.067-0.307%; central native0.07196%. Central linear-field
linear-column closure70.83%; log-layer/log-column0.5447%, held-out column
error0.7163%, T-profile error19.8K and Pe5.00%. Integrated-log-top column
closes by construction, so independent native error is required:0.67-3.52%
across5targets. No best-case selection or interpolation/EOS accuracy claim.
Current production interpolator still withholds columns, Fe packet still
refuses changed T/ne/nH/depth. Final external audit column-interpolation-v2
retained with v1;3byte-identical registered artifacts in column_interpolation_v1.
11focused tests pass. New characterized model_domain limitation ledger row;
no new software crash/adopted-product defect. Statev190 and SEQUENCE updated.
Next structured-cell/simplex comparison and additional independent native cells,
then joint thermodynamic/column validation before canonical-Sun/10K background
regeneration/capped runs. Current O controls and all prior gates remain HOLD.

Checkpoint verification: 11 focused tests pass; register pointer/mirror and
history/PR freshness checks pass (existing tracker/gold provenance notice).
SEQUENCE reconciles against fetched main f8c723a8. Archived three artifacts
are byte-identical to external v2. Generator diagnostics remain exactly2436,
zero added/removed against previous checkpoint. Standing library audit now
finds190 documents/126 accounted for (64 unindexed; one newly present unrelated
PDF since prior63); no literature acquisition or science citation from it.
No global-clean or CI-green claim.
