# RYA-1281: FeII feature comparison and native broadening source audit

Predeclared2026-11-06 before new measurements. Start2235bab9, main9a74d843
unchanged, stacked base4635c742/draftPR598. Mac only; no new RT in this phase.
Historical zero/FeII spectra and both source columns remain frozen, no fitted
gf/width/RV/LSF, source/continuum adoption or abundance inference.

Manifest config/experiments/rya1281_feii_profile.json fixes core+/-0.1A and
symmetric shoulder bands0.35-0.55A from each unchanged fixed-RV observer centre.
Retain raw shipped Kitt and independent IAG native/Gaussian-reference views.
Fixed shoulder contrast cancels a multiplicative level; it is shape bookkeeping,
not an isolated EW, new ORCH continuum, profile fit or physical noise budget.
Interpolation0.005A adds no resolution or independent samples. Preserve both
no-gravity/plus-gravity conventions; no post-hoc recentring or shoulder movement.
Companion lines/molecules can affect both core and shoulders; don't assign the
whole observed feature to FeII. Actual observation arrays remain read-only.

Inspect pinned source21269948750e2177fe32c5df38e1d9bb81363374: original
anstee/bsyn/depth/IRWIN/GES inputs and iSpec export. Execute actual FeII pressure
priority via isolated original anstee routine, then exact native width fragments.
Test stage2ionization and provided versus exported-zero/Cowley Stark separately
across3T(4000/5772/10000),3H(1e19/1e21/1e23m^-3),3ne(0/1e18/1e20m^-3).
Positive H avoids native PH denominators; H2 explicitly0. Compare all54total
widths under the existing128*float32eps arithmetic guard, not a physical-error
criterion. Then evaluate every native depth with saved actual H0/He0/ne values,
without changing opacity or spectra. Never silently use a neutral-Fe Cowley
coefficient for FeII or promote a supplied coefficient over a source lookup.

Existing bibliography barklem2005_fe2 found before searching/download. No matching filename or indexed local_file was found in the held library;
unindexed documents were not comprehensively identified. A&A web page403; author arXiv abstract
https://arxiv.org/abs/astro-ph/0502098 accessible. Abstract describes24188FeII
H-collision calculations and a Kurucz-list gf>-5 applicability selection; target
RU-tag gf=-5.33 does not itself prove or disprove original table coverage. Actual
native lookup is required. Abstract only inspected, no full-paper validation or
new document acquisition claimed. DOI10.1051/0004-6361:20042469.

## Atlas comparison: finite metrics, unresolved background

Profilev1 executed after plan commit9689ef23. Final profilev2 adds shoulder
structure/candidate inventory without changing windows, spectra, normalizations
or the original24metric rows. Every old CSV field and the PNG are byte-identical
between v1/v2; no post-hoc band movement or profile tuning. Both fixed RV variants,
zero/FeII models, Kitt, IAG native and its shipped Gaussian reference are retained.
IAG native/Gaussian and repeated cases are not independent observations.

Fixed plus-gravity core contrasts (dimensionless):

| Profile | Canonical model | Lower-g model | Kitt | IAG native | IAG Gaussian |
|---|---:|---:|---:|---:|---:|
|zero FeII| -0.000272720|-0.000271512|0.004856142|0.005090364|0.005091877|
|natural FeII|0.004280312|0.004760173|0.004856142|0.005090364|0.005091877|

For canonical plus-gravity FeII, descriptive shoulder-scaled core RMS is
0.000831591against Kitt and0.000999972against native IAG; zero-line controls are
0.005324986/0.005546235. Lower-g residuals are0.000715650/0.000774101. These
finite residual reductions are not chi-squared, uncertainties, a preferred gravity
or proof of an isolated transition strength. No star parameter or gf chosen.

Raw canonical model/data core RMS remains0.00853718(Kitt)/0.01297796(IAG), so
shape agreement does not remove the broader level discrepancy. Do not uniquely
assign that discrepancy to normalization rather than missing opacity/source physics.
The observed shoulders have fractional ranges0.00677908(Kitt),0.00656757(IAG),
0.00655721(IAG Gaussian), comparable to or larger than the core contrast. Plot
shows structured absorption in both shoulder bands which the sparse model does
not reproduce. This invalidates physical isolated-profile qualification from the
contrast metric alone; the metric is still recorded exactly as predeclared.

The frozen source centre inventory has blue/red shoulder3/1candidates at the
no-gravity centre,5/1at plus-gravity. Those centres are not blend identifications
or predicted strengths. No candidate/shoulder is dropped to improve agreement.
Native IAG quality flags are exact0through the full plotted neighborhood; that
is not proof of zero telluric residual or physical reference accuracy. No
additional telluric correction. Common0.005A interpolation adds no resolution;
observations are never recentred or fitted. Relative-core-area columns are the
same contrast statistic in mA units, not an independent EW or abundance result.

## Native pressure and Stark source result

Accepted sourcev2 executes unchanged original anstee.f plus exact depth.f width
fragments in an isolated generated driver. Actual stage2lookup yields U followed
by the native negative-fdamp W fallback: loggamma6(10000K)=-7.910. This target does
not receive a quantum-table S width. Sigma/alpha0are unused native W metadata,
not invented zero physical broadening. The native branch priority is preserved.
The paper abstract is context only; it does not validate a pressure uncertainty
or establish a full physical-line model.

All54native total-width checks agree within3.97255e-7 relative, below the
pre-existing128*float32eps arithmetic guard. This is default-REAL/source-port
verification, not physical accuracy. Supplied Stark uses log coefficient-6.530;
IRWIN FeII ionization16.188eV is used for the stage2Cowley control. The source
ionic branch is5.42e-7*Eudiff*ne/(ijon+1)^2 with
Eudiff=13.598*(ijon/(eion-chiu))^2 and ijon=2. No neutral-Fe coefficient substituted.
The inspected iSpec writer does not export the supplied Stark field on this
atomic path; bsyn initializes gamst=0, invoking Cowley. Supplied and export-Cowley
are separate diagnostics, neither selected for agreement. Original files unchanged.

Actual full-column widths use saved H0/He0/ne and all56depths; H2 omitted, no
EOS or population feedback. Near logtau500=0 (native index41):

| Frozen case | Gamma natural | Gamma H/He pressure | Stark supplied | Stark export-Cowley | Total supplied/natural | Total Cowley/natural |
|---|---:|---:|---:|---:|---:|---:|
|canonical|2.95e8|1.472174e9|1.906395e7|1.753413e6|6.05504|5.99636|
|lower-g|2.95e8|1.287048e9|1.768386e7|1.626518e6|5.42282|5.36839|

Units s^-1; actual logtau500-0.009913/-0.008953. Across the full canonical column,
total/natural ranges1.01456-13.31092(supplied) and1.01453-7.41045(Cowley);
lower-g1.01280-11.91627/1.01277-6.48981. These factors multiply the local damping
width, not the emergent flux, EW or oxygen response. Pressure dominates near the
reference row. Correct implementation is depth-dependent Voigt coefficients in
the same transfer, not smoothing a final spectrum with an averaged width.
No new transfer or corrected oxygen result has been computed in this phase.

## Rejected native v1 and concrete bug disposition

Initial v1 completes the native54widths but fails strict recipe parsing because
Fortran list-directed adjacent character values emit FEII_RECIPEW without a
separator. No v1 source report is accepted. Explicit1X record formatting and a
malformed/duplicate/unsupported-record refusal test fix the interface at
842ede414e49b25051710a0ec7e894deda8f26b2, before accepted v2. All54native width
records are exactly identical in v1/v2; only the recipe record is delimited. No
precision/physical threshold changed or recipe inferred to bypass refusal.

Original v1 driver/input/stdout remain immutable externally, with the input archived as rejected and byte receipts for raw outputs,
with source code9689ef23 for reproduction. Accepted v2 outputs separate. Concrete
bug-ledger row rya1281-feii-native-recipe-delimiter is fixed, formatting-only;
previous clipping/fractional-flag and other ledger dispositions remain unchanged.

## Reproduction, QA and next task

```sh
source config/hosts/mac.env
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/prepare_bride_feii_broadening.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-width-source-v2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/compare_bride_feii_profile.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-feii-profile-v2
```

Both final commands exit0, Mac only; original compiler/source files unmodified.
Existing gfortran compiles the isolated driver (60s cap) and runs it (30s cap);
no installs/downloads, benchmark, RT/SE, Sirius or CI jobs. Lightweaver default
config warning retained. Compiler/version/command, input/native stdout/binary/
driver and pinned source fingerprints stored. Binary/data copies remain external;
no binary rename or size-guard bypass. Reproduction requires fresh directories.
Registered original v1/v2 resources preserve both rejected and accepted outputs.

7results/1visually inspected figure copied byte-for-byte under
 data/results/rya1281/feii_profile_source_v1 and
 data/output/rya1281/feii_profile_source_v1. Largest result24787bytes. No raw
external spectral flux arrays copied. Accepted native widths include every case,
mode and all112depth rows; profile summary all24comparisons, no favorable selection.
Actual writers/invocations and historical rejected code version registered.

[Atlas comparison](../../../data/output/rya1281/feii_profile_source_v1/profile_comparison.png).

81tests pass in8.36s,292existing backend warnings. Exact expanded command: the
previous77-test suite in feii_transfer_review.md plus
 tests/test_bride_feii_profile_rya1281.py, OMP/OPENBLAS threads1. Four new tests
cover fixed contrast/scale/symmetric linear baseline, truncated-window refusal,
ionic supplied/Cowley behavior and strict native record parsing. Initial focused
3passes; after parser regression4passes. All54native arithmetic cases additionally
execute on the actual compiler. No new scientific result adopted.

Self-review: both raw and derived views retained; fixed windows/RVs unchanged,
no missing shoulders replaced or data cropped for agreement. Actual source
lookup outranks assumptions from a paper/species label. Stage2Cowley distinct from
neutral Fe and provided coefficient. Entire-column density/T fields retained;
source widths are not presented as flux or physical error bars. Original source,
failed v1, old transfer/continuum/O arrays and historical3percent source rule stay
unchanged; current main15percent model-guided rule is not exercised here.

Next bounded task: replay the natural-only transfer exactly, then insert the
source W+export-Cowley and W+provided-Stark depth widths in the same transfer with
unchanged O/FeI/background, both references/RVs and all refusals. Validate broader
neighbor opacity/continuum rather than moving the visibly structured shoulders.
Only after source/profile qualification assess coupled-O feedback with fresh
independent physical rates. O atom/EOS/reference, missing wings/molecules, broader
parameters/interactions, additional O diagnostics, C/N then Fe/reference star and
adaptive3D/all bands remain OPEN. Scientific HOLD, RYA-1281 In Progress /
RYA-1278 In Review; no native/default/production, abundance/adoption/merge or
publication change. No owned jobs remain.

Native raw driver/stdout copies fail repository whitespace checks because source/Fortran output includes trailing padding. Originals are not trimmed, overwritten or given a guard exemption. Archive source input, full accepted report/widths, and native_receipts/report.json instead; receipt script verifies all54width records and raw inputs equal and fingerprints rejected/accepted raw driver/stdout/binary. The external registered v1/v2 stores retain those exact files. No whitespace rule or size guard changed. Receipt command: python scripts/receipt_bride_feii_broadening.py --first <registered-v1> --final <registered-v2> --output <fresh-registered-receipts>.

Final QA: MIRROR targets/v217 pointer/whitespace/main landing reconciliation pass. Generator diagnostics2436unchanged, zero added/removed. Fresh library195documents/127accounted/68unindexed/12off-Mac registered paths. Historical RYA-587 ledger/gold/tracker debts remain visible; no green CI claim. Parser formatter bug fixed/tested/recorded; no source/physics change, new RT or owned job outstanding.
