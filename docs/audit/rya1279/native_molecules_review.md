# Matched native LTE molecular-line controls

2026-11-07, RYA1281/1278, HOLD. Starting Bride evidence200a2e1a,
worktree `/Users/ryanschmitt/codex/rya1279`, stacked base4635c742.
Actual current ORCH checkout de96ae51 remains isolated and supplies canonical
gf resolution. Existing source iSpec21269948750e2177fe32c5df38e1d9bb81363374.

## Predeclared experiment

`config/experiments/rya1281_native_molecules.json`: independently build unmodified
native babsma/bsyn from hash-checked original source in a new Mac directory.
Existing compiler only; module dependency order before users, original source
equations/constants unchanged. Compile caps90s per target. Record all source,
compiler, command and binary receipts. Installed binaries never replaced.

Use the already-held native MARCS5750/logg4.5/xi1reference file, fingerprint-checked
against the qualified point-producer record. Source has56depths; record native
opacity depth count rather than assume LW/native parity. Header rounding gives
parsed logg4.500003068, retained, not overwritten. Same92element MARCS abundance
vector, O fixed8.70in every case. Fixed atmosphere, no new equilibrium structure.
This reference differs from canonical solar5772/4.438 and from Bride NLTE.

Resolve the full atomic source through actual current-ORCH canonical gf owner
before selecting a +/-10A atomic margin. Record explicit partial-coverage report
and unsupported export rows; no silent fallback. Same exported atomic file for
all legs. Explicit native isotope file using inspected wrapper's3decimal fraction
serialization, unchanged between legs. Shared continuum-opacity file computed
once on0.1Agrid; original native chemistry/scattering remains active. Only molecular
bound-bound line-file lists change, not molecule formation on/off.

Native AIR output7763-7783A,0.005Aline step: off/all and each CH/CN/C2/NH/OH
family separately, using complete9routed files with native wing handling. CO.dat
not silently added or converted into this route. Refine off/all only to0.002A.
Nine serial programs,180s cap each; no measured runtime forecast yet. Retain
native EPS1e-4weak-line cutoff and five-decimal normalized-flux serialization.
Report truncation/rounding limits; do not interpret zero printed response as
zero physical opacity. Family-only effects can overlap/nonlinearly saturate;
do not add them as independent uncertainty terms or multiply spectra.

All programs get private working directories/DATAread-only links, preserved
inputs/stdout/stderr/raw/CSV. DATA source hashes before/after, same continuum
receipt between all legs. No source overwrite or destructive temporary cleanup.

## Scope and verification

Compile/source/count/grid/finite-positive guards precede interpretation. Exact
programs/flags/input bytes retained. Test existing gf/canonical/opacity readers,
then source-matched native outputs, line-step control and identity checks. This
is legacy-source LTE sensitivity, not new molecular-primary provenance, full
chemical/EOS parity, a Bride mixed-opacity/NLTE result, observed abundance or
full3D/all-band capability. No source tuning, adopted continuum/gf/abundance,
default/product/production, merge/publication, Sirius or CI. Mac-only serial
control at one atmosphere, not a production grid. Source/count inventory and
remaining physical qualification holds remain intact.

Precommit plan/build/runner before native execution. Commands use project Mac
interpreter with OMP_NUM_THREADS=1/OPENBLAS_NUM_THREADS=1 after host profile:
`scripts/build_bride_native_lte.py --output <fresh-build>` then
`scripts/run_bride_native_molecules.py --build <build> --orch-worktree <pinned-current-orch>
--output <fresh-run>`.

## Results and interpretation

All nine v2 native cases completed on the Mac. Original source model composition
is C8.39/N7.78/Fe7.45, O fixed8.70; 56depths, xi1km/s. Canonical gf resolution
applies to the full source before margin selection: 922selected atomic rows,
276unsupported export rows,646exported. This is explicit partial coverage.
No primary molecular provenance gate is cleared by running these legacy files.

Correction to the predeclared opacity request above: 0.1A was requested, but
unmodified babsma.f enforces DEL=max(del,1.d0) and adds100points. The actual
shared file has121points,7763–7883A,1Aspacing and56depths. This is native
behavior, not evidence of0.1A continuum convergence.

| O aperture centre AIR A | all minus off intrinsic area mA | CN only minus off mA | refined all minus off mA |
|---|---:|---:|---:|
|7771.944|0.084692|0.084642|0.084460|
|7774.166|0.317250|0.317250|0.316820|
|7775.388|0.465505|0.463555|0.465440|

The fixed apertures span centre+/-0.5A. All/CN maximum normalized-flux change
across20A is0.01081. C2 changes the third area0.001850mA; CH/OH print zero
aperture changes, each with whole-region maximum flux change1e-5; NH reproduces
off exactly in this delivered region. Native EPS1e-4 and five-decimal flux remain.
Paired1A-area rounding alone permits0.01mA worst-case error. Tiny family,
refinement and nonadditivity residuals are below this floor; zeros do not prove
zero physical opacity. No physical accuracy bound follows from this rounding bound.

Actual current RYA1315 ORCH continuum source de96ae51 is imported and fingerprinted:
model-guided bound0.15, envelope0.03, Q90. Both frozen RV views and both references
are retained:108fits,27intrinsic records. With gravitational RV, all-minus-off
synthesis-guided observed areas change +6.630306/-0.710741/+2.320351mA;
without gravity +6.964982/-0.954580/+1.946295mA. IAG-guided observed areas remain
identical. These observed-area changes are estimator sensitivities, not intrinsic
molecular absorption or corrections to Bride NLTE abundances.

Exact replay of108native records, continuum arrays and areas passes. At7771,
actual coefficient support changes86to85pixels, remains entirely blue-sided,
and extrapolation to centre changes2.772to2.582A. Centre weight L1 changes
6.923to5.266. At7774 support72to61changes from72/0blue/red to45/16;
centre extrapolation4.104to0A and L1 14.850to1. At7775 support90to77,
43/47to37/40blue/red, L1 remains1. Both selector and reference denominator can
change; a selector-only attribution is not established here. Historical fifth
clipping survivor counts remain distinct from actual coefficient supports.

The four-panel plot shows shipped Kitt Peak, fixed IAG helper, and broadened
native off/all spectra. No observed continuum refit in that figure. CN adds
structure near FeII; O cores remain weaker than atlases. Fixed macro3.8km/s,
vsini1.8km/s,R500000,limb0.6. AIR comparison; Kitt vacuum is converted once by
the established loader. No additional telluric correction, tuned gf, abundance,
width, RV or selected favorable mask. This is VIS7763–7783A only.

## Engineering, reproducibility and verification

First build v1 failed because recipe parsing truncated clock-gf.f; fixed in
1cea1429 with two tests. First spectrum v1 failed because native input requires
`NLTE :` with the space; blank labels also needed explicit LTE-only export
sentinels. Fixed3746be8e with two tests, including detecting a bad-input STOP
that returns zero. Both original attempts are preserved without accepted reports.
Original source equations, binaries installed in iSpec, inputs and observations
remain unchanged.136source/build/derived-input hash checks pass.

Accepted native run v2 has4001points per nominal case and10001per refinement.
Native calls took0.605–0.689s nominal,1.320/1.368s refined, opacity0.178s.
Reported7.527s excludes Python atomic loading and canonical preflight, so is
not end-to-end wall time. Repeat v1 ran but macOS time was blocked by sandbox
sysctl permissions; preserved, no accepted memory report. Authorized repeat v2
outside sandbox succeeds: off/all native spectrum files byte-identical to original,
0.696/0.637s, maximum resident sizes26247168/26198016bytes. Those memory readings
apply only to these two child programs, not the nine original calls or Python setup.

Focused final suite35passed in2.04s; earlier current-ORCH gf/molecular suite16passed.
Library audit196documents/127accounted/69unindexed/12off-Mac remains incomplete.
Results archive38files, maximum1917269bytes; no binaries or native source copies
committed. Source/build/data and failure logs remain in registered Mac resources.
Every archived file is byte-checked and linked to its actual generator.

## Reproduction and next steps

Use `source config/hosts/mac.env`, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and
`/Users/ryanschmitt/codex/venv/bin/python`. Commands and source fingerprints are
in archived reports and data/results/GENERATORS.yaml; all output arguments require
fresh directories. Build script precedes run; analysis precedes trace and plot;
receipt checks retained failed and accepted directories; repeat uses accepted build/run.

Next qualify native CN chemistry, partition, isotope and line-coefficient interface
against independently sourced semantics before putting molecular opacity into the
Bride transfer. Then matched same-atmosphere/EOS molecular controls and fresh O
NLTE solves, retaining both reference/RV views and support diagnostics. No adding
native LTE flux differences to LW spectra. Continuum leverage, primary source
provenance, atmosphere mismatch, independent transfer/MULTI and3D/all-band gates
remain open. C/N and later Fe abundance work follow those checks; scientific HOLD.

Final repository guards: mirror, version pointer and sequence reconciliation pass.
Generator check retains exactly2436pre-existing errors, zero new or removed.
Global bug-ledger suite remains1failed/1passed on the previously documented
RYA587blank science_impact; both new rows contain required fields. No global-green claim.
