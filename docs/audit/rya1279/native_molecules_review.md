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
