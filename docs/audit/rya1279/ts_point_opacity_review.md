# Exact-input Turbospectrum opacity adapter

RYA-1281 / parent RYA-1278, 2026-11-01. Mac-only opacity/EOS diagnostics.
Base468bf9c54a86956206d6e05a99775b194b15cf0f; fetched main f8c723a8 unchanged.
No oxygen RT/SE, abundance fit, atmosphere-equilibrium solve or production output.

## Implementation and source protection

The previously verified isolated TS baseline/reporting builds initialize the
native abundance, partition and opacity-data machinery. An input-only derivative
loads a complete T/Pe point table just before NEWT=2 and the native opacity calls
in BABSMA. Temperatures are K, electron pressures serialized as dyn/cm2. No ABSKO,
JON, DETABS, cross-section, EOS, partition, abundance or rate formula is changed.
Every supplied point is required. Existing compiler recipe/flags are reused.

Both derivatives are built from hash-verified baseline/reporting source copies,
with every other copied Fortran/include file checked against the source receipt.
Original iSpec sources, DATA and binaries are untouched. Input data are copied
from the existing receipted diagnostic run and verified before execution. Mac
only, serial, 60 s maximum per opacity subprocess; no new dependency or Sirius job.

The native MARCS file is a read-only driver initialization scaffold. Its geometry
and tau are **placeholders for the changed-input opacity evaluation**, not a
claim of a new MARCS equilibrium atmosphere. The resulting .opa files must not
be used as replacement atmosphere products. The original 56-row structure is
retained; previous solar oxygen experiments dropped one bottom row.

## Native identity before changed input

At identical native5750/4.5 T/Pe, adapter baseline continuum.opa is **byte-identical**
to the original native baseline over all56 depths and1601 wavelengths. Adapter
baseline and reporting outputs are byte-identical. Reporting state, component
records and post-DETABS records are array-identical to the original reporting run.

The extracted FeI mass opacity, neutral Fe per kg and effective opacity per neutral
Fe are also **array-identical to the original native Fe packet**. This verifies
the adapter's input path and accounting. It does not independently validate TS
cross-sections/EOS or the physical accuracy of the original atmosphere.

## Five solar input tables

The five existing 5752/5762/5772/5782/5792 K targets at logg4.438 use the nearest
5750-6000/g4-4.5 cell, existing log-positive Pe/rho/log-layer interpolation and the
full56 rows. Plan, source model, binary, adapter, interpolation and data hashes
are retained. Both builds agree byte-for-byte for every changed-input case.
Each produces1601 coefficients across the same finite NUV support as the native
packet (about200.065-360.103nm vacuum); wavelengths originate in native AIR output
and are converted once with the shared reader. No wavelength frame is inferred
from the unverified MARCS OPA companion.

Native default-REAL/output quantization is retained explicitly. Maximum T
quantization across all cases is0.0005K; maximum relative Pe change8.07e-8.
All six cases, including native identity, finish in31.57s summed elapsed time
(includes parsing/checks). These are opacity jobs, not oxygen-solver timings.

| Requested Teff | Max TS EOS vs local input rho difference | Effective FeI coefficient per neutral at276.98178nm, top row (m2) |
|---|---:|---:|
|5752|0.1882%|4.37803e-25|
|5762|0.3898%|4.46682e-25|
|5772|0.5812%|4.55710e-25|
|5782|0.7518%|4.64900e-25|
|5792|0.9011%|4.74246e-25|

For comparison, the native input's maximum density difference is0.06660%.
These differences characterize two EOS/input constructions, not a measurement
of which density is physically correct. Producer density is not substituted for
local density. Raw mass opacity remains explicitly tied to producer EOS density.

Across adjacent solar10K points, the effective coefficient per neutral Fe at the
selected wavelength changes2.01-2.03% at row1,1.519-1.533% at row29. Row29 is a
previously declared diagnostic sample, not a calculated oxygen contribution
maximum. All56 selected point rows/case are retained. The coefficient includes
native LTE level weighting/partition/cross-section and net absorption conventions;
it is **not a newly measured single-level cross-section** or an Fe NLTE result.
This response supports generating exact-input coefficients rather than reusing
the old native packet for changed T/Pe.

## Refusal and regression evidence

Twenty-nine focused tests pass: seven adapter serialization/source controls plus
the22 prior input/column/interpolation controls. Initial development testing found
that reapplying the source transform was not refused because its insertion anchors
remain present after one patch. The explicit already-patched check now refuses;
initial buildv1 is retained, final buildv2 is used. No scientific output was made
with a doubly patched source. This engineering guard defect is recorded fixed
in the run-bug ledger with the regression; no opacity formula was changed.

Actual native executable probes reject wrong point count, truncated table and
extra row with exact stop codes72/73/76. Partial failed .opa files are retained as
failed diagnostic artifacts and are not consumed. Refusal report is archived.

## Evidence and next bounded step

Original resources under experiments.bride_mac_grid:
- rya1281-ts-point-build-v2 (buildv1 retained)
- rya1281-ts-points-v1 (six complete paired runs and full coefficient arrays)
- rya1281-ts-points-review-v1 (compact checked review)
- rya1281-ts-point-refusals-v1 (negative probes)

Six byte-identical compact review artifacts are archived under
data/results/rya1281/ts_points_v1. Full arrays/binaries remain hash-referenced in
the registered Mac experiment store. Builders, runner, analysis and refusal
scripts are committed; no opaque external-only harness is required.

The adapter now supplies exact-input opacity evidence. It does **not** enable a
new provider packet schema or change native-only Fe packet refusal. Next bind
new coefficients to complete immutable local T/ne/nH/depth inputs, keep
per-neutral and producer-mass normalization separate, and qualify joint
pressure/EOS/geometry/column support before any capped interpolated O control.
Retain native-node control, fixed abundance and wavelength-sampling checks.
Visible Kitt Peak AIR residuals and producer-contract gaps remain open. No input
or spectrum is adopted; scientific HOLD remains. O then C/N, Fe, another star;
full3D/all-band validation remains separate.
