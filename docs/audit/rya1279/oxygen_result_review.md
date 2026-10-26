# First Bride oxygen spectra: numerical experiment, scientific admission HOLD

Three real-atom experiments converged on Sirius. They use the full unchanged held
41-level O atom and the exact MARCS 5750 K / cgs logg4.5 / A(O)=8.70 reference
node. **8.70 is an input, not a fitted solar abundance.** These are static 1D
vertical intensities, not solar STAGGER, 1.5D, flux or full-3D results.

| Control | Iterations | Wall seconds | Physical rate residual | Max relative NLTE intensity difference from first run |
|---|---:|---:|---:|---:|
| 3 rays, native line grid, joint tolerance1e-6 | 453 | 136.19 | 1.81e-14 | baseline |
| 5 rays, doubled line sampling, tolerance1e-7 | 832 | 639.46 | 5.54e-14 | 1.91e-5 |
| 3 rays, tilted initial populations, tolerance1e-7 | 1001 | 276.47 | 1.86e-14 | 4.20e-7 |

All use native pure-lambda iteration with third-order Ng extrapolation, explicitly
recorded controls, fixed source structure/ne/abundances, one thread and nice19.
Particle conservation errors at termination are below9e-16. Positive finite
populations and intensities, population changes, J changes and independent physical
R+C residuals are checked. Simultaneous grid/tolerance/initialization changes do
not isolate individual numerical errors; these differences are diagnostics, not
an adopted abundance uncertainty or universal certification budget.

The all-level maximum absolute b differences are0.00576 for refinement and0.00209
for the tilted/tighter run. Small changes per iteration alone do not prove a bound
on total population error. Emergent-intensity stability is reported independently.

## Evidence and reproduction

- [Machine-readable comparison](../../../data/results/rya1281/oxygen_validation_v1/comparison.json)
- [First profiles](../../../data/results/rya1281/oxygen_validation_v1/q1_rays3/triplet_profiles.csv)
- [Refined profiles](../../../data/results/rya1281/oxygen_validation_v1/q2_rays5/triplet_profiles.csv)
- [Different-initial-population profiles](../../../data/results/rya1281/oxygen_validation_v1/q1_rays3_tilted/triplet_profiles.csv)

Each profile has air wavelength in Angstrom, vacuum wavelength in nm, absolute
NLTE and LTE vertical intensities in W m^-2 sr^-1 Hz^-1, and their ratio. Manifests
and complete iteration histories sit alongside them. All original arrays,
including populations and physical matrices, remain in isolated immutable Sirius
directories named in the comparison. Their hashes are recorded there; raw atom,
departure-grid binary, atmosphere, spectra/FITS and frozen gold were not changed.

The original producer is committed at df866b4e; the comparison harness at9bb75f8b.
Generators and exact invocations are registered in data/results/GENERATORS.yaml.
Use a new output directory when reproducing: overwriting evidence is rejected.
The source reference receipt is
`/mnt/codex-data/codex/rya1281/reference-5750-4.5-A8.70-v1.json` on Sirius.

## Why this is not yet a solar oxygen result

The independently produced MULTI departures disagree with all three runs. In the
diagnostic reference logtau500 window[-3,0], upper triplet departures differ by
up to about0.18. This window is not a calculated line contribution function.
The refined NLTE/LTE intensity ratio is close to unity; the mismatch must not be
resolved by fitting rates, gf, atmosphere or an abundance toward the reference.
Default LW background/EOS, missing molecular opacity and exact reduction lineage
of the held41-level atom remain open relative to the reference and the cited
122-level published model. Matched opacity/EOS/rate-convention evidence is needed
before assigning the discrepancy to a particular cause or admitting this model.

The earlier MALI experiment is also a material failed control: population/J
changes settled near1e-13 while physical balance stalled at0.005316. Its failed
outputs/matrices remain preserved, with a scaled-operator residual near8e-16.
No upstream bug has been assigned. Pure-lambda plus Ng passes physical balance;
an unextrapolated pure-lambda control remained above the population/J threshold
at its3000-iteration cap (2.06e-6 after816.47s), and correctly FAILED.

RYA-1279 model/reference admission and RYA-1280 real-model numerical findings stay
open. RYA-1281 now has actual timings and three converged experimental spectra,
but no production or STAGGER/full-3D sign-off. Ryan must review these scientific
outputs before production, per DEV_CYCLE.md sections4–5. Nothing was merged or
published and no adopted abundance, model token or element verdict was changed.

## Verification

Sirius focused suite:111passed /666 upstream deprecation warnings /12.52s.
Local held LW environment, all five Bride/converter modules:144passed /666warnings.
Project environment without LW:38passed8skipped. Source-page generation matches
the211-row bibliography. Original MULTI explicit-continuum semantics and the
initial failed assumptions are preserved in [progress.md](progress.md).

The global result-generator/invocation guard is RED:2436 violations. Comparing
the same guard against origin/main's manifest and tracked artifact tree gives
the same2436 violations, with zero new violations from this branch. All ten new
result artifacts have real producers/invocations. All nine copied experiment
files match the hashes in the Sirius comparison receipt. Register freshness and
whitespace checks pass; no CI completion is claimed.
