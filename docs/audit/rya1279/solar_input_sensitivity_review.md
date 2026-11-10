# Solar 10 K input sensitivity

RYA-1281 / parent RYA-1278, 2026-10-31. Mac input-only diagnostic, no RT/SE.
Base c8f75032bb3e330d7029fb40b864478b771a6e29; fetched main remains f8c723a8.

## Design and provenance

The existing fine-temperature experiment already uses the nearest local native
cell: Teff 5750-6000 K, logg 4-4.5. This audit does not claim the earlier solves
used a coarse cell. Its five requested targets come from the existing plan:
5752/5762/5772/5782/5792 K at canonical solar logg 4.438. Center and gravity are
verified against get_star_params('solar'); configuration and plan hashes retained.

Compare three cells: nearest 250 K by 0.5 dex; wider-temperature 500 K by 0.5 dex
(5500-6000, g4-4.5); wider-both 500 K by 1 dex (5500-6000, g4-5). Retain both
simplex/bilinear methods and all three field conventions: 90 field candidates.
Eight native source models are hash-receipted. No native model exists at these
exact solar targets in this audit, so there is **no independent accuracy estimate**.

All 56 native depth rows are retained (prior fine-temperature solves dropped one
bottom row). Profiles contain T, electron pressure, rho and outward height;
354 comparisons retain operator spread, cell spread and adjacent 10 K response.
Each has both all-depth and fixed native logtauRoss[-3,0] scopes. The latter is
a declared 31-row diagnostic subset, not an oxygen contribution function.
Fractional differences use the explicitly named baseline denominator;
electron-density differences use Pe/T with kB canceling. Three column constructions
are retained separately, with no column metadata admitted to the interpolator.

## Results at the solar point

The following uses the previously proposed log Pe/rho and log-layer-thickness
convention; linear-field negative controls are also archived.

| Comparison at 5772/4.438 | Max T-profile difference, all depths | Max T difference, logtau[-3,0] | Max Pe difference, all depths |
|---|---:|---:|---:|
|Nearest cell, bilinear vs simplex|0.1419 K|0.01746 K|0.04352%|
|Wider T cell, bilinear vs simplex|1.7673 K|0.1754 K|0.54351%|
|Wider T+g cell, bilinear vs simplex|9.5270 K|1.1185 K|3.13837%|
|Wider T vs nearest, simplex|5.6572 K|1.6376 K|3.11012%|
|Wider T vs nearest, bilinear|5.4642 K|1.5741 K|3.62872%|

Nearest-cell method spread in rho is 0.02336%, ne 0.04374%. Across each of the
four adjacent 10 K steps in that cell, the maximum kinetic-temperature profile
response is 12.856 K (simplex) or 12.92048 K (bilinear). The subset response is
10.424 or 10.420032 K. A change in labeled Teff does not require an identical
kinetic-temperature change at every atmospheric depth; these are responses,
not temperature errors. Maximum Pe response per step is about 1.82%, ne 1.64%.

The small operator spread at this particular nearest-cell point is reassuring
about local method sensitivity. It does not contradict the larger independent
held-out errors in the previous coarser-cell audit, and does not establish the
physical correctness of either structure. Wider-cell pressure sensitivity remains
several percent even where the temperature difference is only a few kelvin.

At the solar point, nearest-cell log-column rho*dz closure discrepancy is 0.1536%
(simplex) or 0.1646% (bilinear); subset values 0.04227%/0.04642%. These are internal
consistency checks, **not held-out column errors**. Density-integrated candidates
close by construction. No input is certified from these numbers.

## Corrected diagnostic subset bug

Initial immutable run v1 built subset masks from interpolated logtau values.
Floating roundoff excluded an exact boundary row in six of 90 field candidates:
30 rather than 31 rows. This was discovered by counting selected rows, before
adopting results. The corrected implementation selects fixed native source row
indices once and records them, after verifying all source grids agree and
candidate coordinate differences are roundoff only. Profiles display native
coordinates. v2 is the final run; v1 is retained and not used as the review copy.

Six of 354 comparison rows change, all in the subset; maximum rho-spread correction
is 0.0010254 in fractional units (0.10254 percentage point). Temperature, Pe, ne and
all-depth comparison metrics are unchanged. No source/candidate thermodynamic
field or solve is changed. Regression includes exact boundaries and a deliberate
nextafter counterexample. Bug ledger records the fixed diagnostic row-selection
defect separately from the existing open model-domain limitation.

## Verification and next action

Twenty-two focused tests pass: six input sensitivity/mask controls plus the 16
existing rectangular/column/parameter controls. Original source models, native-only
Fe packet checks, production interpolation and column withholding remain unchanged.
No observation, abundance, atomic rate, opacity scaling or production output changed.
Scientific HOLD remains; child In Progress and parent In Review are preserved.

Reproduce with the Bride environment and Mac host profile:
`python scripts/audit_bride_solar_input_sensitivity.py --output <new-immutable-directory>`.
Optional --plan must satisfy the declared existing fine-temperature contract.
Final original: experiments.bride_mac_grid/rya1281-solar-input-sensitivity-v2.
Four byte-identical review copies: data/results/rya1281/solar_input_sensitivity_v1.

The next useful check is joint pressure/EOS/geometry support for the nearest-cell
inputs, followed by exact-structure background packet generation. No packet from
5750/4.5 may be reused at changed T/ne/nH. The earlier native-node oxygen/Fe
background comparison and visible Kitt Peak AIR residuals remain independent
evidence; do not tune solar parameters to erase them. Local 10 K parameter
sensitivity is feasible, but admission of interpolated caps still needs joint
input support. Oxygen then C/N, Fe, another star; full 3D/all-band gates stay open.

A source-level route for the next opacity check is identified: held iSpec pin
21269948750e2177fe32c5df38e1d9bb81363374, Turbospectrum source absko.f lines
2-23 and139-154. ABSKO accepts T/Pe point arrays (K and dyn/cm2), invokes JON
for gas pressure/density, then DETABS for absorption factors. The held source
matches that pin byte-for-byte (SHA256
f6bf92575be6b938725e25ef01f9143ab1779599fbbfb5a4d88a8e7300c501c7). This suggests an
isolated point-evaluation adapter could generate opacity at exact local T/Pe
without pretending an interpolant is a new MARCS equilibrium model. It is only
source-interface feasibility: initialization, abundance/partition/data identity,
full component accounting, native-node replay identity and per-neutral versus
per-mass normalization must be verified before any changed-input packet use.
No adapter, packet or capped solve has been executed or admitted by this audit.
