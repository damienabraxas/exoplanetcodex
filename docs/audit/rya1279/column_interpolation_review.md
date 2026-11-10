# Held-out column/geometry interpolation checks

RYA-1281 / parent RYA-1278, 2026-10-30. Mac input audit only; no RT/SE,
opacity-packet replacement, production interpolation or scientific admission.

Nine held native MARCS models span5500/5750/6000K and logg4/4.5/5 at the
same composition. Four outer corners form the existing SciPy simplex input;
five distinct native edge/central models are held out. All source files are
hash-receipted. Every existing thermodynamic convention is tested, without
choosing a favorable case: linear fields, log-positive Pe/rho, and log-positive
Pe/rho plus interpolated log layer thickness. Native tauRoss rows remain common.

Column candidates are convex interpolation of native RHOX (linear or log), and
a third construction using an interpolated log top column plus the density
integral through the interpolated geometry. Candidate values are confined to
the audit. The production interpolator still **withholds column metadata**;
existing cap/Fe-packet node refusal remains unchanged.

## Independent reference check

Integrating rho*dz from each native model's printed top column gives maximum
fractional differences from its native RHOX of0.0669-0.3073% across the nine
models. At5750/4.5 it is0.07196%. This is the printed-grid trapezoidal diagnostic
baseline, not an assertion of exact hydrostatic or radiative equilibrium.

For the held-out5750/4.5 model:

| Thermodynamic/geometry interpolation | Column interpolation | Max rho*dz / candidate-column discrepancy |
|---|---|---:|
|Linear fields|Linear RHOX|70.83%|
|Linear fields|Log RHOX|100.84%|
|Log Pe/rho; linear height|Log RHOX|74.00%|
|Log Pe/rho and layer thickness|Log RHOX|0.5447%|

Independent interpolation of density, height and column can therefore create
a large inconsistency even when all individual fields are positive/monotonic.
The existing log-layer option substantially reduces this particular inconsistency.
The held-out column error itself is0.7163% with log-column interpolation at the
central model. Its temperature profile error is still19.8K and Pe error5.00%.

## Density-integrated candidate: construction versus validation

Integrating the candidate density/geometry enforces mass closure by definition.
Its zero closure error is **not independent validation**. The held-out native
column is the independent comparison:

| Held-out Teff/logg | Max temperature-profile error | Max column error, log-layer geometry + integrated column |
|---|---:|---:|
|5500/4.5|8.80K|0.6681%|
|5750/4.0|8.05K|2.2783%|
|5750/4.5|19.80K|0.6770%|
|5750/5.0|4.20K|3.5157%|
|6000/4.5|10.90K|0.7670%|

All45 candidate combinations are retained. The simplex vertices and weights are
recorded; the central target uses a diagonal of the rectangular source cell.
This audit does not compare every possible interpolation operator or grid-cell
shape, nor does it price an atmospheric-model uncertainty from five targets.

The19.8K number is a kinetic-temperature **profile interpolation error** at a
known labeled Teff node, not an inferred stellar Teff error. A10K parameter
sampling interval is possible but cannot by itself establish10K physical accuracy.
No interpolated column or background packet is admitted from these results.

## Reproduction and next work

Run the Bride venv with the Mac host profile:
`scripts/audit_bride_column_interpolation.py --output <new-immutable-directory>`.
Final external audit: experiments.bride_mac_grid/rya1281-column-interpolation-v2.
Initial two-column-convention v1 is retained. Three review artifacts are copied
under data/results/rya1281/column_interpolation_v1 and registered in GENERATORS.

Eleven focused tests pass: analytic SI mass integration, node identity,
monotonic convex columns, no extrapolation, invalid-geometry refusal, plus
existing parameter checks. Original models, existing interpolation/solve behavior,
source atom, adopted/gold products and frozen observations are unchanged.
No new verified software defect: these are explicitly experimental interpolation
limits already protected by the missing-column refusal. Boundary-contract work
remains open.

Next compare structured-cell versus simplex interpolation and additional native
cells, retaining independent column/geometry/pressure checks. Establish joint
thermodynamic/column support before enabling a cap on interpolants, then regenerate
source opacity packets for those exact fields. The existing Fe packet correctly
rejects a changed10K node. Keep the remaining visible residuals, blend/LSF/velocity
and producer contracts separate; no observed renormalization or fitted opacity
factor follows from this audit. O, then C/N, Fe, another star; full3D/all-band gates
remain open.
