# Official column kernel and the actual boundary input

RYA-1281 / parent RYA-1278, 2026-11-01. Mac-only input diagnostics.
Basea29c894faa87cb7f0e5947db3e5a497e9d00f8f5; refreshed main f8c723a8 unchanged.

## What was executed

Extracted the unchanged calcrhox/rinteg/parcoe fragment from the receipted official
Masseron source. A separate guarded driver supplies log10 Rosseland tau and
log10 cgs opacity to the original default-REAL routines; results convert g/cm2
to kg/m2. Original source/manual/test assets are unchanged and not vendored.
Only this isolated mass kernel is compiled/executed, **not** the full empirical
interpolator or any upstream shell script. Existing gfortran, checked build,
bounds checking and binary/compiler/source/fragment receipts are retained.

Six analytic cases (4/9/56 rows, constant opacity and linear inverse opacity)
have maximum fractional errors4.71e-8 to2.45e-7. The predeclared conservative
float32 operation budget is an engineering guard, not an atmospheric accuracy
or admission tolerance. Independent Python piecewise-power-law integration is
also tested against closed forms for opacity proportional to tau^0/tau^1/tau^2.
Unit/order/finite/native-REAL overflow/collapsed-row checks refuse invalid inputs.

## Native replay does not follow from analytic success

Nine native models, five distinct held-out comparisons and five solar10K inputs
are checked with four numerical methods and two top anchors:152 cases.
The numerical methods are original official kernel,2x/4x quadrature refinement
of the **same** declared log-opacity interpolant, and its independent exact
piecewise-power-law integral. Extra quadrature points are not new physical
opacities. Original source grid samples are preserved.

The original top anchor tau0/kappa0 undershoots native top column by40.3-45.7%.
It assumes an opacity continuation above the first sample; it is not a measured
native atmospheric top mass. An alternate anchor uses the known native top for
native replay, or **source-corner log-interpolated top** for held-out/solar cases.
Held-out target top masses are not leaked into the interpolation.

At native5750/4.5, with the native top restored:

| Method | Max fractional native-column error |
|---|---:|
|Official kernel, original56 rows|4.1344%|
|Same interpolant,111 quadrature rows|3.5990%|
|Same interpolant,221 quadrature rows|3.5368%|
|Independent exact piecewise-log-opacity integral|3.5285%|

Across all nine native cases, the original anchored kernel differs4.119-4.164%;
the exact piecewise-log-opacity route differs3.510-3.548%. Refinement removes
part of the integration discrepancy, but does not reproduce native mass. This
is not a float32-unit error or grounds to replace the native density/column.
The exact integral validates a declared interpolant's quadrature, not that
interpolant or opacity table as a complete physical model.

At central held-out5750/4.5, anchored column error changes3.9692%->3.3566% between
original kernel and exact log-opacity integral. Across five heldouts, exact
log-opacity columns miss native RHOX by2.970-3.436%. At actual solar5772/4.438,
the corresponding mass-closure discrepancies are4.3250% and3.6705%; neither is
an independent solar-column accuracy estimate. Full negative results are retained.

## Source-history clarification changes the prior pressure interpretation

The official [April2018 MARCS note](https://marcs.oreme.org/news/) states that
RHOX was recomputed by integrating density above unit area, replacing the old
total-pressure/gravity proxy, including for plane-parallel models. It also says
that this exported mass scale is not used by MARCS or BSYN/Turbospectrum.
The inspected Masseron utility source predates that update (latest bundled
changes2012) and computes mass from tau/kappa instead.

This clarifies the export contract. The earlier3.8-4.2% pressure-difference
diagnostic is not an independent native-column pass/fail criterion and must not
be forced away. Its empirical values remain valid diagnostics; the exact
producer/discretization origin of every residual is not assigned. The existing
ledger anomaly is updated to **characterized**, with the format-history evidence,
rather than claiming a defect was fixed in MARCS. No native file is repaired.

The source-supported route for an explicit diagnostic interpolant is therefore
density/geometry integration with a separately declared native-derived top anchor.
Its zero internal closure is a construction, not physical certification; prior
independent held-out column errors and temperature/pressure limitations still
apply. Full official interpolation is not validated by this kernel experiment.

## Test the quantity that actually drives the cap

The current finite-column cap uses only **m_top/rho_top** as its effective optical
length. Deeper mass entries are validity metadata, not a formal-transfer scale.
An actual cap-constructor regression verifies different deeper profiles with the
same top mass/density yield identical cap length. Geometry/T/ne/nH continue to
drive transfer and populations independently; nothing is inferred about the true
vertical extension above the atmosphere from this effective length alone.

Using the source-corner log-native top and retained log-layer density, held-out
effective length errors against independent native models are:

| Held-out Teff/logg | Fractional cap-length error |
|---|---:|
|5500/4.5|+0.12148%|
|5750/4.0|-0.22266%|
|5750/4.5|-0.24800%|
|5750/5.0|-0.15562%|
|6000/4.5|+0.11875%|

These limited holdouts support a **bounded experimental control**, not universal
accuracy or scientific admission. Actual solar effective length is97762.94m;
the five10K inputs span97393.39-98133.89m with no native solar truth at those
targets. All quadrature methods using the same top anchor yield the same cap
length, despite their different deeper mass curves. A larger deep-column residual
must not be presented as a measured boundary-length error.

## Verification, evidence and concrete next experiment

62 focused tests pass;3 known LW warnings remain. All152 cases and six analytic
controls finish with positive ordered columns. No new code defect, source edit,
opacity generation, O RT/SE, observation/abundance/adopted change or provider
admission occurs. Existing source certificates still explicitly lack solar
columns and reject added columns. No production model or source policy changed.

Final experiment: experiments.bride_mac_grid/rya1281-rosseland-columns-v4;
v1-v3 retained as successive immutable development audits. Three byte-identical
compact artifacts under data/results/rya1281/rosseland_columns_v1; isolated source,
driver, compiler output and binary remain externally hash-referenced and fully
reproducible from the committed harness and receipted official source. Source
intake's execution=false remains a historical intake fact; current source-use
state in the path/holdings/register now records the isolated kernel execution.

Next implement a clearly tagged **experimental** density-integrated column with
log-native top, issue fresh exact-input bindings, and connect a separate opt-in
bound-coefficient provider. Preserve the legacy native-only loader. Prove native
identity and changed-input refusal at the actual backend, then run one matched
solar5772/4.438/A(O)=8.70 control with the existing iteration/time budgets,
source collision convention, per-neutral Fe and explicit wavelength sampling.
Retain rate balance, populations, convergence history and spectra; no fitted
boundary/opacity/abundance factor. This is an authorized diagnostic step, not a
new atmospheric-equilibrium grid, full3D result or abundance admission.

Scientific HOLD remains; O then C/N, Fe, another star. No merge/publication.
