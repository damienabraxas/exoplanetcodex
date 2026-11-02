# Joint pressure, geometry and column audit

RYA-1281 / parent RYA-1278, 2026-11-01. Mac-only input diagnostics.
Base506d53315a9d36dddc6f57bd886daadfd0a90b3b; refreshed main f8c723a8 unchanged.

## Reference and method

Read held Gustafsson2008, printed pp953-954 (PDF3-4), Eqs1-5 and the discussion
of turbulent support. For plane-parallel geometry, the total-pressure equation
suggests d(Pg+Prad+Pturb)/dm=g. Differences remove an unknown pressure intercept;
they do not recover top column mass from absolute pressure. The paper explicitly
warns that the radiation-force boundary is nontrivial. Its chosen grid support
velocity is zero; a spectral microturbulence label is not automatic evidence for
a turbulent-support pressure term.

The official [MARCS reader](https://marcs.astro.uu.se/documents/auxiliary/readmarcs.f)
documents the gas/radiation/turbulent pressures in dyn/cm2 and RHOX as mass above
the point in g/cm2. The supplemental parser checks row/T/Pe/depth identity against
the existing reader and converts to SI. It does not alter that reader or a source.

Nine native models (5500/5750/6000 K by g4/4.5/5), five held-out comparisons and
six exact producer-point cases give40 diagnostic cases. Held-out/solar fields use
the existing log-positive Pe/rho/log-layer simplex convention. Supplemental Pg
is interpolated logarithmically, Prad/Pturb linearly with identical source weights;
this pressure convention is declared, not admitted. Three columns are retained:
log-native RHOX, density-integrated with log-native top anchor, and pressure-
difference construction with the same anchor. The latter two enforce different
identities by construction and are **not** independent validation.

## Native reference limitation

Across the nine native files, maximum adjacent pressure/column discrepancies
are3.814-4.222%; top-anchored differences3.792-4.103%. Their native rho*dz closure
is much smaller,0.0669-0.3073%. At native5750/4.5, adjacent pressure discrepancy
is3.868%; top-anchored3.810%; native mass closure0.07196%.

This is a real **diagnostic/native-reference contract gap**, not proof that MARCS
violates equilibrium. Units and table identities are checked; the relationship
between these printed pressure fields, producer radiation-force treatment and
native discretization remains unresolved. Do not assign the whole discrepancy
to interpolation or silently calibrate it away. The run-bug ledger records an
open model-domain anomaly, not a verified defect in the MARCS producer.

A nonfitted beta1 rho*xi^2 support control is retained, using the header spectral
microturbulence only as a hypothesis. Across native models its maximum adjacent
discrepancy spans2.282-5.229%: it is not a universal fix and sometimes worsens the
baseline. No such term is adopted, no xi/support identity is inferred, and the
printed zero Pturb is preserved.

## Held-out and actual solar comparisons

At central held-out5750/4.5:

| Column candidate | Max mass closure discrepancy | Max native column error | Max adjacent pressure discrepancy |
|---|---:|---:|---:|
|Log-native interpolation|0.5447%|0.7163%|3.9086%|
|Density integral|0 by construction|0.6770%|4.3921%|
|Pressure difference|3.6186%|3.3439%|roundoff by construction|

Pressure-derived columns miss independent native RHOX by2.97-3.42% across all
five held-out targets. Forcing this hydrostatic identity would therefore worsen
this independent comparison; zero constructed residual is not an admission test.
The held-out central Pg error is0.6857%, a different quantity from column error.

At actual solar5772/4.438, log-column mass closure is0.1536%, adjacent pressure
discrepancy3.8467%, compared with density-integral0%/4.1151% and pressure-derived
3.6571%/roundoff. There is **no independent native truth** at that solar target.
Its native pressure discrepancy baseline prevents claiming that this3.85% is
caused by interpolation.

Across five solar10K points, log-column mass closure spans0.1233-0.1766%; adjacent
pressure discrepancy3.837-3.857%. At identical T/Pe, producer TS gas pressure
differs from supplemental local Pg by0.1946-0.9285%, rho by0.1882-0.9011%.
At solar these are0.6039% and0.5812%; at native reference0.05309% and0.06660%.
EOS outputs remain separate; no producer rho/Pg replaces local fields.

## Official interpolation source discovered and receipted

After holdings/local-source checks found no utility copy, retrieved the archive
linked from the official [software page](https://marcs.astro.uu.se/software.html):
https://marcs.astro.uu.se/documents/auxiliary/interpol_marcs.tar.gz.
SHA2568a3d0fbd15a378830ece79c8cdc3c8ee60d74f2fd385f03b1bb3fbecb762f189.
It is safely extracted by a committed intake script; **no upstream script/program
is executed**. Full source/manual/test assets stay outside git. Redistribution
license is unverified; there is no upstream-code/manual vendoring or publication.

Manual text extraction is garbled; rendered PDF2 confirms Thomas Masseron,
title and30May2008 date. PDF3 introduction/limitations was visually inspected.
The bundled Fortran header lists later updates through2012. The existing source
is not an atmosphere-equilibrium solver: its header describes common-optical-
depth resampling, component interpolation and empirical weighting. It does not
directly provide every Bride input such as a jointly qualified retained density.

Importantly, inspected interpol_modeles.f lines1316 onward reconstruct RHOX via
Rosseland optical depth and opacity, integrating1/kappa over tau with native
parabolic quadrature and top anchor tau0/kappa0. It does **not** derive column
from gas pressure. This supplies a concrete independent next comparison, without
pretending pressure-forcing repairs the model. It still requires kernel/units/
top-anchor/native replay tests and then held-out checks; official availability
or empirical calibration alone is not physical validation.

Resource experiments.bride_marcs_interpolator points to the immutable intake.
Holdings registry and bibliography include the source/manual. The library manual
copy is byte-identical to the archive document;191 documents/127 accounted,64
unindexed remains existing debt. Source page regenerated with213 bibliography
rows; not a publication or clean-library/link claim.

## Verification, disposition and next task

51 focused tests pass, including analytic pressure slope with nonzero intercept,
wrong-gravity counterexample, invalid mass/gravity, SI conversion and native row
identity. Three known LW warnings remain. No O solve or new opacity generation.
Final pressure audit experiments.bride_mac_grid/rya1281-pressure-columns-v2;
initialv1 retained. Four compact byte-identical review files are registered under
data/results/rya1281/pressure_columns_v1, including source-intake receipt.
Source intake is at experiments.bride_marcs_interpolator. An initial attempt to
inspect the archive before curl completed hit EOF; it was not accepted/extracted.
The completed archive was hash-checked before the successful immutable intake.

Next isolate the official calcrhox kernel, verify an analytic constant-opacity
case and native mass replay, then compare its tau/kappa columns against density
geometry and held-out native RHOX. Keep the pressure-contract anomaly open and
retain the independent top-anchor caveat. Only after an explicit experimental
column choice and fresh bindings connect a provider/capped solar O control.
No threshold is fitted from these comparisons. Native-only provider refusal,
current certificates and scientific HOLD remain. No merge/publication/adopted
change; O then C/N, Fe, another star; full3D/all-band gates separate.
