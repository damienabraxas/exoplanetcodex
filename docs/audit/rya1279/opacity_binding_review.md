# Exact local-input coefficient bindings

RYA-1281 / parent RYA-1278, 2026-11-01. Mac-only certificate checks, no RT/SE.
Base3fcc5defefdb918d3748df08eeafc734ae8b2447; fetched main f8c723a8 unchanged.

## What is bound

The six previously verified TS coefficient sets now have immutable local-input
certificates. They contain the complete static 1D structure: temperature,
electron pressure, rho, outward height, both depth coordinates, native abundance
vector, hydrogen mass fraction, microturbulence, gravity, depth count and model
identity. The active diagnostic oxygen abundance is separate from the producer's
native mixture. The native certificate uses the same declared diagnostic oxygen
value as the existing solar plan; no producer abundance is changed.

Derived ne and total hydrogen density use the existing shared pinned-engine
unit/constants boundary and are included in the binding. Static zero vlos is
explicit; unmodeled structure fields (including a supplied velocity field) are
refused rather than silently ignored. Canonical JSON and SHA256 retain exact
numeric identity. No comparison tolerance or float rounding is introduced.

Native column metadata is bound when present. Its absence in the five solar
interpolants is explicit. Adding a constructed column later invalidates that
certificate: an opacity certificate cannot become cap approval by adding a field.

Each certificate also binds the completed producer case, full coefficient-array
hash and producer-report hash. Requested T/Pe and local rho in that array must
match the local structure exactly; producer native composition must match its
underlying local mixture. Nonfinite evidence and invalid per-neutral coefficient
shape/wavelength order are refused. Producer density remains separate evidence;
there is no local-density replacement or EOS repair.

## Verification

All six exact identities pass. Ninety-four controlled mismatches are refused:
- 6 array fields changed by one representable float step per case;
- microturbulence, hydrogen fraction and gravity changed by one float step;
- changed Fe composition and active O abundance;
- changed/added column, depth cut, unmodeled velocity;
- legacy provider schema refusal for each certificate;
- four neighboring 10 K inputs offered another temperature's certificate.

The new format is deliberately **not** a background-provider packet. The existing
native Fe packet loader rejects it by schema. Only local-neutral normalization
is declared for future use; raw producer mass coefficients remain diagnostic
evidence. Provider/cap admission flags must remain false and cannot enable either
path. The existing provider, original native packet, interpolation and top-cap
implementation are unchanged.

Seventeen focused binding tests plus the29 preceding input/adapter tests pass
(46 total). These include appended bytes that a ZIP/NPZ reader could otherwise
ignore, a modified source report, a changed binding digest and an attempted
provider-admission flag change. All are rejected. Three actual native-provider
backend regressions also pass separately. No new verified run defect or
scientifically material anomaly was discovered; no fabricated bug-ledger row.

The final immutable audit is
experiments.bride_mac_grid/rya1281-bound-point-coefficients-v3. Earlier v1/v2 are
retained development checkpoints; finalv3 adds normalization of accepted input
array representation, explicit admission-flag refusal, canonical plan-derived
native diagnostic O and producer/local composition checking. It does not change
any producer coefficient, thermodynamic field or scientific result.

Eight byte-identical compact review artifacts (report, checks and six manifests)
are archived at data/results/rya1281/opacity_binding_v1. Full source arrays and
source reports remain in the immutable certificate directories in the registered
Mac store. The builder and checking logic are committed.

## Remaining boundary gate

An exact match establishes provenance, **not physical accuracy**. Native replay
identity from the preceding phase and small nearest-cell method spread remain
useful controls, but do not certify interpolated equilibrium or boundary columns.
Next check pressure/EOS/geometry/column jointly against native MARCS tables and
retained producer EOS outputs before selecting any explicit capped diagnostic
interpolant. Rebinding a changed column must be a visible new operation, not a
fallback. Only then connect an explicitly experimental provider and run matched
native/solar O controls with fixed abundance and wavelength sampling.

Scientific HOLD remains. No new O solve, abundance fit, observation change,
gold/adopted change, merge or publication. O then C/N, Fe, another star; full3D
and all-band validation remain separate gates.
