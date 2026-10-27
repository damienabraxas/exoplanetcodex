# Oxygen population audit and additional solar diagnostics

Ryan authorized the matched LTE/NLTE audit, collision/background investigation,
blend checks and further oxygen diagnostics before testing other solar elements
and only then moving to another star. This phase is a Mac-only, fixed-input
experiment. Scientific admission remains **HOLD**; no source atom, observed spectrum,
gold, adopted abundance, production treatment or published product changed.

## What was actually calculated

Two existing, converged atmospheres/population sets were replayed: native MARCS
5750 K / logg 4.5, and the experimental full-field solar interpolation at
5772 K / logg 4.438. Both retain xi = 1 km/s and A(O) = 8.70. For each, form
disk-integrated flux from five angular rays using the same LW 0.16.1 context,
full held O atom, fixed electron density, geometry, boundaries and background.
Compare local LTE and local NLTE populations with no SE update. Each gets its
own bound-bound-free continuum retaining that mode's bound-free populations.
Radiation converges to dJ < 1e-8; population identity is asserted after replay.

At the exact native reference node only, also apply the independent MULTI
departure pattern to the local LTE populations. Explicitly renormalize each
depth to conserve local total oxygen. Factors span 0.999999919 to 1.000000073.
This preserves relative departure ratios, but is **not** an identical-EOS MULTI
synthesis or an independent absolute-population validation. Do not transport
those reference departures to the interpolated solar atmosphere.

The first development replay is preserved externally as population-audit-v1.
The final reproducible package is population-audit-v3, generated from committed
code, with reviewed population-file hashes, atmosphere hashes, atom hash,
environment, canonical line inventory and observation identity verified/recorded.

## Triplet result: population physics is a substantial part of the discrepancy

Intrinsic absorption areas in the same fixed +/-0.5 Angstrom apertures, mA:

| Native 5750/4.5 control | 7771.944 | 7774.166 | 7775.388 |
|---|---:|---:|---:|
| Kitt Peak corrected observed windows | 77.56 | 67.20 | 53.79 |
| Local LTE | 60.58 | 52.53 | 41.49 |
| Local NLTE | 61.31 | 53.08 | 41.84 |
| Independent b pattern on local LTE | 69.12 | 59.50 | 46.47 |
| CE rates x0.1, fresh SE solve | 62.19 | 53.91 | 42.61 |
| CE rates x10, fresh SE solve | 61.19 | 52.97 | 41.73 |
| CH rates x0.1, fresh SE solve | 64.92 | 55.77 | 43.36 |
| CH rates x10, fresh SE solve | 60.82 | 52.74 | 41.66 |

The local NLTE enhancement over LTE is only 0.84-1.21%. The independent departure
pattern increases absorption by 11.1-12.7% relative to local NLTE, while retaining
the same line-formation setup. It still does not fill the observed windows.
This isolates a significant population-pattern effect; it does not uniquely assign
the original cause to collisions, UV background, photoionization or atom reduction.

For example, at native logtau500 = -1.2621, the first triplet lower/upper b are
1.00884 / 0.99283 locally, versus 1.03018 / 0.84559 in the independent reference.
The upper/lower ratio differs substantially even where the lower population is
quite similar. The depth plots are population diagnostics, **not** contribution
functions or claims about a uniquely determined formation depth.

Four fresh full-atom SE solves tested CE and CH separately by predeclared,
symmetric factors 0.1 and 10. All converged with the unchanged joint 1e-6
population/J/physical-rate criterion and 1e-10 conservation gate. Maximum final
physical-rate residual is 1.83e-14. Summed case wall time is 252.37 s. These are
causal interventions, not error bars or proposed rate calibrations. The larger
CH response does not justify changing source collision data to improve a fit.
CI and CH0 were not varied. Source rate arrays are copied, never modified, and
altered-rate cases are rejected from the fixed-atom parameter-grid analysis.

At interpolated solar parameters, local LTE gives 62.72 / 54.55 / 43.31 mA and
local NLTE gives 63.57 / 55.19 / 43.71 mA. These are intrinsic profiles; the
previous 63.54 / 55.17 / 43.70 numbers included canonical broadening. The tiny
fixed-aperture difference is expected and is not a changed source solution.
Replayed NLTE profiles match the saved intrinsic curves to approximately 1e-13
in normalized flux.

## Additional oxygen diagnostics: calculated, with blends still owed

The held atom contains [O I] near 6300 Angstrom and the permitted 6156-6158 group.
No new model atom was assembled or truncated. Same full-atom populations were
used for these emergent spectra; the atom's internal radiation grid already
contains these transitions. This extends diagnostic wavelength coverage, not
full-band validation.

| Solar 5772/4.438 diagnostic | LTE area, mA | NLTE area, mA | Observed window, mA |
|---|---:|---:|---:|
| 6300.304 +/-0.3 Angstrom | 4.14727 | 4.14733 | 8.95460 |
| 6158.19 +/-0.5 Angstrom | 3.97654 | 5.98163 | 63.23047 |

These are **oxygen-only model versus blended observed windows**, not recovered
oxygen EWs or abundance discrepancies. The 6300 line is almost insensitive to
these local LTE/NLTE choices, unlike the permitted group. The 6158 aperture
includes strong neighbouring absorption; treating its entire 63 mA as oxygen
would be scientifically wrong.

The committed canonical atomic inventory was inspected by wavelength. The 6300
aperture contains 20 candidate rows, including canonical Ni I 6300.342 at
loggf = -2.110 and O I 6300.304 at -9.717. The 6158 aperture contains 46 rows,
including Fe I 6157.728; the three disjoint triplet apertures are inventoried
individually in report.json. Counts alone do not establish line strength,
reliability or blend completeness; no molecular inventory is included.
The source O gf values and nearby canonical O records are listed separately
with their differences, without overriding either source. See also the existing
[6300 provenance audit](../rya367_oi6300.md).

Only oxygen is active/detailed static in the shared constructor. The native
passive species contribute LW background; their presence does not install a
complete neighbouring atomic bound-bound or molecular spectrum. Installed
RadiativeSet.compute_wavelength_grid constructs explicit transitions from its
active/detailed-static sets. A full blend treatment therefore needs genuine
mixed opacity, consistent LTE partitions/EOS and canonical atomic/molecular
line data. Multiplying independently normalized spectra is not an adequate
substitute. No blend-subtracted oxygen abundance is reported here.

## Observation, reproducibility and review

Same corrected Kitt Peak 2005 holding and canonical reader as the previous
phase, observation SHA256 75dd93479b5d3fe63c2e7828f5dbb3ce29f4c672f4f1147e0d30ce41bc8f5adb.
Native vacuum nm converts once to standard-air Angstrom; engine frequencies
remain vacuum. Preserve delivered gravitational shift and provider normalization/
telluric correction. No continuum, RV, broadening or abundance optimization.
These new intrinsic overlays are not fitted line-shape comparisons. No pixel
uncertainty/covariance model, chi-square acceptance or calibrated uncertainty.

Commands from the review worktree, using registered Mac paths:

```sh
source config/hosts/mac.env
nice -n 15 "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/run_bride_mac_checks.py \
  --plan config/experiments/rya1281_mac_collision_controls.json \
  --output "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-collisions-v1"
nice -n 15 "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/audit_bride_oxygen_populations.py \
  --reference "$CODEX_WORK_ROOT/bride_inputs/rya1281/reference-5750-4.5-A8.70-v1.json" \
  --collision-runs "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-collisions-v1" \
  --output "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-population-audit-v3"
```

Output directories must be new; choose new suffixes for reproduction. Four SE
controls ran serially at reduced priority with one LW thread. The preliminary
replay's sandbox denied nice; it completed at normal priority. The final replay
and SE controls used approved reduced-priority execution. The only Sirius work
was read-only location/transfer of the existing 72 KB reference manifest, whose
hash matches the earlier experiment receipt; no Sirius science job was started.

Validation: 167 focused LW/parser/numerical/reference/observable/control tests
passed, with 666 known dependency deprecation warnings; 3 installed-iSpec
policy/effects tests passed with 4 warnings. The new replay explicitly verifies
cached nstar, unchanged populations, source identities, finite flux and observed
coverage. The synthetic population-shape test verifies conservation and relative
ratios; the collision test verifies class isolation, immutability and invalid-input
rejection. An initial roundoff assertion used less than one floating-point epsilon;
it was corrected to a four-epsilon arithmetic bound without changing science.
The real parameter-grid assembler was also invoked on the collision-only run set:
it rejected altered rates with the expected exit 1 before emitting a grid. All
27 archived result files have generator/invocation entries; replay copies match
their originals byte for byte, and five figures were packaged. The global
generator/invocation guard remains at 2436 violations, matching the previous
phase and pinned main baseline; this phase adds zero violations.

Standing bibliography audit is red: 186 documents, 125 accounted for, with
unindexed/duplicate documents reported. No new literature claim relies on those
unindexed files; library cleanup is outside this experiment. Linear's advertised
get_agent_context endpoint returned "Tool ... not found"; issue/comments and
committed sources supplied task context. Previously documented global generator
and stale six-band failures remain open. No full-CI green claim; final evidence
push skips Sirius CI as in the preceding phase. New run-bug-ledger disposition:
none discovered; root cause of the existing population discrepancy remains open.

## Next concrete work

1. Compare per-process physical rate matrices and driving radiation for the
   triplet/cascade levels against native MULTI conventions. Trace CE/CH/CI/CH0,
   photoionization table semantics, LTE normalization, UV background and boundary
   assumptions independently; the current result does not assign a winner.
2. Implement qualified mixed-opacity LTE blends at fixed canonical abundances,
   beginning with 6300 O+Ni and the 6158 neighbourhood; keep an oxygen-only
   switch and demonstrate additive-opacity behavior, units and partition consistency.
3. Assess local telluric residuals, observed continuum provenance and instrumental
   line shape for each new diagnostic using ORCH holding policy. Any more complete
   comparison must distinguish blends, population effects and atmospheric effects.
4. Validate a native/canonical second solar element and its LTE/NLTE controls
   before selecting a different star. Additional temperature/gravity holdouts,
   independent solar atlases and genuine 3D remain separate evidence gates.

The next highest-value question is why the local upper/lower departure pattern
is much closer to LTE than the independent reference under the nominally same
node. More temperature tuning would not answer it.
