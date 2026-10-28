# Oxygen rate conventions and first mixed-opacity Ni test

> 2026-10-28 follow-up: the [native charge audit](native_charge_review.md) finds the held lowercase RCOLL/GENCL implementation missed by the earlier search, tests its compiled spline and a reversible-charge hypothesis, and preserves the unresolved outer-layer/producer gates.

Continuation of RYA-1281 at e75e2b94, on the existing stacked draft PR #598.
Ryan authorized continuing the oxygen investigation on the Mac. Refreshed main
remains 0ade2170. All scientific admission gates remain **HOLD**.

## Finding: a collision label is not a complete physics contract

The earlier experiments used the generic converter's RH/Lightweaver conventions.
The new per-process audit exposed a source-convention hazard, now recorded as an
**open collision_convention defect** in data/audit/run_bug_ledger.csv:

- LW 0.16.1 TemperatureInterpolationRates.setup sorts endpoints. Its CH0 class
  applies a downward-only rate. The held atom's CH0 rows are ordered, e.g. 6,40.
  Sorting loses that raw source/destination information.
- The inspected [legacy GENCOL source](https://github.com/sjbruce3301/FORWARD-coronal-modeling/blob/15145239ffdbb62d19dca0749b11c5920f9fdf13/idl/FORCOMP/NOTUSEDBYFORWARD/EXTRAS/gencol.f)
  explicitly makes CH0 directional in input order, uses linear interpolation of
  log(rate) in temperature, and defines CH coefficients in the downward direction.
- [RH collision.c](https://github.com/han-uitenbroek/RH/blob/c95bdf66ff6bc1c6cea0ed04dc688100170e8ee6/collision.c)
  defines CH coefficients upward; its
  [MULTI converter](https://github.com/han-uitenbroek/RH/blob/c95bdf66ff6bc1c6cea0ed04dc688100170e8ee6/tools/modelatom.c)
  also reorders CH0. The current LW mapping follows those RH conventions.

These sources establish different conventions. They do **not** identify the exact
modified MULTI executable that produced the held departure grid. The inspected
GENCOL copy has a 1999 adaptation header and is not claimed to be that executable.
No new convention is adopted merely because one spectrum is closer to an atlas.
The default converter remains unchanged; explicit diagnostic adapters live in
pipeline/bride_gencol_controls.py. A complete intake contract must declare the
producer's coefficient direction, collider, reverse process, interpolation and
units before using a table. This matters directly to Bride's arbitrary-grid goal.

The new controls preserve the full atom and source values, changing named rate
interpretations only. CH still uses the existing WENO interpolation, so even the
CH-down control is not an exact GENCOL numerical reproduction. Neither a paper
attribution nor an atom filename supplies the missing executable identity.

## Physical rate/radiation audit

scripts/audit_bride_oxygen_rates.py replays local NLTE, local LTE and the normalized
independent departure pattern at native MARCS 5750/4.5/A(O)8.70, with no SE update.
Three radiation solves finish in about 14 s summed mode time on the Mac.

Each native collision class is independently accumulated and compared with the
backend's C matrix. Individual radiative transition rates are accumulated separately;
the full physical destination/source R+C operator is independently reconstructed
and compared with physical_rate_matrix. Signed and gross flows have units m^-3 s^-1.
The accelerated Gamma matrix is never used as a physical rate budget.
Native Jnu, Bnu, populations and full rate arrays remain in the immutable external
Mac store, hash-referenced by the committed report. Selected depth budgets and
process identities are in the review package.

For the original local NLTE populations, radiation replay leaves maximum physical
imbalance 4.77e-7, consistent with the source's joint 1e-6 stopping criterion; this
is a fresh-radiation check, not the lagged-field residual after the original SE step.
LTE and external populations are not equilibria of this local operator: their
maximum relative imbalance is nearly one. That is a diagnostic, not a failed
claim that the foreign populations should equilibrate under different physics.

At logtau500=-0.9737, upper triplet level 9 (zero-based), the original local NLTE
flow budget normalized by total gross is:

| Process | Signed net / total gross |
|---|---:|
| 777.408 nm vacuum radiative transition | -0.00323935 |
| CH class | +0.00283028 |
| 926.855 nm cascade | +0.00026144 |
| 1130.547 nm cascade | +0.00007382 |
| CE class | +0.00000533 |
| Total over every process | +4.46e-12 |

This identifies channels that balance this specific local solution. It does not
rank uncertainty or establish the correct source interpretation. In the external
population replay, the original CH0 interpretation instead dominates the net
flow at this depth; inspecting that discrepancy led to the convention audit.
No native MULTI rate/J dump was available for a direct process-by-process parity test.

## Convention trials, including every failure

Both trials retain all levels/lines/continua and use the native reference node.
The first changes CH0 source direction and log-rate interpolation; the second
also changes CH to downward coefficients. No abundances or coefficients are fitted.

1. **v1, original Ng settings:** both stop after 18 valid history records when an
   accelerated update produces invalid populations. No spectra emitted.
2. **v2, no Ng:** both remain physical but fail the 1500-iteration joint-convergence
   cap. Each lasts about 100.6-100.7 s. No spectra emitted or mislabeled converged.
3. **v3, safeguarded Ng:** reject an invalid accelerated proposal, restore the
   previous populations, clear Ng history, and solve the same equation once more.
   The retry must be finite/positive; otherwise stop. No clipping, floors, physical
   rate change or tolerance relaxation. Both converge: 887 and 392 iterations,
   60.73 and 26.84 s, each rejecting one invalid proposal. Joint population/J/rate
   threshold remains 1e-6, conservation 1e-10; terminal physical residuals ~3e-15.

The safeguard is opt-in. The old solver behavior and old failed outputs remain.
The fixed-atom parameter-grid assembler refuses both altered-rate scales and
experimental convention variants. Synthetic tests discriminate endpoint reversal,
upward/downward CH interpretation, out-of-domain interpolation and invalid retry.

Intrinsic fixed-aperture triplet areas, mA:

| Native 5750/4.5 control | 7771.944 | 7774.166 | 7775.388 |
|---|---:|---:|---:|
| Original RH-convention local NLTE | 61.31 | 53.08 | 41.84 |
| Independent b shape on local LTE | 69.12 | 59.50 | 46.47 |
| Directed CH0 + log-rate interpolation | 59.86 | 51.66 | 40.58 |
| Above plus downward CH | 71.83 | 62.27 | 49.24 |
| Kitt Peak observed blended windows | 77.56 | 67.20 | 53.79 |

The final row of model controls is closer in area but still fails population parity.
At logtau500=-1.2621, its lower/upper b are 0.96177/0.77597 versus reference
1.03018/0.84559. Its upper-level maximum absolute difference is 0.129 over
-3<=logtau500<=0 and 1.25 over all depths; maximum over every level/depth is 12.21.
The CH0-only variant likewise does not reproduce the reference. These checks
explicitly prevent selecting a convention from improved observed agreement.

## First O+Ni calculation with actual mixed opacity

pipeline/bride_lte_blend.py adds trace LTE Ni absorption/emission to native LW
background coefficients. Oxygen line opacity and Ni line opacity enter the same
formal-transfer calculation. All components normalize to one true continuum with
both bound-bound contributions removed. No multiplication of separate residual
profiles is used to produce the combined spectrum.

Ni ionization and partition functions come from the installed Kurucz table at
fixed MARCS ne/nH. A(Ni)=6.20 is read from SOLAR_ASPLUND2021; the canonical
combined Ni6300 line supplies wavelength, excitation 4.266 eV and loggf=-2.110.
Natural damping 1.86e8 s^-1 comes from the explicitly identified installed-iSpec
GES record. The registered data.ges_v6_atomic_lines path is absent on the Mac;
the held installed-iSpec asset under engines.ispec is used explicitly and separately
fingerprinted, with no silent path fallback or source edit.

**This is a natural-Voigt prototype.** Pressure broadening, isotope splitting,
other atomic/molecular blends and a new atmosphere/EOS solution are not implemented.
The original RH-convention oxygen populations remain frozen, so these blend tests
do not resolve the collision-convention issue or constitute a new SE solution.

Areas in the fixed 6300.004-6300.604 Angstrom aperture, mA:

| Model/observation | Native 5750/4.5 | Solar interpolation 5772/4.438 |
|---|---:|---:|
| Oxygen only | 3.94170 | 4.14733 |
| Nickel only | 2.00188 | 1.98792 |
| Mixed oxygen + nickel | 5.83687 | 6.02392 |
| Kitt Peak blended window | 8.95460 | 8.95460 |

The mixed calculation is not the sum of independent absorption areas. Its
normalized spectrum differs from the product of separate normalized spectra by
up to 0.00113/0.00118. That is direct evidence for treating shared opacity and
formation depths rather than combining residual curves. Remaining observed area
is not assigned uniquely to oxygen or used to fit an abundance.

The opacity-kernel tests verify oscillator-strength integrated extinction with a
finite-wing error bound, Kirchhoff emission, and zero-strength behavior. Real
replays assert immutable O populations, convergence and positive finite flux.
These checks validate the prototype's mechanics, not complete Ni/O solar physics.

## Reproduction and scope

All jobs are finished, serial one-thread reduced-priority Mac work. No new Sirius
science job or dependency installation. Raw grids/spectra, gold, adopted values,
production tokens and published surfaces are unchanged. No merge/publication.

```sh
source config/hosts/mac.env
nice -n 15 "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/audit_bride_oxygen_rates.py \
  --reference "$CODEX_WORK_ROOT/bride_inputs/rya1281/reference-5750-4.5-A8.70-v1.json" \
  --output "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-rate-audit-v1"
nice -n 15 "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/run_bride_mac_checks.py \
  --plan config/experiments/rya1281_mac_gencol_controls.json \
  --output "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-gencol-v3-safe-ng"
nice -n 15 "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/check_bride_ni6300_blend.py \
  --output "$CODEX_WORK_ROOT/bride_mac_checks/rya1281-ni6300-v2"
```

Use new suffixes; existing directories cannot be overwritten. Earlier failed plans
are recoverable from source commits named in their summary: v1 used default Ng,
v2 used no Ng; current config is the safeguarded v3. The summary generator records
all three sets, not just successes. Full external arrays/history retain source hashes.

The corrected Kitt Peak standard-air/provider-normalized holding is unchanged.
No telluric/continuum/RV/width or abundance fit. New overlays are intrinsic and
neither a full observed-window model nor a covariance/uncertainty budget.
Standing bibliography audit remains red (186 documents, 125 accounted for);
source-code evidence uses the inspected pinned URLs above. No unindexed paper
is used to choose a convention. Existing global generator/stale-six-band and
generated-surface aging issues remain visible. Evidence push skips Sirius CI.

## Validation and review

172 focused LW/parser/numerical/reference/observable/control tests pass, with
666 known upstream deprecation warnings; 3 installed-iSpec policy/effects tests
pass with 4 warnings. Kernel sum-rule/Kirchhoff/zero-opacity controls and actual
convention convergence/positivity/conservation checks pass. Both initial invalid
Ng trials and both capped unaccelerated trials remain failed, not silently retried
in place. Review copies include 27 registered result artifacts and six figures;
full 15 MB rate/radiation arrays remain externally hash-referenced. Plots were
visually reviewed; the final Ni rerender reproduces v1 numeric results exactly.
The selected-flow CSV review copy normalizes CRLF to LF for repository whitespace;
all parsed rows equal the immutable original. Future generator output uses LF.
Run-bug ledger has one new open convention-contract row. No global CI-green claim.
The global generator/invocation guard remains at 2436 pre-existing violations,
matching the prior phase/pinned main baseline: zero added. Register v180, its
SEQUENCE entry, artifact checks and whitespace checks pass.

## Next work

1. Identify the actual grid producer's modified GENCOL/collision implementation,
   coefficient direction and reverse-process handling, interpolation, background
   and atom version. Obtain a small native rate/population/J dump at this exact node.
   That evidence, not a better atlas fit, decides the validated intake convention.
2. Turn declared source conventions into typed atom-adapter contracts, with thermal,
   direction and per-process parity tests. Keep unsupported producer semantics out
   of production. The new bug-ledger row stays open until that is established.
3. Extend Ni from this natural-Voigt prototype to justified pressure/isotope handling
   and other qualified 6300/6158 blends. Recheck observational residual tellurics,
   continuum provenance and LSF under ORCH holding policy before inference.
4. Additional solar elements precede another star. Adaptive 10 K interpolation,
   other bands and true full-3D transfer retain their separate validation gates.
