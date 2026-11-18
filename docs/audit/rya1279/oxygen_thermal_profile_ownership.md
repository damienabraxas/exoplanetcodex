# Oxygen thermal profile: measured ownership and implementation contract

RYA-1371, 2026-11-15, v293. Read-only qualification evidence; no adopted profile.
The held case is 5750 K/logg 4.5/A(O) 8.70, not canonical solar validation.

| Component | Current owner and convention | Required coherent-profile control |
|---|---|---|
| Background extinction and colliders | Pinned HOP EOS/opacity via `hop_packet.compute_packet`; active O excluded; source abundances/options and prescribed T/ne/rho | Retain forward opacity/EOS inputs; declare EOS, partition, ionization lowering and collider conventions before claiming a new equilibrium state |
| Background thermal emission | `eta=(chi-sca)*B_HOP`, historical h/k | Match the selected Planck profile; test absorption separately from total extinction |
| Background scattering | `sca` added by backend as `sca*J`, exactly once | At J=B, test `eta+sca*B=chi*B`; never interpret scattering as thermal absorption |
| Active BF | Compiled `LwTransition.hpp` shared alpha/U/V; held LW h/k, imported native nStar; optional native Hermite alpha | Thermal inverse coefficient, spontaneous emission and rate integration must change together; alpha alone is already isolated |
| Active BB | LW atomic Einstein relations and line-centre CRD inverse ratios; source width prescription and imported nStar | Declare finite-profile thermodynamic behavior; a line-centre identity does not imply exact frequency-dependent Planck balance over a finite CRD profile |
| Oxygen LTE populations | Native `atom_mod.compute_lte_populations/lte_pop`; consumed energy table, FP32 thermal quantities, source partition/Saha conventions | Preserve energy/table-rounding as distinct factors; same-stage Boltzmann reconstruction does not qualify full Saha/partition/EOS |
| Collisions | Pinned native forward oracle and imported-LTE reverse rates; fixed native collider fields | Preserve source forward rates and collider ownership; regenerate inverse relations only under a qualified common LTE convention |
| Thermal lower boundary | Compiled `FormalScalar.cpp` calls held-LW `LwMisc.hpp::planck_nu` on bottom two depths, then optical-depth gradient correction | Boundary Planck profile must match the declared emission profile; boundary intensity is not simply B and needs an actual formal-solve control |
| Upper boundary | Explicit source-matched zero ghost control | Preserve geometry and vacuum convention; do not fit boundary parameters |

## Measurements on the actual rebuilt packet

`thermal_profile_v1/rya1371-thermal-profile-contract-v2/report.json` covers
19,963 wavelengths × 56 depths = **1,117,928 cells**. It rebuilds the actual
context/background but performs no radiation or SE iteration. Original one-thread,
600-second and 20,000-node controls apply. Rebuilt chi/eta/scattering hashes are
recorded; equality to arrays from an earlier run is not asserted because those
arrays were not saved then.

- Native-Planck background absorption defect is zero by construction; total
  including scattering balances to 1.36e-16. This verifies packet wiring, not an
  independent derivation of HOP opacity physics.
- With held-LW Planck, maximum normalized absorption defect is 2.34811e-4.
  The maximum relative Planck difference is 4.69511e-4 at 90.0064 nm vacuum,
  top depth, T=4045.89892578125 K. Bottom-two-depth Planck differences reach
  1.88736e-4. These are local equilibrium input comparisons, not a measured
  boundary-solve or observational error.
- Maximum scattering fraction is 0.41537. Both component and total defects are
  reported, avoiding suppression of absorption errors in scattering-rich cells.
- Maximum exponent is 39.5101, below HOP's 80 Wien switch and LW's 150 cutoff.
  The tested discrepancy therefore does not arise from their different branches.
- Own and native collision matrices balance at imported nStar to 2.44e-16 and
  2.21e-16 over 43,680 unordered pair/depth edges, with no one-way edges.
  **Level 41 has no collision edges in either source-qualified matrix.** This is
  explicit missing collision coverage, not proof that collisions are physically
  negligible, and not authorization to invent rates. Its only radiative connection
  remains continuum 4->41.
- The same audit independently reconstructs all 31,136 positive radiative
  directional/depth rates to 2.86e-14. Integrated BB and BF canaries remain
  separate: BB/native-Planck/imported-nStar defect 9.35e-6, BB/LW/imported-nStar
  2.25e-4, BF/native-Planck/held-LW factors 2.33e-4. These component maxima are
  not additive, and ordinary NLTE residuals need not vanish at the actual field.

## Implementation order and gates

The next bounded experiment should isolate a **declared emission-profile
control** at fixed source populations and coefficients: background thermal
emission and lower-boundary thermal inputs must use the same declared profile,
with scattering, chi, collider fields and original default reproduction held.
It is a diagnostic axis, not a complete coherent model. Its fresh RT and rates
must be independently checked; a Planck substitution in one component cannot be
described as thermal closure. Source-native and held-LW profiles remain explicit
alternatives, with no selection by smaller abundance/residual.

Before a coherent coupled profile, jointly qualify active BB/BF emission and
inverse-rate relations, energy/partition/Saha/LTE ownership, collision reciprocity
and boundary behavior. Native split BB canaries balance by construction at J=B;
matching that split in rates alone does not qualify active opacity/emission.
Any backend implementation needs source/build fingerprints and default replay
reproduction, followed by fixed-state physical rates, equilibrium canaries,
particle conservation, numerical refinement and coupled closure/restart.

No constants monkeypatch, synthetic collision repair, source-table edits,
level removal, fitted thermal factors, gold/default/production change or
scientific admission is authorized by these diagnostics. UV/VIS/IR observational
and canonical-solar/ORCH continuum/telluric/geometry gates remain open. Fe is next
after Oxygen foundation; C/N and RYA-1374 retain their existing gates.
